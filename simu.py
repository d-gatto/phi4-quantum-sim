from qat.lang.AQASM import Program, QRoutine, X, H, PH, RX, RY, RZ, SWAP
from qat.lang.AQASM.qint import QInt
import qat.lang.AQASM.qftarith as qftarith
import mpmath as mp
from qat.qpus import get_default_qpu

# Arithmetic subroutines. 
# Quantum integer square calculator
def square(bits):
    rout = QRoutine()
    qint = rout.new_wires(bits, QInt, reverse_bit_order=True)
    dummy = rout.new_wires(bits, QInt, reverse_bit_order=True)
    rout.set_ancillae(dummy)
    out = rout.new_wires(bits, QInt, reverse_bit_order=True)
    dummy += qint
    out += qint*dummy
    dummy -= qint
    return rout
# Quantum integer subtractor 
def subtract(bits):
    rout = QRoutine()
    min = rout.new_wires(bits, QInt, reverse_bit_order=True)
    sub = rout.new_wires(bits, QInt, reverse_bit_order=True)
    diff = rout.new_wires(bits, QInt, reverse_bit_order=True)
    diff += min - sub
    return rout
# Multiplication of a quantum integer by a constant
def multiply(inbits, outbits, multiplier):
    rout = QRoutine()
    qint = rout.new_wires(inbits, QInt, reverse_bit_order=True)
    qout = rout.new_wires(outbits, QInt, reverse_bit_order=True)
    for i in range(multiplier):
        qout += qint
    return rout
# This subroutine calculates the absolute value of the difference betwee two quantum integers
def AbsDiff(bits):
    rout = QRoutine()
    a = rout.new_wires(bits, QInt, reverse_bit_order=True)
    b = rout.new_wires(bits, QInt, reverse_bit_order=True)
    c = rout.new_wires(1, QInt)
    rout.set_ancillae(c)
    diff = rout.new_wires(bits, QInt, reverse_bit_order=True)
    (a<b).evaluate(c)
    rout.apply(subtract(bits).ctrl(), c, b, a, diff)
    rout.apply(X, c)
    rout.apply(subtract(bits).ctrl(), c, a, b, diff)
    rout.apply(X, c)
    (a<b).evaluate(c)
    return rout
# This subroutine calculates the integer part of a quantum integer multiplied by a constant real number
def IntPartMult(bits, multiplier):
    rout = QRoutine()
    qint = rout.new_wires(bits, QInt, reverse_bit_order=True)
    short = rout.new_wires(bits, QInt, reverse_bit_order=True)
    long = rout.new_wires(2*bits, QInt, reverse_bit_order=True)
    rout.set_ancillae(short)
    rout.set_ancillae(long)
    qout = rout.new_wires(bits, QInt, reverse_bit_order=True)
    if multiplier>0:
        rout.apply(multiply(bits,bits,int(mp.floor(multiplier))), qint, qout)
        rout.apply(multiply(bits,2*bits,int(mp.floor(2**bits*(multiplier-mp.floor(multiplier))))), qint, long)
        for i in range(bits):
            rout.apply(X.ctrl(), long[i+bits], short[i])
        qout += short
        for i in range(bits):
            rout.apply(X.ctrl(), long[i+bits], short[i])
        rout.apply(multiply(bits,2*bits,int(mp.floor(2**bits*(multiplier-mp.floor(multiplier))))).dag(), qint, long)
    if multiplier<0:
        rout.apply(multiply(bits,bits,int(mp.floor(-multiplier))).dag(), qint, qout)
        rout.apply(multiply(bits,2*bits,int(mp.floor(2**bits*(-multiplier-mp.floor(-multiplier))))), qint, long)
        for i in range(bits):
            rout.apply(X.ctrl(), long[i+bits], short[i])
        qout -= short
        for i in range(bits):
            rout.apply(X.ctrl(), long[i+bits], short[i])
        rout.apply(multiply(bits,2*bits,int(mp.floor(2**bits*(-multiplier-mp.floor(-multiplier))))).dag(), qint, long)
    return rout
# Gaussian state generator subroutine.
# Classical computation of a specific function of the Gaussian parameters is required
def JacobiTheta(mean, var):
    return mp.re(mp.jtheta(3,mean/(mp.j*var**2),mp.exp(-1/var**2))) 
# The function is used to determine the qubit rotation angle employed in the algorithm
def KitaevAngle(mean, var):
    return mp.acos(mp.sqrt(JacobiTheta(mean/2,var/2)/JacobiTheta(mean,var)))
# The actual state preparation subroutine is defined recursively
def discr_gaussian(qbits, mean, var):
    rout = QRoutine()
    wires = rout.new_wires(qbits)
    rout.apply(RY(2*float(KitaevAngle(mean,var))), wires[0]) #Rotate leftmost qubit. The angle must be converted to float for the quantum gate to accept it
    if qbits>1:
        rout.apply(X, wires[0]) #Switch leftmost bit logical values
        rout.apply(discr_gaussian(qbits-1, mean/2, var/2).ctrl(), wires[0], wires[1:qbits]) #Complete one half of the state
        rout.apply(X, wires[0]) #Unswitch leftmost bit logical values
        rout.apply(discr_gaussian(qbits-1, (mean-1)/2, var/2).ctrl(), wires[0], wires[1:qbits]) #Complete the other half of the state
    return rout

# Upon initialisation all the simulation parameters are stored and all the necessary qubits are allocated 
class JLPsim(Program):
    def __init__(self, length, dimension, qbits, mass, coupling, spacing, wavefunction, cutoff):
        super().__init__()
# Store parameters
        self.length = length
        self.dim = dimension
        self.bits = qbits
        self.mass = mass
        self.coupl = coupling
        self.spac = spacing
        self.cutoff = cutoff
# Create a cubic matrix with one register per lattice site. If dimension=2 the "cubic" matrix has only a single square slice.
        self.lattice = []
        for i in range(self.length):
            qmatrix = []
            if i==0 or self.dim==3:
                for j in range(self.length):
                    qrow = []
                    for k in range(self.length):
                        qrow.append(self.qalloc(self.bits, QInt, reverse_bit_order=True))
                    qmatrix.append(qrow)
                self.lattice.append(qmatrix)
# Store the wavefunction flattened in lexicographic order
        self.wvfunc = []
        for i in range(len(self.lattice)):
            for j in range(self.length):
                for k in range(self.length):
                    self.wvfunc.append(wavefunction[i][j][k])
# Allocate ancillas
        self.copy_ancilla = self.qalloc(self.bits, QInt, reverse_bit_order=True)
        self.temp_ancilla = self.qalloc(self.bits, QInt, reverse_bit_order=True)
        self.prep_ancilla = self.qalloc(1)
        self.time_ancilla = self.qalloc(self.bits, QInt, reverse_bit_order=True)
        self.estm_ancilla = self.qalloc(self.bits)
# Dispersion relation function definition
    def energy(self, i, j, k):
        grad = mp.sin(i*mp.pi/self.length)**2 + mp.sin(j*mp.pi/self.length)**2 + mp.sin(k*mp.pi/self.length)**2
        return mp.sqrt(self.mass**2 + 4*grad/self.spac**2)
# Create vacuum matrix
    def build_VacuumMatrix(self):
        self.VacuumMatrix = []  
        for x1 in range(len(self.lattice)):
            for x2 in range(self.length):
                for x3 in range(self.length):
                    vrow = []
                    for y1 in range(len(self.lattice)):
                        for y2 in range(self.length):
                            for y3 in range(self.length):
                                Sigma = 0
                                for p1 in range(len(self.lattice)):
                                    for p2 in range(self.length):
                                        for p3 in range(self.length):
                                            Sigma += mp.cos(2*mp.pi*(p1*(x1-y1) + p2*(x2-y2) + p3*(x3-y3))/self.length) * self.energy(p1,p2,p3,)
                                vrow.append( (self.spac/self.length)**self.dim * Sigma )
                    self.VacuumMatrix.append(vrow)
# Calculate LDL decomposition
    def computeLDL(self):
        self.CholeskyDiagonal = []
        CholeskyTriangular = mp.eye(self.length**self.dim)
        for i in range(self.length**self.dim):
            diag = self.VacuumMatrix[i][i]
            for k in range(i):
                diag -= CholeskyTriangular[k,i]**2 * self.CholeskyDiagonal[k]
            self.CholeskyDiagonal.append(diag)
            for j in range(i+1,self.length**self.dim):
                CholeskyTriangular[i,j] = self.VacuumMatrix[i][j]
                for k in range(i):
                    CholeskyTriangular[i,j] -= CholeskyTriangular[k,i]*CholeskyTriangular[k,j]*self.CholeskyDiagonal[k]
                CholeskyTriangular[i,j] = CholeskyTriangular[i,j]/self.CholeskyDiagonal[i]
        self.InverseTriangular = CholeskyTriangular**-1
# Calculate correlation matrix
    def build_CorrelationFunction(self):
        self.CorrelationFunction = (1/2)*mp.matrix(self.VacuumMatrix)**-1
# Calculate normalisation factor
    def compute_normalisation(self):
        vec = self.CorrelationFunction*mp.matrix(self.wvfunc)
        eta = 0
        for i in range(self.length**self.dim):
            eta += mp.conj(self.wvfunc[i])*vec[i]
        self.norm = 1/mp.sqrt(mp.re(eta))
# Put lattice in discretised vacuum state
    def build_VacuumState(self):
        self.build_VacuumMatrix()
        self.computeLDL()
        for i in range(len(self.lattice)):
            for j in range(self.length):
                for k in range(self.length):
                    var = 2**(self.bits-1)/(self.cutoff*mp.sqrt(self.CholeskyDiagonal[i*self.length**2+j*self.length+k]))
                    self.apply(discr_gaussian(self.bits,0,float(var)), self.lattice[i][j][k])
                    self.lattice[i][j][k] += 2**(self.bits-1)
        for i in range(len(self.lattice)):
            for j in range(self.length):
                for k in range(self.length):
                    for h in range(i*self.bits**2+j*self.bits+k+1,self.length**self.dim):
                        matrix_el = self.InverseTriangular[i*self.bits**2+j*self.bits+k,h]
                        lattice_adj = self.lattice[h//self.length**2][(h//self.length)%self.length][h%self.length]
                        self.apply(IntPartMult(self.bits,matrix_el), lattice_adj, self.lattice[i][j][k])
# Phi-sector wavepacket preparation unitary
    def FieldWaveUnitary(self, labels, time):
        gamma = 0
        for i in range(len(self.lattice)):
            for j in range(self.length):
                for k in range(self.length):
                    gamma += self.wvfunc[i*self.length**2+j*self.length+k] * labels[i][j][k]
                    for h in range(self.bits):
                        if (labels[i][j][k]//2**h)%2 == 0:                   
                            self.apply(X, self.lattice[i][j][k][h])
        gamma *= self.norm/2
        self.apply(RZ(float(mp.arg(gamma))).ctrl(self.bits*self.length**self.dim), self.lattice, self.prep_ancilla)
        self.apply(RX(2*float(mp.fabs(gamma)*time)).ctrl(self.bits*self.length**self.dim), self.lattice, self.prep_ancilla)
        self.apply(RZ(float(-mp.arg(gamma))).ctrl(self.bits*self.length**self.dim), self.lattice, self.prep_ancilla)
        for i in range(len(self.lattice)):
            for j in range(self.length):
                for k in range(self.length):
                    for h in range(self.bits):
                        if (labels[i][j][k]//2**h)%2 == 0:
                            self.apply(X, self.lattice[i][j][k][h])
# Phi-sector wavepacket preparation overall subroutine
    def FieldWaveEvolution(self, time):
        for n in range(2**(self.bits*self.length**self.dim)):
            labels = []
            for i in range(self.length):
                intmatrix = []
                if i==0 or self.dim==3:
                    for j in range(self.length):
                        introw = []
                        for k in range(self.length):
                            introw.append((n//2**(self.bits*(self.length**self.dim-3*i-2*j-k-1)))%2**self.bits)
                        intmatrix.append(introw)
                    labels.append(intmatrix)
            self.FieldWaveUnitary(labels, time)
# Pi-sector wavepacket preparation unitary
    def ConjugateWaveUnitary(self, labels, time):
        tau = 0
        vec = self.CorrelationFunction*mp.matrix(self.wvfunc)
        for i in range(len(self.lattice)):
            for j in range(self.length):
                for k in range(self.length):
                    tau += self.spac**self.dim * vec[i*self.length**2+j*self.length+k] * labels[i][j][k]
                    for h in range(self.bits):
                        if (labels[i][j][k]//2**h)%2 == 0:
                            self.apply(X, self.lattice[i][j][k][h])
        tau *= self.norm/mp.j
        self.apply(RZ(float(mp.arg(tau))).ctrl(self.bits*self.length**self.dim), self.lattice, self.prep_ancilla)
        self.apply(RX(2*float(mp.fabs(tau)*time)).ctrl(self.bits*self.length**self.dim), self.lattice, self.prep_ancilla)
        self.apply(RZ(float(-mp.arg(tau))).ctrl(self.bits*self.length**self.dim), self.lattice, self.prep_ancilla)
        for i in range(len(self.lattice)):
            for j in range(self.length):
                for k in range(self.length):
                    for h in range(self.bits):
                        if (labels[i][j][k]//2**h)%2 == 0:
                            self.apply(X, self.lattice[i][j][k][h])
# Pi-sector wavepacket preparation overall subroutine
    def ConjugateWaveEvolution(self, time):
        for i in range(len(self.lattice)):
            for j in range(self.length):
                for k in range(self.length):
                    self.apply(qftarith.QFT(self.bits), self.lattice[i][j][k])
        for n in range(2**(self.bits*self.length**self.dim)):
            labels = []
            for i in range(self.length):
                intmatrix = []
                if i==0 or self.dim==3:
                    for j in range(self.length):
                        introw = []
                        for k in range(self.length):
                            introw.append((n//2**(self.bits*(self.length**self.dim-3*i-2*j-k-1)))%2**self.bits)
                        intmatrix.append(introw)
                    labels.append(intmatrix)
            self.ConjugateWaveUnitary(labels, time)
        for i in range(len(self.lattice)):
            for j in range(self.length):
                for k in range(self.length):
                    self.apply(qftarith.QFT(self.bits).dag(), self.lattice[i][j][k])
# Build wavepacket
    def MakeWavepacket(self, steps):
        self.build_CorrelationFunction()
        self.compute_normalisation()
        for h in range(steps):
            self.FieldWaveEvolution(mp.pi/(2*h+2))
            for i in range(len(self.lattice)):
                for j in range(self.length):
                    for k in range(self.length):
                        self.apply(qftarith.QFT(self.bits), self.lattice[i][j][k])
            self.ConjugateWaveEvolution(mp.pi/(2*h+2))
            for i in range(len(self.lattice)):
                for j in range(self.length):
                    for k in range(self.length):
                        self.apply(qftarith.QFT(self.bits).dag(), self.lattice[i][j][k])
# Phi-sector phase function
    def compute_FieldPhase(self, time):
        prefactor1 = self.spac**self.dim*self.mass**2*self.cutoff**2*time/(2**self.bits*mp.pi)
        prefactor2 = self.spac**self.dim*self.coupl*self.cutoff**4*time/(8**self.bits*3*mp.pi)
        prefactor3 = self.spac**self.dim*self.cutoff**2*time/(2**self.bits*mp.pi)
        for i in range(len(self.lattice)):
            for j in range(self.length):
                for k in range(self.length):
                    self.copy_ancilla.set_value(2**(self.bits-1))
                    self.apply(AbsDiff(self.bits), self.lattice[i][j][k], self.copy_ancilla, self.temp_ancilla)
                    self.reset(self.copy_ancilla)
                    self.copy_ancilla += self.temp_ancilla
                    self.reset(self.temp_ancilla)
                    self.apply(square(self.bits), self.copy_ancilla, self.temp_ancilla)
                    self.apply(IntPartMult(self.bits,prefactor1), self.temp_ancilla, self.time_ancilla)
                    self.reset(self.copy_ancilla)
                    self.copy_ancilla += self.temp_ancilla
                    self.reset(self.temp_ancilla)
                    self.apply(square(self.bits), self.copy_ancilla, self.temp_ancilla)
                    self.apply(IntPartMult(self.bits,prefactor2), self.temp_ancilla, self.time_ancilla)
                    self.reset(self.copy_ancilla)
                    self.reset(self.temp_ancilla)
                    if self.dim==3:
                        self.apply(AbsDiff(self.bits), self.lattice[i][j][k], self.lattice[(i+1)%self.length][j][k], self.copy_ancilla)
                        self.apply(square(self.bits), self.copy_ancilla, self.temp_ancilla)
                        self.reset(self.copy_ancilla)
                    self.apply(AbsDiff(self.bits), self.lattice[i][j][k], self.lattice[i][(j+1)%self.length][k], self.copy_ancilla)
                    self.apply(square(self.bits), self.copy_ancilla, self.temp_ancilla)
                    self.reset(self.copy_ancilla)
                    self.apply(AbsDiff(self.bits), self.lattice[i][j][k], self.lattice[i][j][(k+1)%self.length], self.copy_ancilla)
                    self.apply(square(self.bits), self.copy_ancilla, self.temp_ancilla)
                    self.apply(IntPartMult(self.bits,prefactor3), self.temp_ancilla, self.time_ancilla)
                    self.reset(self.copy_ancilla)
                    self.reset(self.temp_ancilla)
# Pi-sector phase function
    def compute_ConjugatePhase(self, time):
        prefactor = self.spac**self.dim*self.cutoff**2*time/(2**self.bits*mp.pi)
        for i in range(len(self.lattice)):
            for j in range(self.length):
                for k in range(self.length):
                    self.copy_ancilla.set_value(2**(self.bits-1))
                    self.apply(AbsDiff(self.bits), self.lattice[i][j][k], self.copy_ancilla, self.temp_ancilla)
                    self.reset(self.copy_ancilla)
                    self.copy_ancilla += self.temp_ancilla
                    self.reset(self.temp_ancilla)
                    self.apply(square(self.bits), self.copy_ancilla, self.temp_ancilla)
                    self.apply(IntPartMult(self.bits,prefactor), self.temp_ancilla, self.time_ancilla)
                    self.reset(self.copy_ancilla)
                    self.reset(self.temp_ancilla)
# Simulate time evolution
    def TimeEvolve(self, time, steps):
        self.time_ancilla.set_value(1)
        self.apply(qftarith.QFT(self.bits).dag(), self.time_ancilla)
        for h in range(steps):
            self.compute_FieldPhase(time/(h+1))
            for i in range(len(self.lattice)):
                for j in range(self.length):
                    for k in range(self.length):
                        self.apply(qftarith.QFT(self.bits), self.lattice[i][j][k])
            self.compute_ConjugatePhase(time/(h+1))
            for i in range(len(self.lattice)):
                for j in range(self.length):
                    for k in range(self.length):
                        self.apply(qftarith.QFT(self.bits).dag(), self.lattice[i][j][k])
# Phi-sector measurement function
    def MomentumFieldPhase(self, i, j, k, h, time):
        prefactor = self.spac**self.dim*self.energy(i,j,k)*self.cutoff**2*time/(self.length**self.dim*4**(self.bits-1))
        for x1 in range(len(self.lattice)):
            for x2 in range(self.length):
                for x3 in range(self.length):
                    for y1 in range(len(self.lattice)):
                        for y2 in range(self.length):
                            for y3 in range(self.length):
                                cosine = mp.cos(2*mp.pi*(i*(x1-y1) + j*(x2-y2) + k*(x3-y3))/self.length)
                                self.copy_ancilla += self.lattice[y1][y2][y3]
                                self.temp_ancilla += self.lattice[x1][x2][x3]*self.copy_ancilla
                                self.apply(IntPartMult(self.bits,prefactor*cosine).ctrl(), self.estm_ancilla[h], self.temp_ancilla, self.time_ancilla)
                                self.reset(self.copy_ancilla)
                                self.reset(self.temp_ancilla)
# Pi-sector measurement function
    def MomentumConjugatePhase(self, i, j, k, h, time):
        prefactor = self.spac**self.dim*self.cutoff**2*time/(self.length**self.dim*self.energy(i,j,k)*4**(self.bits-1))
        for x1 in range(len(self.lattice)):
            for x2 in range(self.length):
                for x3 in range(self.length):
                    for y1 in range(len(self.lattice)):
                        for y2 in range(self.length):
                            for y3 in range(self.length):
                                cosine = mp.cos(2*mp.pi*(i*(x1-y1) + j*(x2-y2) + k*(x3-y3))/self.length)
                                self.copy_ancilla += self.lattice[y1][y2][y3]
                                self.temp_ancilla += self.lattice[x1][x2][x3]*self.copy_ancilla
                                self.apply(IntPartMult(self.bits,prefactor*cosine).ctrl(), self.estm_ancilla[h], self.temp_ancilla, self.time_ancilla)
                                self.reset(self.copy_ancilla)
                                self.reset(self.temp_ancilla)
# Quantum phase estimation
    def measure_momentum(self, p1, p2, p3, steps):
        for h in range(self.bits):
            self.apply(H, self.estm_ancilla[h])
            for s in range(steps):
                self.MomentumFieldPhase(p1,p2,p3,h,2**h/(s+1))
                for i in range(len(self.lattice)):
                    for j in range(self.length):
                        for k in range(self.length):
                            self.apply(qftarith.QFT(self.bits).ctrl(), self.estm_ancilla[h], self.lattice[i][j][k])
                self.MomentumConjugatePhase(p1,p2,p3,h,2**h/(s+1))
                for i in range(len(self.lattice)):
                    for j in range(self.length):
                        for k in range(self.length):
                            self.apply(qftarith.QFT(self.bits).dag().ctrl(), self.estm_ancilla[h], self.lattice[i][j][k])
        self.apply(qftarith.QFT(self.bits).dag(), self.estm_ancilla)


simu = JLPsim(2, 2, 2, 1, 1, 0.1, [[[1,1],[1,1]]], 2)
simu.build_VacuumState()
simu.computeLDL()
print(mp.matrix(simu.VacuumMatrix))
print(mp.matrix(simu.CholeskyDiagonal))
print(simu.InverseTriangular)
# simu.build_VacuumState()
# simu.MakeWavepacket(1)
# simu.TimeEvolve(100, 1)
# simu.measure_momentum(0,1,0, 1)

qpu = get_default_qpu()
circ = simu.to_circ(link=[qftarith])
job = circ.to_job()
# job = circ.to_job(qubits=simu.estm_ancilla)
result = qpu.submit(job)
for sample in result:
    print("State {} Amplitude: {}".format(sample.state,sample.amplitude))
    # print("State %s: probability %s +/- %s" % (sample.state, sample.probability, sample.err))