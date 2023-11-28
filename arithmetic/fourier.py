from qat.lang.AQASM import QRoutine, H, PH, SWAP
from mpmath import pi

# Quantum Fourier transform subroutine
def QFT(nqbits):
    rout = QRoutine()
    wires = rout.new_wires(nqbits)  
    for j in range(nqbits):
        rout.apply(H, wires[j])
        for k in range(1, nqbits-j):
            rout.apply(PH(-float(pi/2**k)).ctrl(), wires[k+j], wires[j])
    if nqbits>1:
        for h in range(nqbits//2):
            rout.apply(SWAP, wires[h], wires[nqbits-h-1])
    return rout

def genCZ(nqbits):
    rout = QRoutine()
    wires = rout.new_wires(2*nqbits)  
    for j in range(nqbits):
        for k in range(j, nqbits):
            rout.apply(PH(float(pi/2**k)).ctrl(), wires[k], wires[j+nqbits])
    return rout

# # Inverse transform subroutine
# def qFtInv(nqbits):
#     rout = QRoutine()
#     wires = rout.new_wires(nqbits)  
#     for j in range(nqbits):
#         rout.apply(H, wires[j])
#         for k in range(1, nqbits-j):
#             rout.apply(PH(float(pi/2**k)).ctrl(), wires[k+j], wires[j])
#     if(nqbits>1):
#         for h in range(nqbits//2):
#             rout.apply(SWAP, wires[h], wires[nqbits-h-1])
#     return rout