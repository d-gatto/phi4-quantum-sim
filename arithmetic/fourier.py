from qat.lang.AQASM import QRoutine, H, PH, SWAP
from mpmath import pi

# Quantum Fourier transform subroutine
def QFT(qbits):
    rout = QRoutine()
    wires = rout.new_wires(qbits)  
    for j in range(qbits):
        rout.apply(H, wires[j])
        for k in range(1, qbits-j):
            rout.apply(PH(-float(pi/2**k)).ctrl(), wires[k+j], wires[j])
    if(qbits>1):
        for h in range(qbits//2):
            rout.apply(SWAP, wires[h], wires[qbits-h-1])
    return rout

def genCZ(qbits):
    rout = QRoutine()
    wires = rout.new_wires(2*qbits)  
    for j in range(qbits):
        for k in range(j, qbits):
            rout.apply(PH(float(pi/2**k)).ctrl(), wires[k], wires[j+qbits])
    return rout

# # Inverse transform subroutine
# def qFtInv(qbits):
#     rout = QRoutine()
#     wires = rout.new_wires(qbits)  
#     for j in range(qbits):
#         rout.apply(H, wires[j])
#         for k in range(1, qbits-j):
#             rout.apply(PH(float(pi/2**k)).ctrl(), wires[k+j], wires[j])
#     if(qbits>1):
#         for h in range(qbits//2):
#             rout.apply(SWAP, wires[h], wires[qbits-h-1])
#     return rout