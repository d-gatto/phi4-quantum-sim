from mpmath import re, exp, acos, sqrt, jtheta
from mpmath import j as i


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