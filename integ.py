from qat.lang.AQASM import Program, QRoutine, X, H, PH, RX, RY, RZ, SWAP
from qat.lang.AQASM.qint import QInt
import qat.lang.AQASM.qftarith as qftarith
import mpmath as mp
from qat.qpus import get_default_qpu

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
def subtract(bits):
    rout = QRoutine()
    min = rout.new_wires(bits, QInt, reverse_bit_order=True)
    sub = rout.new_wires(bits, QInt, reverse_bit_order=True)
    diff = rout.new_wires(bits, QInt, reverse_bit_order=True)
    diff += min - sub
    return rout
def multiply(inbits, outbits, multiplier):
    rout = QRoutine()
    qint = rout.new_wires(inbits, QInt, reverse_bit_order=True)
    qout = rout.new_wires(outbits, QInt, reverse_bit_order=True)
    for i in range(multiplier):
        qout += qint
    return rout
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

simu = Program()
x = simu.qalloc(3, QInt, reverse_bit_order=True)
y = simu.qalloc(3, QInt, reverse_bit_order=True)
z = simu.qalloc(1)
x.set_value(2)
y.set_value(6)
simu.apply(AbsDiff(3), x, y, z)

qpu = get_default_qpu()
circ = simu.to_circ(link=[qftarith])
job = circ.to_job()
result = qpu.submit(job)
for sample in result:
    print("State {} Amplitude: {}".format(sample.state,sample.amplitude))