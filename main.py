from qat.qpus import PyLinalg
from qat.lang.AQASM import Program
from qat.lang.AQASM.qint import QInt
from fourier import QFT

def main():
    qbits = 2    
    
    QFT = QFT(1)

    #Testing on vector |0> of the computational basis
    prog = Program()
    x = prog.qalloc(1, QInt)
    y = prog.qalloc(1, QInt)
    prog.apply(qFt_1, reg)
    circ = prog.to_circ()

    qpu = PyLinalg()

    #%qatdisplay circ
    job = circ.to_job()
    result = qpu.submit(job)

    for sample in result:
        print("State {} Amplitude: {}".format(sample.state,sample.amplitude))

if __name__=='__main__':
    main()