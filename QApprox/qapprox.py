#! /usr/bin/python3

import numpy as np
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit import execute
from qiskit import BasicAer
from qiskit import IBMQ

# access to IBM Account
#IBMQ.save_account('MY_TOKEN')
#IBMQ.load_accounts()

backend = BasicAer.get_backend('statevector_simulator')
#backend = BasicAer.get_backend('unitary_simulator')
#backend = BasicAer.get_backend('qasm_simulator')

def qc_approx_sim(x, t1, t2):
    theta1 = x - t1;
    theta2 = x - t2;

    q = QuantumRegister(2, 'q')
    c = ClassicalRegister(2, 'c')
    qc = QuantumCircuit(q, c)

    qc.h( q[0] )
    qc.h( q[1] )

    qc.u3(t1, 0.0, 0.0, q[0]);
    qc.u3(t2, 0.0, 0.0, q[1]);

    qc.barrier( q )
    #qc.measure(q,c)
    qc.measure( q[0], c[0] )
    qc.measure( q[1], c[1] )

    job = execute(qc, backend, shots=1024)

    rslt = job.result()
    #counts = rslt.get_counts(qc)
    #print(counts)

    outputstate = rslt.get_statevector( qc, decimals=13 )
    #print(outputstate)

    qval = outputstate;

    return qval;


E = np.array(  [[1.0, 0.0, 0.0, 0.0],
		[0.0, 2.0, 0.0, 0.0],
		[0.0, 0.0, 3.0, 0.0],
		[0.0, 0.0, 0.0, 4.0]] )

q = qc_approx_sim(0.0, 0.1, 0.3);
q = np.matrix( q ).reshape( (4, 1) )

print( np.transpose(q) * (E * q) );

