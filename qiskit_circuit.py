import matplotlib.pyplot as plt
from qiskit import *
import numpy as np 
import os
import sys
from Compiler_1qubit.Compiler_mcas_qiskit.compiler_qiskit_mcas import transpile_ciruit_for_diamond, construct_full_mcas_file


########## Circuit part ############
qc = QuantumCircuit(3, 3)
qc.x(0)
qc.y(0)

########## Compiler part ############
directory = os.path.dirname(os.path.abspath(sys.argv[0]))
result = construct_full_mcas_file(qc, 'tom', directory)

########## show circuit ############
qc.draw(output='mpl')
plt.show()
plt.savefig('circuit.png')
tranpiled_circuit = transpile_ciruit_for_diamond(qc)
tranpiled_circuit.draw(output='mpl')
plt.show()
