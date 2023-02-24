import matplotlib.pyplot as plt
from qiskit import *
import numpy as np 
import qiskit
import pickle
from json import JSONEncoder
import os
import sys
import random

from qiskit.circuit import Parameter
from qiskit.circuit import measure
from Compiler_1qubit.Compiler_mcas_qiskit.compiler_qiskit_mcas import transpile_ciruit_for_diamond, construct_full_mcas_file
qc = QuantumCircuit(1, 1)

########## Circuit part ############


qc.x(0)
qc.y(0)



#directory = os.path.dirname(os.path.abspath(sys.argv[0]))
#result = construct_full_mcas_file(qc, 'tom', directory)
qc.draw(output='mpl')
plt.show()
#plt.savefig('circuit.png')
tranpiled_circuit = transpile_ciruit_for_diamond(qc)
tranpiled_circuit.draw(output='mpl')
plt.show()
