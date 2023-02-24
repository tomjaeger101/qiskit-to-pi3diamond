# Compiler-mcas-qiskit
The repository for the qiskit QuantumCircuit object/circuit to mcas compiler. <br />
Code by Max Keller <br />

The file builder has the basic stettings and structures, that are needed to run pi3diamond code. Inital state, readout, thresholds, etc. can be changed here <br />
The compiler_qiskit_mcas.py executes the compalation process and appends the neacessary operations depending on the qiskit circuit. <br />
The mcas_file_max_1617364943.py is an exmple of the output. <br />
The operations.py file asigns the qiskit building blocks pulse sequences. <br />
qiskit_circuit.py shows how the compiler is used.  <br />
<br />
NV_backend and First_optimal_control_pulse are not finished
