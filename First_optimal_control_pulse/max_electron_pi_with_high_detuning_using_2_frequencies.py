import matlab.engine
import qutip as qt
import numpy as np
import OptimalControl as OC; reload(OC)
from qutip_enhanced import *


with open(__file__, 'r') as f:
    script_code = f.read()

##############################
# Extra rotation operators for qubit rotations in a qutrit
##############################
def ry(phi, target):
    if target == '+1':
        qobj = np.pi * qt.Qobj([[np.cos(phi / 2.), -1.j * np.sin(phi / 2.), 0.], # Why do I multiply these Gates with PI ?
                                [1.j * np.sin(phi / 2.), np.cos(phi / 2.), 0.],
                                [0., 0., 1.]])
    if target == '-1':
        qobj = np.pi * qt.Qobj([[0., 0., 0.],
                                [0., np.cos(phi / 2.), -1.j * np.sin(phi / 2.)],
                                [0., 1.j * np.sin(phi / 2.), np.cos(phi / 2.)]])
    if target == '+-':
        qobj = np.pi * qt.Qobj([[np.cos(phi / 2.), 0., -1.j * np.sin(phi / 2.)],
                                [0., 0., 0.],
                                [1.j * np.sin(phi / 2.), 0., np.cos(phi / 2.)]]
                               )

    return qobj.data.toarray()


def sigma_x(target):
    if target == '+1':
        qobj = np.pi * qt.Qobj([[0., 1., 0.],
                                [1., 0., 0.],
                                [0., 0., 1.]])
    if target == '-1':
        qobj = np.pi * qt.Qobj([[1., 0., 0.],
                                [0., 0., 1.],
                                [0., 1., 0.]])
    return qobj.data.toarray()


def rx(phi, target):
    if target == '+1':
        qobj = np.pi * qt.Qobj([[np.cos(phi / 2.), np.sin(phi / 2.), 0.],
                                [np.sin(phi / 2.), np.cos(phi / 2.), 0.],
                                [0., 0., 1.]])
    if target == '-1':
        qobj = np.pi * qt.Qobj([[0., 0., 0.],
                                [0., np.cos(phi / 2.), np.sin(phi / 2.)],
                                [0., np.sin(phi / 2.), np.cos(phi / 2.)]])
    if target == '+-':
        qobj = np.pi * qt.Qobj([[np.cos(phi / 2.), 0., np.sin(phi / 2.)],
                                [0., 0., 0.],
                                [np.sin(phi / 2.), 0., np.cos(phi / 2.)]]
                               )

    return qobj.data.toarray()



##############################
# Time and bin settings
##############################

n_bins = 50 # Number of time slices
T = 200. # Overall pulse time 
dims = [3] 
nf = 2 # Number of frequencies

##############################
# Slices and fields
##############################

maxrabis = np.array([1.0 / np.sqrt(2), 1.0 / np.sqrt(2), 1.0 / np.sqrt(2), 1.0 / np.sqrt(2)]) #, 1.0 / np.sqrt(2)]) # Sqrt(2) comes from the normalisation!
# take a look for this into Sebastians thesis

# field_list = [0.01*np.random.rand(1, n_bins+add_slices)[0], 0.01*np.random.rand(1, n_bins+add_slices)[0]]
# List for the initial amplitude lists!
field_list = [0.0001 * np.ones(n_bins), 0.0001 * np.ones(n_bins), 0.0001 * np.ones(n_bins), 0.0001 * np.ones(n_bins)] # 0.0001 * np.ones(n_bins)]  # [1x50, 1x50, 1x50, 1x50]
fields = np.dstack(field_list)[0] # changes the dimension to 50 x 4
D = np.prod(dims) # Product of the elements over a given axis
control_type = nf * 'mm' # the number of m's gives the number of rotation axie: for example x,y & would mean 'mm' 
control_par = [matlab.double([-i, 2 * i]) for i in maxrabis] # denotes the range of the Amplitudes (how much they can vary ): In this cas [-1/np.sqrt(2), 2/np.sqrt(2)]

##############################
# Controls. Here we put together our control Hamiltonians. Each sepearate Entry has his own control Hamiltonian! 
##############################

L_Bc_xy = [
    rx(np.pi, '+1'),
    ry(np.pi, '+1'),
    rx(np.pi, '-1'),
    ry(np.pi, '-1')
]

L_Bc = L_Bc_xy

##############################
# Setting up the ensemble and weighting
##############################

p_scales_base = [np.array(len(L_Bc_xy) * [i]) for i in [0.999, 1., 1.001]]  # [0.98, 0.99, 1., 1.01, 1.02]] # Account for the frequency detuning  
spectral_width = 0.001
freqspan = np.arange(-1.31, 1.31, 0.1) 
detunings_base = freqspan
detunings = np.repeat(detunings_base, len(p_scales_base))
p_scales = np.repeat(p_scales_base, len(detunings_base), axis=0)

weight = np.ones([len(detunings)])
weight = weight / np.sum(weight)
weight = matlab.double(weight.tolist())


##############################
# Hamiltonian
##############################

# hfl = OC.misc.mfl({'14N': [1, 0, -1]})
def H_drift(det):
    out = 2*np.pi * np.diag([-det, -det, -det])
    return matlab.double(out.tolist(), is_complex=True)


H_drift_list = [H_drift(det) for det in detunings]
H_ctrl_list = [[matlab.double((p[i] * Bc).tolist(), is_complex=True) for i, Bc in enumerate(L_Bc)] for p in p_scales]


#pi_operator = get_rot_operator_all_spins(dims=dims, angle=np.pi, rotated_spin=0, rotation_axis={'y': 0})


##############################
# Setting up the gates
##############################

initial_gate = qt.Qobj(np.eye(D))

# gate = qt.Qobj(np.array([[0,1,0],[1,0,0],[0,0,1]]))
gate = qt.Qobj([[0, 1, 0], [1, 0, 0], [0, 0, 1]])

final_gate = gate

##############################
# Dynamo and optimization parameters
##############################

task = 'closed gate'
desc = 'double quantum'
c_labels = (['x+', 'y+ '] + ['x-', 'y- ']) if nf == 2 else (['x', 'y'])#(['x']) #, 'y-'])#(['x+', 'y+ '] + ['x-', 'y- ']) if nf == 2 else (['x', 'y'])

TL = [[T / n_bins, 0]] * n_bins if n_bins > 0 else []
TL = matlab.double(TL)

##############################
# Fields
##############################

# fields[:n_bins, 2:] = 0
# fields[n_bins:, :2] = 0

##############################
# Set up Dynamo
##############################

dyn_path = '/Users/maxkeller/Documents/Uni/Hiwi/Dynamo_smite/'

dp = OC.dynamo_helpers.DynPython(dyn_path, initial=initial_gate, final=final_gate, dims=dims, add_slices=0)

dp.set_dyn(task, H_drift_list, H_ctrl_list, weight)
dp.set_labels(desc, c_labels)
dp.seq_init(n_bins=int(n_bins), add_slices=0, TL=TL, control_type=control_type, control_par=control_par)
dp.set_controls(matlab.double(fields.tolist()))

##############################
# Fields
##############################

mask = np.array(dp.get_mask(optimize_times=False))

##############################
# Start optimaization
##############################

dp.open_ui()
dp.search(mask, options=dict(max_walltime=100000., plot_interval=5, StepTolerance=1e-13)) # the search function searches then after the optimal control pulse!

##############################
# Save
##############################

dp.subset_names = ['+1', '-1']
subset = [slice(1,3,1), slice(3,5,1)]

dp.save(script_path=__file__, script_code=script_code, substr='chrestenson', subsets=subset)