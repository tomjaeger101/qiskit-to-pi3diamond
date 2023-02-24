import matlab.engine
import qutip as qt
import numpy as np
import OptimalControl as OC; reload(OC)
from qutip_enhanced import *


with open(__file__, 'r') as f:
    script_code = f.read()


############################
#NV_Hamiltonian (not needed for this pulse ... I just played with it)
############################
dims = [2, 3, 2, 2, 2] #NV_el, N14, C13414,13C90,13C13, 13C6

spin_num2name = {1: '14n', 2: '13c414', 3: '13c90', 4: '13c13', 5: '13c6', 6: '13c-5', 7: '13c-6'}
spin_name2num = {v: k for k, v in spin_num2name.items()}
hf_para_n = \
    {'13c-5': -0.0059054,
     '13c-6': -0.0065554,
     '13c13': 0.0123119136673,
     '13c414': 0.412969173247,
     '13c6': 0.00588444108978,
     '13c90': 0.0887466495476,
     '14n': -2.16418839453}
Azz = hf_para_n
gamma = {'13c': 10.70620445, '14n': 3.07522462025, 'e': -28031.67956263535}
B = 0.66556876294903522
qp = {'14n': -4.94572259602}
nvham = nv_hamilton.NVHam(
    magnet_field={'z': B},
    n_type='14n',
    nitrogen_levels=[0, 1, 2],
    electron_levels=[0, 1, 2], 
    gamma=gamma,
    qp=qp,
    hf_para_n=hf_para_n)



for c13 in ['13c414', '13c90', '13c13', '13c6', '13c-6', '13Cm3'][:len(dims) - 2]:
    nvham.add_spin(np.diag([0, 0, Azz[c13]]), nvham.h_13c(), [0, 1])

nv_center_rot_frame_hamiltonian = nvham.h_nv

erotframe_nodet_after_carbon_spins = nvham.h_e
nrotframe14n_nodet = nvham.h_n
nrotframe_13cnodet = nvham.h_13c()




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
T = 200. # Overall pulse time ns 
dims = [3] 
nf = 1 # Number of frequencies

##############################
# Slices and fields
##############################

maxrabis = np.array([1.0/np.sqrt(2), 1.0 / np.sqrt(2)]) # 1.0 / np.sqrt(2), 1.0 / np.sqrt(2)]) #, 1.0 / np.sqrt(2)]) # Sqrt(2) comes from the normalisation!
# take a look for this into Sebastians thesis

# field_list = [0.01*np.random.rand(1, n_bins+add_slices)[0], 0.01*np.random.rand(1, n_bins+add_slices)[0]]
# List for the initial amplitude lists!
field_list = [0.0001 * np.ones(n_bins), 0.0001 * np.ones(n_bins)] #, 0.0001 * np.ones(n_bins)] #0.0001 * np.ones(n_bins), 0.0001 * np.ones(n_bins)] # 0.0001 * np.ones(n_bins)]  # [1x50, 1x50, 1x50, 1x50]
fields = np.dstack(field_list)[0] # changes the dimension to 50 x 4
D = np.prod(dims) # Product of the elements over a given axis
control_type = nf * 'mm'
control_par = [matlab.double([-i, 2 * i]) for i in maxrabis] # denotes the range of the Amplitudes (how much they can vary ): In this cas [-1/np.sqrt(2), 2/np.sqrt(2)]

##############################
# Controls. Here we put together our control Hamiltonians. Each sepearate Entry has his own control Hamiltonian! 
##############################

L_Bc_xy = [
    rx(np.pi, '+1'),
    ry(np.pi, '+1'),
    #rx(np.pi, '-1'),
    #ry(np.pi, '-1')
]

L_Bc = L_Bc_xy

##############################
# Setting up the ensemble and weighting
##############################

p_scales_base = [np.array(len(L_Bc_xy) * [i]) for i in [0.999, 1., 1.001]]  # [0.98, 0.99, 1., 1.01, 1.02]] # Account for the frequency detuning  
spectral_width = 0.001
freqspan = np.arange(-0.3, 0.3, 0.03) 
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
c_labels = (['x+', 'y+'])# + ['x-', 'y- ']) if nf == 2 else (['x', 'y'])#(['x']) #, 'y-'])#(['x+', 'y+ '] + ['x-', 'y- ']) if nf == 2 else (['x', 'y'])

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

dp.subset_names = ['+1']
subset = [slice(0,2,1)]

dp.save(script_path=__file__, script_code=script_code, substr='electron_pi_pulse', subsets=subset)