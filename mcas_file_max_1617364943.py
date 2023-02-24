# coding=utf-8 
from pi3diamond import pi3d 
import numpy as np 
import os 
import UserScripts.helpers.sequence_creation_helpers as sch; reload(sch) 
import multi_channel_awg_seq as MCAS; reload(MCAS) 
import UserScripts.helpers.operations as op; reload(op) 
from collections import OrderedDict 
 
seq_name = os.path.basename(__file__).split('.')[0] 
nuclear = sch.create_nuclear(__file__) 
with open(os.path.abspath(__file__).split('.')[0] + '.py', 'r') as f: 
	meas_code = f.read() 
 
def ret_ret_mcas(pdc):
	def ret_mcas(current_iterator_df):
		mcas = MCAS.MultiChSeq(seq_name=seq_name, ch_dict={'2g': [1, 2], '128m': [1, 2]})
		for idx, _I_ in current_iterator_df.iterrows():

			op.full_initialisation(mcas, '+++')
			op.readout_nuclear_spin_state(mcas, '+++')
			
			pi3d.gated_counter.set_n_values(mcas)
		return mcas
	return ret_mcas
 
def settings(pdc={}):
	ana_seq=[
	['init', '<', 0, 0, 10, 2],
	['init', '<', 0, 0, 10, 1],
	['result', '<', 0, 0, 10, 2],
	]

	sch.settings(
		nuclear=nuclear,
		ret_mcas=ret_ret_mcas(pdc),
		analyze_sequence=ana_seq,
		pdc=pdc,
		meas_code=meas_code)

	nuclear.x_axis_title = 'tau_half [mus]'
	nuclear.analyze_type = 'standard'

	pi3d.gated_counter.trace.analyze_type = 'standard'
	pi3d.gated_counter.trace.consecutive_valid_result_numbers = [0]
	pi3d.gated_counter.trace.average_results = True

	nuclear.number_of_simultaneous_measurements = 1

	nuclear.parameters = OrderedDict(
		(
			('sweeps', range(30)),
		)
	)
def run_fun(abort, **kwargs):
	pi3d.readout_duration = 150e6
	nuclear.debug_mode = False
	settings()
	nuclear.run(abort)