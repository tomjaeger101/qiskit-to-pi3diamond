import os 
import textwrap


operations_file_name = "operations"
indent_level = "\t\t\t"
alias_for_operations = "op"

class Mcas_file: 

    def __init__(self, sequences, parameters=None):
        self.__content = ""
        self.add_header()
        self.add_sequence_part(sequences)
        self.add_settings()
        self.add_nuclear_params(parameters)
        self.add_run_fun()
    
    def get_content(self): 
        return self.__content

    def add_header(self): 
        self.__content += mcas_header

    def add_nuclear_params(self, parameters): 
        self.__content += ("\tnuclear.parameters = OrderedDict(\n"
        "\t\t(\n")
        if isinstance(parameters, list):
            for param in parameters: 
                self.__content += "\t\t\t{}\n".format(param)
        self.__content += "\t\t)\n\t)\n"
    
    def add_sequence_part(self, sequences): 
        if isinstance(sequences, list):
            for line in sequences:
                line = ".".join((alias_for_operations, line))
                indented_line = textwrap.indent(text=line, prefix=indent_level) 
                self.__content += indented_line + '\n'
        else:
            raise TypeError("The sequences must be in a list!")

    def add_settings(self): 
        self.__content += mcas_settings
    
    def add_run_fun(self): 
        self.__content += run_fun


mcas_header = ("# coding=utf-8 \n" 
"from pi3diamond import pi3d \n"
"import numpy as np \n"
"import os \n" 
"import UserScripts.helpers.sequence_creation_helpers as sch; reload(sch) \n"
"import multi_channel_awg_seq as MCAS; reload(MCAS) \n"
"import UserScripts.helpers."+ operations_file_name + " as op; reload(op) \n"
"from collections import OrderedDict \n \n"
"seq_name = os.path.basename(__file__).split('.')[0] \n"
"nuclear = sch.create_nuclear(__file__) \n"
"with open(os.path.abspath(__file__).split('.')[0] + '.py', 'r') as f: \n"
"\tmeas_code = f.read() \n \n"
"def ret_ret_mcas(pdc):\n"
"\tdef ret_mcas(current_iterator_df):\n"
"\t\tmcas = MCAS.MultiChSeq(seq_name=seq_name, ch_dict={'2g': [1, 2], '128m': [1, 2]})\n"
"\t\tfor idx, _I_ in current_iterator_df.iterrows():\n\n")



mcas_settings = ("\n\t\t\tpi3d.gated_counter.set_n_values(mcas)\n"
"\t\treturn mcas\n"
"\treturn ret_mcas\n \n"
"def settings(pdc={}):\n"
"\tana_seq=[\n"
"\t['init', '<', 0, 0, 10, 2],\n"
"\t['init', '>', 3, 0, 0, 1],\n"
"\t['result', '<', 0, 0, 10, 2],\n"
"\t]\n\n"
"\tsch.settings(\n"
"\t\tnuclear=nuclear,\n"
"\t\tret_mcas=ret_ret_mcas(pdc),\n"
"\t\tanalyze_sequence=ana_seq,\n"
"\t\tpdc=pdc,\n"
"\t\tmeas_code=meas_code)\n\n"
"\tnuclear.x_axis_title = 'tau_half [mus]'\n"
"\tnuclear.analyze_type = 'standard'\n\n"
"\tpi3d.gated_counter.trace.analyze_type = 'standard'\n"
"\t#pi3d.gated_counter.trace.consecutive_valid_result_numbers = [0]\n"
"\tpi3d.gated_counter.trace.average_results = False\n"
"\tnuclear.odmr_interval = 1\n"
"\tnuclear.refocus_interval = 1\n"
"\tnuclear.maximum_odmr_drift = 0.015\n"
"\tnuclear.number_of_simultaneous_measurements = 1\n")





run_fun = ("def run_fun(abort, **kwargs):\n"
"#\tpi3d.readout_duration = 150e6\n"
"\tpi3d.gated_counter.readout_duration = 150e6\n"
"\tnuclear.debug_mode = False\n"
"\tsettings()\n"
"\tnuclear.run(abort)")

