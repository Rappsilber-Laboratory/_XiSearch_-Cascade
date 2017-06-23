# initiate object for each experiment
list_of_experiments = []

Fr14 = {
    'list_of_peak_files': [
        r'../../Data/Input/170323_iBAQ_based_opt/Fr14/peak_files/B160805_05_Lumos_ML_IN_190_FOMix_Tryp_SEC_Fr14.HCD.FTMS.peak.apl',
        r'../../Data/Input/170323_iBAQ_based_opt/Fr14/peak_files/B160805_05_Lumos_ML_IN_190_FOMix_Tryp_SEC_Fr14.HCD.FTMS.sil0.apl'
    ],
    'maxquant_result_file': r'../../Data/Input/170323_iBAQ_based_opt/Fr14/MaxQuant_results_Fr14_proteinGroups.txt',
    'exp_name': "Fr14"}
list_of_experiments += [Fr14]

Fr15 = {
    'list_of_peak_files': [
        r'../../Data/Input/170323_iBAQ_based_opt/Fr15/peak_files/B160805_06_Lumos_ML_IN_190_FOMix_Tryp_SEC_Fr15.HCD.FTMS.peak.apl',
        r'../../Data/Input/170323_iBAQ_based_opt/Fr15/peak_files/B160805_06_Lumos_ML_IN_190_FOMix_Tryp_SEC_Fr15.HCD.FTMS.sil0.apl'
    ],
    'maxquant_result_file': r"../../Data/Input/170323_iBAQ_based_opt/Fr15/MaxQuant_results_Fr15_proteinGroups.txt",
    'exp_name': "Fr15"}
list_of_experiments += [Fr15]

Fr16 = {
    'list_of_peak_files': [
        r'../../Data/Input/170323_iBAQ_based_opt/Fr16/peak_files/B160805_07_Lumos_ML_IN_190_FOMix_Tryp_SEC_Fr16.HCD.FTMS.peak.apl',
        r'../../Data/Input/170323_iBAQ_based_opt/Fr16/peak_files/B160805_07_Lumos_ML_IN_190_FOMix_Tryp_SEC_Fr16.HCD.FTMS.sil0.apl'
    ],
    'maxquant_result_file': r"../../Data/Input/170323_iBAQ_based_opt/Fr16/MaxQuant_results_Fr16_proteinGroups.txt",
    'exp_name': "Fr16"}
list_of_experiments += [Fr16]

fasta_base_file = r'../../Data/Input/170323_iBAQ_based_opt/fasta_HomoSapiens_UP000005640_170306/uniprot.fasta'

xi_xifdr_settings_dict = {
    # same config file for all searches
    'xi_config': r'../../Data/Input/170323_iBAQ_based_opt/Fr14/xisearch_ribosome_10770_altered_BufferOutput-altered_MaxPeptideMass.cfg',
    'xi_memory': '100G',
    'additional_xi_parameters': ["--xiconf=TOPMATCHESONLY:true", "--xiconf=MAXPEAKCANDIDATES:5000"],
    'xifdr_settings': {
        'pepfdr': "10",
        'additional_xifdr_arguments': [],
        'reportfactor': "10000"
    }
}
