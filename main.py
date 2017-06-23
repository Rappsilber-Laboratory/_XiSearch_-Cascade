"""
Execution of XiSearch on several DBs of increasing size
After each DB, Spactra satisfying a certain FDR have to be removed from the peak files for the following XiSearch
"""

import logging
from lib import pipeline
import os
import sys
import re
import pandas as pd
import numpy as np
import subprocess
from lib.XiWrapper import XiSearchOutOfMemoryException
from lib.iBAQ_FASTA_handler import IbaqExtraction, FastaHandler

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


# OBJECT: Experiment: holds all the necessary settings for each Fraction/replicate
class Experiment:
    """
    one object per MS analysis

    unique for each object:
    peak files
    shared between objects:
    fasta ressource file

    """

    def __init__(self, list_of_peak_files, maxquant_result_file, exp_name):
        self.ibaq_object = IbaqExtraction(maxquant_result_file)
        self.peak_files = list_of_peak_files
        self.name = exp_name


def xifdr_result_to_spectra_list(xifdr_results, out_file):
    """
    takes a list of xifdr result files
    picks out the "_false_PSM_" file
    removes all entries containing decoys
    writes column 2 and 3 to new csv-file
    RETURN:
        csv-file path
    """
    input_files = [s for s in xifdr_results if "_false_PSM_" in s]
    if len(input_files) > 1:
        raise AttributeError("More than one candidate file in 'xifdr results' matching '*_false_PSM_*'")
    df = pd.read_csv(input_files[0])
    df = df[df["isDecoy"] == False]
    df = df.loc[:, ("run", "scan")]
    df.to_csv(out_file, index=False)
    return out_file

# # # Test cases
# result_files = ['/home/henning/mnt/xitu/scripts/170323_Ribsome_IBAQ/results/with_MAXCANDIDATES/Fr14/0.9/xifdr_output/FDR_1.000_0.050_1.000_1.000_1.000_10000.000_false_summary_xiFDR1.0.14.csv',
#  '/home/henning/mnt/xitu/scripts/170323_Ribsome_IBAQ/results/with_MAXCANDIDATES/Fr14/0.9/xifdr_output/FDR_1.000_0.050_1.000_1.000_1.000_10000.000_false_Linear_Peptides1.0.14_xiFDR1.0.14.csv',
#  '/home/henning/mnt/xitu/scripts/170323_Ribsome_IBAQ/results/with_MAXCANDIDATES/Fr14/0.9/xifdr_output/FDR_1.000_0.050_1.000_1.000_1.000_10000.000_false_Linear_PSM_xiFDR1.0.14.csv',
#  '/home/henning/mnt/xitu/scripts/170323_Ribsome_IBAQ/results/with_MAXCANDIDATES/Fr14/0.9/xifdr_output/FDR_1.000_0.050_1.000_1.000_1.000_10000.000_false_ppi_xiFDR1.0.14.csv',
#  '/home/henning/mnt/xitu/scripts/170323_Ribsome_IBAQ/results/with_MAXCANDIDATES/Fr14/0.9/xifdr_output/FDR_1.000_0.050_1.000_1.000_1.000_10000.000_false_PSM_xiFDR1.0.14.csv',
#  '/home/henning/mnt/xitu/scripts/170323_Ribsome_IBAQ/results/with_MAXCANDIDATES/Fr14/0.9/xifdr_output/FDR_1.000_0.050_1.000_1.000_1.000_10000.000_false_PeptidePairs_xiFDR1.0.14.csv',
#  '/home/henning/mnt/xitu/scripts/170323_Ribsome_IBAQ/results/with_MAXCANDIDATES/Fr14/0.9/xifdr_output/FDR_1.000_0.050_1.000_1.000_1.000_10000.000_false_Links_xiFDR1.0.14.csv',
#  '/home/henning/mnt/xitu/scripts/170323_Ribsome_IBAQ/results/with_MAXCANDIDATES/Fr14/0.9/xifdr_output/FDR_1.000_0.050_1.000_1.000_1.000_10000.000_false_proteingroups_xiFDR1.0.14.csv']
#
# xifdr_result_to_spectra_list(
#     xifdr_results=result_files,
#     out_file=r'testing/matched_spectra.csv'
# )


def spectra_filter(csv_of_scans, peak_files, out_file, xisearch_path="XiSearch.jar"):
    assert type(peak_files) == list, "peak_files needs to be of type list but is: {}".format(type(peak_files))

    # does out dir exist, create if not
    if not os.path.exists(os.path.split(out_file)[0]):
        os.makedirs(os.path.split(out_file)[0])

    # create *.msmlist text file as input file for filter
    msmfile = os.path.join(os.path.split(out_file)[0], "peaklist.msmlist")
    with open(msmfile, 'w+') as f:
        for peak_f in peak_files:
            f.write(os.path.abspath(peak_f) + '\n')

    # execute XiSearch spectra filter
    # "java -cp XiSearch.jar rappsilber.ms.dataAccess.filter.spectrafilter.ScanFilteredSpectrumAccess [list_of_scans.csv] [peaklist] true"
    xi_arguments = ["java", "-cp",
                    xisearch_path,
                    "rappsilber.ms.dataAccess.filter.spectrafilter.ScanFilteredSpectrumAccess",
                    csv_of_scans]
    # peaks files
    # for peak_file in peak_files:
    #     xi_arguments.append("--peaks=" + os.path.abspath(peak_file))
    xi_arguments += [msmfile]
    # exclude specified scans
    xi_arguments += ["true"]
    # output settings
    # xi_arguments += [">", out_file]
    # runnning of XiSearch
    print("xisearch arguments: {}".format(" ".join(map(str, xi_arguments))))
    with open(out_file, 'w+') as mgf_output:
        process = subprocess.Popen(xi_arguments, stdout=mgf_output, stderr=subprocess.PIPE, universal_newlines=True)
        while True:
            output = process.stderr.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                logger.debug("XiSearch: " + output.strip())
        # ret_code = process.wait()
        mgf_output.flush()
    return out_file

# # # Test cases
# spectra_filter(csv_of_scans=r"testing/matched_spectra.csv",
#                peak_files=[r'../../Data/Input/170323_iBAQ_based_opt/Fr14/peak_files/B160805_05_Lumos_ML_IN_190_FOMix_Tryp_SEC_Fr14.HCD.FTMS.peak.apl',
#                            r'../../Data/Input/170323_iBAQ_based_opt/Fr14/peak_files/B160805_05_Lumos_ML_IN_190_FOMix_Tryp_SEC_Fr14.HCD.FTMS.sil0.apl'],
#                out_file="testing/filtered_peaks.mgf")


def create_db(method, fasta_file, additional_params):
    if method == 'relIbaq':


# METHOD: pipeline(list_of_experiments, xi_xifdr_settings_dict, fasta_file, output_basedir):
def pipeline_execution(list_of_experiments, xi_xifdr_settings_dict, fasta_file, output_basedir):
    for experiment in list_of_experiments:
        peak_files = experiment.peak_files
        xifdr_results = []
        for i, DB in enumerate(list_of_DBs):    # TODO is the list stored in each experiment or can we generalize?
            db_res_dir = os.path.join(output_basedir, experiment.name, "DB_"+str(i))
            if not os.path.exists(db_res_dir):
                os.makedirs(db_res_dir)

            db = create_db(method, fasta_file, additional_params)

            # filter out matched spectra from peak files
            if i != 0:  # do not filter spectra for first iteration
                # assert that xifdr result list is not empty
                assert xifdr_results, "This if condition should not be reached with empty 'xifdr_results'"
                matched_spectra_list = xifdr_result_to_spectra_list(
                    xifdr_results=xifdr_results,
                    out_file=os.path.join(db_res_dir, "excluded_spectra.csv"))
                peak_files = spectra_filter(
                    csv_of_scans=matched_spectra_list,
                    peak_files=peak_files,
                    out_file=os.path.join(db_res_dir, 'peak_file_of_remaining_spectra.mgf'))

            # Carry out XiSearch and FDR Analysis
            # make sure that a failed pipeline run leads to an empty list
            xifdr_results = []
            try:
                xi_result, xifdr_results = pipeline.execute_pipeline(
                    # optional general settings
                    output_basedir=db_res_dir,
                    # xi settings
                    list_of_fasta_dbs=[db],
                    xi_config=xi_xifdr_settings_dict['xi_config'],
                    peak_files=peak_files,
                    # optional xi settings
                    xi_memory=xi_xifdr_settings_dict['xi_memory'],
                    additional_xi_parameters=xi_xifdr_settings_dict['additional_xi_parameters'],
                    # optional xifdr settings
                    pepfdr=xi_xifdr_settings_dict['xifdr_settings']['pepfdr'],
                    reportfactor=xi_xifdr_settings_dict['xifdr_settings']['reportfactor'],
                    additional_xifdr_arguments=xi_xifdr_settings_dict['xifdr_settings']['additional_xifdr_arguments'])
            except XiSearchOutOfMemoryException as e:
                msg1 = "XiSearch produced a Memory Exception for cmd '{}'".format(e.cmd)
                msg2 = "removing output stub '{}'".format(e.out_file)
                logging.error(msg1)
                print msg1
                logging.error(msg2)
                print msg2
                os.remove(e.out_file)
                # TODO flow after exception: search for next DB? search for next experiment? Decision needs to go in log
                raise NotImplementedError("Flow after catching of exception not implemented")

#             # if number_of_satisfying_matches < k:
#                 # break         # this ensures a well controlled FDR. Too few matches lead to bad estimation of FDR.
#                 # TODO is this feasible with the small amount of links identified by crosslink search?
            raise NotImplementedError("Abortion of search for small number of Crosslinks is not decided yet")


# METHOD: main
if __name__ == "__main__":
    # print help message if script is called without argument
    if len(sys.argv) != 1:
        print """Script has to be called with output dir as argument.
            Output dir has to contain config file "config.py"."""
        sys.exit(1)

    # set output dir
    output_basedir = sys.argv[1]

    log_file = os.path.join(output_basedir, __name__ + '.log')

    logging.basicConfig(filename=log_file, level=logging.DEBUG,
                        format='%(asctime)s - %(levelname)s - %(name)s - %(message)s')

    # import of config.py
    sys.path.append(output_basedir)
    import config

    list_of_experiments = config.list_of_experiments
    fasta_base_file = config.fasta_base_file
    xi_xifdr_settings_dict = config.xi_xifdr_settings_dict


    # list_of_experiments = []
    #
    # Fr14 = Experiment(
    #     databases=[],
    #     list_of_peak_files=[r'Data/Fr14/peak_files/B160805_05_Lumos_ML_IN_190_FOMix_Tryp_SEC_Fr14.HCD.FTMS.peak.apl',
    #                         r'Data/Fr14/peak_files/B160805_05_Lumos_ML_IN_190_FOMix_Tryp_SEC_Fr14.HCD.FTMS.sil0.apl'],
    #     maxquant_result_file=r'Data/Fr14/Fr14_proteinGroups.txt',
    #     exp_name="Fr14")
    # list_of_experiments += [Fr14]
    #
    # Fr15 = Experiment(
    #     databases=[],
    #     list_of_peak_files=[r'Data/Fr15/peak_files/B160805_06_Lumos_ML_IN_190_FOMix_Tryp_SEC_Fr15.HCD.FTMS.peak.apl',
    #                         r'Data/Fr15/peak_files/B160805_06_Lumos_ML_IN_190_FOMix_Tryp_SEC_Fr15.HCD.FTMS.sil0.apl'],
    #     maxquant_result_file=r"Data/Fr15/proteinGroups.txt",
    #     exp_name="Fr15")
    # list_of_experiments += [Fr15]
    #
    # Fr16 = Experiment(
    #     databases=[],
    #     list_of_peak_files=[r'Data/Fr16/peak_files/B160805_07_Lumos_ML_IN_190_FOMix_Tryp_SEC_Fr16.HCD.FTMS.peak.apl',
    #                         r'Data/Fr16/peak_files/B160805_07_Lumos_ML_IN_190_FOMix_Tryp_SEC_Fr16.HCD.FTMS.sil0.apl'],
    #     maxquant_result_file=r"Data/Fr16/proteinGroups.txt",
    #     exp_name="Fr16")
    # list_of_experiments += [Fr16]
    #
    # fasta_base_file = r'Data/fasta_HomoSapiens_UP000005640_170306/uniprot.fasta'
    # output_basedir = '../../Data/Results/170609_cascade_search/'
    # xi_xifdr_settings_dict = {
    #     'xi_config': r'Data/Fr14/xisearch_ribosome_10770.cfg',  # same config file for all searches
    #     'xi_memory': '100G',
    #     'additional_xi_parameters': ["--xiconf=TOPMATCHESONLY:true", "--xiconf=MAXPEAKCANDIDATES:5000"],
    #     'xifdr_settings': {
    #         'pepfdr': "10",
    #         'additional_xifdr_arguments': [r"--lenghtgroups=4"],
    #         'reportfactor': "10000"
    #     }
    # }

    pipeline_execution(
        list_of_experiments=list_of_experiments,
        xi_xifdr_settings_dict=xi_xifdr_settings_dict,
        fasta_file=fasta_base_file,
        output_basedir=output_basedir)
    logging.shutdown()

# TODO
# save the databases in an initial step as attribute of each experiment
    # this way adjustment of what databases to use is easier
# initialize each experiment with the necessary parameters
#   parameters include the databases to search on
    # TODO do the DBs differ for each Fraction or can we generalize by just giving parameters for DB creation?
# TODO adjust paths to new data structure
