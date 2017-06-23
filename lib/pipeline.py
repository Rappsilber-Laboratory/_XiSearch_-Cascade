"""
author: 'Henning Schiebenhoefer'

Input:
Peak Lists
fasta files
Xi cfg
Output:
XiFDR results
"""

import os
import subprocess
import pandas as pd
import sys  # Anzahl der Parameter beim Aufruf
import argparse     # parsen von Befehlen aus der Kommandozeile
import logging
import time
import datetime
import XiWrapper


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


# class XiSearchException(subprocess.CalledProcessError):
#     def __init__(self, returncode, cmd, out_file, output=None):
#         subprocess.CalledProcessError.__init__(returncode, cmd, output)
#         self.out_file = out_file
#         pass
#
#
# class XiSearchOutOfMemoryException(XiSearchException):
#     def __init__(self, returncode, cmd, out_file, output=None):
#         XiSearchException.__init__(returncode, cmd, out_file, output)
#         pass


def setup_xi_logger(logger_name, log_file):
    l = logging.getLogger(logger_name)
    formatter = logging.Formatter('%(message)s')
    fileHandler = logging.FileHandler(log_file)
    fileHandler.setFormatter(formatter)
    streamHandler = logging.StreamHandler()
    streamHandler.setFormatter(formatter)

    l.setLevel(logging.DEBUG)
    l.addHandler(fileHandler)
    l.addHandler(streamHandler)

# setup_xi_logger("XiSearch", "XiSearch.log")
# xi_logger = logging.getLogger("XiSearch")


def calculate_elapsed_time(starttime):
    """ returns elapsed time since starttime in minutes """
    seconds = (time.time() - starttime)
    # m, s = divmod(seconds, 60)
    # h, m = divmod(m, 60)
    # return "%d:%02d:%02d hours" % (h, m, s)
    return str(datetime.timedelta(seconds=seconds))

# # test cases
# starttime = 1491209339.95
# print calculate_elapsed_time(starttime)


# def build_xi_arguments(xi_config, peak_files, fasta_files, memory, output, additional_parameters=()):
#     cmd = []
#     cmd.extend(["java", "-Xmx" + memory,
#                 # deprecated native jave memory tracking
#                 # "-XX:NativeMemoryTracking=detail",
#                 # creates a memory dump on out of memory error
#                 # "-XX:+HeapDumpOnOutOfMemoryError",
#                 # kills process on out of memory error and displays message
#                 # '-XX:OnOutOfMemoryError="echo OutOfMemory occured killing;kill %p"',
#                 "-cp", "XiSearch.jar", "rappsilber.applications.Xi"])
#     # config
#     cmd.append("--config=" + xi_config)
#     # additional parameters
#     for par in additional_parameters:
#         cmd.append(par)
#     # peak_files
#     for peak_file in peak_files:
#         cmd.append("--peaks=" + peak_file)
#     # fasta_files
#     for fasta_file in fasta_files:
#         cmd.append("--fasta=" + fasta_file)
#     # output
#     cmd.append("--output=" + output)
#     return cmd
#
# # print build_xi_arguments("das_ist_die_cfg", ["erste_peaks", "zweite_peaks"], ["fasta1", "fasta2"], "hier_gehts_hin")
# # print build_xi_arguments("das_ist_die_cfg", ["zweite_peaks"], ["fasta1", "fasta2"], "hier_gehts_hin",
# #                          additional_parameters=["--xiconfig=TOPMATCHESONLY:true"])
#
#
# def xi_execution(xi_config, peak_files, fasta_files, memory, output_folder, additional_parameters=list()):
#     """
#     Calls Xi and gives back Xi result
#     Xi gives back a single csv file
#
#     Example cmd line:
#     # java -cp XiSearch.jar rappsilber.applications.Xi --conf=[path to config]
#     # --xiconfig=TOPMATCHESONLY:true --peaks=[path to peaklist1] --peaks=[path to peaklist2] --peaks=[path to peaklist3]
#     # --fasta=[path to fasta file1] --fasta=[path to fasta file2] --fasta=[path to fasta file3]
#     # --output=[path to result file]
#
#     :param xi_config: string
#     :param peak_files: list of strings
#     :param fasta_files: list of strings
#     :param output_folder: string
#     :param additional_parameters: list of strings
#     :return: string: output file name
#     """
#     assert type(peak_files) == type(fasta_files) == type(additional_parameters) == list, \
#         """type of the following files needs to be list. It is actually:
#         peak_files: {}
#         fasta_files: {}
#         additional_parameters: {}""".format(type(peak_files), type(fasta_files), type(additional_parameters))
#     # generate output_file name
#     output_file = os.path.join(output_folder, "xi_results.csv")  # filename example: "xi_results.csv"
#     # generate xi commands
#     xi_cmd = build_xi_arguments(xi_config, peak_files, fasta_files, memory, output_file, additional_parameters)
#     # call xi
#     logger.debug("xisearch arguments: {}".format(" ".join(map(str, xi_cmd))))
#     process = subprocess.Popen(xi_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
#     while True:
#         output = process.stdout.readline()
#         exit_code = process.poll()
#         if output == '' and exit_code is not None:
#             break
#         elif "java.lang.OutOfMemoryError" in output:
#             process.kill()
#             raise XiSearchOutOfMemoryException(returncode=1, cmd=xi_cmd, out_file=output_file, output=output)
#         if output:
#             # print output.strip()
#             logger.debug("XiSearch: " + output.strip())
#     if exit_code != 0:     # if process exit code is non zero
#         raise XiSearchException(exit_code, xi_cmd, output_file, 'XiSearch exited with error message!')
#     # try:
#     #     print subprocess.check_output(xi_cmd)
#     # except subprocess.CalledProcessError:
#     #     raise AttributeError('XiSearch exited with error message!')
#     return output_file

# # Test Cases
# xi_execution("Data/xi_config.cfg", ["Data/B161022_OCCM_BS3_Tryp_SECFr18-15_16_42-31_Oct_2016/B161022_07_Lumos_ML_IN_205_OCCM_BS3_Tryp_SECFr18.HCD.FTMS.peak.apl",
#              "Data/B161022_OCCM_BS3_Tryp_SECFr18-15_16_42-31_Oct_2016/B161022_07_Lumos_ML_IN_205_OCCM_BS3_Tryp_SECFr18.HCD.FTMS.sil0.apl"],
#              ["Data/OCCM_Scerevisiae-11_17_09-09_Feb_2016/OCCM_Scerevisiae.fasta"], "test_Xi_results.csv",
#              additional_parameters=["--xiconf=TOPMATCHESONLY:true"])

# xi_cmd = ["java", "-cp", "XiSearch.jar", "-Xmx1G", "rappsilber.applications.Xi", "--help"]
# print subprocess.check_output(xi_cmd)

class XiFdrWrapper:
    def __init__(self):
        pass

    @staticmethod
    def build_xifdr_arguments(fdr_input_csv, fdr_output_dir, pepfdr, memory="1G", reportfactor="10000",
                              additional_xifdr_arguments=list(), xifdr_filename="xiFDRDB-1.0.14.34-jar-with-dependencies.jar"):
        assert type(fdr_input_csv) == list, """type of fdr_input_csv needs to be list but is: {}""".format(type(fdr_input_csv))
        # Example cmd line:
        # java -Xmx1g -cp xiFDRDB-1.0.13.32-jar-with-dependencies.jar org.rappsilber.fdr.CSVinFDR --psmfdr=X --pepfdr=X
        # --proteinfdr=X --reportfactor=X --linkfdr=X --ppifdr=X --csvOutDir=X --csvBaseName=X csv-file1 csv-file2
        cmd = []
        # memory
        cmd.extend(["java", "-Xmx" + memory, "-cp", xifdr_filename,
                    "org.rappsilber.fdr.CSVinFDR", '--reportfactor='+reportfactor])
        # pepfdr
        cmd.append("--pepfdr=" + pepfdr)
        # additional arguments
        for par in additional_xifdr_arguments:
            cmd.append(par)
        # reportfactor
        # taken into default config as it will probably never be changed
        # cmd.append("--reportfactor=" + reportfactor)
        # csvOutDir
        cmd.append("--csvOutDir=" + fdr_output_dir)
        # input csv
        for i in fdr_input_csv:
            cmd.append(i)
        return cmd

    @staticmethod
    def xifdr_execution(xifdr_input_csv, xifdr_output_dir, pepfdr="5", memory="1G", reportfactor="10000",
                        additional_xifdr_arguments=list()):
        """
        takes XiSearch output and gives back fdr csv
        :param xifdr_input_csv:
        :param xifdr_output_dir:
        :param additional_xifdr_arguments:
        :return:
        """
        list_of_results = []
        starttime = time.time()
        # TODO funzt das auch, wenn man die Klasse mit nem Alias laed? Ich haette lieber etwas mit 'self.'
        # TODO check whether @classmethod can do this
        xifdr_cmd = XiFdrWrapper.build_xifdr_arguments(xifdr_input_csv, xifdr_output_dir, pepfdr, memory, reportfactor,
                                                       additional_xifdr_arguments)
        logger.debug("xifdr arguments: {}".format(" ".join(map(str, xifdr_cmd))))
        # # # # #
        process = subprocess.Popen(xifdr_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        # real time output of Xi messages
        while True:
            output = process.stdout.readline()
            exit_code = process.poll()
            if output == '' and exit_code is not None:
                break
            if output:
                # print output.strip()
                logging.debug("XiFdr: " + output.strip())
        if exit_code != 0:  # if process exit code is non zero
            raise subprocess.CalledProcessError(exit_code, xifdr_cmd)
        logging.debug("Search execution took {} for cmd: {}"
                      .format(calculate_elapsed_time(starttime), xifdr_cmd))
        # # # # #
        for rel_dir, sub_dirs, files in os.walk(xifdr_output_dir):
            list_of_results = [os.path.join(rel_dir, f) for f in files]
        return list_of_results


def fun_makedirs(list_of_dirs):
    """ checks whether directory exists for each directory in list.
    Creates each that does not exist."""
    for directory in list_of_dirs:
        if not os.path.exists(directory):
            os.makedirs(directory)


def execute_pipeline(
        # mandatory xi settings
        list_of_fasta_dbs, xi_config, peak_files,
        # optional general settings
        output_basedir=".",
        # optional xi settings
        xi_memory='1G', additional_xi_parameters=list(), xi_path="XiSearch.jar",
        # optional xifdr settings
        pepfdr="5", xifdr_memory="1G", reportfactor="10000",
        additional_xifdr_arguments=list()):
    """
    This pipeline executes:
    xiSeqrch
    and Xifdr
    RETURNS:
        Xi result file: str
        XiFDR result files: list of str
    """
    assert type(list_of_fasta_dbs) in [list, str], \
        "type of list_of_fasta_dbs must be 'list' or 'str' type but is: {}".format(type(list_of_fasta_dbs))
    if type(list_of_fasta_dbs) == str:
        list_of_fasta_dbs = [list_of_fasta_dbs]
    # base dirs for all the files
    pre_list_of_dirs = ["xi_output", "xifdr_output"]
    # generate necessary dirs for this run
    list_of_dirs = []
    for directory in pre_list_of_dirs:
        list_of_dirs.append(os.path.join(output_basedir, directory))
    fun_makedirs(list_of_dirs)
    # call xisearch
    # starttime = time.time()
    xi_result = XiWrapper.XiWrapper.xi_execution(
        xi_path=xi_path,
        xi_config=xi_config,
        peak_files=peak_files,
        fasta_files=list_of_fasta_dbs,
        memory=xi_memory,
        output_file=os.path.join(list_of_dirs[0], "xi_results.csv"),
        additional_parameters=additional_xi_parameters)
    # logger.info("xi search execution for '{}' took {}"
    #             .format(list_of_dirs[0], calculate_elapsed_time(starttime)))
    # xi_result is a string but xifdr needs a list as input
    xifdr_input = [xi_result]
    # call xifdr
    starttime = time.time()
    xifdr_results = XiFdrWrapper.xifdr_execution(
        xifdr_input_csv=xifdr_input,
        xifdr_output_dir=list_of_dirs[1],
        pepfdr=pepfdr,
        memory=xifdr_memory,
        reportfactor=reportfactor,
        additional_xifdr_arguments=additional_xifdr_arguments
    )
    logger.info("xifdr execution for '{}' took {}"
                .format(list_of_dirs[1], calculate_elapsed_time(starttime)))
    return xi_result, xifdr_results


# # Test cases
if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(levelname)s - %(name)s - %(message)s'
    )
    execute_pipeline(
        xi_path=r"../XiSearch.jar",
        list_of_fasta_dbs=r"../Data/fasta_HomoSapiens_UP000005640_170306/uniprot.fasta",
        xi_config=r"../Data/Fr14/xisearch_ribosome_10770.cfg",
        peak_files=[r"../Data/Fr14/peak_files/B160805_05_Lumos_ML_IN_190_FOMix_Tryp_SEC_Fr14.HCD.FTMS.peak.apl",
                    r"../Data/Fr14/peak_files/B160805_05_Lumos_ML_IN_190_FOMix_Tryp_SEC_Fr14.HCD.FTMS.sil0.apl"],
        xi_memory='100M',  # xi settings
        output_basedir="test",  # optional general settings
        pepfdr="5", xifdr_memory="1G", reportfactor="10000", additional_xi_parameters=list(),  # optional xi settings
        additional_xifdr_arguments=list())

    # # test cases
    # print fun_create_dir_for_experiment("testing", "occm", 5, "Ecoli", 7)


# todo Annahme der Argumente ueber die commandline
