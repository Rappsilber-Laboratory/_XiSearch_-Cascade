"""
This script simulates a cascade search on already existing XiSearch results.

This is done by carrying out a xifdr evaluation of the XiResults and removing scans of certain runs that satisfy
FDR-conditions from all successive XiResults PSMs. This process is repeated for every successive XiResult.
"""

import pandas as pd
import timeit
import os
from lib.pipeline import XiFdrWrapper, calculate_elapsed_time
from lib.iBAQ_FASTA_handler import FastaHandler
import logging
import sys
import time
import re


pd.options.display.max_rows = 999


# OBJECT: Experiment: holds all the necessary settings for each Fraction/replicate
class Experiment:
    """
    one object per MS analysis
    """
    def __init__(self, ordered_list_of_xi_results, exp_name):
        self.name = exp_name
        self.ordered_list_of_xi_results = ordered_list_of_xi_results


def guess_dtypes_of_columns(
        csv_to_guess_dtypes=r"/home/henning/ownCloud/masterieren/Data/Results/170323_iBAQ_based_opt/RAW/with_MAXCANDIDATES/Fr14/0.9/xi_output/xi_results.csv"
):
    df = pd.read_csv(csv_to_guess_dtypes)
    # df.dtypes.to_csv(out_csv)
    return df.dtypes
# TODO: make default paths point to actual data, change to relative


def read_dtypes_from_file(
        dtypes_csv_file=r"/home/henning/ownCloud/masterieren/Data/Results/170323_iBAQ_based_opt/RAW/with_MAXCANDIDATES/Fr14/0.9/xi_output/dtypes.csv"
):
    dtypes = pd.read_csv(
        dtypes_csv_file,
        index_col=0,
        header=None
    )
    # print dtypes
    dtypes = pd.Series(dtypes.iloc[:, 0], index=dtypes.index).to_dict()
    return dtypes
# TODO: make default paths point to actual data, change to relative


def xifdr_result_to_spectra_df(xifdr_results):
    """
    takes a list of xifdr result files
    picks out the "_false_PSM_" file
    removes all entries containing decoys
    writes column "run" and "scan" to new DataFrame
    RETURN:
        DataFrame containing run and scan of xifdr hits
    """
    input_files = [s for s in xifdr_results if "_false_PSM_" in s]
    if len(input_files) > 1:
        raise AttributeError("More than one candidate file in 'xifdr results' matching '*_false_PSM_*'")
    elif len(input_files) < 1:
        raise AttributeError("No candidate file in 'xifdr results' matching '*_false_PSM_*'")

    df = pd.read_csv(input_files[0])
    # do not keep decoy hits
    df = df[~df["isDecoy"]]
    df = df.loc[:, ("run", "scan")]
    return df


def rm_xifdr_spectra_from_xi_results(
        df_of_spectra_to_remove,
        lst_dct_ordered_xi,
        xi_result_dir,
        **kwargs
):
    """
    Function that removes spectra specified in df_of_spectra_to_remove from all xi_results in lst_dct_ordered_xi

    :param df_of_spectra_to_remove:
    :param lst_dct_ordered_xi: untreated xi_results
    :param xi_result_dir: output folder for treated xi_results
    :param kwargs:
    :return: list of dicts with cleaned xi_result files
    """
    df_fdr = df_of_spectra_to_remove
    lst_dct_xi_filtered = []

    # read dtypes for xiresult files, either from from standard dtypes.csv, guess from standard xiresult.csv or from
        # provides xiresult.csv
    """dtypes is not respected as this leads to loss of precision compared to using dtype='object'"""
    # if "dtypes_csv_file" in kwargs.keys():
    #     dtypes = read_dtypes_from_file(kwargs["dtypes_csv_file"])
    # elif "csv_to_guess_dtypes" in kwargs.keys():
    #     dtypes = guess_dtypes_of_columns(kwargs["csv_to_guess_dtypes"])
    # else:
    #     dtypes = read_dtypes_from_file()

    # iterate over xi_results and filter out spectra to keep
    for dct_xi_result in lst_dct_ordered_xi:
        xi_result = dct_xi_result['filename']

        # create subdir, named after specific dir, where xi_result is from
        subdir = dct_xi_result['subdir']

        result_dir = os.path.join(xi_result_dir, subdir)
        result_file = os.path.join(result_dir, "xi_results.csv")
        if not os.path.exists(result_dir):
            os.makedirs(result_dir)

        # read xi_result into pandas dataframe, respecting dtypes
        df_xi = pd.read_csv(xi_result, dtype='object', float_precision='high')

        # remove spectra specified by df_of_spectra_to_remove from xi_result_dataframe
        logging.info("Filtering {}...".format(xi_result))
        df_xi_filtered = df_xi[~(df_xi["Run"]+" "+df_xi["Scan"].map(str)).isin(df_fdr["run"]+" "+df_fdr["scan"].map(str))]
        logging.warning("Column 'OpenModWindow' occurs twice, second occurence is therefore renamed to 'OpenModWindow.1'")
        # raise Exception("pandas rounds floats from xi_results! Stop this!")

        # store xi_result_dataframe in csv in created subdir
        df_xi_filtered.to_csv(result_file, index=False)
        lst_dct_xi_filtered.append({'filename': result_file, 'subdir': subdir})

    return lst_dct_xi_filtered

# # # Testing
# df_fdr = xifdr_result_to_spectra_df([r"/home/henning/ownCloud/masterieren/Data/Results/170323_iBAQ_based_opt/RAW/with_MAXCANDIDATES/Fr14/0.9/xifdr_output/FDR_1.000_0.050_1.000_1.000_1.000_10000.000_false_PSM_xiFDR1.0.14.csv"])
# rm_xifdr_spectra_from_xi_results(
#     df_of_spectra_to_remove=df_fdr,
#     ordered_list_of_xi_results=[r"/home/henning/ownCloud/masterieren/Data/Results/170323_iBAQ_based_opt/RAW/with_MAXCANDIDATES/Fr14/0.9/xi_output/xi_results.csv"],
#     xi_result_dir=r"test/xi_filtered",
#     dtypes_csv_file=r"/home/henning/ownCloud/masterieren/Data/Results/170323_iBAQ_based_opt/RAW/with_MAXCANDIDATES/Fr14/0.9/xi_output/dtypes_object.csv"
# )


def rm_links_explainable_by_fasta(xi_result, fasta, result_dir):
    """
    removes peptide links from xi_result, where both precursor proteins are in fasta
    :param xi_result:
    :param fasta:
    :param result_dir: where to write "xi_results.csv" to
    :return: cleaned xi_result
    """
    fasta_object = FastaHandler(fasta_filename=fasta, re_id_pattern=r'^>(.*)$')
    df_xi = pd.read_csv(xi_result, dtype='object', float_precision='high')
    file_result = os.path.join(result_dir, "xi_results.csv")
    file_dropped_scans = os.path.join(result_dir, "dropped_xi_results.csv")

    def are_proteins_in_fasta(prot_a, prot_b):
        if (type(prot_a) is str) and (type(prot_b) is str):
            prot_a = prot_a.replace('REV_', '')
            prot_b = prot_b.replace('REV_', '')
            if (prot_a in fasta_object.dict) and (prot_b in fasta_object.dict):
                return True
        return False
    # lst_rows_to_drop = []
    ser_to_drop = pd.Series(index=df_xi.index, dtype='bool')
    for index, row in df_xi.iterrows():
        # lst_rows_to_drop.append(are_proteins_in_fasta(row['Protein1'], row['Protein2']))
        ser_to_drop[index] = are_proteins_in_fasta(row['Protein1'], row['Protein2'])
    logging.info("Dropped {} columns from xi_result '{}' that were entirely in fasta '{}'"
                 .format(sum(ser_to_drop), xi_result, fasta))
    df_xi_results_cleaned = df_xi[~ser_to_drop]
    df_xi_results_dropped = df_xi[ser_to_drop]

    df_xi_results_cleaned.to_csv(file_result, index=False)
    df_xi_results_dropped.to_csv(file_dropped_scans, index=False)
    return file_result

# # # Testing
# logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
# rm_links_explainable_by_fasta(
#     xi_result=r"/home/henning/mnt/xitu/Data/Results/170323_iBAQ_based_opt/170704-reviewed_fasta-majority_protein_ids-new_iBAQ/Fr16/Fr16/log_norm_-5.0/xi_output/xi_results.csv",
#     fasta=r"/home/henning/mnt/xitu/Data/Results/170323_iBAQ_based_opt/170704-reviewed_fasta-majority_protein_ids-new_iBAQ/Fr16/Fr16/log_norm_-4.0/fasta.fasta",
#     result_dir=r"test"
# )


def simulation_for_single_exp(
        lst_ordered_xi,
        xifdr_settings_dict,
        out_dir
):
    """
    Iterate over lst_ordered_xi, calculate FDR fopr each and remove spectra satisfying FDR from all following xiresults

    """
    lst_dct_ordered_xi = []
    dct_fasta = {}

    # build list of dicts of xi-results with keys:
    #   'file': input xi_result file
    #   'subdir': xi_result-specific subdir
    #   'fasta': xi_result-specific fasta
    while lst_ordered_xi:
        file_dict = {}
        filename = lst_ordered_xi.pop(0)
        file_dict['filename'] = filename

        # split off everything after last "/" twice, afterwards take what is behind last "/" as subdir
        for i in range(2):
            filename = os.path.split(filename)[0]
        subdir = os.path.split(filename)[1]
        file_dict['subdir'] = subdir
        dct_fasta[subdir] = os.path.join(filename, "fasta.fasta")

        lst_dct_ordered_xi.append(file_dict)

    logging.debug("Xi Dictionary: {}".format(lst_dct_ordered_xi))
    while lst_dct_ordered_xi:
        """lst_ordered_xi has to be ordered by increasing DB size"""
        dct_xi_result = lst_dct_ordered_xi.pop(0)

        # read variables from dict
        xi_result = dct_xi_result['filename']
        subdir = dct_xi_result['subdir']
        fasta = dct_fasta[subdir]

        # create dirname for this specific iteration of loop
        result_dir = os.path.join(out_dir, subdir)
        xifdr_out_dir = os.path.join(result_dir, "xifdr_output")
        xifdr_results = XiFdrWrapper.xifdr_execution(
            xifdr_input_csv=xi_result,
            xifdr_output_dir=xifdr_out_dir,
            pepfdr=xifdr_settings_dict['pepfdr'],
            reportfactor=xifdr_settings_dict['reportfactor'],
            additional_xifdr_arguments=xifdr_settings_dict['additional_xifdr_arguments']
        )

        # execute if there are remaining xiresults
        if lst_dct_ordered_xi:
            df_of_spectra_to_remove = xifdr_result_to_spectra_df(xifdr_results)
            xi_result_dir = os.path.join(result_dir, "xi_output")
            logging.info("Filtering out spectra satisfying FDR in '{}'".format(xi_result))
            starttime = time.time()
            lst_dct_ordered_xi = \
                rm_xifdr_spectra_from_xi_results(
                    df_of_spectra_to_remove,
                    lst_dct_ordered_xi,
                    xi_result_dir
                )

            # remove matches that can entirely be explained with the current fasta from the next xi_result
            fasta_filtered_xi_result_dir = os.path.join(xi_result_dir, "fasta_filtered")
            lst_dct_ordered_xi[0]['filename'] = rm_links_explainable_by_fasta(
                xi_result=lst_dct_ordered_xi[0]['filename'],
                fasta=fasta,
                result_dir=fasta_filtered_xi_result_dir
            )

            logging.info("Spectra filtering took {}".format(calculate_elapsed_time(starttime)))


def exp_iterator(list_of_experiments, xifdr_settings_dict, out_dir, **kwargs):
    for exp in list_of_experiments:
        # set output basedir for this experiment
        exp_out_dir = os.path.join(out_dir, exp.name)
        logging.info("Starting simulation for '{}'"
                     .format(exp.name))
        starttime = time.time()
        simulation_for_single_exp(
            lst_ordered_xi=exp.ordered_list_of_xi_results,
            xifdr_settings_dict=xifdr_settings_dict,
            out_dir=exp_out_dir
        )
        logging.info("Simulation for '{}' took {}"
                     .format(exp.name, calculate_elapsed_time(starttime)))
        # TODO integrate hand over of kwargs to simulation


def logging_setup(log_file):

    str_format = '%(asctime)s - %(levelname)s - %(name)s - %(message)s'
    logging.basicConfig(filename=log_file, level=logging.DEBUG,
                        format=str_format)
    logger = logging.getLogger('')

    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter(str_format)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)
    return logger


def main():
    # print help message if script is called without argument
    if len(sys.argv) != 2:
        print """Script has to be called with output dir as argument.
                Output dir has to contain config file "config.py"."""
        sys.exit(1)

    # set output dir
    output_basedir = sys.argv[1]

    # set up logging
    log_file = os.path.join(output_basedir, __name__ + '.log')
    global logger
    logger = logging_setup(log_file=log_file)

    # import of config.py
    sys.path.append(output_basedir)
    import myconfig

    list_of_experiment_dicts = myconfig.list_of_experiments
    xi_xifdr_settings_dict = myconfig.xifdr_settings_dict

    list_of_experiments = []
    for experiment_dict in list_of_experiment_dicts:
        list_of_experiments.append(
            Experiment(
                exp_name=experiment_dict['exp_name'],
                ordered_list_of_xi_results=experiment_dict['ordered_list_of_xi_results']
            )
        )

    starttime = time.time()
    exp_iterator(
        list_of_experiments=list_of_experiments,
        xifdr_settings_dict=xi_xifdr_settings_dict,
        out_dir=output_basedir
    )
    logging.info("Script execution took {}"
                 .format(calculate_elapsed_time(starttime)))
    logging.shutdown()


if __name__ == "__main__":
    main()
