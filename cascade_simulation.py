"""
This script simulates a cascade search on already existing XiSearch results.

This is done by carrying out a xifdr evaluation of the XiResults and removing scans of certain runs that satisfy
FDR-conditions from all successive XiResults PSMs. This process is repeated for every successive XiResult.
"""

import pandas as pd
import timeit
import os
from lib.pipeline import XiFdrWrapper, calculate_elapsed_time
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


def remove_spectra_from_xi_results(
        df_of_spectra_to_remove,
        lst_dct_ordered_xi,
        xi_result_dir,
        **kwargs
):
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
# remove_spectra_from_xi_results(
#     df_of_spectra_to_remove=df_fdr,
#     ordered_list_of_xi_results=[r"/home/henning/ownCloud/masterieren/Data/Results/170323_iBAQ_based_opt/RAW/with_MAXCANDIDATES/Fr14/0.9/xi_output/xi_results.csv"],
#     xi_result_dir=r"test/xi_filtered",
#     dtypes_csv_file=r"/home/henning/ownCloud/masterieren/Data/Results/170323_iBAQ_based_opt/RAW/with_MAXCANDIDATES/Fr14/0.9/xi_output/dtypes_object.csv"
# )


def simulation_for_single_exp(
        lst_ordered_xi,
        xifdr_settings_dict,
        out_dir
):
    # TODO what is this function doing? docstring
    lst_dct_ordered_xi = []

    # build list of dicts of xi-results with keys:
    #   'file': xi_results
    #   'subdir': xi_result-specific subdir
    while lst_ordered_xi:
        file_dict = {}
        filename = lst_ordered_xi.pop(0)
        file_dict['filename'] = filename
        for i in range(2):
            filename = os.path.split(filename)[0]
        file_dict['subdir'] = os.path.split(filename)[1]
        lst_dct_ordered_xi.append(file_dict)

    logging.debug("Xi Dictionary: {}".format(lst_dct_ordered_xi))
    while lst_dct_ordered_xi:
        """lst_ordered_xi has to be ordered by increasing DB size"""
        dct_xi_result = lst_dct_ordered_xi.pop(0)

        # read variables from dict
        xi_result = dct_xi_result['filename']
        subdir = dct_xi_result['subdir']

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
        if lst_dct_ordered_xi:
            df_of_spectra_to_remove = xifdr_result_to_spectra_df(xifdr_results)
            xi_result_dir = os.path.join(result_dir, "xi_output")
            logging.info("Filtering out spectra satisfying FDR in '{}'".format(xi_result))
            starttime = time.time()
            lst_dct_ordered_xi = \
                remove_spectra_from_xi_results(
                    df_of_spectra_to_remove,
                    lst_dct_ordered_xi,
                    xi_result_dir
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
