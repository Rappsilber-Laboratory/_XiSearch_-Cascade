"""
This script simulates a cascade search on already existing XiSearch results.

This is done by carrying out a xifdr evaluation of the XiResults and removing scans of certain runs that satisfy
FDR-conditions from all successive XiResults PSMs. This process is repeated for every successive XiResult.
"""

import pandas as pd
import timeit
import os
from lib.pipeline import XiFdrWrapper


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
        in_csv=r"/home/henning/ownCloud/masterieren/Data/Results/170323_iBAQ_based_opt/RAW/with_MAXCANDIDATES/Fr14/0.9/xi_output/xi_results.csv",
        out_csv=r"/home/henning/ownCloud/masterieren/Data/Results/170323_iBAQ_based_opt/RAW/with_MAXCANDIDATES/Fr14/0.9/xi_output/dtypes.csv"
):
    df = pd.read_csv(in_csv)
    df.dtypes.to_csv(out_csv)
    return out_csv
# TODO: make default paths point to actual data, change to relative


def read_dtypes_from_file(
        csv_file=r"/home/henning/ownCloud/masterieren/Data/Results/170323_iBAQ_based_opt/RAW/with_MAXCANDIDATES/Fr14/0.9/xi_output/dtypes.csv"
):
    dtypes = pd.read_csv(
        csv_file,
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
    df = pd.read_csv(input_files[0])
    # do not keep decoy hits
    df = df[df["isDecoy"] == False]
    df = df.loc[:, ("run", "scan")]
    return df


df = pd.read_csv(filename, dtype=dtypes)
print df.dtypes

def remove_spectra_from_xi_results(df_of_spectra_to_remove, ordered_list_of_xi_results, xi_result_dir, **kwargs):
    # read dtypes for xiresult files, either from from standard dtypes.csv, guess from standard xiresult.csv or from
        # provides xiresult.csv
    # for xi_result in ordered_list_of_xi_results:
        # create subdir, named after specific dir, where xi_result is from
        # read xi_result into pandas dataframe, respecting dtypes
        # remove spectra specified by df_of_spectra_to_remove from xi_result_dataframe
        # store xi_result_dataframe in csv in created subdir
    # return list_of_processed_xi_results
    

def simulation_for_single_exp(ordered_list_of_xi_results, xifdr_settings_dict, out_dir):
    while ordered_list_of_xi_results: # TODO macht "1" Sinn?
        """ordered_list_of_xi_results has to be ordered by increasing DB size"""
        xi_result = ordered_list_of_xi_results[0]
        remaining_ordered_list_of_xi_results = ordered_list_of_xi_results[1:]

        # read subdir out of filepath
        filename = xi_result
        for i in range(2):
            filename = os.path.split(filename)[0]
        subdir = os.path.split(filename)[1]

        # create dirname for this specific iteration of loop
        result_dir = os.path.join(out_dir, subdir)
        xifdr_out_dir = os.path.join(result_dir, "xifdr_output")
        xifdr_results = XiFdrWrapper.xifdr_execution(
            xifdr_input_csv=xi_result,
            xifdr_output_dir=xifdr_out_dir,
            pepfdr=xifdr_settings_dict['xifdr_settings']['pepfdr'],
            reportfactor=xifdr_settings_dict['xifdr_settings']['reportfactor'],
            additional_xifdr_arguments=xifdr_settings_dict['xifdr_settings']['additional_xifdr_arguments']
        )
        if remaining_ordered_list_of_xi_results:
            df_of_spectra_to_remove = xifdr_result_to_spectra_df(xifdr_results)
            xi_result_dir = os.path.join(result_dir, "xi_output")
            remaining_ordered_list_of_xi_results = remove_spectra_from_xi_results(df_of_spectra_to_remove, remaining_ordered_list_of_xi_results, xi_result_dir)
        ordered_list_of_xi_results = remaining_ordered_list_of_xi_results
#

def main_function(list_of_experiments, xifdr_settings_dict, out_dir):
    for exp in list_of_experiments:
        # set output basedir for this experiment
        exp_out_dir = os.path.join(out_dir, exp.name)
        simulation_for_single_exp(
            ordered_list_of_xi_results=exp.ordered_tuple_of_xi_results,
            xifdr_settings_dict=xifdr_settings_dict,
            out_dir=exp_out_dir
        )

if __name__ == "__main__":
    # read command line args
    # initiate logging
    # read config
    # execute main function
    # shutdown logging
