"""
This script simulates a cascade search on already existing XiSearch results.

This is done by carrying out a xifdr evaluation of the XiResults and removing scans of certain runs that satisfy
FDR-conditions from all successive XiResults PSMs. This process is repeated for every successive XiResult.
"""

import configparser
import pandas as pd
import os
from lib.pipeline import XiFdrWrapper, calculate_elapsed_time
from lib.iBAQ_FASTA_handler import FastaHandler
import logging
import sys
import time
import re

g_allowed_modification = ["Mox", "bs3"]


# OBJECT: Experiment: holds all the necessary settings for each Fraction/replicate
class Experiment:
    """
    one object per MS analysis
    """
    def __init__(
            self,
            ordered_list_of_xi_results,
            exp_name,
            ordered_list_of_fastas=None,
            ordered_list_of_subdirs=None,
            ordered_list_of_missing_modifications=None,
            rm_links_explainable_by_prev_db=True,
            unmod_mod_cascade=False
    ):
        self.name = exp_name
        self.rm_links_explainable_by_prev_db = rm_links_explainable_by_prev_db
        self.unmod_mod_cascade = unmod_mod_cascade
        self.ordered_list_of_xi_results = ordered_list_of_xi_results
        self.ordered_list_of_fastas = self.construct_ordered_list_of_fastas(ordered_list_of_fastas)
        self.ordered_list_of_subdirs = self.construct_ordered_list_of_subdirs(ordered_list_of_subdirs)
        self.all_files = self.ordered_list_of_xi_results + self.ordered_list_of_fastas
        self.ordered_list_of_missing_modifications = \
            self.construct_ordered_list_of_missing_modifications(ordered_list_of_missing_modifications)

    def construct_ordered_list_of_missing_modifications(self, ordered_list_of_missing_modifications):
        if self.unmod_mod_cascade:
            if ordered_list_of_missing_modifications is None:
                assert len(self.ordered_list_of_xi_results) == 2, \
                    ("For unmodified versus modified searches without specified missing modifications,\n" +
                     "you must specify 2 xi results:\n" +
                     "\t 1. unmodified xiresult (without [bs3, Mox])\n" +
                     "\t 2. modified xiresult")
                res_ordered_list_of_missing_modifications = [["Mox", "bs3"], None]
            else:
                res_ordered_list_of_missing_modifications = ordered_list_of_missing_modifications
        else:
            res_ordered_list_of_missing_modifications = list()
        for mods in res_ordered_list_of_missing_modifications:
            if mods is not None:
                assert set(mods).issubset(set(g_allowed_modification)), \
                    "allowed modifications are: {}".format(g_allowed_modification)
        return res_ordered_list_of_missing_modifications

    def construct_ordered_list_of_subdirs(self, ordered_list_of_subdirs):
        if ordered_list_of_subdirs is None:
            ordered_list_of_subdirs = list()
            for f in self.ordered_list_of_xi_results:
                # split off everything after last "/" twice, afterwards take what is behind last "/" as subdir
                # for "/Data/Results/170323_iBAQ_based_opt/Ribosome/Ribosome/fasta_0.01/xi_output/xi_results.csv" this would
                #   lead to "fasta_0.01" as subdir
                for i in range(2):
                    f = os.path.split(f)[0]
                subdir = os.path.split(f)[1]
                ordered_list_of_subdirs.append(subdir)
            return ordered_list_of_subdirs
        elif isinstance(ordered_list_of_subdirs, list):
            if len(ordered_list_of_subdirs) == len(self.ordered_list_of_xi_results):
                return ordered_list_of_subdirs
            else:
                raise AttributeError("ordered_list_of_subdirs must be of same length as ordered_list_of_xi_results")
        else:
            raise AttributeError("ordered_list_of_subdirs must be list or None but is type: {}"
                                 .format(type(ordered_list_of_subdirs)))

    def construct_ordered_list_of_fastas(self, ordered_list_of_fastas):
        if (ordered_list_of_fastas is None) and (self.rm_links_explainable_by_prev_db is True):
            ordered_list_of_fastas = list()
            for filename in self.ordered_list_of_xi_results:
                # split off everything after last "/" twice, afterwards take what is behind last "/" as subdir
                # for "/Data/Results/170323_iBAQ_based_opt/Ribosome/Ribosome/fasta_0.01/xi_output/xi_results.csv" this would
                #   lead to "fasta_0.01" as subdir
                for i in range(2):
                    filename = os.path.split(filename)[0]
                subdir = os.path.split(filename)[1]
                fasta = os.path.join(filename, "fasta.fasta")

                ordered_list_of_fastas.append(fasta)

            return ordered_list_of_fastas
        elif isinstance(ordered_list_of_fastas, list) \
                and (len(ordered_list_of_fastas) == len(self.ordered_list_of_xi_results)):
            return ordered_list_of_fastas
        elif self.rm_links_explainable_by_prev_db is False:
            return list()
        else:
            raise AttributeError("ordered_list_of_fastas must be list or None but is type: {}"
                                 .format(type(ordered_list_of_fastas)))

    def check_files_exist(self):
        lst_unfound_files = []
        for f in self.all_files:
            b = os.path.exists(f)
            if not b:
                lst_unfound_files.append(f)
        if len(lst_unfound_files) > 0:
            raise IOError("File(s) do(es) not exist: \n" + "\n".join(lst_unfound_files))
        return True


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
    re_pattern = r"(false|true)_PSM_xiFDR"
    input_files = [s for s in xifdr_results if re.search(re_pattern, s)]
    if len(input_files) > 1:
        raise AttributeError("More than one candidate file in 'xifdr results' matching regex '{}'".format(re_pattern))
    elif len(input_files) < 1:
        raise AttributeError("No candidate file in 'xifdr results' matching regex '{}'".format(re_pattern))

    df = pd.read_csv(input_files[0])
    # do not keep decoy hits
    df = df[~df["isDecoy"]]
    df = df.loc[:, ("run", "scan")]
    return df


# def rm_xifdr_spectra_from_xi_results(
#         df_of_spectra_to_remove,
#         lst_dct_ordered_xi,
#         xi_result_dir,
#         **kwargs
# ):
#     """
#     Function that removes spectra specified in df_of_spectra_to_remove from all xi_results in lst_dct_ordered_xi
#
#     :param df_of_spectra_to_remove: pandas Dataframe with columns ["run", "scan"]. Column-scan-combination of each
#         row is removed from each xi_result in lst_dct_ordered_xi
#     :param lst_dct_ordered_xi: untreated xi_results
#     :param xi_result_dir: output folder for treated xi_results
#     :param kwargs:
#     :return: list of dicts with cleaned xi_result files
#     """
#     df_fdr = df_of_spectra_to_remove
#     lst_dct_xi_filtered = []
#
#     # read dtypes for xiresult files, either from from standard dtypes.csv, guess from standard xiresult.csv or from
#         # provides xiresult.csv
#     """dtypes is not respected as this leads to loss of precision compared to using dtype='object'"""
#     # if "dtypes_csv_file" in kwargs.keys():
#     #     dtypes = read_dtypes_from_file(kwargs["dtypes_csv_file"])
#     # elif "csv_to_guess_dtypes" in kwargs.keys():
#     #     dtypes = guess_dtypes_of_columns(kwargs["csv_to_guess_dtypes"])
#     # else:
#     #     dtypes = read_dtypes_from_file()
#
#     # iterate over xi_results and filter out spectra to keep
#     for dct_xi_result in lst_dct_ordered_xi:
#         xi_result = dct_xi_result['filename']
#
#         # create subdir, named after specific dir, where xi_result is from
#         subdir = dct_xi_result['subdir']
#
#         result_dir = os.path.join(xi_result_dir, subdir)
#         result_file = os.path.join(result_dir, "xi_results.csv")
#         if not os.path.exists(result_dir):
#             os.makedirs(result_dir)
#
#         # read xi_result into pandas dataframe, respecting dtypes
#         df_xi = pd.read_csv(xi_result, dtype='object', float_precision='high')
#
#         # remove spectra specified by df_of_spectra_to_remove from xi_result_dataframe
#         logging.info("Filtering {}...".format(xi_result))
#         df_xi_filtered = df_xi[~(df_xi["Run"]+" "+df_xi["Scan"].map(str)).isin(df_fdr["run"]+" "+df_fdr["scan"].map(str))]
#         logging.warning("Column 'OpenModWindow' occurs twice, second occurence is therefore renamed to 'OpenModWindow.1'")
#         # raise Exception("pandas rounds floats from xi_results! Stop this!")
#
#         # store xi_result_dataframe in csv in created subdir
#         df_xi_filtered.to_csv(result_file, index=False)
#         lst_dct_xi_filtered.append({'filename': result_file, 'subdir': subdir})
#
#     return lst_dct_xi_filtered

# # # Testing
# df_fdr = xifdr_result_to_spectra_df([r"/home/henning/ownCloud/masterieren/Data/Results/170323_iBAQ_based_opt/RAW/with_MAXCANDIDATES/Fr14/0.9/xifdr_output/FDR_1.000_0.050_1.000_1.000_1.000_10000.000_false_PSM_xiFDR1.0.14.csv"])
# rm_xifdr_spectra_from_xi_results(
#     df_of_spectra_to_remove=df_fdr,
#     ordered_list_of_xi_results=[r"/home/henning/ownCloud/masterieren/Data/Results/170323_iBAQ_based_opt/RAW/with_MAXCANDIDATES/Fr14/0.9/xi_output/xi_results.csv"],
#     xi_result_dir=r"test/xi_filtered",
#     dtypes_csv_file=r"/home/henning/ownCloud/masterieren/Data/Results/170323_iBAQ_based_opt/RAW/with_MAXCANDIDATES/Fr14/0.9/xi_output/dtypes_object.csv"
# )


def rm_xifdr_spectra_from_xi_result(
        df_of_spectra_to_remove,
        xi_result,
        xi_result_dir,
        **kwargs
):
    """
    Function that removes spectra specified in df_of_spectra_to_remove from all xi_results in lst_dct_ordered_xi

    :param df_of_spectra_to_remove: pandas Dataframe with columns ["run", "scan"]. Column-scan-combination of each
        row is removed from xi_result
    :param xi_result: path to xi_result.csv which should be filtered
    :param xi_result_dir: output folder for treated xi_results
    :param subdir:
    :param kwargs:
    :return: list of dicts with cleaned xi_result files
    """
    df_fdr = df_of_spectra_to_remove

    # read dtypes for xiresult files, either from from standard dtypes.csv, guess from standard xiresult.csv or from
        # provides xiresult.csv
    """dtypes is not respected as this leads to loss of precision compared to using dtype='object'"""
    # if "dtypes_csv_file" in kwargs.keys():
    #     dtypes = read_dtypes_from_file(kwargs["dtypes_csv_file"])
    # elif "csv_to_guess_dtypes" in kwargs.keys():
    #     dtypes = guess_dtypes_of_columns(kwargs["csv_to_guess_dtypes"])
    # else:
    #     dtypes = read_dtypes_from_file()
    result_file = os.path.join(xi_result_dir, "xi_results.csv")
    if not os.path.exists(xi_result_dir):
        os.makedirs(xi_result_dir)

    # read xi_result into pandas dataframe, respecting dtypes
    df_xi = pd.read_csv(xi_result, dtype='object', float_precision='high')

    # remove spectra specified by df_of_spectra_to_remove from xi_result_dataframe
    logging.info("Filtering {}...".format(xi_result))
    df_xi_filtered = df_xi[~(df_xi["Run"]+" "+df_xi["Scan"].map(str)).isin(df_fdr["run"]+" "+df_fdr["scan"].map(str))]
    logging.warning("Column 'OpenModWindow' occurs twice, second occurence is therefore renamed to 'OpenModWindow.1'")
    # raise Exception("pandas rounds floats from xi_results! Stop this!")

    # store xi_result_dataframe in csv in created subdir
    df_xi_filtered.to_csv(result_file, index=False)

    return result_file


def are_proteins_in_fasta(prot_a, prot_b, fasta_object):
    if (type(prot_a) is str) and (type(prot_b) is str):
        prot_a = prot_a.replace('REV_', '')
        prot_b = prot_b.replace('REV_', '')
        if (prot_a in fasta_object.dict) and (prot_b in fasta_object.dict):
            return True
    return False


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

    # lst_rows_to_drop = []
    ser_to_drop = pd.Series(index=df_xi.index, dtype=bool)
    for index, row in df_xi.iterrows():
        # lst_rows_to_drop.append(are_proteins_in_fasta(row['Protein1'], row['Protein2']))
        ser_to_drop[index] = are_proteins_in_fasta(row['Protein1'], row['Protein2'], fasta_object)
    logging.info("Dropped {} rows from xi_result '{}' that were entirely in fasta '{}'"
                 .format(sum(ser_to_drop), xi_result, fasta))
    file_result, file_dropped_scans = write_cleaned_and_not_cleaned_xi_results(df_xi, result_dir, ser_to_drop)
    return file_result


def write_cleaned_and_not_cleaned_xi_results(df_xi, result_dir, ser_to_drop):
    file_result = os.path.join(result_dir, "xi_results.csv")
    file_dropped_scans = os.path.join(result_dir, "dropped_xi_results.csv")

    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    df_xi_results_cleaned = df_xi[~ser_to_drop]
    df_xi_results_dropped = df_xi[ser_to_drop]
    df_xi_results_cleaned.to_csv(file_result, index=False)
    df_xi_results_dropped.to_csv(file_dropped_scans, index=False)
    return file_result, file_dropped_scans


# # # Testing
# logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
# rm_links_explainable_by_fasta(
#     xi_result=r"/home/henning/mnt/xitu/Data/Results/170323_iBAQ_based_opt/170704-reviewed_fasta-majority_protein_ids-new_iBAQ/Fr16/Fr16/log_norm_-5.0/xi_output/xi_results.csv",
#     fasta=r"/home/henning/mnt/xitu/Data/Results/170323_iBAQ_based_opt/170704-reviewed_fasta-majority_protein_ids-new_iBAQ/Fr16/Fr16/log_norm_-4.0/fasta.fasta",
#     result_dir=r"test"
# )


def rm_links_not_modified_by(xi_result, result_dir, modifications):
    """
    Filter out xi-results that do not contain any of the specified modifications.

    :param xi_result: path to xi_results.csv
    :param result_dir: folder for filtered output xi_results.csv
    :param modifications: list of modifications that should be dropped
    :return:
    """
    assert set(modifications).issubset(set(g_allowed_modification)), \
        "allowed modifications are: {}".format(g_allowed_modification)
    df_xi = pd.read_csv(xi_result, dtype='object', float_precision='high')

    # construct string of both Peptides, second has to be converted as it can be nan. This would lead to sum of strings
    #    being nan, if not converted
    s_pep_str = df_xi["Peptide1"] + df_xi["Peptide2"].astype(str)

    regex_pattern = "|".join(modifications)
    ser_to_drop = ~(s_pep_str.str.contains(regex_pattern))

    file_result, file_dropped_scans = write_cleaned_and_not_cleaned_xi_results(df_xi, result_dir, ser_to_drop)

    logging.info("Dropped {} rows from xi_result '{}' that were entirely unmodified"
                 .format(sum(ser_to_drop), xi_result))

    return file_result


def simulation_for_single_exp(
        exp,
        xifdr_settings_dict,
        out_dir
):
    """
    Iterate over lst_ordered_xi, calculate FDR for each and remove spectra satisfying FDR from all following xiresults

    """
    rm_links_explainable_by_prev_db = exp.rm_links_explainable_by_prev_db
    unmod_mod_cascade = exp.unmod_mod_cascade
    assert isinstance(rm_links_explainable_by_prev_db, bool), "rm_links_explainable_by_prev_db must be of type bool"
    assert isinstance(unmod_mod_cascade, bool), "unmod_mod_cascade must be of type bool"
    assert unmod_mod_cascade & rm_links_explainable_by_prev_db is not True, \
        "unmod_mod_cascade and rm_links_explainable_by_prev_db are exclusive options."

    # build list of dicts of xi-results with keys:
    #   'file': input xi_result file
    #   'subdir': xi_result-specific subdir
    #   'fasta': xi_result-specific fasta
    lst_dct_ordered_xi = build_lst_dct_ordered(exp)

    logging.debug("List of Xi dictionaries: {}".format(lst_dct_ordered_xi))
    while lst_dct_ordered_xi:
        """lst_ordered_xi has to be ordered by increasing DB size"""
        dct_xi_result_this_run = lst_dct_ordered_xi.pop(0)

        # read variables from dict
        xi_result = dct_xi_result_this_run['filename']
        subdir = dct_xi_result_this_run['subdir']

        # create dirname for this specific iteration of loop
        result_dir = os.path.join(out_dir, subdir)
        xifdr_out_dir = os.path.join(result_dir, "xifdr_output")
        xifdr_results = XiFdrWrapper.xifdr_execution(
            xifdr_input_csv=xi_result,
            xifdr_output_dir=xifdr_out_dir,
            pepfdr=xifdr_settings_dict['pepfdr'],
            reportfactor=xifdr_settings_dict['reportfactor'],
            additional_xifdr_arguments=xifdr_settings_dict['additional_xifdr_arguments'],
            xifdr_filename=xifdr_settings_dict['executable']
        )

        # execute if there are remaining xiresults
        if lst_dct_ordered_xi:
            df_of_spectra_to_remove = xifdr_result_to_spectra_df(xifdr_results)
            xi_result_dir = os.path.join(result_dir, "filtered_xi_outputs")
            logging.info("Filtering out spectra from following xi-results that satisfy FDR in '{}'".format(xi_result))
            starttime = time.time()

            # rm spectra that passed fdr
            new_lst_dct_ordered_xi = []
            for xi_dct_unanalyzed in lst_dct_ordered_xi:
                xi_res = xi_dct_unanalyzed['filename']
                subdir = xi_dct_unanalyzed['subdir']
                ith_xi_result_dir = os.path.join(xi_result_dir, "spectra_filtered", subdir)
                xi_res_filtered = rm_xifdr_spectra_from_xi_result(
                    df_of_spectra_to_remove=df_of_spectra_to_remove,
                    xi_result=xi_res,
                    xi_result_dir=ith_xi_result_dir,
                    subdir=subdir
                )
                xi_dct_unanalyzed['filename'] = xi_res_filtered
                new_lst_dct_ordered_xi.append(xi_dct_unanalyzed)

            lst_dct_ordered_xi = new_lst_dct_ordered_xi

            # remove matches that can entirely be explained with the current fasta from the next xi_result
            if rm_links_explainable_by_prev_db:
                fasta_already_searched = dct_xi_result_this_run['fasta']
                fasta_filtered_xi_result_dir = os.path.join(xi_result_dir, "fasta_filtered")
                next_xi_result = lst_dct_ordered_xi[0]['filename']
                lst_dct_ordered_xi[0]['filename'] = rm_links_explainable_by_fasta(
                    xi_result=next_xi_result,
                    fasta=fasta_already_searched,
                    result_dir=fasta_filtered_xi_result_dir
                )

            # remove matches that do not contain modifcations missing in "this run"
            if unmod_mod_cascade:
                missing_mods = dct_xi_result_this_run["missing_mods"]
                modification_filtered_xi_result_dir = os.path.join(xi_result_dir, "modification_filtered")
                next_xi_result = lst_dct_ordered_xi[0]['filename']
                lst_dct_ordered_xi[0]['filename'] = rm_links_not_modified_by(
                    xi_result=next_xi_result,
                    result_dir=modification_filtered_xi_result_dir,
                    modifications=missing_mods
                )

            logging.info("Spectra filtering took {}".format(calculate_elapsed_time(starttime)))


def build_lst_dct_ordered(exp):
    lst_ordered_xi = list(exp.ordered_list_of_xi_results)
    lst_ordered_fasta = list(exp.ordered_list_of_fastas)
    lst_ordered_subdirs = list(exp.ordered_list_of_subdirs)
    lst_ordered_missing_modifications = list(exp.ordered_list_of_missing_modifications)

    lst_dct_ordered_xi = []
    save_fastas = len(lst_ordered_fasta) > 0
    save_missing_mods = len(lst_ordered_missing_modifications) > 0

    while lst_ordered_xi:
        file_dict = {}
        file_dict['filename'] = lst_ordered_xi.pop(0)

        subdir = lst_ordered_subdirs.pop(0)
        file_dict['subdir'] = subdir
        if save_fastas:
            file_dict['fasta'] = lst_ordered_fasta.pop(0)

        if save_missing_mods:
            file_dict["missing_mods"] = lst_ordered_missing_modifications.pop(0)

        lst_dct_ordered_xi.append(file_dict)
    return lst_dct_ordered_xi


def exp_iterator(list_of_experiments, xifdr_settings_dict, out_dir,
                 **kwargs):
    for exp in list_of_experiments:
        # set output basedir for this experiment
        exp_out_dir = os.path.join(out_dir, exp.name)
        logging.info("Starting simulation for '{}'"
                     .format(exp.name))
        starttime = time.time()
        simulation_for_single_exp(
            exp=exp,
            xifdr_settings_dict=xifdr_settings_dict,
            out_dir=exp_out_dir
        )
        logging.info("Simulation for '{f}' took {t}"
                     .format(f=exp.name, t=calculate_elapsed_time(starttime)))
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


def build_experiment_from_dict(experiment_dict, rm_links_explainable_by_prev_db, unmod_mod_cascade):
    ordered_list_of_fastas = None
    ordered_list_of_subdirs = None
    ordered_list_of_missing_modifications = None

    exp_name = experiment_dict['exp_name']
    ordered_list_of_xi_results = experiment_dict['ordered_list_of_xi_results']

    if 'ordered_list_of_fastas' in experiment_dict.keys():
        ordered_list_of_fastas = experiment_dict['ordered_list_of_fastas']
        assert len(ordered_list_of_fastas) == len(ordered_list_of_xi_results), \
            "Xi results and fasta files need to be of same number."

    if "ordered_list_of_subdirs" in experiment_dict.keys():
        ordered_list_of_subdirs = experiment_dict["ordered_list_of_subdirs"]
        assert len(ordered_list_of_subdirs) == len(ordered_list_of_xi_results), \
            "Xi results and subdirs must be of same number."

    if unmod_mod_cascade:
        if "ordered_list_of_missing_modifications" in experiment_dict.keys():
            ordered_list_of_missing_modifications = experiment_dict["ordered_list_of_missing_modifications"]
            assert len(ordered_list_of_missing_modifications) == len(ordered_list_of_xi_results), \
                "Xi results and missing modifications need to be of same number."
        else:
            ordered_list_of_missing_modifications = None

    o_exp = Experiment(
        exp_name=exp_name,
        ordered_list_of_xi_results=ordered_list_of_xi_results,
        ordered_list_of_fastas=ordered_list_of_fastas,
        ordered_list_of_subdirs=ordered_list_of_subdirs,
        ordered_list_of_missing_modifications=ordered_list_of_missing_modifications,
        rm_links_explainable_by_prev_db=rm_links_explainable_by_prev_db,
        unmod_mod_cascade=unmod_mod_cascade
    )
    return o_exp

# # # Test
# exp_dict = {
#     'exp_name': "Ribosome",
#     'ordered_list_of_xi_results': [
#         r"/home/henning/mnt/xitu/Data/Results/170323_iBAQ_based_opt/Ribosome/Ribosome/fasta_0.01/xi_output/xi_results.csv",
#         r"/home/henning/mnt/xitu/Data/Results/170323_iBAQ_based_opt/Ribosome/Ribosome/fasta_0.001/xi_output/xi_results.csv",
#         r"/home/henning/mnt/xitu/Data/Results/170323_iBAQ_based_opt/Ribosome/Ribosome/fasta_0.0001/xi_output/xi_results.csv"
#     ],
#     'ordered_list_of_fastas': [
#         r"/home/henning/mnt/xitu/Data/Input/170323_iBAQ_based_opt/fastas/180319-based_on_iBAQ/0.01.fasta",
#         r"/home/henning/mnt/xitu/Data/Input/170323_iBAQ_based_opt/fastas/180319-based_on_iBAQ/0.001.fasta",
#         r"/home/henning/mnt/xitu/Data/Input/170323_iBAQ_based_opt/fastas/180319-based_on_iBAQ/0.0001.fasta"
#     ]
# }
# exp = build_experiment_from_dict(exp_dict)
# exp.check_files_exist()


def parse_config(config_file):
    list_of_experiments = []
    rm_links_explainable_by_prev_db = True
    unmod_mod_cascade = False

    config = configparser.ConfigParser()
    config.read(config_file)

    assert "xifdr settings" in config, "'xifdr settings' section must be in config file"
    exp_keys = filter(lambda x: "Experiment" in x, config.sections())
    assert len(exp_keys) >= 1, "You must specifiy at least one experiment"

    # general section
    if "general settings" in config:
        if "rm_links_explainable_by_prev_db" in config["general settings"]:
            rm_links_explainable_by_prev_db = config.getboolean("general settings", "rm_links_explainable_by_prev_db")
        if "unmodified_vs_modified" in config["general settings"]:
            unmod_mod_cascade = config.getboolean("general settings", "unmodified_vs_modified")

    # xifdr settings
    xifdr_settings_dict = {
        'pepfdr': config["xifdr settings"]["pepfdr"],
        'reportfactor': config["xifdr settings"]["reportfactor"]
    }
    xifdr_settings_dict["executable"] = config.get("xifdr settings", "executable location")
    if not os.path.exists(xifdr_settings_dict["executable"]):
        raise IOError("xiFDR executable not found under: \n\t{}".format(xifdr_settings_dict["executable"]))

    if "additional xifdr arguments" in config:
        xifdr_settings_dict["additional_xifdr_arguments"] = [
            "--"+k+"="+v for k, v in config["additional xifdr arguments"].items()
        ]

    # Experiments
    for key in exp_keys:
        exp_sec_dct = config[key]
        exp_dict = {
            "exp_name": exp_sec_dct["name"]
        }

        lst_peak_keys = filter(lambda x: "peak_file_" in x, exp_sec_dct.keys())
        lst_peak_files = [exp_sec_dct[f] for f in sorted(lst_peak_keys)]
        exp_dict["ordered_list_of_xi_results"] = lst_peak_files

        lst_fasta_keys = filter(lambda x: "fasta_file_" in x, exp_sec_dct.keys())
        lst_fasta_files = [exp_sec_dct[f] for f in sorted(lst_fasta_keys)]
        if len(lst_fasta_files) >= 1:
            exp_dict["ordered_list_of_fastas"] = lst_fasta_files

        lst_subdir_keys = filter(lambda x: "subdir_" in x, exp_sec_dct.keys())
        lst_subdirs = [exp_sec_dct[f] for f in sorted(lst_subdir_keys)]
        if len(lst_subdirs) >= 1:
            exp_dict["ordered_list_of_subdirs"] = lst_subdirs

        if unmod_mod_cascade:
            lst_mod_keys = filter(lambda x: "missing_modifications_" in x, exp_sec_dct.keys())
            lst_mods = [exp_sec_dct[f].split(", ") if exp_sec_dct[f] != u"None" else None for f in sorted(lst_mod_keys)]
            exp_dict["ordered_list_of_missing_modifications"] = lst_mods

        list_of_experiments.append(exp_dict)

    return list_of_experiments, xifdr_settings_dict, rm_links_explainable_by_prev_db, unmod_mod_cascade


def parse_cmd_line():
    # print help message if script is called without argument
    if len(sys.argv) != 2:
        print \
            """
            Script has to be called with config file as argument.
            The directory of the config file will be the output dir.
            """
        sys.exit(1)
    return sys.argv[1]


def main():
    rel_config_file = parse_cmd_line()

    main_execution(rel_config_file)


def main_execution(rel_config_file):
    config_file = os.path.abspath(rel_config_file)

    list_of_experiment_dicts, xifdr_settings_dict, rm_links_explainable_by_prev_db, unmod_mod_cascade = \
        parse_config(config_file)

    # set output dir
    out_dir = os.path.split(config_file)[0]
    # set up logging
    log_file = os.path.join(out_dir, __name__ + '.log')
    logging.basicConfig(filename=log_file, level=logging.DEBUG,
                        format='%(asctime)s - %(levelname)s - %(name)s - %(message)s')
    logger = logging.getLogger()

    list_of_experiments = []
    for experiment_dict in list_of_experiment_dicts:
        list_of_experiments.append(
            build_experiment_from_dict(experiment_dict, rm_links_explainable_by_prev_db, unmod_mod_cascade)
        )
    for exp in list_of_experiments:
        exp.check_files_exist()

    starttime = time.time()
    exp_iterator(
        list_of_experiments=list_of_experiments,
        xifdr_settings_dict=xifdr_settings_dict,
        out_dir=out_dir
    )
    logging.info("Execution of entire script took {}"
                 .format(calculate_elapsed_time(starttime)))
    logging.shutdown()


if __name__ == "__main__":
    main()

    # testing
    # rel_config_file = r"/home/henning/mnt/xitu/Data/Results/170626_cascade_search/Proteasome/20180328-modified_unmodified/180329-3step_cascade-umod_Mox_bs3/myconfig.ini"
    # main_execution(rel_config_file)
