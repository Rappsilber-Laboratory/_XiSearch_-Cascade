"""
author: Henning Schiebenhoefer

This Script takes the internal TT and in between TT from each xifdr summary file and plots them in a stacked plot.
For values from the same sample set, averages and standard deviation of internal and between are calculated.

"""

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from lib.iBAQ_FASTA_handler import IbaqExtraction
# import main
from decimal import Decimal     # for scientific notation of numbers

try:
    import seaborn  # styling of plots
except ImportError:
    print "seaborn not installed!"
else:
    seaborn.set()


def get_list_of_files(location, file_regex=r"FDR_.*_false_summary_xiFDR(\d+\.)*csv"):
    """generates a list of files, that satisfy specific conditions, such as filename and location
    INPUT: constraints
    RETURNS a list with all the experiment files as values"""
    list_of_files = []
    regex = re.compile(file_regex)
    for rel_dir, sub_dirs, files in os.walk(location):
        for f in files:
            if regex.match(f):
                list_of_files.append(os.path.join(rel_dir, f))
    return list_of_files

# # Test cases
# files = get_list_of_files("../../170317_OCCM_random_fasta_analysis/", r"FDR_.*_false_summary_xiFDR(\d+\.)*csv")
# for f in files:
#     print f
# print os.path.commonprefix(files)
# raw_path = r"../../occm-proteins_14_e-coli-proteins_10/xifdr_output/3/FDR_1.000000_0.000500_1.000000_1.000000_1.000000_10000.000000_false_summary_xiFDR1.0.12.csv"
# result_path = raw_path
# for i in range(3):
#     result_path = os.path.split(result_path)[0]
# print result_path


def convert_list_to_dict_of_files(list_of_files):
    """groups replicates as lists in dict value for experiment key"""
    exp_name_last_run = ""
    dict_of_files = {}
    for f in list_of_files:
        # read out folder of experiment
        result_path = f
        for i in range(2):
            result_path = os.path.split(result_path)[0]
        exp_name_this_run = os.path.split(result_path)[1]
        if exp_name_this_run == exp_name_last_run:
            dict_of_files[exp_name_this_run].append(f)
        else:
            dict_of_files[exp_name_this_run] = [f]
        exp_name_last_run = exp_name_this_run
    return dict_of_files


# # test cases
# file_list = get_list_of_files("../../170317_OCCM_random_fasta_analysis/", r"FDR_.*_false_summary_xiFDR(\d+\.)*csv")
# file_dict = convert_list_to_dict_of_files(file_list)
# print file_dict
# for key in file_dict:
#     print "key: %s , value: %s" % (key, file_dict[key])


class ReadResults:
    def __init__(self):
        pass

    @staticmethod
    def fun_get_xifdr_results(file_path):
        """
        INPUT: filepath of xifdr result 
        (i.e. "FDR_1.000_0.050_1.000_1.000_1.000_10000.000_false_summary_xiFDR1.0.14.csv") 
        OUTPUT: tuple of internal TT and between TT for fdr on peptide pair level
        """
        result_line = 18
        internal_tt_column = 2
        between_tt_column = 5
        csv_file = pd.read_csv(file_path, names=range(10))
        return csv_file.iloc[result_line, internal_tt_column], csv_file.iloc[result_line, between_tt_column]

    @staticmethod
    def read_xifdr_values_out_of_files(dict_of_files):
        """reads 'quantile', 'internal TT' and 'between TT' of each file in dict_of_files into a pandas dataframe
        INPUT: dict of files, specifier for values to read
        OUTPUT: a pandas datafram that contains: experiment, run number, internal TT and between TT"""
        # find out number of runs
        number_of_runs = 0
        for key in dict_of_files:
            if len(dict_of_files[key]) > number_of_runs:
                number_of_runs = len(dict_of_files[key])
        # construct name of columns
        column_names = []
        # regex_exp_specifier = re.compile(r"(.+)_(\d+)_(.+)_(\d+)")
        # exp_specifier = regex_exp_specifier.match(dict_of_files.keys()[0])
        # column_names.extend(exp_specifier.group(1, 3))
        quantile_column = "quantile"
        internal_column = "internal TT"
        between_column = "between TT"
        column_names.extend([quantile_column, internal_column, between_column])
        # construction of data frame to calculate on
        value_data_frame = pd.DataFrame(index=dict_of_files, columns=column_names)
        # value_data_frame.sort_index(inplace=True)
        # print value_data_frame
        for exper_key in dict_of_files:
            f = dict_of_files[exper_key][0]
            # catching of quantile
            quantile_string = f
            quantile = np.NAN
            for i in range(3):
                quantile_string, quantile = os.path.split(quantile_string)
            value_data_frame.loc[exper_key, quantile_column] = quantile
            # call of function to get internal TT and between TT out of file
            internal_tt, between_tt = ReadResults.fun_get_xifdr_results(file_path=f)
            value_data_frame.loc[exper_key, internal_column] = internal_tt
            value_data_frame.loc[exper_key, between_column] = between_tt
            # print i
        # conversion of columns to numeric
        value_data_frame = value_data_frame.apply(pd.to_numeric)
        value_data_frame.sort_values(quantile_column, inplace=True)
        # for i in value_data_frame.columns:
        #     print value_data_frame[i].dtype
        return value_data_frame

    # # test cases
    # file_list = get_list_of_files("../../170317_OCCM_random_fasta_analysis/")
    # file_dict = convert_list_to_dict_of_files(file_list)
    # # print file_dict.keys()[0]
    # values_data_frame = read_xifdr_values_out_of_files(file_dict)

    @staticmethod
    def read_relative_iBAQ_values(values_data_frame, maxquant_result_file):
        """returns values dataframe with new column 'relative iBAQ' containing relative iBAQ value associated 
        with each quantile"""
        relative_iBAQ_dict = {}
        list_of_quantiles = values_data_frame['quantile'].tolist()
        ibaq_object = IbaqExtraction(filename=maxquant_result_file)
        max_iBAQ_value = max(ibaq_object.results['iBAQ'])
        for quantile in list_of_quantiles:
            relative_iBAQ_dict[quantile] = ibaq_object.results.quantile(quantile)['iBAQ'] / max_iBAQ_value
        values_data_frame['relative iBAQ'] = pd.Series(relative_iBAQ_dict)
        return values_data_frame


def calculate_values_of_interest(value_data_frame):
    """takes the object generated by 'read_xifdr_values_out_of_files' and calculates mean and sdev for internal and between TT
    for each experiment.
    """
    # for each row calculate sd and mean of internal and between TT
    internal_tt_df = value_data_frame.filter(regex=re.compile("internal tt", re.I))
    between_tt_df = value_data_frame.filter(regex=re.compile("between tt", re.I))
    value_data_frame['mean internal TT'] = internal_tt_df.mean(1)
    value_data_frame['std internal TT'] = internal_tt_df.std(1)
    value_data_frame['mean between TT'] = between_tt_df.mean(1)
    value_data_frame['std between TT'] = between_tt_df.std(1)
    return value_data_frame


class Plotting:
    def __init__(self):
        pass

    @staticmethod
    def values_for_plotting(df_calculated_values):
        ind = df_calculated_values["quantile"]
        internal_tt = df_calculated_values['internal TT']
        between_tt = df_calculated_values['between TT']
        return ind, internal_tt, between_tt

    @staticmethod
    def plot_stacked_barplot(df_calculated_values, plot_title):
        """
        y-axis: number of internal and between TT crosslinks stacked on each other
        x-axis: sample ID
        :return:
        """
        ind, internal_tt, between_tt = Plotting.values_for_plotting(df_calculated_values)

        # calculation of distance between bars
        minimal_distance = ind[1] - ind[0]
        for i in range(len(ind)-1):
            if ind[i+1] - ind[i] < minimal_distance:
                minimal_distance = ind[i+1] - ind[i]
        width = 0.8 * minimal_distance  # the width of the bars: can also be len(x) sequence
        #
        errorbar_capsize = 2
        # p1 = plt.bar(ind, menMeans, width, color='#d62728', yerr=menStd)
        p1 = plt.bar(ind, internal_tt, width, color='#d62728')
        # p2 = plt.bar(ind, womenMeans, width,
        #              bottom=menMeans, yerr=womenStd)
        p2 = plt.bar(ind, between_tt, width, bottom=internal_tt)
        #
        plt.xlabel('fasta DB made out of proteins with iBAQ >= depicted quantile')
        # plt.ylabel('Scores')
        plt.ylabel('no. of peptide links')
        # plt.title('Scores by group and gender')
        plt.title(plot_title)
        # plt.xticks(ind, ('G1', 'G2', 'G3', 'G4', 'G5'))
        xticks = ind

        # print df_calculated_values.columns
        # for index, row in df_calculated_values.iterrows():
        #     # xticks.append("{} OCCM, {} E.coli".format(int(row['occm-proteins']), int(row['e-coli-proteins'])))
        #     xticks.append("{}\n{}".format(int(row['occm-proteins']), int(row['e-coli-proteins'])))
        plt.xticks(ind, xticks, rotation=0)
        # plt.yticks(np.arange(0, 81, 10))
        # plt.legend((p1[0], p2[0]), ('Men', 'Women'))
        plt.legend((p1[0], p2[0]), ('internal TT', 'between TT'))
        plt.subplots_adjust(bottom=0.2)
        #
        # plt.show()
        plt.show()

    @staticmethod
    def plot_stacked_barplot_relIBAQ(df_calculated_values, plot_title):
        """
        y-axis: number of internal and between TT crosslinks stacked on each other
        x-axis: log-scale, relative iBAQ value
        :return:
        """
        _, internal_tt, between_tt = Plotting.values_for_plotting(df_calculated_values)

        # indices on logarithmic scale
        ind = np.log(df_calculated_values['relative iBAQ'])
        # print min(ind)
        if min(ind) == -np.inf:
            ind['0.0'] = min(ind[ind != -np.inf]) - 0.5
        xticks = ['%.2E' % Decimal(x) for x in df_calculated_values['relative iBAQ']]

        # set bar width: calculation of distance between bars
        minimal_distance = ind[1] - ind[0]
        for i in range(len(ind) - 1):
            if ind[i + 1] - ind[i] < minimal_distance:
                minimal_distance = ind[i + 1] - ind[i]
        width = 0.8 * minimal_distance  # the width of the bars: can also be len(x) sequence

        #
        errorbar_capsize = 2
        # p1 = plt.bar(ind, menMeans, width, color='#d62728', yerr=menStd)
        p1 = plt.bar(ind, internal_tt,
                     width
                     , color='#d62728'
                     )
        # p2 = plt.bar(ind, womenMeans, width,
        #              bottom=menMeans, yerr=womenStd)
        p2 = plt.bar(ind, between_tt,
                     width,
                     bottom=internal_tt)
        #
        plt.xlabel('relative iBAQ value of lowest abundant Protein in fasta DB\n (log scale)')
        # plt.ylabel('Scores')
        plt.ylabel('no. of peptide links')
        # plt.title('Scores by group and gender')
        plt.title(plot_title)

        # print df_calculated_values.columns
        # for index, row in df_calculated_values.iterrows():
        #     # xticks.append("{} OCCM, {} E.coli".format(int(row['occm-proteins']), int(row['e-coli-proteins'])))
        #     xticks.append("{}\n{}".format(int(row['occm-proteins']), int(row['e-coli-proteins'])))

        plt.xticks(ind, xticks, rotation=45)

        # plt.yticks(np.arange(0, 81, 10))
        # plt.legend((p1[0], p2[0]), ('Men', 'Women'))
        plt.legend((p1[0], p2[0]), ('internal TT', 'between TT'))
        plt.subplots_adjust(bottom=0.2)
        #
        # plt.show()
        plt.show()

        # todo
        # make clear, that 0 value on x-axis is not in correct distance to the other values (

    @staticmethod
    def lineplot(df_calculated_values, plot_title):
        ind, internal_tt, between_tt = Plotting.values_for_plotting(df_calculated_values)
        p1 = plt.plot(ind, internal_tt)
        p2 = plt.plot(ind, between_tt)
        plt.legend((p1[0], p2[0]), ('internal TT', 'between TT'))
        plt.show()


def runner(file_location,
           plot_title,
           file_regex=r"FDR_.*_false_summary_xiFDR(\d+\.)*csv",
           maxquant_result_file=None):
    # get list of files
    file_list = get_list_of_files(file_location, file_regex)
    # convert list to dict
    file_dict = convert_list_to_dict_of_files(file_list)
    # read out values from files
    values_data_frame = ReadResults.read_xifdr_values_out_of_files(file_dict)
    if maxquant_result_file:
        values_data_frame = ReadResults.read_relative_iBAQ_values(
            values_data_frame=values_data_frame,
            maxquant_result_file=maxquant_result_file)
    print values_data_frame
    # calculate mean and sd for internal and between TT for each experiment
    df_with_calculated_values = calculate_values_of_interest(values_data_frame)
    # Plotting.plot_stacked_barplot(df_with_calculated_values, plot_title)
    Plotting.plot_stacked_barplot_relIBAQ(df_with_calculated_values, plot_title)
    # lineplot(df_with_calculated_values, plot_title)


# # Test cases
# runner("../../170317_OCCM_random_fasta_analysis/")
# runner("/home/henning/mnt/xitu/xi_and_xifdr_pipeline/results_with_optional_runs_and_larger_decoyDB/")

if __name__ == "__main__":
    runner(
        file_location=r"/home/henning/mnt/xitu/Data/Results/170626_cascade_search/170626_cascade_simulation-ON-iBAQ_with_MAXCANDIDATES/Fr16/",
        plot_title=r"""xiFDR peptide links as a function of fasta DB derived from iBAQ quantification (10% pepfdr)""",
        maxquant_result_file=r'/home/henning/mnt/xitu/Data/Input/170323_iBAQ_based_opt/Fr14/MaxQuant_results_Fr14_proteinGroups.txt'
    )
    # builds pandas data frame for a single experiment
    # Fr14 = main.Experiment(
    #     list_of_peak_files=[],
    #     maxquant_result_file=r'../../Data/Fr14/Fr14_proteinGroups.txt',
    #     initial_quantiles=[],
    #     exp_name="Fr14")
    print "penis"
