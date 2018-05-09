import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib # for setting of figsize
import numpy as np
from lib.iBAQ_FASTA_handler import IbaqExtraction
from lib.list_files import get_list_of_files
from decimal import Decimal     # for scientific notation numbers


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


class ReadResults:
    def __init__(self):
        pass

    @staticmethod
    def fun_get_xifdr_results(file_path, level="Peptide Pairs"):
        """
        INPUT: filepath of xifdr result
        (i.e. "FDR_1.000_0.050_1.000_1.000_1.000_10000.000_false_summary_xiFDR1.0.14.csv")
        OUTPUT: tuple of internal TT and between TT for fdr on peptide pair level
        """
        lst_levels = ["PSMs", "Peptide Pairs", "Link", "Protein Group Pairs"]
        assert level in lst_levels, "Level has to be in '{}'".format(lst_levels)
        level_full = "fdr " + level
        col_int = 2
        col_bet = 5
        df = pd.read_csv(file_path, names=range(10))
        ser_bool = df.iloc[:, 0] == level_full
        assert sum(ser_bool) == 1, "level '{}' was not found/ found more that once in first column.".format(level_full)
        row = df[ser_bool]
        no_int = row.iloc[0, col_int]
        no_bet = row.iloc[0, col_bet]
        return no_int, no_bet

    @staticmethod
    def read_xifdr_values_out_of_files(dict_of_files, level="Peptide Pairs"):
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
        #         quantile_column = "quantile"
        internal_column = "internal TT"
        between_column = "between TT"
        column_names.extend([
            #             quantile_column,
            internal_column,
            between_column])
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
            #             value_data_frame.loc[exper_key, quantile_column] = quantile
            # call of function to get internal TT and between TT out of file
            internal_tt, between_tt = ReadResults.fun_get_xifdr_results(file_path=f, level=level)
            value_data_frame.loc[exper_key, internal_column] = internal_tt
            value_data_frame.loc[exper_key, between_column] = between_tt
            # print i
        # conversion of columns to numeric
        value_data_frame = value_data_frame.apply(pd.to_numeric)
        value_data_frame.sort_index(inplace=True)
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
        exp = main.Experiment(
            list_of_peak_files=[],
            maxquant_result_file=maxquant_result_file,
            initial_quantiles=[],
            exp_name="generic exp")
        max_iBAQ_value = max(exp.ibaq_object.results['iBAQ'])
        for quantile in list_of_quantiles:
            relative_iBAQ_dict[quantile] = exp.ibaq_object.results.quantile(quantile)['iBAQ'] / max_iBAQ_value
        values_data_frame['relative iBAQ'] = pd.Series(relative_iBAQ_dict)
        return values_data_frame


def calculate_values_of_interest(value_data_frame, scaling_method='median'):
    """takes the object generated by 'read_xifdr_values_out_of_files' and calculates mean and sdev for internal and between TT
    for each experiment.
    """
    lst_scaling_methods = ['median', 'mean', 'max']
    assert scaling_method in lst_scaling_methods, 'scaling_method has to be in {}'.format(lst_scaling_methods)
    # for each row calculate sd and mean of internal and between TT
    internal_tt_df = value_data_frame.filter(regex=re.compile("internal tt", re.I))
    between_tt_df = value_data_frame.filter(regex=re.compile("between tt", re.I))
    value_data_frame['mean internal TT'] = internal_tt_df.mean(1)
    value_data_frame['std internal TT'] = internal_tt_df.std(1)
    value_data_frame['mean between TT'] = between_tt_df.mean(1)
    value_data_frame['std between TT'] = between_tt_df.std(1)

    # scaling
    if scaling_method == 'median':
        scaling_factor = value_data_frame['internal TT'].median()
    elif scaling_method == 'mean':
        scaling_factor = value_data_frame['internal TT'].mean()
    elif scaling_method == 'max':
        scaling_factor = value_data_frame['internal TT'].max()

    value_data_frame['scaled internal TT'] = value_data_frame['internal TT'] / scaling_factor

    return value_data_frame


def extract_searchtime_from_log(
        logfile,
        time_regex=r"Search execution took (.*) for cmd.*(?:'rappsilber\.applications\.Xi')",
        identifier_regex=r"\/Fr\d{2}\/([^\/]*)\/xi_output\/"
):
    row = -1
    #     dct_keys = ["id", "time", "line"]
    dct_keys = ["time", "line"]
    dict_single_hit = {}
    df_hits = pd.DataFrame(columns=dct_keys)

    re_time = re.compile(time_regex)
    re_id = re.compile(identifier_regex)
    with open(logfile) as db:
        for line in db.readlines():
            re_hit_time = re_time.search(line)
            if re_hit_time:
                #                 print "time hit"
                re_hit_id = re_id.search(line)
                if re_hit_id:
                    row += 1
                    #                     df_hits.loc[row, dct_keys[0]] = re_hit_id.group(1)
                    #                     df_hits.loc[row, dct_keys[1]] = re_hit_time.group(1)
                    #                     df_hits.loc[row, dct_keys[2]] = line
                    df_hits.loc[re_hit_id.group(1)] = {dct_keys[0]: re_hit_time.group(1),
                                                       dct_keys[1]: line}

    #                     print "id hit"
    #                     dict_single_hit = {
    #                         dct_keys[0]: re_hit_id.group(1),
    #                         dct_keys[1]: re_hit_time.group(1),
    #                         dct_keys[2]: line
    #                     }
    # #                     print dict_single_hit
    #                     df_dict = pd.DataFrame(dict_single_hit, index=[dict_single_hit[dct_keys[0]]])
    #                     df_hits.append(df_dict, ignore_index=True)
    #                     print df_hits
    return df_hits


def convert_strtime_to_hours(series_str_time):
    """
    converts a series of time strings to hours
    INPUT:
    series_str_time: series of strings
    RETURN:
    series: series of floats, hours
    """
    series_hours = pd.Series(index=series_str_time.index, name="hours")
    time_regex = re.compile(r"^(?:([0-9]*)\sdays?,\s)?(?:([0-9]{1,2}):)?([0-9]{1,2}):([0-9\.]{9})$")
    for index in series_str_time.index:
        re_hit = time_regex.match(series_str_time[index])
        # seconds
        if re_hit.group(4) != None:
            minutes = float(re_hit.group(4)) / 60
        else:
            minutes = 0

        # minutes
        if re_hit.group(3) != None:
            minutes += int(re_hit.group(3))

        # hours
        if re_hit.group(2) != None:
            hours = int(re_hit.group(2)) + minutes / 60.
        else:
            hours = minutes / 60.

        # days
        if re_hit.group(1) != None:
            hours += int(re_hit.group(1)) * 24

        #         timedict = {
        #             "days": int(re_hit.group(1)),
        #             "hours": int(re_hit.group(2)),
        #             "minutes": int(re_hit.group(3)),
        #             "seconds": float(re_hit.group(4))
        #         }
        #         minutes = timedict["minutes"] + timedict["seconds"] / 60.
        #         hours = timedict["hours"] + minutes / 60.
        #         hours += timedict["days"] * 24
        #         df.loc[index, "search_time [hours]"] = hours
        series_hours[index] = hours
    return series_hours


class Plotting:
    def __init__(self):
        pass

    @staticmethod
    def values_for_plotting(df_calculated_values, index="quantile", y1='internal TT', y2='between TT'):
        ind = df_calculated_values[index]
        internal_tt = df_calculated_values[y1]
        between_tt = df_calculated_values[y2]
        return ind, internal_tt, between_tt

    @staticmethod
    def calc_minimal_distance(lst):
        ind = lst
        # calculation of distance between bars
        minimal_distance = ind[1] - ind[0]
        for i in range(1, len(ind)):
            if ind[i] - ind[i - 1] < minimal_distance:
                minimal_distance = ind[i] - ind[i - 1]
        return minimal_distance

    @staticmethod
    def plot_stacked_barplot(df_calculated_values, plot_title, **kwargs):
        """
        y-axis: number of internal and between TT crosslinks stacked on each other
        x-axis: sample ID
        :return:
        """
        ind, internal_tt, between_tt = Plotting.values_for_plotting(df_calculated_values)

        # calculation of distance between bars
        minimal_distance = ind[1] - ind[0]
        for i in range(len(ind) - 1):
            if ind[i + 1] - ind[i] < minimal_distance:
                minimal_distance = ind[i + 1] - ind[i]
        width = 0.8 * minimal_distance  # the width of the bars: can also be len(x) sequence
        #
        errorbar_capsize = 2
        # set up new figure
        plt.figure()
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
        if 'xlim' in kwargs:
            plt.xlim(kwargs['xlim'])

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
        minimal_distance = Plotting.calc_minimal_distance(ind)
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

    @staticmethod
    def multilineplot(list_of_dfs, plot_title, **kwargs):
        # setup figure
        plt.figure()
        # read data
        for df in list_of_dfs:
            ind, internal_tt, between_tt = \
                Plotting.values_for_plotting(
                    df, index="relative iBAQ",
                    y1="scaled internal TT")
            plt.plot(ind, internal_tt, label=df.name)
            plt.xscale('log')
        plt.title(plot_title)
        plt.legend()
        plt.show()

    @staticmethod
    def cascade_barplot(df, plot_title, **kwargs):
        """
        dataframe has to have column "sum internal" which gives the cumulative link number if all the DBs
        up to the "row" one are searched.
        """
        # sort df in descending order
        #         df.sort_index(ascending=False, inplace=True)
        if "xvals" in kwargs:
            ind = kwargs["xvals"]
        else:
            ind = df["quantile"]

        # bar width
        minimal_distance = Plotting.calc_minimal_distance(ind)
        if "width" in kwargs:
            width = kwargs["width"]
        else:
            width = 0.8 * minimal_distance
        y = df["sum internal"]

        plt.bar(ind, y, width, color='#d62728')

        plt.title(plot_title)
        plt.xticks(ind, df.index, rotation=45)

        if 'xlim' in kwargs:
            plt.xlim(kwargs['xlim'])


def plot_multi_cascade_plot(dct_of_dfs, title="{}", xvals_column=""):
    xlim = []

    fig = plt.figure(figsize=(9., 12.))
    for i, dct_key in enumerate(sorted(dct_of_dfs)):
        df = dct_of_dfs[dct_key]

        if not xvals_column:
            xvals = range(df.shape[0])
        else:
            xvals = df[xvals_column]

        # fig.SubplotParams(hspace=0.4)
        plt.subplot(3, 1, i + 1)
        # plt.subplots_adjust(hspace=0.4)
        Plotting.cascade_barplot(df=df, plot_title=title.format(df.name), xvals=xvals)
        plt.tick_params(
            axis='x',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            bottom='off',  # ticks along the bottom edge are off
            top='off',  # ticks along the top edge are off
            labelbottom='off')  # labels along the bottom edge are off
        if not xlim:
            xlim = plt.xlim()
        if xlim:
            plt.xlim(xlim)

    # put ticks from subplot with maximum number of datapoint below figure
    n_rows_max = 0
    for key in dct_of_dfs:
        n_rows_df = dct_of_dfs[key].shape[0]
        if n_rows_df > n_rows_max:
            n_rows_max = n_rows_df
            max_key = key
    ind = range(n_rows_max)
    label = dct_of_dfs[key].index
    plt.xticks(ind, label, rotation=45)
    plt.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom='on',  # ticks along the bottom edge
        top='off',  # ticks along the top edge
        labelbottom='on')  # labels along the bottom edge


def plot_cascade_versus_seq(dct_df_cas, dct_df_seq):
    xlim = []
    n_plots = len(dct_df_cas.keys())
    fig = plt.figure()
    axes = []
    for i, sample in enumerate(sorted(dct_df_cas)):
        df_cas = dct_df_cas[sample]
        df_seq = dct_df_seq[sample]

        xvals = np.arange(df_cas.shape[0])

        # fig.SubplotParams(hspace=0.4)
        # ax = fig.add_subplot(n_plots, 1, i + 1)
        ax = fig.add_subplot(1, n_plots, i + 1)
        # plt.subplots_adjust(hspace=0.4)

        width = 0.5 * 0.8 * (xvals[1] - xvals[0])
        xvals_seq = xvals - width / 2
        yvals_seq = df_seq["internal TT"]
        # raise NotImplementedError
        ax.bar(xvals_seq, yvals_seq, width, label="separate searches")

        xvals_cas = xvals + width / 2
        yvals_cas = df_cas["sum internal"]
        ax.bar(xvals_cas, yvals_cas, width, label="cumulative searches")
        ax.set_title(sample)

        # Plotting.cascade_barplot(df=df, plot_title=title.format(df.name), xvals=xvals)
        # ax.tick_params(
        #     axis='x',  # changes apply to the x-axis
        #     which='both',  # both major and minor ticks are affected
        #     bottom='off',  # ticks along the bottom edge are off
        #     top='off',  # ticks along the top edge are off
        #     labelbottom='off'  # labels along the bottom edge are off
        # )
        if not xlim:
            xlim = plt.xlim()
        if xlim:
            plt.xlim(xlim)
        axes.append(ax)

    # put ticks from subplot with maximum number of datapoint below figure
    n_rows_max = 0
    max_key = 0
    for key in dct_df_cas:
        n_rows_df = dct_df_cas[key].shape[0]
        if n_rows_df > n_rows_max:
            n_rows_max = n_rows_df
            max_key = key
    ind = range(n_rows_max)
    labels = dct_df_cas[max_key].index.values
    for ax in axes:
        plt.sca(ax)
        plt.xticks(ind, labels, rotation=45, ha="right")
    axes[0].legend(loc="upper left")
    # axes[-1].tick_params(
    #     axis='x',  # changes apply to the x-axis
    #     which='both',  # both major and minor ticks are affected
    #     bottom='on',  # ticks along the bottom edge
    #     top='off',  # ticks along the top edge
    #     labelbottom='on')  # labels along the bottom edge
    return fig, axes


def build_index_series(dataframe_index, percentage_1_index=None, top_no_200_index=None):
    ind = range(len(dataframe_index))
    ser_ind = pd.Series(ind, index=dataframe_index)
    #     df['index'] = ind
    if percentage_1_index is None or top_no_200_index is None:
        percentage_1_index = ser_ind['percentage_1'] + 0.5
        top_no_200_index = ser_ind['top_no_200'] + 0.5
    ser_ind.loc['percentage_1'] = percentage_1_index
    ser_ind.loc['top_no_200'] = top_no_200_index
    return ser_ind


def calculate_sums(df):
    """
    df has to be in ascending order of DB size
    """
    # calculate sums
    lst_summed_vals = [df["internal TT"][0]]
    for row in df["internal TT"][1:]:
        lst_summed_vals.append(lst_summed_vals[-1] + row)
    df["sum internal"] = lst_summed_vals
    return df


def build_df(
    file_location,
    file_regex=r".*/FDR_.*_false_summary_xiFDR(\d+\.)*csv",
    maxquant_result_file=None,
    level="Peptide Pairs",
    **kwargs
):
    # get list of files
    file_list = get_list_of_files(file_location, file_regex)
    # convert list to dict
    file_dict = convert_list_to_dict_of_files(file_list)
    # read out values from files
    values_data_frame = ReadResults.read_xifdr_values_out_of_files(file_dict, level)
    if maxquant_result_file:
        values_data_frame = ReadResults.read_relative_iBAQ_values(
            values_data_frame=values_data_frame,
            maxquant_result_file=maxquant_result_file
        )
#     print values_data_frame
    # calculate some additional values for each experiment
    df_with_calculated_values = calculate_values_of_interest(values_data_frame)
#     print df_with_calculated_values
#     Plotting.plot_stacked_barplot(df_with_calculated_values, plot_title, **kwargs)
    # Plotting.plot_stacked_barplot_relIBAQ(df_with_calculated_values, plot_title)
    # lineplot(df_with_calculated_values, plot_title)
#     Plotting.multilineplot(list_of_dfs, plot_title, **kwargs)
    return df_with_calculated_values


def single_run_df(identifier, xifdr_basedir, index_order=None, level="Peptide Pairs"):
    df = build_df(
        file_location=xifdr_basedir.format(identifier),
        file_regex=r".*/FDR_.*_false_summary_xiFDR((\d+\.)+|null)\.?csv",
        level=level# ,
        #             maxquant_result_file=maxquant_file.format(identifier)
        )
    df.sort_index(ascending=False, inplace=True)
    if index_order is not None:
        df = df.reindex(index_order)
    df["sum internal"] = df["internal TT"].cumsum()
    df.name = identifier
    return df


def extract_fdr_results(lst_identifiers, xifdr_basedir, index_order=None, level="Peptide Pairs"):
    dct_of_dfs = {}
    for i, identifier in enumerate(lst_identifiers):
        df = single_run_df(identifier, xifdr_basedir, index_order, level)
        dct_of_dfs[identifier] = df
    return dct_of_dfs
