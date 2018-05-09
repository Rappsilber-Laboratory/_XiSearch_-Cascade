"""
Evaluate results of ribosome fraction 14.
To do so, use pipeline for xifdr analysis.
"""

from lib.XiFdrWrapper import XiFdrWrapper
import os
import logging
import sys
import re


def get_list_of_files(location, file_regex=r"xi_results.csv"):
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


def convert_list_to_dict_of_files(list_of_files):
    """groups replicates as lists in dict value for experiment key"""
    exp_name_last_run = ""
    dict_of_files = {}
    for f in list_of_files:
        # read out folder of experiment
        result_path = f
        for i in range(1):
            result_path = os.path.split(result_path)[0]
        exp_name_this_run = os.path.split(result_path)[1]
        if exp_name_this_run == exp_name_last_run:
            dict_of_files[exp_name_this_run].append(f)
        else:
            dict_of_files[exp_name_this_run] = [f]
        exp_name_last_run = exp_name_this_run
    return dict_of_files


# get list of files for input
def generate_dict_of_files(location=r".",
                           file_regex=r"xi_results.csv"):
    list_of_files = get_list_of_files(location, file_regex)
    file_dict = convert_list_to_dict_of_files(list_of_files)
    return file_dict

# generate output dir
# ouput_dir=../../results/Fr14_xifdr_reevaluation
# subdir for each quantile
# subdir "xifdr_output


def generate_input_output_dict(output_base, input_base, lst_input_files):
    input_output_dict = {}
    for input_file in lst_input_files:
        input_dir = os.path.split(input_file)[0]
        # input dir relative to input base dir
        rel_input_dir = os.path.relpath(input_dir, input_base)
        rel_input_dir_wo_lowest_dir, lowest_dir = os.path.split(rel_input_dir)
        if lowest_dir == "xi_output":
            rel_input_dir = rel_input_dir_wo_lowest_dir
        output_dir = os.path.join(output_base, rel_input_dir, "xifdr_ouput")
        input_output_dict[input_file] = output_dir

    return input_output_dict


def execute_xifdr(input_output_dict, xifdr_settings_dict):
    for input_file, out_dir in input_output_dict.items():
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        XiFdrWrapper.xifdr_execution(
            xifdr_input_csv=input_file,
            xifdr_output_dir=out_dir,
            pepfdr=xifdr_settings_dict['pepfdr'],
            additional_xifdr_arguments=xifdr_settings_dict['additional_xifdr_arguments'],
            reportfactor=xifdr_settings_dict['reportfactor'])


if __name__ == "__main__":

    # print help message if script is called without argument
    if len(sys.argv) != 2:
        print \
"""
Script has to be called with config file as argument.
The directory of the config file will be the output dir.
"""
        sys.exit(1)

    rel_config_file = sys.argv[1]

    # # testing
    # print "testing func active"
    # rel_config_file = r"../../Data/Results/170823_Chaetomium/171002/xifdr_reevaluation/ppifdr_10/xifdr_reevaluation_config.py"

    config_file = os.path.abspath(rel_config_file)
    execfile(config_file)
    # set output dir
    out_dir = os.path.split(config_file)[0]

    xifdr_settings_dict = xifdr_settings
    ressource_basedir = os.path.abspath(ressource_basedir)

    logging.basicConfig(filename=os.path.join(out_dir, 'log.txt'), level=logging.DEBUG,
                        format='%(asctime)s - %(levelname)s - %(name)s - %(message)s')
    lst_input_files = get_list_of_files(
        location=ressource_basedir,
        file_regex=r"xi_results.csv"
    )
    input_output_dict = generate_input_output_dict(
        output_base=out_dir,
        input_base=ressource_basedir,
        lst_input_files=lst_input_files
    )

    execute_xifdr(
        input_output_dict=input_output_dict,
        xifdr_settings_dict=xifdr_settings_dict
    )
