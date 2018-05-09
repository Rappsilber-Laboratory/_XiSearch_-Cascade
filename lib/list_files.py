import re
import os


def get_list_of_files(location, file_regex=r".*/FDR_.*_false_summary_xiFDR(\d+\.)*csv"):
    """
    searches location recursively for files matching regex.
    regex has to match full file path.
    INPUT:
    location: base dir from where to search recursively
    file_regex: regex that files have to match

    RETURNS
    sorted list of files
    """
    list_of_files = []
    regex = re.compile(file_regex)
    if not os.path.exists(location):
        raise IOError("Specified location '{loc}' does not exist in filesystem.".format(loc=location))
    for rel_dir, sub_dirs, files in os.walk(location):
        for f in files:
            f_full = os.path.join(rel_dir, f)
            # print f_full
            if regex.match(f_full):
                list_of_files.append(f_full)
    if len(list_of_files) == 0:
        raise Warning("No files were found in '{loc}' that match '{re}'"\
                      .format(loc=location, re=file_regex))
    return sorted(list_of_files)


# # # testing
if __name__ == "__main__":
    list_files = \
        get_list_of_files(location=r"/home/henning/mnt/xitu/Data/Results/170823_Chaetomium/171020-random_decoys/",
                          file_regex=r".*/xi_output/xi_results\.csv")
    print "\n".join(list_files)
