import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
import re

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class IbaqExtraction:
    """takes a MaxQuant results file and
    returns a sorted list of tuples,
    each tuple containing the protein identifier and the respective iBAQ value of the protein"""
    def __init__(self, filename):
        self.filename = filename
        # pandas data frame that holds protein IDs and their iBAQ intensities, absolute and normalized
        self.results = self.read_file()

    def read_file(self):
        # read MaxQuant result file with protein ID as index
        raw_file = pd.read_table(self.filename, index_col=0)
        raw_file[['Potential contaminant', "Reverse"]] = raw_file[['Potential contaminant', "Reverse"]].astype(str)
        raw_file_wo_contaminants = raw_file[raw_file["Potential contaminant"] != "+"]
        raw_file_wo_reverse = raw_file_wo_contaminants[raw_file_wo_contaminants["Reverse"] != "+"]
        unsorted_results = pd.DataFrame(raw_file_wo_reverse["iBAQ"])
        sorted_results = unsorted_results.sort_values(by='iBAQ', ascending=False)
        # add column with logarithmic iBAQ values, natural logarithm
        sorted_results['logIBAQ'] = np.log(sorted_results['iBAQ'])
        # add columns with normalized iBAQ values
        sorted_results['normIBAQ'] = sorted_results['iBAQ'] / sorted_results['iBAQ'].max()
        sorted_results['normlogIBAQ'] = sorted_results['logIBAQ'] / sorted_results['logIBAQ'].max()
        return sorted_results

    def split_up_protein_groups(self, list_to_clean):
        cleaned_list_of_proteins = []
        for element in list_to_clean:
            cleaned_list_of_proteins.extend(element.split(';'))
        return cleaned_list_of_proteins

    def get_rel_log_higher_than(self, fraction):
        """returns a list of proteins of at least 'percentage' intensity of the maximum iBAQ value"""
        if fraction == 0:   # to include '-inf' log values
            list_of_proteins = self.results.index.tolist()
        else:
            top_frac = self.results[self.results['normlogIBAQ'] >= fraction]
            list_of_proteins = top_frac.index.tolist()
        # split single elements containing multiple proteins into multiple elements
        cleaned_list_of_proteins = self.split_up_protein_groups(list_of_proteins)
        return cleaned_list_of_proteins

    def get_top_quant(self, quantile):
        """returns a list of proteins that have an intensity >= the specified quantile"""
        top_quant = self.results[self.results['iBAQ'] >= self.results.quantile(quantile)['iBAQ']]
        top_quant = top_quant.index.tolist()
        top_quant_list = self.split_up_protein_groups(top_quant)
        return top_quant_list

    def visualization(self):
        """histogram comparison of the different preprocessed iBAQ intensities"""
        bins = 100
        vals = self.results['normIBAQ'].dropna()
        logvals = self.results[self.results['normlogIBAQ'] >= 0]['normlogIBAQ']

        plt.subplot(211)    # the first subplot in the first figure (figure consists of two rows)
        # x-axis ends at the end of the values
        # plt.xlim(0, 1)
        plt.hist(vals, bins=bins)
        # plt.xlabel('normalized intensity')
        plt.ylabel('Frequency')
        plt.title('normalized intensity')
        plt.grid(True)

        # get 2 rows, 1 column and fignum
        plt.subplot(212)    # the second subplot in the first figure (figure consists of two rows)
        # x-axis ends at the end of the values
        # plt.xlim(0, 1)
        # the histogram of the data
        plt.hist(logvals, bins=bins)
        plt.xlabel('normalized log intensity')
        plt.ylabel('Frequency')
        # plt.title('normalized logarithmic values')
        plt.grid(True)

        plt.show()


# # test cases
# ibaq_list = IbaqExtraction("../Data/Fr14/Fr14_proteinGroups.txt")
# print ibaq_list.results
# protein_list = ibaq_list.get_intensity_top_log_frac(0.9)
# print ibaq_list.get_intensity_top_log_frac(0)
# print ibaq_list.normalized_results
# dfgui.show(ibaq_list.results)
# ibaq_list.visualization()
# print ibaq_list.get_top_quant(0.5)


class FastaHandler:
    def __init__(self, fasta_filename):
        self.filename = fasta_filename
        self.dict = self.read_fasta()

    def read_fasta(self):
        """read self.filename fasta file into a dict with protein id as key and its sequencce as value"""
        mydict = {}
        dict_key = ""
        list_of_non_unique_ids = []
        protein_id_regex = re.compile(r'^>.*\|(.*)\|.*')
        protein_seq_regex = re.compile(r'^[A-Za-z]\B')
        with open(self.filename) as db:
            for line in db.readlines():
                protein_id_regex_hit = protein_id_regex.match(line)
                protein_seq_regex_hit = protein_seq_regex.match(line)
                if protein_id_regex_hit:
                    dict_key = protein_id_regex_hit.group(1)
                    if dict_key not in mydict.keys():
                        mydict[dict_key] = ''
                        protein_not_in_mydict = True
                    else:
                        list_of_non_unique_ids.append(dict_key)
                        protein_not_in_mydict = False
                elif protein_seq_regex_hit and protein_not_in_mydict:
                    mydict[dict_key] += line.strip()
        if not list_of_non_unique_ids:
            logging.debug("no duplicates in fasta file '{}'".format(self.filename))
        else:
            msg = "Duplicates in fasta file '{}': \n{}".format(self.filename, list_of_non_unique_ids)
            print msg
            logging.warning(msg)
        return mydict
        # TODO
        # put duplicate statement in log if no duplicates were found

    def build_fasta(self, protein_id_list, filename):
        with open(filename, 'w+') as f:
            for protein_id in protein_id_list:
                if protein_id in self.dict.keys():
                    f.write('>' + protein_id + '\n')
                    f.write(self.dict[protein_id] + '\n')
                else:
                    msg = "protein ID '{}' not in fasta file '{}'".format(protein_id, self.filename)
                    print msg
                    logging.warning(msg)
