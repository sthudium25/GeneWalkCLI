from os import mkdir, chdir, listdir
from os.path import splitext, basename, isdir
from typing import Dict, Any
import pandas as pd
from datetime import date
import re
import seaborn as sns
import matplotlib.pyplot as plt


class GenewalkObj:
    user_dir: str = ""
    analysis_name: str = ""
    go_cat: str = ""
    go_type: str = ""

    def __init__(self, project_name: str = f"Project_{date.today()}"):
        print("Enter the directory containing your Genewalk results: ")
        self.user_dir = input()
        self.analysis_name = project_name
        self.directory_setup()

    def directory_setup(self) -> None:
        chdir(self.user_dir)
        try:
            mkdir(f'../results_{self.analysis_name}')
        except FileExistsError:
            print("Warning: there is already an analysis with this name")

    def load_results(self, significance: float = 0.1):
        """pull in the gene walk results to dictionary"""
        files = [f for f in listdir(self.user_dir) if f.endswith(".csv")]
        fns = [splitext(basename(x))[0] for x in files]
        raw_gw: dict = {}
        for i in range(len(fns)):
            raw_gw[fns[i]] = pd.read_csv(files[i])
        self.filter_results(input_dfs=raw_gw, significance=significance)

    def filter_results(self, input_dfs: dict, significance: float):
        """Filter the raw results of Genewalk for only those gene:term pairs that fall below the supplied
            significance threshold.
            :param input_dfs: A dictionary containing raw Genewalk results dataframes from each provided Target
            :type significance: float """
        out = {}
        for name, df in input_dfs.items():
            out[name] = df[(df["global_padj"] < significance) & (df["gene_padj"] < significance)]
        print("select GO type ('go_id' or 'go_name'")
        go_type = input()
        print("select GO cat ('biological process', 'molecular function', 'cellular component', or 'all'")
        go_cat = input()
        self.make_go_set(set_dict=out, go_type=go_type, go_cat=go_cat)

    def make_go_set(self, set_dict: dict, go_type: str,
                    go_cat: str = "all"):
        """Takes a dictionary of Genewalk results (those provided in the input directory) and
            produces a new dictionary containing the set of all significant GO terms of the specified
            arguments (i.e. GO IDs in the biological process category. Default behavior is to
            produce the set containing all GO categories (CC, BP, MF).
            :type go_cat: str
            :type go_type: str
            :type set_dict: Dict(Set)"""
        self.go_cat = go_cat
        self.go_type = go_type
        self.validate_go_parameters(self.go_type, self.go_cat)
        if go_cat == "all":
            self.make_go_set_all_terms(set_dict, go_type)
        else:
            out_sets: dict = {}
            for name, df in set_dict.items():
                subset = df[df["go_domain"] == go_cat]
                out_sets[name] = list(set(subset[go_type]))
            self.make_membership_matrix(set_dictionary=out_sets)

    def make_go_set_all_terms(self, set_dict: dict, go_type: str):
        out_sets = {}
        for name, df in set_dict.items():
            out_sets[name] = list(set(df[go_type]))
        self.make_membership_matrix(set_dictionary=out_sets)

    def make_membership_matrix(self, set_dictionary: dict):
        """Takes the output of make_go_set and produces a dataframe containing a binary representation of GO term (rows)
            membership in each dataset (columns). Additionally, the last column describes how many of the input datasets
            contain a given term
            :type set_dictionary: dict"""
        all_names = [y for x in set_dictionary.values() for y in x]
        name_affiliation_matrix = pd.DataFrame(index=list(set(all_names)),
                                               columns=[name for name, _ in set_dictionary.items()])
        for index, row in name_affiliation_matrix.iterrows():
            for dataset, go_set in set_dictionary.items():
                name_affiliation_matrix.loc[index, dataset] = 1 if index in go_set else 0
        name_affiliation_matrix["row_sums"] = name_affiliation_matrix.sum(axis=1)
        self.truncate_go_values(self.go_type, self.go_cat)
        name_affiliation_matrix.to_csv(
            path_or_buf=f"../results_{self.analysis_name}/"
            f"{self.analysis_name}_GO_{self.go_type}-{self.go_cat}_membership_matrix.csv")
        self.make_colored_matrix(membership_mat=name_affiliation_matrix)

    def validate_go_parameters(self, go_type, go_cat):
        if go_cat not in ["biological process", "molecular function", "cellular component", "all"]:
            raise ValueError("GO category must be one of: biological process, "
                             "molecular function, cellular component, all")
        if go_type not in ["go_name", "go_id"]:
            raise ValueError("GO type must be either 'go_name' or 'go_id'")

    def truncate_go_values(self, go_type, go_cat):
        self.go_type = re.sub("^go_*", "", go_type)
        if go_cat == "all":
            return
        else:
            print(type(self.go_cat))
            temp = go_cat.split()
            self.go_cat = temp[0][0].upper() + temp[1][0].upper()

    def min_x_shared_processes(self, affiliation_matrix: pd.DataFrame, threshold: float, go_type: str = "name"):
        """creates a pandas series of the GO terms that are shared among a minimum of x sets, where x is the 
        threshold provided"""
        index_list = affiliation_matrix[affiliation_matrix.row_sums >= threshold].index
        return index_list.to_series(name=go_type)

    def get_common_processes_at_threshold(self, signif_gw_results: pd.DataFrame,
                                          affiliation_matrix: pd.DataFrame, threshold: float, go_type: str = "name"):
        x_common_processes = {k: pd.merge(df, GenewalkObj.min_x_shared_processes(affiliation_matrix, threshold),
                                          left_on="go_name", right_on=go_type)
                                        for k, df in signif_gw_results.items()}
        return x_common_processes

    def make_colored_matrix(self, membership_mat):
        mat_sorted = membership_mat.sort_values('row_sums', ascending=False)
        print(any(mat_sorted.isna()))
        mat_sorted = mat_sorted.astype(float)
        fig, ax = plt.subplots(figsize=(10,11))
        heatmap = sns.heatmap(mat_sorted, ax=ax)
        heatmap.axes.yaxis.set_visible(False)
        plt.savefig(f"../results_{self.analysis_name}/"
                    f"{self.analysis_name}_{self.go_cat}-{self.go_type}_heatmap.png")
