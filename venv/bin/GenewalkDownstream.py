from os import mkdir, chdir, listdir, getcwd
from os.path import splitext, basename, isdir
from typing import Dict, Any
import pandas as pd
from collections import defaultdict
from datetime import date
import logging
import re
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import warnings
import sys


def read_arguments():
    parser = argparse.ArgumentParser(description="Downstream analysis of Genewalk (Churchman lab)", add_help=True)
    required_args = parser.add_argument_group("Required arguements")
    optional_args = parser.add_argument_group("Optional arguements")
    input_arguments = parser.add_argument_group('Optional input file and data handling arguments')
    heatmap_arguments = parser.add_argument_group('Optional scatterplot arguments')
    output_arguments = parser.add_argument_group('Optional outfile arguments')
    
    required_args.add_argument('-i', '--input', help="This is the input directory containing the Genewalk output files "
                                                     "that you would like to use. The directory could contain a set "
                                                     "of .csv files that to be compared.",
                               required=True, action='store')
    required_args.add_argument('-p', '--project', help="The name of the project. This will be used in naming of output "
                                                       "files and figures.", required=True, action='store')
    required_args.add_argument('-o', '--outdir', help="This is the output directory in which results will be stored.",
                               required=True, action='store')
    
    # Input arguments 
    input_arguments.add_argument('-n', '--ontology',
                                 help="Which ontology to use: biological process ('bp'), molecular function ('mf'), "
                                      "cellular component ('cc'), or all ontologies ('all'). Default is all.",
                                 required=False, action='store', default='all', choices=['bp', 'mf', 'cc', 'all'])
    input_arguments.add_argument('-t', '--identifier',
                                 help="Which GO identifier to use. description ('name') or ID ('id')",
                                 required=False, action='store', default='name', choices=['name', 'id'])
    input_arguments.add_argument('-s', '--significance',
                                 help="The significance level to use as a filtering cutoff. Default = 0.1",
                                 required=False, action='store', default=0.1, type=float)
    input_arguments.add_argument('-m', '--number', help="The minimum number of data sets in which a term must be "
                                                        "shared in order to be included in the output 'shared genes' "
                                                        "list.",
                                 required=False, action='store', default=5, type=int)
    return parser


# Creating logfiles. Logger_std is for writing to the console, logger is for writing to a file.
def create_logger(output_dir_path, project_name):
    try:
        mkdir(f"{output_dir_path}/results_{project_name}")
    except FileExistsError:
        print("Warning: there is already an analysis with this name")

    logger = logging.getLogger(__name__)
    logger_std = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    logger_std.setLevel(logging.INFO)

    error_file_handler = logging.FileHandler(f"{output_dir_path}/results_{project_name}/logfile_{project_name}.txt",
                                             mode='w')
    error_console_handler = logging.StreamHandler()
    error_console_handler.setLevel(logging.CRITICAL)
    logger.addHandler(error_file_handler)
    logger.addHandler(error_console_handler)

    console_handler = logging.StreamHandler()

    logger_std.addHandler(console_handler)
    logger_std.propagate = False
    logger.propagate=False
    return logger, logger_std


# Further handling of the command line input arguments
def argparse_error_handling(parser):
    required = parser.parse_args()
    optional = parser.parse_args()
    input_file_path = required.input
    output_file_path = required.outdir
    proj_name = required.project
    ontology = optional.ontology
    go_identifier = optional.identifier
    threshold = optional.significance
    min_shared_sets = optional.number

    # Check output directory and initialize logger instances
    if not isdir(output_file_path):
        exit()
    else:
        logger, logger_std = create_logger(output_dir_path=output_file_path, project_name=proj_name)

    # Check file path and create output directory
    if not isdir(input_file_path):
        logger.critical(f"GenewalkDownstream.py: The provided path does not exist. You entered {input_file_path}")
        exit()

    # Check the ontology
    if ontology.lower() not in ['bp', 'mf', 'cc', 'all']:
        logger.critical(f"GenewalkDownstream.py: the provided ontology, {ontology}, is not valid. It must be one of"
                        " ['bp', 'mf', 'cc', 'all']")
        exit()
    if ontology.lower() == 'bp':
        ontology = 'biological process'
    elif ontology.lower() == 'mf':
        ontology = 'molecular function'
    elif ontology.lower() == 'cc':
        ontology = 'cellular component'

    # Check identifier
    if go_identifier.lower() not in ['name', 'id']:
        logger.critical(f"GenewalkDownstream.py: the provided identifier, {go_identifier}, is not valid. "
                        "It must be one of ['name', 'id']")
        exit()
    if go_identifier == 'name':
        go_identifier = 'go_name'
    elif go_identifier == 'id':
        go_identifier = 'go_id'

    if not isinstance(threshold, float):
        logger.critical(f"GenewalkDownstream.py: the significance threshold must be a float. {threshold} is invalid.")
        exit()
    elif threshold <= 0 or threshold >= 1:
        logger.critical("GenewalkDownstream.py: the significance threshold should between 0 and 1. "
                        f"{threshold} is invalid.")
        exit()

    if not isinstance(min_shared_sets, int):
        logger.critical(f"GenewalkDownstream.py: the shared sets threshold must be an integer. {min_shared_sets} "
                        "is invalid.")
        exit()

    return input_file_path, output_file_path, proj_name, ontology, go_identifier, \
        threshold, min_shared_sets, logger, logger_std


def process_input(input_file_path):
    """pull in the gene walk results to dictionary"""
    files = [f for f in listdir(input_file_path) if f.endswith(".csv")]
    fns = [splitext(basename(x))[0] for x in files]
    raw_gw: dict = {}
    for i in range(len(fns)):
        raw_gw[fns[i]] = pd.read_csv(f"{input_file_path}/{files[i]}")
    return raw_gw


def filter_results(input_dfs: dict, significance: float):
    """Filter the raw results of Genewalk for only those gene:term pairs that fall below the supplied
        significance threshold.
        :param input_dfs: A dictionary containing raw Genewalk results dataframes from each provided Target
        :type significance: float """
    out = {}
    for name, df in input_dfs.items():
        out[name] = df[(df["global_padj"] < significance) & (df["gene_padj"] < significance)]
    return out


def make_go_set(set_dict: dict, go_type: str,
                go_cat: str):
    """Takes a dictionary of Genewalk results (those provided in the input directory) and
        produces a new dictionary containing the set of all significant GO terms of the specified
        arguments (i.e. GO IDs in the biological process category. Default behavior is to
        produce the set containing all GO categories (CC, BP, MF).
        :type go_cat: str
        :type go_type: str
        :type set_dict: Dict(Set)"""
    out_sets: dict = {}
    for name, df in set_dict.items():
        subset = df[df["go_domain"] == go_cat]
        out_sets[name] = list(set(subset[go_type]))
    return out_sets


def make_go_set_all_terms(set_dict: dict, go_type: str):
    out_sets: dict = {name: set(df[go_type]) for name, df in set_dict.items()}
    return out_sets


def make_membership_matrix(set_dictionary: dict):
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
    return name_affiliation_matrix


def min_x_shared_processes(affiliation_matrix: pd.DataFrame, threshold: float, go_type: str = "go_name"):
    """creates a pandas series of the GO terms that are shared among a minimum of x sets, where x is the
    threshold provided"""
    index_list = affiliation_matrix[affiliation_matrix.row_sums >= threshold].index
    return index_list.to_series(name=go_type)


def get_common_processes_at_threshold(signif_gw_results: pd.DataFrame,
                                      affiliation_matrix: pd.DataFrame, threshold: float, go_type: str = "go_name"):
    x_common_processes = {k: pd.merge(df2, min_x_shared_processes(affiliation_matrix, threshold),
                                      left_on="go_name", right_on=go_type) for k, df2 in signif_gw_results.items()}
    return x_common_processes


def make_colored_matrix(membership_mat):
    mat_sorted = membership_mat.sort_values(by=[*reversed(membership_mat.columns)],
                                            ascending=[False]*len(membership_mat.columns))
    mat_sorted = mat_sorted.astype(float)
    fig, ax = plt.subplots(figsize=(10, 11))
    cmap = sns.diverging_palette(250, 10, l=40, center="light", as_cmap=True)
    heatmap = sns.heatmap(mat_sorted.iloc[:, 0:len(mat_sorted.columns)-1], cmap=cmap, yticklabels=False)
    plt.xticks(rotation=45)
    return heatmap


parser_entrance = read_arguments()
args = parser_entrance.parse_args()

input_path, output_path, project_name, go_category, go_id, signif, \
    num_shared, log, log_std = argparse_error_handling(parser_entrance)

log.info('GenewalkDownStream Analysis')
log.info('Date run: '+str(date.today())+'\n')
lof.info('Run from dir:' + str(os.getcwd()))
log.info('Command line input: '+str(sys.argv)+'\n')


print(listdir(input_path))
raw_gw_res: dict = process_input(input_path)
filtered_gw_res: dict = filter_results(raw_gw_res, signif)

if go_category == 'all':
    gw_go_set = make_go_set_all_terms(set_dict=filtered_gw_res, go_type=go_id)
else:
    gw_go_set: dict = make_go_set(set_dict=filtered_gw_res, go_type=go_id, go_cat=go_category)

membership_matrix: pd.DataFrame = make_membership_matrix(set_dictionary=gw_go_set)
membership_matrix.to_csv(
    path_or_buf=f"{output_path}/results_{project_name}/"
    f"{project_name}_GO_{go_id}-{go_category}_membership_matrix.csv")
go_heatmap = make_colored_matrix(membership_matrix)
go_heatmap.set_title(f"{project_name}_{go_category}_{go_id}")
plt.savefig(f"{output_path}/results_{project_name}/"
            f"{project_name}_{go_category}-{go_id}_heatmap_{date.today()}.png")

num_shared_processes = {}
for i in range(1, len(raw_gw_res.keys())+1):
    num_shared_processes[i] = get_common_processes_at_threshold(signif_gw_results=filtered_gw_res,
                                                                affiliation_matrix=membership_matrix,
                                                                threshold=i, go_type=go_id)

for i, i_dict in num_shared_processes.items():
    mkdir(f"{output_path}/results_{project_name}/min_{i}_shared_results/")
    for k, v in i_dict.items():
        v.to_csv(path_or_buf=f"{output_path}/results_{project_name}/min_{i}_shared_results/"
                             f"{project_name}_{k}_min_{i}_shared_terms.csv")

