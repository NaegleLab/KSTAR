#!/usr/bin/env python3
#%%
import pandas as pd 
import argparse 
import subprocess
import config

def main():
    parser = argparse.ArgumentParser(description='Parse KSTAR Activity Arguments')
    parser.add_argument('--exp_file', '--evidence_file', action='store', dest= 'exp_file', help='Experiment file location. csv or tsv file', required=True)
    parser.add_argument("--net_dir","--network_directory", action='store', dest= 'network_directory', help='Network directory of individual networks', required=True)
    parser.add_argument('--pevent','--phosphorylation_events', action = 'store', dest='pevent', help ='phosphorylation event type', choices=['Y','ST'], default='Y')
    parser.add_argument('--name', action = 'store', dest='name', help = 'experiment name', default='Experiment')
    parser.add_argument("--max_cpus", action="store", dest="max_cpus", default=1, type=int)

    results = parser.parse_args()

    experiment = pd.read_table(results.exp_file, nrows=5)
    data_columns = list(experiment.columns) 
    data_columns.remove(config.KSTAR_ACCESSION)
    data_columns.remove(config.KSTAR_SITE)
    # data_columns = " ".join(data_columns)

    arguments = ['hypergeometric_activity_binary.py',
        "--evidence_file", results.exp_file, 
        "--network_directory", results.network_directory,
        "--pevent", results.pevent,
        "--name", results.name,
        "--max_cpus", str(results.max_cpus)]
    data_column_arguments = ["--data_columns"] + data_columns
    arguments = arguments + data_column_arguments
    subprocess.call(
        arguments)
#%%
if __name__=="__main__":
    main()


#%%
# experiment_file = "/Users/bj8th/Documents/GitHub/KSTAR/nextflow/results/pdx/Y/individual_experiments/WHIM30.PDXBC03/random_experiments/WHIM30.PDXBC03_random_experiments.tsv"
# network_directory = "/Users/bj8th/Documents/GitHub/KSTAR/RESOURCE_FILES/NETWORKS/NetworKIN/Y/INDIVIDUAL_NETWORKS"
# pevent = "Y"
# name = "test"
# max_cpus = 8
# experiment = pd.read_table(experiment_file)
# columns = list(experiment.columns) 
# columns.remove(config.KSTAR_ACCESSION)
# columns.remove(config.KSTAR_SITE)
# data_columns = " ".join(columns)

# subprocess.call(
#         ['hypergeometric_activity_binary.py',
#         "--experiment_file", experiment_file, 
#         "--network_directory", network_directory,
#         "--pevent", pevent,
#         "--name", name,
#         "--data_columns", data_columns,
#         "--max_cpus", str(max_cpus)
#         ])
# %%
