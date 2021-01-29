#!/usr/bin/env python3
import pandas as pd
import multiprocessing
import itertools
import argparse
import config
from os import path
import logging

#%%
#%%
def build_filtered_experiment(experiment, compendia, filtered_compendia, num_random_experiments, name ,selection_type='KSTAR_NUM_COMPENDIA_CLASS'):
    rand_experiments = compendia[[config.KSTAR_ACCESSION, config.KSTAR_SITE]]
    if len(experiment) == 0:
        empty_columns = [f"{name}:{i}" for i in range(num_random_experiments)]
        rand_experiments = pd.concat([rand_experiments,pd.DataFrame(columns=empty_columns)])
        return rand_experiments


    sizes = experiment.groupby(selection_type).size()
    for i in range(num_random_experiments):
        rand_experiment_list = []
        for num, size in sizes.iteritems():
            filtered = filtered_compendia[int(num)]
            filtered_random = filtered.sample(size)
            filtered_random[f"{name}:{i}"] = 1
            rand_experiment_list.append(filtered_random)
        rand_experiment = pd.concat(rand_experiment_list)
        rand_experiments = pd.merge(rand_experiments, rand_experiment, how = 'left', on = [config.KSTAR_ACCESSION, config.KSTAR_SITE])
    # name = name.replace(":","_")
    rand_experiments = rand_experiments.set_index([config.KSTAR_ACCESSION, config.KSTAR_SITE])
    rand_experiments = rand_experiments.dropna(how='all', axis = 0)
    rand_experiments = rand_experiments.reset_index()


    rand_experiments.to_csv(f"{name}_random_experiments.tsv", sep="\t",index=False)
    return rand_experiments

def build_random_experiments_single_data_column(binary_evidence, compendia, num_random_experiments, phosphorylation_event, data_column, selection_type='KSTAR_NUM_COMPENDIA_CLASS'):
    """
    Given an experimental dataframe and the human phospho compendia, build random experiments such that each random experiment takes on the same
    distribution with respect to the study bias defined as either NUM_COMPENDIA (total number of compendia a site is annotated in) or 
    NUM_COMPENDIA_CLASS (whether it has low < 1, medium (1-3), or high study bias(>3)).
    
    Parameters
    ----------
    binary_evidence: df
        KSTAR mapped experimental dataframe that has been binarized by kstar_activity generation
    greater: Boolean
        whether to keep sites that have a numerical value >=threshold (TRUE, default) or <=threshold (FALSE)
    num_random_experiments: int
        Number of random experiments to generate for each data_column
    phosphorylation_event: {'Y', 'ST'}
        Which substrate/kinaset-type to generate random experiments for
    data_columns : list
        columns that represent experimental result
    selection_type: {'KSTAR_NUM_COMPENDIA', 'KSTAR_NUM_COMPENDIA_CLASS'}
        Whether to sample according to the absolute number of compendia or low, medium, or high study bias groupings

    Returns
    -------
    rand_experiments: pandas.DataFrame
        Dataframe of random experiments with NaN where phosphorylation sites are not selected, and 1 if they are for that experiment
    """
    #check parameters
    if selection_type != 'KSTAR_NUM_COMPENDIA':
        if selection_type != 'KSTAR_NUM_COMPENDIA_CLASS':
            raise ValueError('selection_type must be either KSTAR_NUM_COMPENDIA or KSTAR_NUM_COMPENDIA_CLASS')

    experiment = binary_evidence
    if phosphorylation_event == 'ST':
        compendia = compendia[(compendia.KSTAR_SITE.str.contains('S')) | (compendia.KSTAR_SITE.str.contains('T'))]
        experiment = experiment[(experiment.KSTAR_SITE.str.contains('S')) | (experiment.KSTAR_SITE.str.contains('T'))]
            
    elif phosphorylation_event == 'Y':
        compendia = compendia[(compendia.KSTAR_SITE.str.contains('Y'))]
        experiment = experiment[(experiment.KSTAR_SITE.str.contains('Y'))]
    else:
        raise ValueError('phosphorylation_event must be Y or ST')


    compendia = compendia[[config.KSTAR_ACCESSION, config.KSTAR_SITE, selection_type]]
    compendia = compendia.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE]).max().reset_index() #uniquify the compendia by KSTAR_ACCESSION and KSTAR_SITE
    

    #experiment = experiment.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE]).agg(agg).reset_index()
    
    sizes = compendia[selection_type].unique()
    filtered_compendia = {}
    for s in sizes:
        filtered_compendia[s] = compendia[compendia[selection_type] == s][[config.KSTAR_ACCESSION, config.KSTAR_SITE]]
    
    
    filtered_experiment = experiment[experiment[data_column] ==1 ]
    build_filtered_experiment(experiment, compendia, filtered_compendia, num_random_experiments, data_column)
     
def parse_args():
    parser = argparse.ArgumentParser(description='Parse Mapping Inference Arguments')
    parser.add_argument('-e', '--exp_file', '--experiment_file', action='store', dest= 'exp_file', help='Experiment file location. csv or tsv file', required=True)
    parser.add_argument('-p','--pevent','--phosphorylation_events', action = 'store', dest='pevent', help ='phosphorylation event type', choices=['Y','S','T','ST', 'STY'], default='STY')
    parser.add_argument('-d', '--data', '--data_column', action='store', dest='data_column', help = 'data_column to use')
    parser.add_argument('-num', '--num_random_experiments',  action='store', dest='num', help = 'Number of random experiments to generate', type = int, default=150)
    parser.add_argument("-comp", "--compendia", action="store", dest="compendia", help = "human reference compnedia", )
    parser.add_argument('-s', '--selection_type', action='store', dest='selection_type', help='KSTAR_NUM_COMPENDIA or KSTAR_NUM_COMPENDIA_CLASS', type=str, default='KSTAR_NUM_COMPENDIA_CLASS')
    results = parser.parse_args()
    return results


def main():
    results = parse_args()
    experiment = pd.read_table(results.exp_file)
    compendia = pd.read_csv(results.compendia)
    data_column = results.data_column
    build_random_experiments_single_data_column(binary_evidence = experiment, compendia = compendia, num_random_experiments=results.num, phosphorylation_event=results.pevent, data_column=data_column)
    
if __name__ == "__main__":
    main()

