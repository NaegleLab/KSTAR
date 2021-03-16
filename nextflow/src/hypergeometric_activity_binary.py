#!/usr/bin/env python3
#%%
import pandas as pd
import scipy.stats as stats
import pickle
from os import path
import os
import argparse

from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Pool
from multiprocessing import cpu_count
from itertools import repeat

import helpers
import config
import summarize_activities

logger = helpers.get_logger("hypergeometric", "hypergeometric.log")

def chunk_data_columns(data_columns, chunk_size):
    """Yields chunks of data columns to use in hypergeometric activity

    Args:
        data_columns (list): list of data columns that are being processed
        chunk_size (int): size of each chunk

    Yields:
        list: list of lists of data column chunks
    """
    for i in range(0, len(data_columns), chunk_size):
        yield data_columns[i:i + chunk_size]

def calculate_hypergeometric_activities(evidence, data_columns, network_directory, max_cpus, chunk_size, network_size):
    """Calculates the hypergeometric activity of evidence for all data columns on all networks

    Args:
        evidence (pandas DataFrame): binary evidence of phosphorylation
        data_columns (string): data columns in evidence to use in calculating hypergeometric activity
        network_directory (string): network directory of all site-kinase networks
        max_cpus (int): number of cpus to use in multiprocessing
        chunk_size (int): chunk size to break data columns into - save memory in multiprocessing
        network_size (int): size of network, if <=0 then network size calculated based on number of accession sites in network

    Returns:
        pandas DataFrame: hypergeometric activites list
    """
    logger.info("Calculating hypergeometric activities on all experiments")
    network_files = []

    for file in os.listdir(network_directory):
        if file.endswith(".tsv"):
            network_files.append(file )
    
    activities_list = []
    if max_cpus > 1:
        chunked_data_columns = chunk_data_columns(data_columns, chunk_size)

        with ProcessPoolExecutor(max_workers=max_cpus) as executor:
            for chunk in chunked_data_columns:
                chunked_evidence = evidence[[config.KSTAR_ACCESSION, config.KSTAR_SITE] + chunk]
                chunked_evidence = chunked_evidence[(chunked_evidence == 1).any(axis=1)].copy()

                for single_network_activity in executor.map(calculate_hypergeometric_activities_single_network, repeat(chunked_evidence), repeat(chunk), repeat(network_directory), network_files, repeat(network_size)):
                    activities_list = activities_list + single_network_activity
    else:
        for network_file in network_files:
            single_network_activity = calculate_hypergeometric_activities_single_network(evidence, data_columns, network_directory, network_file)
            activities_list = activities_list + single_network_activity


    
    # zipped_arguments = zip(repeat(evidence), repeat(data_columns), repeat(network_directory), network_files)
    # activities_list = pool.starmap(calculate_hypergeometric_activities_single_network, zipped_arguments)


    activities_list = pd.concat(activities_list)
    return activities_list

def calculate_hypergeometric_activities_single_network(evidence, data_columns, network_directory, network_file, network_size):
    """Calculates the hypergeometric activity of all data columns of a single network after loading the network file

    Args:
        evidence (pandas DataFrame): binary evidence of phosphorylation
        data_columns (list): data columns to calculate hypergeometric activity on
        network_directory (string): directory of network file
        network_file (string): filename of network
        network_size (int): size of network, if <=0 then network size calculated based on number of accession sites in network

    Returns:
        pandas DataFrame: hypergeometric activity list
    """
    network = pd.read_table(os.path.join(network_directory, network_file))
    if network_size <= 0:
        network_size = len(network.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE]).size())

    # intersect = pd.merge(network, evidence, how='inner',
    #         on=[config.KSTAR_ACCESSION, config.KSTAR_SITE])

    activities_list =[]
    for col in data_columns:
        filtered_evidence = evidence[evidence[col] == 1][[config.KSTAR_ACCESSION, config.KSTAR_SITE, col]]
        activity = calculate_hypergeometric_activities_single_data_column_single_network(filtered_evidence, network, network_size )
        activity['data'] = col
        activity['network'] = network_file
        activity = activity.reset_index()
        activities_list.append(activity)

    return activities_list

def calculate_hypergeometric_activities_single_data_column_single_network(evidence, network, network_size):
    """
    Hypergeometric Cumulative Distribution Function calculated for each kinase given evidence
        k : number of times kinase seen in evidence
        M : number of unique sites in network
        n : number of times kinase seen in network
        N : size of evidence
    
    Args:
        evidence (pandas df): subset of kstar evidence that has been filtered to only include evidence associated with experiment
        network (pandas df): network to use for analysis
        network_size (int): size of network, if <=0 then network size calculated based on number of accession sites in network
    
    Returns:
        pandas DataFrame: Hypergeometric results of evidence for given network
            index : kinase_id
            columns
                frequency : number of times kinase seen in network
                kinase_activity : activity derived from hypergometric cdf 
    """
    intersect = pd.merge(network, evidence, how='inner',
        on=[config.KSTAR_ACCESSION, config.KSTAR_SITE])
    counts = intersect.groupby(config.KSTAR_KINASE).size()
    N = len(intersect.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE]).size())
    results = pd.DataFrame(counts, columns = ['frequency'])

    results['kinase_activity'] = 1.0

    K = network.groupby(config.KSTAR_KINASE).size()
    
    kinases = counts.index
    for kin in kinases:
        k = 0
        if counts.loc[kin] > 0:
            k = counts.loc[kin] - 1
        prb = stats.hypergeom.sf(
            k = int(k), 
            M = int(network_size), 
            n = int(K.loc[kin]), 
            N = int(N)
            )
        results.at[kin, 'kinase_activity'] = prb

    kinases = network[config.KSTAR_KINASE].unique()
    for kin in kinases:
        if kin not in results.index:
            results.at[kin,'frequency'] = 0
            results.at[kin,'kinase_activity'] = 1.0
    return results
    
def aggregate_activities(activities):
    """
    Aggregate network activity using median for all activities

    Args:
        activities (pandas DataFrame): activity list from hypergeometric activity calculations
    
    Returns:
        pandas DataFrame : activities aggregated to calculate median p-value over all networks for each data column and kinase
    """
    
    agg_activities = activities.groupby(['data', config.KSTAR_KINASE ]).agg(
        median_activity = ('kinase_activity', 'median'),
    ).reset_index()
    return agg_activities

def run_kstar_analysis(experiment, network_directory, phospho_type, data_columns, max_cpus, chunk_size, network_size):
    """Runs the kstar hypergeometric analysis including aggregating and summarizing results

    Args:
        experiment (pandas DataFrame): experiment to use for analysis
        network_directory (string): directory of networks
        phospho_type (string): {Y,ST}
        data_columns (list): data columns to use in analysis
        max_cpus (int): number of cpus to use in parallelization
        chunk_size (int): chunk size to use to reduce memory
        network_size (int): size of network, if <=0 then network size calculated based on number of accession sites in network

    Raises:
        TypeError: raised if phospho type is not correct. Phospho typee must be Y or ST

    Returns:
        activities_list (pandas DataFrame): activities list of all data columns
        agg_activities (pandas DataFrame): aggregate median activites
        activites (pandas DataFrame): transformed median activities. columns : data columns, index : kinase
    """

    
    if phospho_type == 'ST':
        experiment_sub = experiment[(experiment.KSTAR_SITE.str.contains('S')) | (experiment.KSTAR_SITE.str.contains('T'))]
        logger.info("Running Serine/Threonine Kinase Activity Analysis")
    elif phospho_type == 'Y':
        experiment_sub = experiment[(experiment.KSTAR_SITE.str.contains('Y'))]
        logger.info("Running Tyrosine Kinase Activity Analysis")

    else:
        logger.error(f"ERROR: Did not recognize phosphoType {phospho_type}, which should only include 'Y' or 'ST'")
        raise TypeError(f"ERROR: Did not recognize phosphoType {phospho_type}, which should only include 'Y' or 'ST'") 
    
    activities_list = calculate_hypergeometric_activities(experiment_sub, data_columns, network_directory, max_cpus, chunk_size, network_size)

    agg_activities = aggregate_activities(activities_list)
    activities = summarize_activities.summarize_activities(agg_activities, 'median_activity')
    return activities_list, agg_activities, activities

def save_kstar_results(activities_list, agg_activities, activities, name):
    """Saves results of hypergeometric activity calculations

    Args:
        activities_list (pandas DataFrame): activities list of all data columns
        agg_activities (pandas DataFrame): aggregate median activites
        activites (pandas DataFrame): transformed median activities. columns : data columns, index : kinase
        name (string): name of experiment
    """
    activities_list.to_csv(f"{name}_activities_list.tsv", sep = '\t', index=False)
    agg_activities.to_csv(f"{name}_aggregated_activities.tsv", sep = '\t', index = False)
    activities.to_csv(f"{name}_activities.tsv", sep = '\t', index = True)

def main():
    """Reads in parser arguments, loads experiment, and runs analysis
    """
    parser = argparse.ArgumentParser(description='Parse KSTAR Activity Arguments')
    parser.add_argument('--exp_file', '--experiment_file', action='store', dest= 'exp_file', help='Experiment file location. csv or tsv file', required=True)
    parser.add_argument("--net_dir","--network_directory", action='store', dest= 'network_directory', help='Network directory of individual networks', required=True)
    parser.add_argument('--pevent','--phosphorylation_events', action = 'store', dest='pevent', help ='phosphorylation event type', choices=['Y','ST'], default='Y')
    parser.add_argument('--name', action = 'store', dest='name', help = 'experiment name', default='Experiment')
    parser.add_argument('--cols', '--data_columns', action='store', dest='data_columns', help = 'data_columns to use', nargs='*', default=None)
    parser.add_argument("--max_cpus", action="store", dest="max_cpus", default=1, type=int)
    parser.add_argument("--chunk_size", action="store", dest="chunk_size", default=5, type=int)
    parser.add_argument("--network_size", action="store", dest="network_size", default=-1, type=int)
    results = parser.parse_args()

    # pool = Pool(results.max_cpus)

    if path.exists(results.exp_file) and path.isfile(results.exp_file):
        filetype = results.exp_file.split('.')[-1]
        if filetype == 'csv':
            experiment = pd.read_csv(results.exp_file)
        elif filetype == 'tsv':
            experiment = pd.read_csv(results.exp_file, sep = '\t')
        else:
            logger.error("Unrecognized experiment filetype. Please use a csv or tsv file")
            exit()
    else:
        logger.error("Please provide a valid experiment file")
        exit()
    
    data_columns = helpers.set_data_columns(experiment, results.data_columns)

    # networks = pickle.load(open(results.networks, "rb"))
   
    activities_list, agg_activities, activities = run_kstar_analysis(experiment, results.network_directory, results.pevent, data_columns, results.max_cpus, results.chunk_size, results.network_size)
    
    save_kstar_results(activities_list, agg_activities, activities, results.name)




#%%
    
if __name__ == "__main__":
    main()

