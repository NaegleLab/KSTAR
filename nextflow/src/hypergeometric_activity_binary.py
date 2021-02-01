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



def read_evidence_file(evidence_file, phospho_type):
    # read evidence 
    if path.exists(evidence_file) and path.isfile(evidence_file):
        filetype = evidence_file.split('.')[-1]
        if filetype == 'csv':
            evidence = pd.read_csv(evidence_file)
        elif filetype == 'tsv':
            evidence = pd.read_csv(evidence_file, sep = '\t')
        else:
            logger.error("Unrecognized evidence filetype. Please use a csv or tsv file")
            exit()
    else:
        logger.error("Please provide a valid evidence file")
        exit()
    # filter to only contain specified phospho types
    if phospho_type == 'ST':
        evidence = evidence[(evidence.KSTAR_SITE.str.contains('S')) | (evidence.KSTAR_SITE.str.contains('T'))]
        logger.info("Running Serine/Threonine Kinase Activity Analysis")
    elif phospho_type == 'Y':
        evidence = evidence[(evidence.KSTAR_SITE.str.contains('Y'))]
        logger.info("Running Tyrosine Kinase Activity Analysis")
    else:
        logger.error(f"ERROR: Did not recognize phosphoType {phospho_type}, which should only include 'Y' or 'ST'")
        raise TypeError(f"ERROR: Did not recognize phosphoType {phospho_type}, which should only include 'Y' or 'ST'") 
    return evidence


# def chunk_data_columns(data_columns, chunk_size):
#     """Yield successive n-sized chunks from lst."""
#     for i in range(0, len(data_columns), chunk_size):
#         yield data_columns[i:i + chunk_size]

def calculate_hypergeometric_activities(evidence_file, phospho_type, data_columns, network_directory, max_cpus):
    logger.info("Calculating hypergeometric activities on all evidences")
    network_files = []

    for file in os.listdir(network_directory):
        if file.endswith(".tsv"):
            network_files.append(file)
    
    activities_list = []
    if max_cpus > 1:
        with ProcessPoolExecutor(max_workers=max_cpus) as executor:

            for single_network_activity in executor.map(calculate_hypergeometric_activities_single_network, repeat(evidence_file), repeat(phospho_type), repeat(data_columns), repeat(network_directory), network_files):
                activities_list = activities_list + single_network_activity
    else:
        evidence = read_evidence_file(evidence_file, phospho_type)
        for network_file in network_files:
            single_network_activity = calculate_hypergeometric_activities_single_network(evidence, data_columns, network_directory, network_file)
            activities_list = activities_list + single_network_activity


    
    # zipped_arguments = zip(repeat(evidence), repeat(data_columns), repeat(network_directory), network_files)
    # activities_list = pool.starmap(calculate_hypergeometric_activities_single_network, zipped_arguments)


    activities_list = pd.concat(activities_list)
    return activities_list

def calculate_hypergeometric_activities_single_network(evidence_file, phospho_type, data_columns, network_directory, network_file):
    evidence = read_evidence_file(evidence_file, phospho_type)

    network = pd.read_table(os.path.join(network_directory, network_file))
    network_size = len(network)

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
        
        Parameters
        ----------
        evidence : pandas df
            subset of kstar evidence that has been filtered to only include evidence associated with evidence
        network_id : str
            network to use for analysis
        
        Returns
        -------
        results : pandas df
            Hypergeometric results of evidence for given network
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

    Parameters
    ---------
    activities : dict
        key : evidence
        value : kinase activity result
    
    Returns 
    ---------
    summaries : dict
        key : evidence
        value : summarized kinase activities accross networks
    """
    
    agg_activities = activities.groupby(['data', config.KSTAR_KINASE ]).agg(
        median_activity = ('kinase_activity', 'median'),
    ).reset_index()
    return agg_activities

def run_kstar_analysis(evidence_file, network_directory, phospho_type, data_columns, max_cpus):

    activities_list = calculate_hypergeometric_activities(evidence_file, phospho_type, data_columns, network_directory, max_cpus)

    agg_activities = aggregate_activities(activities_list)
    activities = summarize_activities.summarize_activities(agg_activities, 'median_activity')
    return activities_list, agg_activities, activities

def save_kstar_results(activities_list, agg_activities, activities, name):
    activities_list.to_csv(f"{name}_activities_list.tsv", sep = '\t', index=False)
    agg_activities.to_csv(f"{name}_aggregated_activities.tsv", sep = '\t', index = False)
    activities.to_csv(f"{name}_activities.tsv", sep = '\t', index = True)

def main():
    
    parser = argparse.ArgumentParser(description='Parse KSTAR Activity Arguments')
    parser.add_argument('--exp_file', '--evidence_file', action='store', dest= 'exp_file', help='evidence file location. csv or tsv file', required=True)
    parser.add_argument("--net_dir","--network_directory", action='store', dest= 'network_directory', help='Network directory of individual networks', required=True)
    parser.add_argument('--pevent','--phosphorylation_events', action = 'store', dest='pevent', help ='phosphorylation event type', choices=['Y','ST'], default='Y')
    parser.add_argument('--name', action = 'store', dest='name', help = 'evidence name', default='evidence')
    parser.add_argument('--cols', '--data_columns', action='store', dest='data_columns', help = 'data_columns to use', nargs='*', default=None)
    parser.add_argument("--max_cpus", action="store", dest="max_cpus", default=1, type=int)
    parser.add_argument("--chunk_size", action="store", dest="chunk_size", default=5, type=int)
    results = parser.parse_args()

    # pool = Pool(results.max_cpus)

    evidence = read_evidence_file(results.exp_file, results.pevent)
    
    data_columns = helpers.set_data_columns(evidence, results.data_columns)
    del evidence
    # networks = pickle.load(open(results.networks, "rb"))
   
    activities_list, agg_activities, activities = run_kstar_analysis(results.exp_file, results.network_directory, results.pevent, data_columns, results.max_cpus)
    
    save_kstar_results(activities_list, agg_activities, activities, results.name)




#%%
    
if __name__ == "__main__":
    main()

#%%

# experiment = pd.read_table("/Users/bj8th/Documents/GitHub/KSTAR/nextflow/results/pdx/Y/individual_experiments/WHIM30.PDXBC03/random_experiments/WHIM30.PDXBC03_random_experiments.tsv")
# network_directory = "/Users/bj8th/Documents/GitHub/KSTAR/RESOURCE_FILES/NETWORKS/NetworKIN/Y/INDIVIDUAL_NETWORKS"
# data_columns = helpers.set_data_columns(experiment, None)
# %%
