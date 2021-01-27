#!/usr/bin/env python3

#%%
import pandas as pd
import scipy.stats as stats
import pickle
from os import path
import argparse

from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count
from itertools import repeat

import helpers
import config
import summarize_activities

logger = helpers.get_logger("hypergeometric", "hypergeometric.log")




def calculate_hypergeometric_activities(evidence, data_columns, networks):
    logger.info("Calculating hypergeometric activities on all experiments")
    activities_list =[]
    for col in data_columns:
        filtered_evidence = evidence[evidence[col] == 1]
        act = calculate_hypergeometric_activities_single_data_column(filtered_evidence, networks, col)
        act['data'] = col
        activities_list.append(act)
    

    activities_list = pd.concat(activities_list)
    return activities_list

def calculate_hypergeometric_activities_single_data_column(evidence, networks, name):
        """
        Perform hypergeometric kinase activity analysis given evidence on all networks
        
        Parameters
        ----------
        evidence : pandas df
            subset of class evidence variable where data is filtered based on experiment
        name : str
            name of experiment being performed
            
        Returns
        ---------
        fdr_act : pd DataFrame
            network : network name, from networks key
            frequency : number of times kinase was seen in subgraph of evidence and network
            kinase_activity : hypergeometric kinase activity
            fdr_corrected_kinase_activity : kinase activity after fdr correction
            significant : whether kinase activity is significant based on fdr alpha
        combined : pd DataFrame
            significant : number of networks where kinase was found to be significant
            fraction_significant : fraction of networks kinase was found to be significant through FDR correction
            avg_p_value : combined p-values of kinase using mean
            median_p_value : combined p-values of kinase using median
        """
        
        logger.info(f"Running hypergeometric analysis on {name}")
        hyp_act_list = []
        with ProcessPoolExecutor(max_workers=cpu_count()) as executor:
            for hyp_act in executor.map(calculate_hypergeometric_activities_single_data_column_single_network, repeat(evidence), networks.items()):
                hyp_act_list.append(hyp_act)

        
        # for network_id in networks.keys(): # calculate kinase activity within each network 
        #     result =calculate_hypergeometric_single_network(evidence, networks[network_id], network_sizes[network_id], network_id) 
        #     results.append(result)

        # combine results into single dataframe
        logger.info(f"hypergeometric analysis run on {name}")
        activity = pd.concat(hyp_act_list)
        activity = activity.reset_index()
        activity['data'] = name
        return activity


def calculate_hypergeometric_activities_single_data_column_single_network(evidence, network_info ):
        """
        Hypergeometric Cumulative Distribution Function calculated for each kinase given evidence
            k : number of times kinase seen in evidence
            M : number of unique sites in network
            n : number of times kinase seen in network
            N : size of evidence
        
        Parameters
        ----------
        evidence : pandas df
            subset of kstar evidence that has been filtered to only include evidence associated with experiment
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
        network_id, network = network_info
        network_size = len(network)
        
        
        intersect = pd.merge(network, evidence, how='inner',
            on=[config.KSTAR_ACCESSION, config.KSTAR_SITE])
        
        logger.info(network.columns)
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
        results['network'] = network_id
        return results
        
def aggregate_activities(activities):
    """
    Aggregate network activity using median for all activities

    Parameters
    ---------
    activities : dict
        key : Experiment
        value : kinase activity result
    
    Returns 
    ---------
    summaries : dict
        key : experiment
        value : summarized kinase activities accross networks
    """
    
    agg_activities = activities.groupby(['data', config.KSTAR_KINASE ]).agg(
        median_activity = ('kinase_activity', 'median'),
    ).reset_index()
    return agg_activities

def run_kstar_analysis(experiment, networks, phospho_type, data_columns):

    experiment_sub = experiment[(experiment.KSTAR_SITE.str.contains(phospho_type))]
    activities_list = calculate_hypergeometric_activities(experiment_sub, data_columns, networks)

    agg_activities = aggregate_activities(activities_list)
    activities = summarize_activities.summarize_activities(agg_activities, 'median_activity')
    return activities_list, agg_activities, activities

def save_kstar_results(activities_list, agg_activities, activities, name):
    activities_list.to_csv(f"{name}_activities_list.tsv", sep = '\t', index=False)
    agg_activities.to_csv(f"{name}_aggregated_activities.tsv", sep = '\t', index = False)
    activities.to_csv(f"{name}_activities.tsv", sep = '\t', index = True)




def check_data_columns(evidence, data_columns):
        """
        Checks data columns to make sure column is in evidence and that evidence filtered on that data column 
        has at least one point of evidence. Removes all columns that do not meet criteria
        """
        new_data_columns = []
        for col in data_columns:
            if col in evidence.columns:
                if len(evidence[evidence[col] >= 1]) > 0:
                    new_data_columns.append(col)
                else:
                    logger.warning(f"{col} does not have any evidence")

            else:
                logger.warning(f"{col} not in evidence")
        data_columns = new_data_columns
        return data_columns

def set_data_columns(evidence, data_columns = None):
    """
    Sets the data columns to use in the kinase activity calculation
    If data_columns is None or an empty list then set data_columns to 
    be all columns that start with data:

    Checks all set columns to make sure columns are vaild after filtering evidence
    """
    if data_columns is None or data_columns == []:
        data_columns = []
        for col in evidence.columns:
            if col.startswith('data:'):
                data_columns.append(col)  
    else:
        data_columns = data_columns
    
    data_columns = check_data_columns(evidence, data_columns)
    return data_columns

def main():
    
    parser = argparse.ArgumentParser(description='Parse KSTAR Activity Arguments')
    parser.add_argument('--exp_file', '--experiment_file', action='store', dest= 'exp_file', help='Experiment file location. csv or tsv file', required=True)
    parser.add_argument("--networks", action='store', dest= 'networks', help='Network pickle file location', required=True)
    parser.add_argument('--pevent','--phosphorylation_events', action = 'store', dest='pevent', help ='phosphorylation event type', choices=['Y','ST'], default='Y')
    parser.add_argument('--name', action = 'store', dest='name', help = 'experiment name', default='Experiment')
    parser.add_argument('--cols', '--data_columns', action='store', dest='data_columns', help = 'data_columns to use', nargs='*', default=None)
    parser.add_argument("--cpus", action="store", dest="cpus", default=1, type=int)
    results = parser.parse_args()

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
    
    data_columns = set_data_columns(experiment, results.data_columns)

    networks = pickle.load(open(results.networks, "rb"))
   
    activities_list, agg_activities, activities = run_kstar_analysis(experiment, networks, results.pevent, data_columns)
    
    save_kstar_results(activities_list, agg_activities, activities, results.name)




#%%
    
if __name__ == "__main__":
    main()

# #%%
# experiment = pd.read_table("/Users/bj8th/Documents/GitHub/KSTAR/nextflow/results/Y/binary_experiment/test_binarized_experiment.tsv")
# networks = pickle.load(open("/Users/bj8th/Documents/GitHub/KSTAR/RESOURCE_FILES/NETWORKS/NetworKIN/Y/NetworKIN_Y_compendia_2000_limit_10.p", "rb"))
# data_columns = ['data:PRE', 'data:EOE']
# #%%
# activities_list, agg_activities, activities = run_kstar_analysis(experiment, networks, "Y", data_columns)

