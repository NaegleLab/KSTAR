#!/usr/bin/env python3
import pandas as pd
import numpy as np
import scipy.stats as stats
from statsmodels.stats import multitest
import multiprocessing
from collections import defaultdict
import pickle
import os
from datetime import datetime
import itertools
from itertools import repeat
import argparse
import helpers
from os import path

import config
import summarize_activities

logger = helpers.get_logger("hypergeometric", "hypergeometric.log")

class KinaseActivity:
    def __init__(self, evidence, logger, phospho_type = 'Y'):
        """
        Kinase Activity calculates the estimated activity of kinases given an experiment using hypergeometric distribution.
        Hypergeometric distribution examines the number of protein sites found to be active in evidence compared to the 
        number of protein sites attributed to a kinase on a provided network.

        Parameters
        ----------
        evidence : pandas df
            a dataframe that contains (at minimum, but can have more) data columms as evidence to use in analysis and KSTAR_ACCESSION and KSTAR_SITE
        evidence_columns : dict
            columns corresponding to substrate and site
            required keys : substrate, site
        logger : Logger object
        """
        
        self.phospho_type = phospho_type
        self.set_evidence(evidence)
        

        self.networks = defaultdict()
        self.network_sizes = defaultdict()
        self.logger = logger
        self.normalizers = defaultdict()
        

        self.num_networks = None

        self.activities_list = None
        self.agg_activities = None
        self.activities = None


        self.aggregate = 'count'
        self.threshold = 1.0
        self.greater = True

        self.run_date = datetime.now().strftime("%d/%m/%Y %H:%M:%S")


    def check_data_columns(self):
        """
        Checks data columns to make sure column is in evidence and that evidence filtered on that data column 
        has at least one point of evidence. Removes all columns that do not meet criteria
        """
        new_data_columns = []
        evidence = self.evidence.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE]).agg(self.aggregate).reset_index()
        for col in self.data_columns:
            if col in self.evidence.columns:
                if self.greater:
                    if len(evidence[evidence[col] >= self.threshold]) > 0:
                        new_data_columns.append(col)
                    else:
                        self.logger.warning(f"{col} does not have any evidence")
                else:
                    if len(evidence[evidence[col] <= self.threshold]) > 0:
                        new_data_columns.append(col)
                    else:
                        self.logger.warning(f"{col} does not have any evidence")
            else:
                self.logger.warning(f"{col} not in evidence")
        self.data_columns = new_data_columns

    def set_data_columns(self, data_columns = None):
        """
        Sets the data columns to use in the kinase activity calculation
        If data_columns is None or an empty list then set data_columns to 
        be all columns that start with data:

        Checks all set columns to make sure columns are vaild after filtering evidence
        """
        if data_columns is None or data_columns == []:
            self.data_columns = []
            for col in self.evidence.columns:
                if col.startswith('data:'):
                    self.data_columns.append(col)  
        else:
            self.data_columns = data_columns
        
        self.check_data_columns()      

    def add_networks_batch(self,networks):
        for nid, network in networks.items():
            self.add_network(nid, network)
        self.num_networks = len(networks)

    def add_network(self, network_id, network, network_size = None):
        """
        Add network to be analyzed

        Parameters
        ----------
        network_id : str
            name of the network
        network : pandas DataFrame
            network with columns substrate_id, site, kinase_id
        """
        self.networks[network_id] = network
        if network_size is not None and isinstance(network_size, int):
            self.network_sizes[network_id] = network_size
        else:
            self.network_sizes[network_id] = len(network.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE]).size())
            self.logger.info(f'ADD NETWORK : Number of Accession Sites : {self.network_sizes[network_id]}')
    
    def get_run_date(self):
        return self.run_date

    def set_evidence(self, evidence):
        """
        Evidence to use in analysis

        Parameters
        ----------
        evidence : pandas DataFrame
            substrate sites with activity seen. 
            columns : dict for column mapping
                substrate : Uniprot ID (P12345)
                site : phosphorylation site (Y123)
        """
        columns = list(evidence.columns)
        if config.KSTAR_ACCESSION in columns and config.KSTAR_SITE in columns:
            phospho_type = tuple(self.phospho_type)
            self.evidence = evidence[evidence[config.KSTAR_SITE].str.startswith(phospho_type)].copy()
        else:
            self.logger.warning(f"Evidence not set. Evidence columns must include '{config.KSTAR_ACCESSION}' and '{config.KSTAR_SITE}' keys")
    

    def create_binary_evidence(self, data_columns = None, agg = 'count', threshold = 1.0,  greater = True):
        """
        Returns a binary evidence data frame according to the parameters passed in for method for aggregating
        duplicates and considering whether a site is included as evidence or not

        Parameters
        ----------
        data_columns : list
            columns that represent experimental result, if None, takes the columns that start with `data:'' in experiment. 
            Pass this value in as a list, if seeking to calculate on fewer than all available data columns
        threshold : float
            threshold value used to filter rows 
        agg : {'count', 'mean'}
            method to use when aggregating duplicate substrate-sites. 
            'count' combines multiple representations and adds if values are non-NaN
            'mean' uses the mean value of numerical data from multiple representations of the same peptide.
                NA values are droped from consideration.
        greater: Boolean
            whether to keep sites that have a numerical value >=threshold (TRUE, default) or <=threshold (FALSE)
        
        Returns
        -------
        evidence_binary : pd.DataFrame
            Matches the evidence dataframe of the kinact object, but with 0 or 1 if a site is included or not.
            This is uniquified and rows that are never used are removed.
        
        """
        evidence = self.evidence.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE]).agg(agg).reset_index()
        
        #set the binary evidence for whether a site is included
        evidence_binary = evidence.copy()
        for col in self.data_columns:
            if greater:
                evidence_binary[col].mask(evidence[col] >= threshold, 1, inplace=True)
                evidence_binary[col].mask(evidence[col] < threshold, 0, inplace=True)
            else:
                evidence_binary[col].mask(evidence[col] <= threshold, 1, inplace=True)
                evidence_binary[col].mask(evidence[col] > threshold, 0, inplace=True)

        #remove phosphorylation sites that were not selected in any experiment (useful for very large experiments where removing the need to copy data reduces time)
        evidence_binary.drop(evidence_binary[evidence_binary[self.data_columns].sum(axis=1) == 0].index, inplace = True) 
        return evidence_binary


    def calculate_kinase_activities(self, data_columns = None, agg = 'count', threshold = 1.0,  greater = True, cpus = 1, pool=None):
        """
        Calculates combined activity of experiments based that uses a threshold value to determine if an experiment sees a site or not
        To use values use 'mean' as agg
            mean aggregation drops NA values from consideration
        To use count use 'count' as agg - present if not na

        Parameters
        ----------
        data_columns : list
            columns that represent experimental result, if None, takes the columns that start with `data:'' in experiment. 
            Pass this value in as a list, if seeking to calculate on fewer than all available data columns
        threshold : float
            threshold value used to filter rows 
        agg : {'count', 'mean'}
            method to use when aggregating duplicate substrate-sites. 
            'count' combines multiple representations and adds if values are non-NaN
            'mean' uses the mean value of numerical data from multiple representations of the same peptide.
                NA values are droped from consideration.
        greater: Boolean
            whether to keep sites that have a numerical value >=threshold (TRUE, default) or <=threshold (FALSE)
        
        Returns
        -------
        activities : dict
            key : experiment
            value : pd DataFrame
                network : network name, from networks key
                kinase : kinase examined
                frequency : number of times kinase was seen in subgraph of evidence and network
                kinase_activity : hypergeometric kinase activity
        
        """

        #for each phosphoType, check that the 
        self.aggregate = agg
        self.threshold = threshold
        self.greater = greater
        self.set_data_columns(data_columns)

        self.evidence_binary = self.create_binary_evidence(data_columns = self.data_columns, agg = self.aggregate, threshold = self.threshold,  greater = self.greater)
        
        # if no data columns are provided use all columns that start with data:
        # data columns that filtered have no evidence are removed
        self.logger.info(f"Kinase Activity will be run on the following data columns: {','.join(self.data_columns)}")

        # MULTIPROCESSING
        if cpus > 1:
            filtered_evidence_list  = [self.evidence_binary[self.evidence_binary[col] == 1 ] for col in self.data_columns] 
            networks = itertools.repeat(self.networks)  
            network_sizes = itertools.repeat(self.network_sizes)
            iterable = zip(filtered_evidence_list, networks, network_sizes, self.data_columns)
            activities_list = pool.starmap(calculate_hypergeometric_activities, iterable)
        
        # SINGLE CORE PROCESSING
        else:
            activities_list =[]
            for col in self.data_columns:
                filtered_evidence = self.evidence_binary[self.evidence_binary[col] == 1]
                act = calculate_hypergeometric_activities(filtered_evidence, self.networks, self.network_sizes, col)
                act['data'] = col
                activities_list.append(act)
        self.num_networks = len(self.network_sizes)

        self.activities_list = pd.concat(activities_list)
        return self.activities_list


    def aggregate_activities(self, activities = None):
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
        if activities is None:
            activities = self.activities_list
        self.agg_activities = activities.groupby(['data', config.KSTAR_KINASE ]).agg(
            median_activity = ('kinase_activity', 'median'),
        ).reset_index()
        return self.agg_activities

    

def calculate_hypergeometric_single_network(evidence, network, network_size, network_id ):
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
        
def calculate_hypergeometric_activities(evidence, networks, network_sizes, name):
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
        
        results = []
        for network_id in networks.keys(): # calculate kinase activity within each network 
            result =calculate_hypergeometric_single_network(evidence, networks[network_id], network_sizes[network_id], network_id) 
            results.append(result)

        # combine results into single dataframe
        logger.info(f"hypergeometric analysis run on {name}")
        hyp_act = pd.concat(results)
        hyp_act = hyp_act.reset_index()
        hyp_act['data'] = name
        return hyp_act


"""
****************************************
Methods for running KSTAR pipeline
****************************************
"""
# def run_kstar_analysis(experiment, log, networks, phospho_types =['Y', 'ST'], data_columns = None, agg = 'count', threshold = 1.0,  greater = True):
#     """
#     A super method to establish a kstar KinaseActivity object from an experiment with an activity log
#     add the networks, calculate, aggregate, and summarize the activities into a final activity object

#     Parameters
#     ----------
#     experiment: pandas df
#         experiment dataframe that has been mapped, includes KSTAR_SITE, KSTAR_ACCESSION, etc.
#     log: logger object
#         Log to write activity log error and update to
#     networks: dictionary of dictionaries
#         Outer dictionary keys are 'Y' and 'ST'.
#         Establish a network by loading a pickle of desired networks. See the helpers and config file for this.
#         If downloaded from FigShare, then the GLOBAL network pickles in config file can be loaded
#         For example: networks['Y'] = pickle.load(open(config.NETWORK_Y_PICKLE, "rb" ))
#     phospho_types: {['Y', 'ST'], ['Y'], ['ST']}
#         Which substrate/kinaset-type to run activity for: Both ['Y, 'ST'] (default), Tyrosine ['Y'], or Serine/Threonine ['ST']
#     data_columns : list
#         columns that represent experimental result, if None, takes the columns that start with `data:'' in experiment. 
#         Pass this value in as a list, if seeking to calculate on fewer than all available data columns
#     agg : {'count', 'mean'}
#         method to use when aggregating duplicate substrate-sites. 
#         'count' combines multiple representations and adds if values are non-NaN
#         'mean' uses the mean value of numerical data from multiple representations of the same peptide.
#             NA values are droped from consideration.
#     threshold : float
#         threshold value used to filter rows
#     greater: Boolean
#         whether to keep sites that have a numerical value >=threshold (TRUE, default) or <=threshold (FALSE)
    
#     Returns
#     -------
#     kinactDict: dictionary of Kinase Activity Objects
#         Outer keys are phosphoTypes run 'Y' and 'ST'
#         Includes the activities dictionary (see calculate_kinase_activities)
#         aggregation of activities across networks (see aggregate activities)
#         activity summary (see summarize_activities)

#     """


#     kinact_dict = {}

#     # For each phosphoType of interest, establish a kinase activity object on a filtered dataset and run, aggregate, and summarize activity
#     for phospho_type in phospho_types:
#         #first check that networks for the phosphotypes were passed in
#         if phospho_type not in networks:
#             print("ERROR: Please pass networks as dictionary with phosphotype key")
#         #filter the experiment (log how many are of that type)
#         if phospho_type == 'ST':
#             experiment_sub = experiment[(experiment.KSTAR_SITE.str.contains('S')) | (experiment.KSTAR_SITE.str.contains('T'))]
#             log.info("Running Serine/Threonine Kinase Activity Analysis")
#         elif phospho_type == 'Y':
#             experiment_sub = experiment[(experiment.KSTAR_SITE.str.contains('Y'))]
#             log.info("Running Tyrosine Kinase Activity Analysis")

#         else:
#             print("ERROR: Did not recognize phosphoType %s, which should only include 'Y' or 'ST' "%(phospho_type))
#             return
#         kinact = KinaseActivity(experiment_sub, log, phospho_type=phospho_type)
#         kinact.add_networks_batch(networks[phospho_type])
#         kinact.calculate_kinase_activities(data_columns, agg=agg, threshold=threshold, greater=greater)
#         kinact.aggregate_activities()
#         kinact.activities = kinact.summarize_activities()
#         kinact_dict[phospho_type] = kinact
#     return kinact_dict



def run_kstar_analyis(experiment, log, networks, phospho_type ='Y', data_columns = None, agg = 'count', threshold = 1.0,  greater = True, cpus = 1, pool=None):

    experiment_sub = experiment[(experiment.KSTAR_SITE.str.contains(phospho_type))]
    kinact = KinaseActivity(experiment_sub, log, phospho_type=phospho_type)
    kinact.add_networks_batch(networks)
    kinact.calculate_kinase_activities(data_columns=data_columns, agg=agg, threshold=threshold, greater=greater, cpus=cpus, pool=pool)
    kinact.aggregate_activities()
    kinact.activities = summarize_activities.summarize_activities(kinact.agg_activities, 'median_activity')
    return kinact

def save_kstar_results(kinact, name):
    kinact.activities_list.to_csv(f"{name}_activities_list.tsv", sep = '\t', index=False)
    kinact.agg_activities.to_csv(f"{name}_aggregated_activities.tsv", sep = '\t', index = False)
    kinact.activities.to_csv(f"{name}_activities.tsv", sep = '\t', index = True)
    kinact.evidence_binary.to_csv(f"{name}_binarized_experiment.tsv", sep='\t', index=False)

def main():
    
    parser = argparse.ArgumentParser(description='Parse KSTAR Activity Arguments')
    parser.add_argument('--exp_file', '--experiment_file', action='store', dest= 'exp_file', help='Experiment file location. csv or tsv file', required=True)
    parser.add_argument("--networks", action='store', dest= 'networks', help='Network pickle file location', required=True)
    parser.add_argument('--pevent','--phosphorylation_events', action = 'store', dest='pevent', help ='phosphorylation event type', choices=['Y','ST'], default='Y')
    parser.add_argument('--name', action = 'store', dest='name', help = 'experiment name', default='Experiment')
    parser.add_argument('--cols', '--data_columns', action='store', dest='data_columns', help = 'data_columns to use', nargs='*', default=None)
    parser.add_argument('--agg', '--activity_agg', action='store', dest='activity_agg', help = 'activity agg to use', default='count', choices =['count','mean'])
    parser.add_argument('--thresh', '--threshold',  action='store', dest='threshold', help = 'threshold to use for analysis', type = float, default=1.0)
    parser.add_argument("--cpus", action="store", dest="cpus", default=1, type=int)
    results = parser.parse_args()

    if path.exists(results.exp_file) and path.isfile(results.exp_file):
        filetype = results.exp_file.split('.')[-1]
        if filetype == 'csv':
            experiment = pd.read_csv(results.exp_file)
        elif filetype == 'tsv':
            experiment = pd.read_csv(results.exp_file, sep = '\t')
        else:
            log.error("Unrecognized experiment filetype. Please use a csv or tsv file")
            exit()
    else:
        log.error("Please provide a valid experiment file")
        exit()
    
    log = helpers.get_logger("hypergeometric", "hypergeometric.log")
    networks = pickle.load(open(results.networks, "rb"))

    if results.cpus > 1:
        pool = multiprocessing.Pool(processes = results.cpus)
    else:
        pool = None

   
    kinact = run_kstar_analyis(experiment, log, networks, results.pevent, results.data_columns, results.activity_agg, threshold=results.threshold, cpus = results.cpus, pool=pool)
    
    save_kstar_results(kinact, results.name)




    
    
if __name__ == "__main__":
    main()







