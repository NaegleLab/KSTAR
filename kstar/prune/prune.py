import pandas as pd 
import logging
import random 
import numpy as np
import multiprocessing
import itertools
import argparse
from collections import defaultdict
import pickle
import os
from os import path
from kstar import helpers, config


class Prune:
    def __init__(self, network, logger, phospho_type = 'Y',
        columns = {'accession' : 'substrate_acc','site' : 'site', 'kinase' : 'Kinase Name', 'score' : 'score' }):

        """
        Pruning Algorithm used for Kinase Activity Algorithm. 
        Parameters
        ----------
        network : pandas df
            kinase-site prediction network where there is an accession, site, kinase, and score column
        logger : 
            logger used for pruning
        phospho_type : str
            phospho_type(s) to use when building pruned networks
        columns : dict
            relevant columns in network
        """
        self.phospho_type = tuple(phospho_type)
        self.logger = logger
        self.columns = columns
        self.network = network[network[self.columns['site']].str.startswith(self.phospho_type)]
        self.network = self.network[self.network[self.columns['score']]>0]
        self.kinases = self.network[self.columns['kinase']].unique()

        self.compendia = config.HUMAN_REF_COMPENDIA
        self.compendia = self.compendia[self.compendia[config.KSTAR_SITE].str.startswith(self.phospho_type)]
        self.compendia = self.compendia[[config.KSTAR_ACCESSION, config.KSTAR_SITE, 'KSTAR_NUM_COMPENDIA']]
        self.compendia = self.compendia.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE]).max().reset_index()

        self.network = pd.merge(self.network, self.compendia, how = 'inner', left_on = [self.columns['accession'], self.columns['site']], right_on=[config.KSTAR_ACCESSION, config.KSTAR_SITE] )
        self.network = self.network[[config.KSTAR_ACCESSION, config.KSTAR_SITE, self.columns['kinase'], self.columns['score'], 'KSTAR_NUM_COMPENDIA']]
        self.logger.info("Pruning initialization complete")
    
    def build_pruned_network(self, network, kinase_size, site_limit):
        """
        Builds a heuristic pruned network where each kinase has a specified number of connected sites and
        each site has an upper limit to the number of kinases it can connect to

        Parameters
        ----------
        network : pandas DataFrame
            network to build pruned network on
        kinase_size: int
            number of sites each kinase should connect to
        site_limit :int
            upper limit of number of kinases a site can connect to
        
        Returns
        ---------
        pruned network : pandas DataFrame
            subset of network that has been pruned
        """
        network = network.copy()
        kinases = network[self.columns['kinase']].unique()
        pruned_network = pd.DataFrame(columns = network.columns)
        for i in range(kinase_size):
            self.logger.info(f"Pruning {i}")
            random.shuffle(kinases)
            for kinase in kinases:
                kinase_network = network[network[self.columns['kinase']] == kinase]
                
                sample = kinase_network.sample(weights=self.columns['score'])           # Heuristically select kinase-substrate edge
                network.drop(sample.index, inplace = True)                              # Remove edge from remaining network
                pruned_network = pd.concat([pruned_network, sample])                    # add edge to pruned network
                
                # Remove substrate site from network if 
                sample = sample.iloc[0]
                pruned_site = pruned_network[config.KSTAR_SITE] == sample[config.KSTAR_SITE]
                pruned_accession = pruned_network[config.KSTAR_ACCESSION] == sample[config.KSTAR_ACCESSION]
                pruned_match = pruned_network[pruned_site & pruned_accession]
                if len(pruned_match) >= site_limit:
                    network_site = network[config.KSTAR_SITE] == sample[config.KSTAR_SITE]
                    network_accession = network[config.KSTAR_ACCESSION] == sample[config.KSTAR_ACCESSION]
                    network = network[~(network_site & network_accession)]
                    # network.drop(network_match.index, inplace=True)
        return pruned_network
    
    def compendia_pruned_network(self, compendia_sizes, site_limit):
        """
        Builds a compendia-pruned network that takes into account compendia size limits per kinase
        
        Parameters
        ----------
        compendia_sizes : dict
            key : compendia size
            value : number of sites to connect to kinase 
        site_limit : int 
            upper limit of number of kinases a site can connect to
        
        Returns
        -------
        pruned_network : pandas DataFrame
            subset of network that has been pruned according to compendia ratios
        """
        pruned_network = []
        for comp_size, kinase_size in compendia_sizes.items():
            self.logger.info(f"Building pruned network for compendia {comp_size} with {kinase_size} sites per kinase")
            network = self.network[self.network['KSTAR_NUM_COMPENDIA'] == comp_size]
            pruned_network.append(self.build_pruned_network(network, kinase_size, site_limit))
        pruned_network = pd.concat(pruned_network)
        return pruned_network
        
    def calculate_compendia_sizes(self, kinase_size):
        """
        Calculates the number of sites per compendia size that a kinase should connect to using same ratios 
        of compendia sizes as found in compendia

        Parameters
        ----------
        kinase_size: int
            number of sites each kinase should connect to
        Returns
        --------
        sizes : dict
            key : compendia size
            value : number of sites each kinase should pull from given compendia size
        """
        self.logger.info(f"Calculating compendia kinase sizes when kinase size is {kinase_size}")
        sizes = [-np.inf]
        multiplier = kinase_size
        while sum(sizes) < kinase_size:
            sizes = self.compendia.groupby('KSTAR_NUM_COMPENDIA').size() / len(self.compendia) * multiplier
            sizes = sizes.astype(int)
            multiplier = multiplier + 1
        sizes = sizes.to_dict()
        self.logger.info(f"Compendia kinase sizes : {sizes}")
        return sizes

    def build_multiple_compendia_networks(self, kinase_size, site_limit, num_networks, network_id):
        """
        Builds multiple compendia-limited networks

        kinase_size: int
            number of sites each kinase should connect to
        site_limit :int
            upper limit of number of kinases a site can connect to
        
        num_networks: int
            number of networks to build
        network_id : str
            id to use for each network in dictionary
        
        Returns
        -------
        pruned_networks : dict
            key : <network_id>_<i>
            value : pruned network
        """

        compendia_sizes = self.calculate_compendia_sizes(kinase_size)

        # MULTIPROCESSING
        if config.PROCESSES > 1:
            pool = multiprocessing.Pool(processes = config.PROCESSES)
            comp_sizes = itertools.repeat(compendia_sizes, num_networks)
            limits = itertools.repeat(site_limit, num_networks)
            iterable = zip(comp_sizes, limits)
            pruned_networks = pool.starmap(self.compendia_pruned_network, iterable)
        
        # SINGLE CORE PROCESSING
        else:
            pruned_networks = []
            for i in range(num_networks):
                net = self.compendia_pruned_network(compendia_sizes, site_limit)
                pruned_networks.append(net)

        pruned_dict = {}
        for i in range(len(pruned_networks)):
            pruned_dict[f"{network_id}_{i}"] = pruned_networks[i]
        return pruned_dict

    def build_multiple_networks(self, kinase_size, site_limit, num_networks, network_id):
        # MULTIPROCESSING
        if config.PROCESSES > 1:
            pool = multiprocessing.Pool(processes = config.PROCESSES)
            network_iter = itertools.repeat(self.network, num_networks)
            size_iter = itertools.repeat(kinase_size, num_networks)
            limit_iter = itertools.repeat(site_limit, num_networks)
            iterable = zip(network_iter, size_iter, limit_iter)
            pruned_networks = pool.starmap(self.build_pruned_network, iterable)
        
        # SINGLE CORE PROCESSING
        else:
            pruned_networks = []
            for i in range(num_networks):
                net = self.build_pruned_network(self.network, kinase_size, site_limit)
                pruned_networks.append(net)

        pruned_dict = {}
        for i in range(len(pruned_networks)):
            pruned_dict[f"{network_id}_{i}"] = pruned_networks[i]
        return pruned_dict


"""
********************** RUNNING PRUNING ALGORITHM FROM COMMAND LINE ****************************
Arguments
---------
--network_file          Experiment file location. csv or tsv file       (required)
--output_directory      output directory for results (required)
--accession             Accession Column
--site                  site column
--score                 score column
--kinase                kinase column
--phospho_type          phospho type
--kinase_size           number of sites a kinase connects to (required)
--site_limit            upper limit of number of kinases can connect to (required)
--num_networks          number of networks to generate (required)
--network_id            name of network to use in building dictionary
--use_compendia         whether to use compendia ratios to build network    

"""
def parse_args():
    parser = argparse.ArgumentParser(description='Parse Network Pruning Arguments')
    parser.add_argument('--network_file', action='store', dest= 'network_file', help='Experiment file location. csv or tsv file', required=True)
    parser.add_argument('--odir','--output_directory', action = 'store', dest='odir', help = 'output directory for results', required=True)
    parser.add_argument('--accession', action = 'store', dest = 'accession', help = 'Accession Column', default = 'substrate_acc', required = False)
    parser.add_argument('--site', action = 'store', dest = 'site', help='site column', default='site',required = False)
    parser.add_argument('--score', action = 'store', dest = 'score', help='score column', default='score', required = False)
    parser.add_argument('--kinase', action = 'store', dest = 'kinase', help='kinase column', default='Kinase Name', required = False)
    parser.add_argument('--phospho_type', action='store', dest = 'phospho_type', help = 'phospho type', default='Y', choices = ['Y','ST'], required = False)
    parser.add_argument('--kinase_size', action='store', dest = 'kinase_size', type = int,help='number of sites a kinase connects to',required=True)
    parser.add_argument('--site_limit', action='store', dest='site_limit', type=int,help='upper limit of number of kinases can connect to',required=True)
    parser.add_argument('--num_networks', action='store', dest='num_networks', help='number of networks to generate',type=int, required=True)
    parser.add_argument('--network_id', action='store', dest='network_id',help='name of network to use in building dictionary', default='network', required=False)
    parser.add_argument('--use_compendia', action='store', dest='use_compendia', help = 'whether to use compendia ratios to build netwokr', default = 'yes')
    results = parser.parse_args()
    return results

def process_args(results):
    # get logger
    if results.odir is None or not (path.exists(results.odir) and path.isdir(results.odir)):
        log = logging.getLogger(results.network_id)
    else:
        log = helpers.get_logger(results.network_id, f"{results.odir}/{results.network_id}_pruning.log")
    
    #check if output directory exists
    if not (path.exists(results.odir) and path.isdir(results.odir)):
        log.error("Please provide a valid output directory")
        exit()
    # check if experiment file exists and is either csv or tsv file. 
    # Load experiment if valid
    if path.exists(results.network_file) and path.isfile(results.network_file):
        log.info("Loading Network File")
        filetype = results.network_file.split('.')[-1]
        if filetype == 'csv':
            network = pd.read_csv(results.network_file)
        elif filetype == 'tsv':
            network = pd.read_csv(results.network_file, sep = '\t')
        else:
            log.error("Unrecognized network filetype. Please use a csv or tsv file")
            exit()
    else:
        log.error("Please provide a valid network file")
        exit()
    
    # Map accession, peptide, site column if valid
    log.info("Mapping Network Columns")
    map_columns = defaultdict()
    columns = list(network.columns)
    if results.accession in columns:
        map_columns['accession'] = results.accession
    else:
        log.error(f"{results.accession} not found in network columns. Please provide a valid accession column")
        exit() 
    if results.site in columns:
        map_columns['site'] = results.site
    else:
        log.error(f"{results.site} not found in network columns. Please provide a valid site column")
        exit() 
    
    if results.score in columns:
        map_columns['score'] = results.score
    else:
        log.error(f"{results.score} not found in network columns. Please provide a valid score column")
        exit() 
    if results.kinase in columns:
        map_columns['kinase'] = results.kinase
    else:
        log.error(f"{results.kinase} not found in experiment columns. Please provide a valid site column")
        exit()

    use_compendia = helpers.string_to_boolean(results.use_compendia)
    return network, log, map_columns, use_compendia


def run_pruning(network, log, map_columns, use_compendia, phospho_type, kinase_size, site_limit, num_networks, network_id):
    log.info("Running pruning algorithm")
    pruner = Prune(network, log, phospho_type, map_columns)
    if use_compendia:
        log.info("Pruning using compendia ratios")
        pruned_networks = pruner.build_multiple_compendia_networks(kinase_size, site_limit, num_networks, network_id)
    else:
        log.info("Pruning without using compendia")
        pruned_networks = pruner.build_multiple_networks(kinase_size, site_limit, num_networks, network_id)
    return pruned_networks

def save_pruning(network_map, phospho_type, network_id, kinase_size, site_limit, use_compendia, odir):
    log.info("Saving pruning results")
    if not os.path.exists(f"{odir}/INDIVIDUAL_NETWORKS"): 
            os.mkdir(f"{odir}/INDIVIDUAL_NETWORKS") 
    
    if use_compendia:
        suffix = f"{phospho_type}_compendia_kinase_sites:{kinase_size}_limit:{site_limit}"
    else:
        suffix = f"{phospho_type}_kinase_sites:{kinase_size}_limit:{site_limit}"
    for nid, network in network_map.items():
        networkd.to_csv(f"{odir}/INDIVIDUAL_NETWORKS/{nid}_{suffix}.tsv", sep = '\t', index=False)
    pickle.dump( kinact_dict, open( f"{odir}/{network_id}_{suffix}.p", "wb" ) )
    
    
def main():
    results = parse_args()
    network, log, map_columns, use_compendia = process_args(results)
    phospho_type = results.phospho_type
    kinase_size = results.kinase_size
    site_limit = results.site_limit
    num_networks = results.num_networks
    network_id = results.network_id
    odir = results.odir
    log.info("Beginning to build pruning networks")
    pruned_networks = run_pruning(network, log, map_columns, use_compendia, phospho_type, kinase_size, site_limit, num_networks, network_id)
    save_pruning(pruned_networks, phospho_type, network_id, kinase_size, site_limit, use_compendia, odir)
    
if __name__ == "__main__":
    main()



