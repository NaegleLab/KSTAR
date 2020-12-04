import pandas as pd 
import logging
import random 
import numpy as np
import multiprocessing
import itertools

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
        self.compendia = self.compendia[self.compendia['KSTAR_SITE'].str.startswith(self.phospho_type)]
        self.compendia = self.compendia[['KSTAR_ACCESSION', 'KSTAR_SITE', 'KSTAR_NUM_COMPENDIA']]
        self.compendia = self.compendia.groupby(['KSTAR_ACCESSION', 'KSTAR_SITE']).max().reset_index()

        self.network = pd.merge(self.network, self.compendia, how = 'left', left_on = [self.columns['accession'], self.columns['site']], right_on=['KSTAR_ACCESSION', 'KSTAR_SITE'] )
        self.network = self.network[['KSTAR_ACCESSION', 'KSTAR_SITE', self.columns['kinase'], self.columns['score'], 'KSTAR_NUM_COMPENDIA']]
    
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
        pruned_network = pd.DataFrame(columns = network.columns)
        for i in range(kinase_size):
            random.shuffle(self.kinases)
            for kinase in self.kinases:
                kinase_network = network[network[self.columns['kinase']] == kinase]
                
                sample = kinase_network.sample(weights=self.columns['score'])           # Heuristically select kinase-substrate edge
                network.drop(sample.index, inplace = True)                              # Remove edge from remaining network
                pruned_network = pd.concat([pruned_network, sample])                    # add edge to pruned network
                
                # Remove substrate site from network if 
                sample = sample.iloc[0]
                pruned_site = pruned_network['KSTAR_SITE'] == sample['KSTAR_SITE']
                pruned_accession = pruned_network['KSTAR_ACCESSION'] == sample['KSTAR_ACCESSION']
                pruned_match = pruned_network[pruned_site & pruned_accession]
                if len(pruned_match) >= site_limit:
                    network_site = network['KSTAR_SITE'] == sample['KSTAR_SITE']
                    network_accession = network['KSTAR_ACCESSION'] == sample['KSTAR_ACCESSION']
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





