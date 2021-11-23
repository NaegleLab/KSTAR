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
from datetime import datetime
from kstar import helpers, config
import shutil
import random
import string

#%%
#%%

class Pruner:
    """
    Pruning Algorithm used for KSTAR. 
        
    Parameters
    -----------
    network : pandas df
        kinase-site prediction network where there is an accession, site, kinase, and score column
    logger : 
        logger used for pruning
    phospho_type : str
        phospho_type(s) to use when building pruned networks
    columns : dict
        relevant columns in network
    """
    def __init__(self, network, logger, phospho_type = 'Y',):
        

        self.phospho_type = tuple(phospho_type)
        self.logger = logger
        self.network = network[network['site'].str.startswith(self.phospho_type)]
        

        self.compendia = config.HUMAN_REF_COMPENDIA
        self.compendia = self.compendia[self.compendia[config.KSTAR_SITE].str.startswith(self.phospho_type)]
        self.compendia = self.compendia[[config.KSTAR_ACCESSION, config.KSTAR_SITE, 'KSTAR_NUM_COMPENDIA']]
        self.compendia = self.compendia.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE]).max().reset_index()

        self.network = pd.merge(self.network, self.compendia, how = 'inner', left_on = ['substrate_acc', 'site'], right_on=[config.KSTAR_ACCESSION, config.KSTAR_SITE] )
        self.network = self.network.drop(columns=['substrate_acc','site','substratefor_id', 'substrate_name', 'pep'])
        self.network = self.network.set_index([config.KSTAR_ACCESSION, config.KSTAR_SITE])

        self.network = self.network.replace('-',np.nan)
        self.network = self.network.dropna(axis='columns', how ='all')

        self.kinases = list(self.network.columns)
        self.kinases.remove('KSTAR_NUM_COMPENDIA')

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
        
        site_sizes = defaultdict()
        for i in network.index:
            site_sizes[i] = 0

        pruned_network = pd.DataFrame(index = network.index, columns = network.columns, data=np.nan)
        for i in range(kinase_size):
            random.shuffle(self.kinases)
            for kinase in self.kinases:
                sample = network.sample(weights=kinase).iloc[0]
                # return sample
                network.at[sample.name, kinase] = np.nan
                pruned_network.at[sample.name, kinase] = 1
                site_sizes[sample.name] = site_sizes[sample.name] + 1
                if site_sizes[sample.name] >= site_limit:
                    network.at[sample.name,:] = np.nan 
        
        pruned_network = pruned_network.reset_index().melt(id_vars=[config.KSTAR_ACCESSION, config.KSTAR_SITE]).dropna()
        pruned_network = pruned_network.rename(columns={'variable':config.KSTAR_KINASE})
        pruned_network = pruned_network.drop(columns=['value'])   
        return pruned_network
        
    
    def compendia_pruned_network(self, compendia_sizes, site_limit, odir):
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
        name = os.path.join(odir,''.join(random.choices(string.ascii_uppercase + string.digits, k=10)))
        for comp_size, kinase_size in compendia_sizes.items():
            self.logger.info(f"Building pruned network for compendia {comp_size} with {kinase_size} sites per kinase")
            network = self.network[self.network['KSTAR_NUM_COMPENDIA'] == comp_size]
            pruned_network.append(self.build_pruned_network(network, kinase_size, site_limit))
        pruned_network = pd.concat(pruned_network)

        pruned_network.to_csv(f"{name}.tsv", sep = "\t", index=False)
        # return pruned_network
        
    def calculate_compendia_sizes(self, kinase_size):
        """
        Calculates the number of sites per compendia size that a kinase should connect to using same ratios 
        of compendia sizes as found in compendia

        Parameters
        ----------
        kinase_size: int
            number of sites each kinase should connect to
            
        Returns
        -------
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

    def build_multiple_compendia_networks(self, kinase_size, site_limit, num_networks, network_id, odir, PROCESSES = 1):
        """
        Builds multiple compendia-limited networks
        
        Parameters
        ----------
        
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
        
        odir = os.path.join(odir, "work")
        if not path.exists(odir):
            os.mkdir(odir)
        num_existing_networks = len(os.listdir(odir))
        num_networks = num_networks - num_existing_networks
        if num_existing_networks > 0:
            self.logger.info(f"Found {num_existing_networks} networks... building {num_networks} additional networks")
        # MULTIPROCESSING
        if PROCESSES > 1:
            pool = multiprocessing.Pool(processes = PROCESSES)
            comp_sizes = itertools.repeat(compendia_sizes, num_networks)
            limits = itertools.repeat(site_limit, num_networks)
            # names = [os.path.join(odir,f"{network_id}_{i}") for i in range(num_networks)] 
            odirs = itertools.repeat(odir, num_networks)
            iterable = zip(comp_sizes, limits, odirs)

            pool.starmap(self.compendia_pruned_network, iterable)
        
        # SINGLE CORE PROCESSING
        else:
            pruned_networks = []
            for i in range(num_networks):
                self.compendia_pruned_network(compendia_sizes, site_limit, os.path.join(odir,f"{network_id}_{i}") )
                
        
        # pruned_dict = {}
        # for i in range(num_networks):
        #     name = os.path.join(odir,f"{network_id}_{i}")
        #     pruned_dict[f"{network_id}_{i}"] = pd.read_table(f"{name}.tsv")
        

        # return pruned_dict

    # def build_multiple_networks(self, kinase_size, site_limit, num_networks, network_id):
    #     """
    #     Basic Network Generation - only takes into account score when determining sites a kinase 
    #     connects to
    #     """
    #     # MULTIPROCESSING
    #     if config.PROCESSES > 1:
    #         pool = multiprocessing.Pool(processes = config.PROCESSES)
    #         network_iter = itertools.repeat(self.network, num_networks)
    #         size_iter = itertools.repeat(kinase_size, num_networks)
    #         limit_iter = itertools.repeat(site_limit, num_networks)
    #         iterable = zip(network_iter, size_iter, limit_iter)
    #         pruned_networks = pool.starmap(self.build_pruned_network, iterable)
        
    #     # SINGLE CORE PROCESSING
    #     else:
    #         pruned_networks = []
    #         for i in range(num_networks):
    #             net = self.build_pruned_network(self.network, kinase_size, site_limit)
    #             pruned_networks.append(net)

    #     pruned_dict = {}
    #     for i in range(len(pruned_networks)):
    #         pruned_dict[f"{network_id}_{i}"] = pruned_networks[i]
    #     return pruned_dict


"""
********************** RUNNING PRUNING ALGORITHM FROM COMMAND LINE ****************************
Arguments
---------
--network_file          Experiment file location. csv or tsv file           (required)
--output_directory      output directory for results                        (required)
--phospho_type          phospho type (Y, ST, ...)
--kinase_size           number of sites a kinase connects to (required)
--site_limit            upper limit of number of kinases can connect to     (required)
--num_networks          number of networks to generate                      (required)
--network_id            name of network to use in building dictionary
--use_compendia         whether to use compendia ratios to build network  



"""
def parse_args():
    parser = argparse.ArgumentParser(description='Parse Network Pruning Arguments')
    parser.add_argument('--network_file', action='store', dest= 'network_file', help='Experiment file location. csv or tsv file', required=True)
    parser.add_argument('--odir','--output_directory', action = 'store', dest='odir', help = 'output directory for results', required=True)
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


    use_compendia = helpers.string_to_boolean(results.use_compendia)
    return network, log, use_compendia


def run_pruning(network, log, use_compendia, phospho_type, kinase_size, site_limit, num_networks, network_id, odir, PROCESSES = 1):
    """
    Generate pruned networks from a weighted kinase-substrate graph and log run information
    
    Parameters
    ----------
    network: pandas dataframe
        kinase substrate network matrix, with values indicating weight of kinase-substrate relationship
    log: logger
        logger to document the pruning process from start to finish
    use_compendia: string
        whether to use compendia ratios to build network
    phospho_type: string
        phospho type ('Y', 'ST', ...)
    kinase_size: int
        number of sites a kinase connects to
    site_limit: int
        upper limit of number of kinases can connect to
    num_networks: int
        number of networks to generate
    network_id: string
        name of network to use in building dictionary
    odir: string
        output directory for results
        
    Returns
    -------
    pruner: Prune object
        prune object that contains the number of pruned networks indicated by the num_networks paramater
    """
    log.info("Running pruning algorithm")
    pruner = Prune(network, log, phospho_type)
    if use_compendia:
        log.info("Pruning using compendia ratios")
        pruner.build_multiple_compendia_networks(kinase_size, site_limit, num_networks, network_id, odir, PROCESSES = PROCESSES)
    # else:
    #     log.info("Pruning without using compendia")
    #     pruned_networks = pruner.build_multiple_networks(kinase_size, site_limit, num_networks, network_id)
    return pruner

def save_pruning(phospho_type, network_id, kinase_size, site_limit, use_compendia, odir, log):
    """
    Save the pruned networks generated by run_pruning function as a pickle to be loaded by KSTAR
    
    Parameters
    ----------
    phosho_type: string
        type of phosphomodification to networks were generated for (either 'Y' or 'ST')
    network_id: string
        name of network used to build dictionary
    kinase_size: int
        number of sites a kinase connects to
    site_limit: int
        upper limit of number of kinases can connect to
    use_compendia: string 
        whether compendia was used for ratios to build networks
    odir: string
        output directory for results
    log: logger
        logger to document pruning process from start to finish
        
    Returns
    -------
    Nothing
    """
    log.info("Saving pruning results")
    # if not os.path.exists(f"{odir}/INDIVIDUAL_NETWORKS"): 
    #         os.mkdir(f"{odir}/INDIVIDUAL_NETWORKS") 
    
    if use_compendia:
        suffix = f"{phospho_type}_compendia_{kinase_size}_limit_{site_limit}"
    else:
        suffix = f"{phospho_type}_{kinase_size}_limit_{site_limit}"

    # TODO: rename work directory and rename files
    os.rename(os.path.join(odir, "work"), os.path.join(odir, "INDIVIDUAL_NETWORKS"))
    odir_ind = os.path.join(odir, "INDIVIDUAL_NETWORKS")
    # rename files and make network map
    temp_filenames = os.listdir(odir_ind)
    network_map = {}
    for i in range(len(temp_filenames)):
        network_map[f"{network_id}_{i}"] = pd.read_table(os.path.join(odir_ind,temp_filenames[i]))
        os.rename(os.path.join(odir_ind, temp_filenames[i]), os.path.join(odir_ind, f"{network_id}_{i}_{suffix}.tsv"))

    # for nid, network in network_map.items():
    #     network.to_csv(f"{odir}/INDIVIDUAL_NETWORKS/{nid}_{suffix}.tsv", sep = '\t', index=False)
    pickle.dump( network_map, open( f"{odir}/{network_id}_{suffix}.p", "wb" ) )

    
def save_run_information(results, use_compendia, pruner):
    """
    Save information about the generation of networks during run_pruning, 
    including the parameters used for generation. Primarily used when running bash script.
    
    Parameters
    ----------
    results: 
        object that stores all parameters used in the pruning process
    use_compendia: string
        whether compendia was used for ratios to build network
    pruner: Prune object
        output of the run_pruning() function
        
    Returns
    -------
    Nothing
    """
    with open(f"{results.odir}/RUN_INFORMATION.txt", "w") as info_file:
        info_file.write("*************************************************\n")
        info_file.write(f"Pruning Information for {results.network_id}\n")
        info_file.write("*************************************************\n")
        info_file.write(f"Date Run\t\t{datetime.now()}\n")
        info_file.write(f"Network Used\t{results.network_file}\n")
        info_file.write(f"Phospho Type\t{results.phospho_type}\n")
        info_file.write(f"Kinase Size\t\t{results.kinase_size}\n")
        info_file.write(f"Site Limit\t\t{results.site_limit}\n")
        info_file.write(f"# of Networks\t{results.num_networks}\n")
        info_file.write(f"Use Compendia\t{results.use_compendia}\n")
        if use_compendia:
            compendia_sizes = pruner.calculate_compendia_sizes(results.kinase_size)
            for comp, size in compendia_sizes.items():
                info_file.write(f"\tCompendia {comp}\t{size}\n")




def main():
    results = parse_args()
    network, log, use_compendia = process_args(results)
    phospho_type = results.phospho_type
    kinase_size = results.kinase_size
    site_limit = results.site_limit
    num_networks = results.num_networks
    network_id = results.network_id
    odir = results.odir
    log.info("Beginning to build pruning networks")
    pruner = run_pruning(network, log, use_compendia, phospho_type, kinase_size, site_limit, num_networks, network_id, odir)
    save_pruning( phospho_type, network_id, kinase_size, site_limit, use_compendia, odir, log)
    save_run_information(results, use_compendia, pruner)
    
if __name__ == "__main__":
    main()



