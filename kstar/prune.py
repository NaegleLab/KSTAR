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
import random
import string
import warnings
import math
#add import required to generate unique hash for the run_information.txt
import hashlib

#%%
#%%

class Pruner:
    """
    Pruning Algorithm used for KSTAR. 
        
    Parameters
    -----------
    network : pandas df
        weighted kinase-site prediction network where there is an accession, site, kinase, and score column
    logger : 
        logger used for pruning
    phospho_type : str
        phospho_type(s) to use when building pruned networks
    acc_col : str
        the name of the column containing Uniprot Accession IDs for each substrate in the weighted network
    site_col : str
        the name of the column containing the residue type and location of each substrate in the weighted network (Y1268, S44, etc.)
    nonweight_cols : list
        indicates the non-weight containing columns in the network (these will be removed in the final processed network, as they are not needed). If None, will automatically look
            for any non-numeric columns and removes them.
    """
    def __init__(self, network, logger, phospho_type = 'Y',acc_col = 'substrate_acc', site_col = 'site', nonweight_cols = ['substrate_acc','site','substrate_id', 'substrate_name', 'pep']):
        

        self.phospho_type = phospho_type
        self.logger = logger
        self.network = network[network[site_col].str.startswith(tuple(self.phospho_type))]
        
        #load human reference phosphoproteome, obtain the correct phosphosites (Y or ST), and combine duplicate entries
        self.compendia = config.HUMAN_REF_COMPENDIA
        self.compendia = self.compendia[self.compendia[config.KSTAR_SITE].str.startswith(tuple(self.phospho_type))]
        self.compendia = self.compendia[[config.KSTAR_ACCESSION, config.KSTAR_SITE, 'KSTAR_NUM_COMPENDIA']]
        self.compendia = self.compendia.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE]).max().reset_index()
        #merge the human phosphoproteome compendia, which contains study bias information for all phosphosites, with the network
        self.network = pd.merge(self.network, self.compendia, how = 'inner', left_on = [acc_col, site_col], right_on=[config.KSTAR_ACCESSION, config.KSTAR_SITE] )
        
        print(f'There are {self.network.shape[0]} unique substrates found in the weighted network that mapped to KinPred reference phosphoproteome\n')
        
        
        #update index to indicate site associated with each row
        self.network = self.network.set_index([config.KSTAR_ACCESSION, config.KSTAR_SITE])
        #drop unneeded columns (or maybe drop any non number columns)
        if nonweight_cols is not None:
            nonweight_cols = [col for col in nonweight_cols if col in self.network.columns]
            self.network = self.network.drop(columns=nonweight_cols)
            
        #check to make sure all remaining columns contain only numeric values
        non_numeric_cols = self.network.select_dtypes(exclude = np.number).columns
        if len(non_numeric_cols) > 0:
            print('The following non-numeric columns were identified and removed:', ', '.join(non_numeric_cols))
            print('It is assumed that these columns do not contain kinase-substrate weights. If this is incorrect, please make sure weight containing columns only contain numeric values or indicate non-weight containing columns via "nonweight_cols" parameter\n')
            self.network = self.network.select_dtypes(include = np.number)

        self.network = self.network.replace('-',np.nan)
        self.network = self.network.dropna(axis='columns', how ='all')
        
        #check for duplicate entries/indices
        if len(self.network.index.duplicated()) > 0:
            tmp = self.network[self.network.index.duplicated(keep = False)]
            #check to see if duplicate entries contain different information. If they do not, drop one row. If they do, return error.
            if any([i and v for i,v in zip(self.network.index.duplicated(keep = False), ~self.network.duplicated(keep = False).values)]):
                raise ValueError('Duplicate entries found in network that contain conflicting weights, please fix before running pruning')

            else:
                print('Duplicate entries found in network, removing repeat rows\n')
                self.network = self.network[~self.network.index.duplicated()]

        self.kinases = list(self.network.columns)
        self.kinases.remove('KSTAR_NUM_COMPENDIA')
        print(f'There are {len(self.kinases)} {phospho_type} kinases found in the weighted network')
        print(f'There are {self.network.shape[0]} unique substrates found in the weighted network that mapped to KinPred reference phosphoproteome\n')

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
        num_edges = {}
        pruned_network = pd.DataFrame(index = network.index, columns = network.columns, data=np.nan)
        for i in range(kinase_size):
            random.shuffle(self.kinases)
            for kinase in self.kinases:
                if network[kinase].isna().all():
                    if kinase not in num_edges:
                        num_edges[kinase] = i+1
                        logger.info(f'{kinase} has no more available edges. Total edge count for this network = {num_edges[kinase]}')
                else:
                    sample = network.sample(weights=kinase).iloc[0]
                    network.at[sample.name, kinase] = np.nan
                    pruned_network.at[sample.name, kinase] = 1
                    site_sizes[sample.name] = site_sizes[sample.name] + 1
                    if site_sizes[sample.name] >= site_limit:
                        network.loc[sample.name,:] = np.nan 


        
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
        self.kinase_size = kinase_size
        self.site_limit = site_limit
        self.num_networks = num_networks
        self.network_id = network_id
        self.network_dir = odir
        self.use_compendia = True
        
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
            odirs = itertools.repeat(odir, num_networks)
            iterable = zip(comp_sizes, limits, odirs)

            pool.starmap(self.compendia_pruned_network, iterable)
        
        # SINGLE CORE PROCESSING
        else:
            pruned_networks = []
            for i in range(num_networks):
                self.compendia_pruned_network(compendia_sizes, site_limit, odir )
        
  

    def build_multiple_networks(self, kinase_size, site_limit, num_networks, network_id, odir, PROCESSES = 1):
        """
        Basic Network Generation - only takes into account score when determining sites a kinase 
        connects to
        """
        self.kinase_size = kinase_size
        self.site_limit = site_limit
        self.num_networks = num_networks
        self.network_id = network_id
        self.network_dir = odir
        self.use_compendia = False
        
        odir = os.path.join(odir, "work")
        if not path.exists(odir):
            os.mkdir(odir)
        # MULTIPROCESSING
        if PROCESSES > 1:
            pool = multiprocessing.Pool(processes = PROCESSES)
            network_iter = itertools.repeat(self.network, num_networks)
            size_iter = itertools.repeat(kinase_size, num_networks)
            limit_iter = itertools.repeat(site_limit, num_networks)
            iterable = zip(network_iter, size_iter, limit_iter)
            pruned_networks = pool.starmap(self.build_pruned_network, iterable)
            #save each network
            for i in range(len(pruned_networks)):
                net = pruned_networks[i]
                name = os.path.join(odir,''.join(random.choices(string.ascii_uppercase + string.digits, k=10)))
                net.to_csv(f"{name}.tsv", sep = "\t", index=False) 
        # SINGLE CORE PROCESSING
        else:
            pruned_networks = []
            for i in range(num_networks):
                name = os.path.join(odir,''.join(random.choices(string.ascii_uppercase + string.digits, k=10)))
                net = self.build_pruned_network(self.network, kinase_size, site_limit)
                net.to_csv(f"{name}.tsv", sep = "\t", index=False)
                
                
    def getMaximumKinaseSize(self, site_limit):
        """
        Given a network and site_limit (maximum number of kinases a phosphorylation site can provide evidence to), will calculate the theoretical maximum number of connections each kinase can have (kinase_size parameter)
        
        Theoretical maximum exists when each substrate hits the maximum site_limit
        
        Parameters
        ----------
        site_limit: int
            Parameter used in pruning: indicates the maximum number of kinases a phosphorylation site can be connected to in the final pruned network
        
        Returns
        -------
        theoretical_max_ksize: int
            largest possible value that 'kinase_size' parameter can have without throwing any errors
        """
        
        num_substrates = self.network.shape[0]
        
        theoretical_max_ksize = (num_substrates*site_limit)/len(self.kinases)
        return theoretical_max_ksize
        
        
    def getRecommendedKinaseSize(self, site_limit):
        """
        Given a network and site_limit (maximum number of kinases a phosphorylation site can provide evidence to), will calculate the theoretical maximum number of connections each kinase can have (kinase_size parameter) and recommend a range of values for kinase_size
        
        Theoretical maximum exists when each substrate hits the maximum site_limit
        
        Parameters
        ----------
        site_limit: int
            Parameter used in pruning: indicates the maximum number of kinases a phosphorylation site can be connected to in the final pruned network
        
        Returns
        -------
        Nothing, prints theoretical maximum of kinase size and the recommened values for the parameter given the site_limit
        """
        #get the total number of substrates in the network
        theoretical_max_ksize = self.getMaximumKinaseSize(site_limit)
        print(f'Theoretical maximum value of kinase_size parameter is {theoretical_max_ksize}, given site_limit of {site_limit}')
        print(f'Recommended value of kinase size parameter is between 20-40% of theoretical maximum: {theoretical_max_ksize*0.2} - {theoretical_max_ksize*0.4}')
        
    def checkParameters(self, kinase_size, site_limit):
        """
        Given the site_limit and kinase_size parameters to be used during pruning, raise errors if not feasible, and raise warnings if value is higher than we would recommend (>40% of the maximum kinase_size value)
        
        Parameters
        ----------
        kinase_size: int
            Parameter used in pruning: indicates the number of substrates each kinase will be connected to
        site_limit: int
            Parameter used in pruning: indicates the maximum number of kinases a phosphorylation site can be connected to in the final pruned network
            
        Returns
        -------
        Nothing, will only raise errors/warnings if parameters are not feasible
        """
        
        theoretical_max_ksize = self.getMaximumKinaseSize(site_limit)
        
        if kinase_size > theoretical_max_ksize:
            raise ValueError(f"Value of 'kinase_size' is not feasible. Please decrease 'kinase_size' to {math.floor(theoretical_max_ksize)} or increase the value of site_limit")
        elif kinase_size > theoretical_max_ksize*0.4:
            warnings.warn("Value of 'kinase_size' may be higher than desired. Check the log file for potential errors after completion")
            
        
                
    def save_networks(self, network_file):
        """
        Save the pruned networks generated by the 'build_multiple_networks' or 'build_multiple_compendia_networks' as a pickle to be loaded by KSTAR
        
        """
        self.logger.info("Saving pruning results")
        self.network_file_used = network_file
        print(network_file)
        if self.use_compendia:
            suffix = f"{self.phospho_type}_compendia_{self.kinase_size}_limit_{self.site_limit}"
        else:
            suffix = f"{self.phospho_type}_{self.kinase_size}_limit_{self.site_limit}"

        # TODO: rename work directory and rename files
        os.rename(os.path.join(self.network_dir, "work"), os.path.join(self.network_dir, "INDIVIDUAL_NETWORKS"))
        odir_ind = os.path.join(self.network_dir, "INDIVIDUAL_NETWORKS")
        # rename files and make network map
        temp_filenames = os.listdir(odir_ind)
        network_map = {}
        for i in range(len(temp_filenames)):
            network_map[f"{self.network_id}_{i}"] = pd.read_table(os.path.join(odir_ind,temp_filenames[i]))
            os.rename(os.path.join(odir_ind, temp_filenames[i]), os.path.join(odir_ind, f"{self.network_id}_{i}_{suffix}.tsv"))

        pickle.dump(network_map, open( f"{self.network_dir}/{self.network_id}_{suffix}.p", "wb" ))

        # Generate a unique hash
        self.unique_id = hashlib.sha256(f"{self.network_id}_{datetime.now()}".encode()).hexdigest()

        #Call save_run_information method with the unique hash
        self.save_run_information()
        
        

        
    def save_run_information(self):
        """
        Save information about the generation of networks during run_pruning, 
        including the parameters used for generation. Primarily used when running bash script.
        
        """
        with open(f"{self.network_dir}/RUN_INFORMATION.txt", "w") as info_file:
            info_file.write("*************************************************\n")
            info_file.write(f"Pruning Information for {self.network_id}\n")
            info_file.write("*************************************************\n")
            info_file.write(f"Unique ID: {self.unique_id}\n")
            info_file.write(f"Date Run\t\t{datetime.now()}\n")
            info_file.write(f"Network Used\t{self.network_file_used}\n")
            info_file.write(f"Phospho Type\t{self.phospho_type}\n")
            info_file.write(f"Kinase Size\t\t{self.kinase_size}\n")
            info_file.write(f"Site Limit\t\t{self.site_limit}\n")
            info_file.write(f"# of Networks\t{self.num_networks}\n")
            info_file.write(f"Use Compendia\t{self.use_compendia}\n")
            if self.use_compendia:
                compendia_sizes = self.calculate_compendia_sizes(self.kinase_size)
                for comp, size in compendia_sizes.items():
                    info_file.write(f"\tCompendia {comp}\t{size}\n")


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
    parser.add_argument('--PROCESSES', action='store', dest='PROCESSES', help = 'number of processes to run in parallel', default = 1)
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


def run_pruning(network, log, phospho_type, kinase_size, site_limit, num_networks, network_id, odir, use_compendia = True,
                acc_col = 'substrate_acc', site_col = 'site', netcols_todrop = ['substrate_acc','site','substrate_id', 'substrate_name', 'pep'], PROCESSES = 1):
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
    pruner = Pruner(network, log, phospho_type, acc_col = acc_col, site_col = site_col, nonweight_cols = netcols_todrop)
    pruner.checkParameters(kinase_size, site_limit)
    if use_compendia:
        log.info("Pruning using compendia ratios")
        pruner.build_multiple_compendia_networks(kinase_size, site_limit, num_networks, network_id, odir, PROCESSES = PROCESSES)
    else:
        log.info("Pruning without using compendia")
        pruned_networks = pruner.build_multiple_networks(kinase_size, site_limit, num_networks, network_id, odir, PROCESSES = PROCESSES)
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
    pickle.dump( network_map, open( f"{odir}/{network_id}_{suffix}.p", "wb" ))
    # Generate a unique hash
    unique_id = hashlib.sha256(f"{network_id}_{datetime.now()}".encode()).hexdigest()


    
def save_run_information(results, use_compendia, pruner, unique_id):
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
        info_file.write(f"Unique ID:{unique_id}\n")
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
    try:
        PROCESSES = int(results.PROCESSES)
    except:
        log.info('PROCESSES parameter must be an integer')
        exit()
        
    log.info("Beginning to build pruning networks")
    pruner = run_pruning(network, log, phospho_type, kinase_size, site_limit, num_networks, network_id, odir)
    save_pruning(phospho_type, network_id, kinase_size, site_limit, use_compendia, odir, log)
    save_run_information(results, use_compendia, pruner)

if __name__ == "__main__":
    main()



