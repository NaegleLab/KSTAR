import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt
import seaborn as sns

from kstar import config

def numUniqueSubstrates(networks, acc_col = 'KSTAR_ACCESSION', site_col = 'KSTAR_SITE'):
    """
    Given a KSTAR network(s), return the number of unique substrates within the network (across all kinases). If a dictionary of multiple pruned networks is provided, will calculate the total number of unique substrates across ALL networks.
    
    Parameters
    ----------
    network: pandas dataframe or dict of pandas dataframes
        pruned KSTAR network, or dictionary containing multiple pruned networks
    acc_col: str
        name of column in network dataframe which indicates UniProt ID of substrates
    site_col: str
        name of column in network dataframe which indicates residue and site number (i.e. Y1197)
    
    Returns
    -------
    Number of unique substrates within network(s)
    """
    if isinstance(networks, dict):
        net_substrates = []
        for nname in networks.keys():
            net = networks[nname]
            net_substrates = net_substrates + list(net[acc_col]+'_'+net[site_col])
        net_substrates = np.unique(net_substrates)
    else:
        net_substrates = np.unique(networks[acc_col]+'_'+networks[site_col])
    return len(net_substrates)
    
def averageUniqueSubstrates_KSTAR(mod_types = ['Y','ST']):
    """
    Calculate the average number of unique substrates covered by each KSTAR pruned network
    
    Parameters
    ----------
    mod_types: list
        list containing which networks to calculate average for. Either ['Y'], ['ST'], or ['Y','ST']
    
    Returns
    -------
    averageSub: dict
        indicates the average number of substrates across all pruned networks for indicated modification types
        
    """
    averageSub = {}
    for mod in mod_types:
        numSubstrates = []
        if mod == 'Y':
            networks = pickle.load(open(config.NETWORK_Y_PICKLE, 'rb'))
        else:
            networks = pickle.load(open(config.NETWORK_ST_PICKLE, 'rb'))
        
        for nname in networks.keys():
            net = networks[nname]
            numSubstrates.append(numUniqueSubstrates(net))
        
        averageSub[mod] = np.mean(numSubstrates)
        
    return averageSub
    
def getStudyBiasDistribution_InPhosphoproteome(mod_type = 'Y', ax = None, figsize = (4,3), return_dist = False):
    """
    Plot the distribution of study bias across the reference phosphoproteome
    
    Parameters
    ----------
    mod_type: str
        indicates which modification type, tyrosine ('Y') or serine/threonine ('ST'), you would like to plot. Default is 'Y'
    ax: matplotlib axes object
        axis to plot the distribution on. If none, will create subplot
    figsize: tuple
        size of matplotlib figure. Default is (4,3)
    return_dist: bool
        whether you would like to also return the distribution values. Default is False.
        
    Returns
    -------
    Histogram plotting the distribution of study bias found in overall phosphoproteome, as defined by the number of compendia a phosphorylation site is recorded in. If return_dist = True, will also return a series object containing the same data as the histogram.
    """
    #get compendia info specific to modification type
    if mod_type == 'Y':
        comp = config.HUMAN_REF_COMPENDIA[config.HUMAN_REF_COMPENDIA['type'] == 'Y']
    elif mod_type == 'ST':
        comp = config.HUMAN_REF_COMPENDIA[config.HUMAN_REF_COMPENDIA['type'] != 'Y']
    else:
        raise ValueError("mod_type must be either 'Y' or 'ST'")
    #make histogram, return corresponding data if requested
    if ax is None:
        fig, ax = plt.subplots(figsize = figsize)
    sns.histplot(data = comp, x = 'KSTAR_NUM_COMPENDIA', bins = 5, ax = ax)
    ax.set_xlabel('Number of Compendia')
    ax.set_ylabel('Number of Sites')
    if return_dist:
        return comp.groupby('KSTAR_NUM_COMPENDIA').size()
    
def getStudyBiasDistribution_InExperiment(binary_experiment, ax = None, figsize = (4,3), return_dist = False):
    """
    Plot the distribution of study bias within a single phosphoproteomic experiment
    
    Parameters
    ----------
    mapped_experiment: pandas dataframe
        phosphoproteomic experiment that has been mapped by KSTAR (contains 'KSTAR_SITE','KSTAR_ACCESSION', and 'KSTAR_NUM_COMPENDIA' columns)
    ax: matplotlib axes object
        axis to plot the distribution on. If none, will create subplot
    figsize: tuple
        size of matplotlib figure. Default is (4,3)
    return_dist: bool
        whether you would like to also return the distribution values. Default is False.
        
    Returns
    -------
    Histogram plotting the distribution of study bias found in the provided experiment, as defined by the number of compendia a phosphorylation site is recorded in. If return_dist = True, will also return a series object containing the same data as the histogram.
    """
    #aggregate sites into unique rows (no repeat sites)
    plt_data = binary_experiment.groupby(['KSTAR_ACCESSION','KSTAR_SITE']).mean()
    #create histogram, return corresponding data if requested
    if ax is None:
        fig, ax = plt.subplots(figsize = figsize)
    sns.histplot(data = plt_data, x = 'KSTAR_NUM_COMPENDIA', bins = 5, ax = ax)
    ax.set_xlabel('Number of Compendia')
    ax.set_ylabel('Number of Sites')
    if return_dist:
        return plt_data.groupby('KSTAR_NUM_COMPENDIA').size()
        
def getStudyBiasDistribution_InSample(binary_experiment, data_column, ax = None, figsize = (4,3), return_dist = False):
    """
    Plot the distribution of study bias within a single phosphoproteomic experiment
    
    Parameters
    ----------
    mapped_experiment: pandas dataframe
        phosphoproteomic experiment that has been mapped by KSTAR (contains 'KSTAR_SITE','KSTAR_ACCESSION', and 'KSTAR_NUM_COMPENDIA' columns)
    ax: matplotlib axes object
        axis to plot the distribution on. If none, will create subplot
    figsize: tuple
        size of matplotlib figure. Default is (4,3)
    return_dist: bool
        whether you would like to also return the distribution values. Default is False.
        
    Returns
    -------
    Histogram plotting the distribution of study bias found in the provided experiment, as defined by the number of compendia a phosphorylation site is recorded in. If return_dist = True, will also return a series object containing the same data as the histogram.
    """
    #aggregate sites into unique rows (no repeat sites)
    plt_data = binary_experiment[binary_experiment[data_column] == 1]
    plt_data = plt_data.groupby(['KSTAR_ACCESSION','KSTAR_SITE']).mean()
    #create histogram, return corresponding data if requested
    if ax is None:
        fig, ax = plt.subplots(figsize = figsize)
    sns.histplot(data = plt_data, x = 'KSTAR_NUM_COMPENDIA', bins = 5, ax = ax)
    ax.set_xlabel('Number of Compendia')
    ax.set_ylabel('Number of Sites')
    if return_dist:
        return plt_data.groupby('KSTAR_NUM_COMPENDIA').size()
    
    
def experimentCoverage(experiment, networks, mod = 'Y', exp_cols = ['KSTAR_ACCESSION', 'KSTAR_SITE'], net_cols = ['KSTAR_ACCESSION', 'KSTAR_SITE']):
    """
    Given an experiment, determine how many of the sites observed in the experiment can be captured by a kinase-substrate network (function was designed for KSTAR pruned networks, but should work with any ks-network that indicates UniProt ID and site number)
    
    Parameters
    ----------
    experiment: pandas dataframe
        phosphoproteomic experiment, ideally that has been mapped to KinPred by KSTAR already
    network: pandas dataframe
        binarized kinase-substrate network (unweighted), ideally having been mapped to KinPred/KSTAR already
    exp_cols: list
        list indicating the columns in experiment dataframe that contain uniprot id and site number
    net_cols: list
        list indicating the columns in network dataframe that contain the uniprot id and site number
        
    Returns
    -------
    fraction_of_sites_covered: dict
        indicates the fraction of phosphorylation sites observed in experiment that are also found within the kinase-substrate network, for each modification type (tyrosine, serine/threonine).
    """
    #extract unique substrates found in network(s)
    if isinstance(networks, dict):
        net_substrates = []
        for nname in networks.keys():
            net = networks[nname]
            net_substrates = net_substrates + list(net[net_cols[0]]+'_'+net[net_cols[1]])
        net_substrates = np.unique(net_substrates)
    else:
        net_substrates = np.unique(networks[net_cols[0]]+'_'+networks[net_cols[1]])
    #extract sites found in experiment
    exp_substrates = np.unique(experiment[exp_cols[0]]+'_'+experiment[exp_cols[1]])
    #restrict analysis to the correct 
    if mod == 'ST':
        mod_sub = [sub for sub in exp_substrates if '_S' in sub or '_T' in sub]
        total_substrates = len(mod_sub)
    
    else:
        mod_sub = [sub for sub in exp_substrates if '_Y' in sub]
        total_substrates = len(mod_sub)
        
    
    if total_substrates > 0:
        overlap = len(set(mod_sub).intersection(set(net_substrates)))
        fraction_of_sites_covered = overlap/total_substrates
        return fraction_of_sites_covered
    else:
        print(f"No phosphorylation sites of type '{mod}' found in experiment")
        return np.nan
    return fraction_of_sites_covered
            
def sampleCoverage(binary_experiment, data_col, networks, mod = 'Y', exp_cols = ['KSTAR_ACCESSION', 'KSTAR_SITE'], net_cols = ['KSTAR_ACCESSION', 'KSTAR_SITE']):
    """
    Given a sample within an experiment, determine how many of the sites observed in the experiment can be captured by KSTAR pruned networks. Essentially the same as experimentCoverage(), but restricts experiment sites to those used as evidence for a given sample
    
        Parameters
    ----------
    binary_experiment: pandas dataframe
        binarized phosphoproteomic experiment, with each 1 indicating that site was observed in sample. Ideally has been mapped to KinPred by KSTAR already
    data_col: str
        column name of the sample of interest
    network: pandas dataframe
        binarized kinase-substrate network (unweighted), ideally having been mapped to KinPred/KSTAR already
    exp_cols: list
        list indicating the columns in experiment dataframe that contain uniprot id and site number
    net_cols: list
        list indicating the columns in network dataframe that contain the uniprot id and site number
        
    Returns
    -------
    fraction_of_sites_covered: dict
        indicates the fraction of phosphorylation sites observed in sample that are also found within the kinase-substrate network, for each modification type (tyrosine, serine/threonine).
    """
    if isinstance(networks, dict):
        net_substrates = []
        for nname in networks.keys():
            net = networks[nname]
            net_substrates = net_substrates + list(net[net_cols[0]]+'_'+net[net_cols[1]])
        net_substrates = np.unique(net_substrates)
    else:
        net_substrates = np.unique(networks[net_cols[0]]+'_'+networks[net_cols[1]])
    binary_experiment = binary_experiment[binary_experiment[data_col] == 1]
    exp_substrates = np.unique(binary_experiment['KSTAR_ACCESSION']+'_'+binary_experiment['KSTAR_SITE'])
    if mod == 'ST':
        mod_sub = [sub for sub in exp_substrates if '_S' in sub or '_T' in sub]
        total_substrates = len(mod_sub)
    
    else:
        mod_sub = [sub for sub in exp_substrates if '_Y' in sub]
        total_substrates = len(mod_sub)
        
    
    if total_substrates > 0:
        overlap = len(set(mod_sub).intersection(set(net_substrates)))
        fraction_of_sites_covered = overlap/total_substrates
        return fraction_of_sites_covered
    else:
        print(f"No phosphorylation sites of type '{mod}' found in experiment")
        return np.nan
    return fraction_of_sites_covered

    