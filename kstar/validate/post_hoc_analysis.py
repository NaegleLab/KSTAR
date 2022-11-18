import pandas as pd 
import numpy as np
from collections import defaultdict
from kstar import config

def getSubstrateInfluence(networks, kinase, substrate_subset = None):
    """
    Given the pruned networks and kinase of interest, return the number of networks each substrate is connected to that kinase in (the 'substrate influence' on that kinase's activity prediction). If subset of substrates is provided, will only do this for the given subset
    
    Parameters
    ----------
    networks: dict
        dictionary storing all pruned networks used in activity calculation
    kinase: str
        name of the kinase of interest: should match the name found in provided networks
    substrate_subset: list
        subset of substrates to analyze, indicated by '{KSTAR_ACCESSION}_{KSTAR_SITE}'. If none, will return a series containing info on all substrates with at least one connection to the given kinase
    
    Returns
    -------
    Pandas series indicating the number of networks each substrate is connected the indicated kinase, sorted from the most connections (highest influence) to the least (lowest influence). Sites with no connection will not be included.
    """
    num_networks = defaultdict(int)
    for nname in networks:
        net = networks[nname]
        net_trim = net[net['KSTAR_KINASE'] == kinase]
        net_trim['KSTAR_SUBSTRATE'] = net_trim['KSTAR_ACCESSION']+'_'+net_trim['KSTAR_SITE']
        
        #if subset of substrates provided, reduce to this size
        if substrate_subset is not None:
            keep = [True if sub in substrate_subset else False for sub in net_trim['KSTAR_SUBSTRATE'].values]
            net_trim = net_trim[keep]
        
        #loop through each substrate in network and add to num_network dict if present
        for sub in net_trim['KSTAR_SUBSTRATE']:
            num_networks[sub] += 1
        
    return pd.Series(num_networks, name = f'Number of Networks connected to {kinase}').sort_values(ascending = False)
    
def getSubstrateInfluence_inExperiment(networks, experiment, threshold, kinase, data_cols = None, greater = True):
    if data_cols is None
        data_cols = [col for col in experiment.columns if 'data:' in col]
    
    experiment['KSTAR_SUBSTRATE'] = experiment['KSTAR_ACCESSION'] + '_' + experiment['KSTAR_SITE']
    experiment_influence = {}
    for col in data_cols:
        if greater:
            trim_experiment = experiment[experiment[col] >= threshold]
        else:
            trim_experiment = experiment[experiment[col] <= threshold]
        substrate_subset =  trim_experiment['KSTAR_SUBSTRATE'].values
        experiment_influence[col] = getSubstrateInfluence(networks, kinase, substrate_subset = substrate_subset)
    return experiment_influence