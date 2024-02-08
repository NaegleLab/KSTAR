import pandas as pd 
import numpy as np
from collections import defaultdict


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
        net_trim = net[net['KSTAR_KINASE'] == kinase].copy()
        net_trim['KSTAR_SUBSTRATE'] = net_trim['KSTAR_ACCESSION']+'_'+net_trim['KSTAR_SITE']
        
        #if subset of substrates provided, reduce to this size
        if substrate_subset is not None:
            keep = [True if sub in substrate_subset else False for sub in net_trim['KSTAR_SUBSTRATE'].values]
            net_trim = net_trim[keep]
        
        #loop through each substrate in network and add to num_network dict if present
        for sub in net_trim['KSTAR_SUBSTRATE']:
            num_networks[sub] += 1
        
    return pd.Series(num_networks, dtype = float, name = f'Number of Networks connected to {kinase}').sort_values(ascending = False)
    
def getSubstrateInfluence_inExperiment(networks, binary_evidence, kinase, data_cols = None):
    """
    Given the binary evidence used for activity prediction, identify which sites are found across the most networks for a given kinase and each sample. 
    
    Parameters
    ----------
    networks: dictionary
        dictionary containing all 50 pruned networks used for activity prediction
    binary_evidence: pandas dataframe
        binarized dataset (using the same threshold/criteria as the one used for activity prediction)
    kinase: str
        name of the kinase to probe
    data_cols: list or None
        name of the data columns in binary_evidence to probe. If None, will analyze all columns with 'data:' at the start of the column name.
    """
    experiment = binary_evidence.copy()
    if data_cols is None:
        data_cols = [col for col in experiment.columns if 'data:' in col]
    
    experiment['KSTAR_SUBSTRATE'] = experiment['KSTAR_ACCESSION'] + '_' + experiment['KSTAR_SITE']
    experiment_influence = {}
    for col in data_cols:
        substrate_subset =  experiment.loc[experiment[col] == 1,'KSTAR_SUBSTRATE'].values
        experiment_influence[col] = getSubstrateInfluence(networks, kinase, substrate_subset = substrate_subset)
    return experiment_influence
    
    
    
"""
def plotMannWhitney(real_enrichment, random_enrichment, kinase, data_col, mann_whitney_activities = None, mann_whitney_fpr = None):
    #transform enrichment values and plot
    transformed_real = -np.log10(real_enrichment.loc[(real_enrichment['kinase'] == kinase) & (real_enrichment['data'] == data_col), 'kinase_activity'].values)
    plt.hist(transformed_real, alpha = 0.5, bins = 50, label = 'Real')
    transformed_random = -np.log10(random_enrichment['kinase'] == kinase) & (random_enrichment['data'] == data_col) 'kinase_activity'].values)
    plt.hist(transformed_random, alpha = 0.5, bins = 50, label = 'Random')
    plt.legend()
    
    #annotate plot with mann whitney results if provided
    activity = -np.log10(mann_whitney_activities.loc[kinase, data_col])
    plt.annotate(
"""