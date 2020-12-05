import pandas as pd
import numpy as np
import networkx as nx
import copy
import seaborn as sns
import matplotlib.pyplot as plt

    
def plot_kinase_heatmap(heatmap, use_mask = True, annotate = True):
    """
    Plots Kinase network heatmap
    
    Parameters
    ----------
    heatmap : pandas dataframe
        Network Heatmap to plot
    use_mask : bool
        If true a mask is applied to the heatmap 
    annotate : bool
        If true then numbers are annotated into each heatmap square
    """
    mask = np.zeros_like(heatmap)
    if use_mask:
        mask[np.triu_indices_from(mask)] = True
    with sns.axes_style("white"):
        sns.heatmap(heatmap, mask=mask, square=True, annot=annotate)


def kinase_network_similarity(network1, network2, kinase_column, substrate_column):
    """
    Compares two networks for shared substrates between kinases
    Similarity based on Jaccard Index. Returned dataframes set index to be network1 kinases
    and columns to be network2 kinases

    Parameters
    ----------
    network1 : pandas dataframe
        Network 1 to compare
    network2 : pandas dataframe
        Network 2 to compare

    Returns
    --------
    heatmap : pandas dataframe
        Number of substrates that overlap between kinases of two networks
    normalized : pandas dataframe
        Normalized mutual information into Jaccard Index. 
        size of intersection of two kinase networks / size of union of two kinase networks.
    heatlist : pandas dataframe
        intersction of kinase networks
    """
    n1_kinases = list(network1[kinase_column].unique())
    n2_kinases = list(network2[kinase_column].unique())
    n1_num_kinases=len(n1_kinases)
    n2_num_kinases=len(n2_kinases)

    n1_kinasemap = {}
    n2_kinasemap = {}

    for kin in n1_kinases:
        n1_kinasemap[kin] = set(network1[network1[kinase_column] == kin][substrate_column])
    for kin in n2_kinases:
        n2_kinasemap[kin] = set(network2[network2[kinase_column] == kin][substrate_column])

    heatlist = [[set() for i in range(n2_num_kinases)] for j in range(n1_num_kinases)]                                                           
    heatmap = np.zeros((n1_num_kinases, n2_num_kinases))
    normalized = np.zeros((n1_num_kinases, n2_num_kinases))

    for i in range(n1_num_kinases):
        n1_kin = n1_kinases[i]
        for j in range(n2_num_kinases):
            n2_kin = n2_kinases[j]
            heatlist[i][j] = n1_kinasemap[n1_kin].intersection(n2_kinasemap[n2_kin])
            heatmap[i][j] = len(heatlist[i][j])
            normalized[i][j] = len(heatlist[i][j]) / len(n1_kinasemap[n1_kin].union(n2_kinasemap[n2_kin]))
    
    heatmap = pd.DataFrame(heatmap, index = n1_kinases, columns = n2_kinases)
    normalized = pd.DataFrame(normalized, index = n1_kinases, columns = n2_kinases)
    heatlist = pd.DataFrame(heatlist, index = n1_kinases, columns = n2_kinases)

    return heatmap, normalized, heatlist
    

def kinase_mutual_information(network, kinase_column, substrate_column):
    """
    Finds mutual information shared between kinases based on the substrate phosphorylated
    Mutual Information is defined as the intersection substrates between two kinases
    A substrate is defined as the substrate accession and site, i.e. P54760_Y596.
    Normalization is performed by comparing intersection of kinases vs union of the two kinases
    This the the Jaccard Index. Jaccard Distance can be calcualted by taking 1 - JI

    Parameters
    ----------
    network : pandas dataframe
        The network to analyze for mutual kinase information
    kinase_column : str
        Column in network that contiains kinase information
    substrate_column : str
        Column in network that contains substrate information
    
    Returns
    --------
    heatmap : pandas dataframe
        Number of substrates that overlap between kinases
    normalized : pandas dataframe
        Normalized mutual information into Jaccard Index. 
        size of intersection of two kinase networks / size of union of two kinase networks.
    heatlist : pandas dataframe
        intersction of kinase networks
    """
    kinases = list(network[kinase_column].unique())
    num_kinases = len(kinases)
    
    heatlist = [[set() for i in range(num_kinases)] for j in range(num_kinases)]                                                           
    heatmap = np.zeros((num_kinases, num_kinases))
    normalized = np.zeros((num_kinases, num_kinases))
    
    
    for i in range(num_kinases):                                                                                    # Get kinase-substrate network of each kinase
        substrates = set(network[network[kinase_column] == kinases[i]][substrate_column])
        heatlist[i][i] = substrates
    
    for i in range(num_kinases):                                                                                    # Iterate through each row in heatmaps
        for j in range(num_kinases):                                                                                # Iterate through each column in heatmaps
            heatlist[i][j] = heatlist[i][i].intersection(heatlist[j][j])                                            # Find intersection of kinase networks
            heatmap[i][j] = len(heatlist[i][j])                                                                     # Get size of intersection
            normalized[i][j] = len(heatlist[i][j]) / len(heatlist[i][i].union(heatlist[j][j]))                      # Normalize via size of the union of two kinase networks
    
    heatmap = pd.DataFrame(heatmap, index = kinases, columns = kinases)
    normalized = pd.DataFrame(normalized, index = kinases, columns = kinases)
    heatlist = pd.DataFrame(heatlist, index = kinases, columns = kinases)
                                                    
    return heatmap, normalized, heatlist


