import pandas as pd
import numpy as np
import seaborn as sns

    
def plot_kinase_heatmap(heatmap, use_mask = True, annotate = False):
    """
    Plots Kinase network heatmap
    
    Parameters
    ----------
    heatmap : pandas dataframe
        Network Heatmap to plot (must be square matrix)
    info_type: str
        Indicates what type of informatin is included in heatmap variable. Default is mutual information, equivalent to the normalized matrix obtained from kinase_mutual_information function
    use_mask : bool
        If true a mask is applied to the heatmap 
    annotate : bool
        If true then numbers are annotated into each heatmap square
    """
    mask = np.zeros_like(heatmap)
    if use_mask:
        mask[np.triu_indices_from(mask)] = True
    with sns.axes_style("white"):
        sns.heatmap(heatmap, mask=mask, square=True, annot=annotate, xticklabels=1, yticklabels=1)
    


def kinase_mutual_information(network, kinase_column = 'KSTAR_KINASE', accession_column = 'KSTAR_ACCESSION', site_column ='KSTAR_SITE', substrate_list=None):
    """
    Finds mutual information shared between kinases based on the substrate phosphorylated
    Mutual Information is defined as the intersection substrates between two kinases
    A substrate is defined as the substrate accession and site, i.e. P54760_Y596.
    Normalization is performed by comparing intersection of kinases vs union of the two kinases
    This the the Jaccard Index. Jaccard Distance can be calcualted by taking 1 - JI

    Parameters
    ----------
    network : pandas dataframe or dictionary of pandas dataframe
        The network to analyze for mutual kinase information. Can send a dictionary of multiple pandas dataframes and this will average the MI across all networks in dictionary
    kinase_column : str
        Column in network that contiains kinase information
    substrate_column : str
        Column in network that contains substrate information
    substrate_list : list
        Optional and default is no subset list to use. You can calculate the MI within network(s) for only the evidence given in a substrate_evidence_list (must matche substrate_column of network passed in)
    
    Returns
    --------
    heatmap : pandas dataframe
        Number of substrates that overlap between kinases
    normalized : pandas dataframe
        Normalized mutual information into Jaccard Index. 
        size of intersection of two kinase networks / size of union of two kinase networks.
    heatlist or heatdict: list or dictionary of lists
        intersection of kinase networks. If a single network it is a list. If multiple networks it is a dict of lists with keys the same as the network name

    """
    AVG = 0
    if isinstance(network, dict):
        heatlist = {}
        AVG = 1
        networkNames = list(network.keys())
        #get a network to start and then march through the rest, averaging
        net = network[networkNames[0]]
        #create substrate columns
        net['Substrate'] = net[accession_column]+'_'+net[site_column]
        if substrate_list:
            net = net[net['Substrate'].isin(substrate_list)]
            if net.empty:
                raise ValueError("Mapping of substrate_list to the network failed")
                return -1, -1, -1
        heatmap_all, normalized_all, heatlist[networkNames[0]] = kinase_mutual_information_singleNetwork(net, kinase_column, 'Substrate')
        for i in range(1,len(networkNames)):
            net = network[networkNames[i]]
            #create substrate columns
            net['Substrate'] = net[accession_column]+'_'+net[site_column]
            if substrate_list:
                net = net[net['Substrate'].isin(substrate_list)]
                if net.empty:
                    raise ValueError("Mapping of substrate_list to the network failed")
                    return -1, -1, -1
            heatmap, normalized, heatlist[networkNames[i]] = kinase_mutual_information_singleNetwork(net, kinase_column, 'Substrate')
            normalized_all = pd.concat([normalized_all, normalized], sort=True)
            heatmap_all = pd.concat([heatmap_all, heatmap], sort=True)
        average_MI = normalized_all.groupby(level=0).mean()
        average_heatmap = heatmap_all.groupby(level=0).mean()
        return average_heatmap, average_MI, heatlist
    else:
        net = network
        #create substrate columns
        net['Substrate'] = net[accession_column]+'_'+net[site_column]
        if substrate_list:
            net = net[net['Substrate'].isin(substrate_list)]
            if net.empty:
                raise ValueError("Mapping of substrate_list to the network failed")
                return -1, -1, -1
        heatmap, normalized, heatlist = kinase_mutual_information_singleNetwork(net, kinase_column, 'Substrate')
        return heatmap, normalized, heatlist



def kinase_mutual_information_singleNetwork(network, kinase_column, substrate_column):
    """
    Finds mutual information shared between kinases based on the substrate phosphorylated
    Mutual Information is defined as the intersection substrates between two kinases
    A substrate is defined as the substrate accession and site, i.e. P54760_Y596.
    Normalization is performed by comparing intersection of kinases vs union of the two kinases
    This the the Jaccard Index. Jaccard Distance can be calcualted by taking 1 - JI

    Parameters
    ----------
    network : pandas dataframe or dictionary of pandas dataframe
        The network to analyze for mutual kinase information. Can send a dictionary of multiple pandas dataframes and this will average the MI across tehm
    kinase_column : str
        Column in network that contiains kinase information
    substrate_column : str
        Column in network that contains substrate information
    substrate_list : array
        Optional and default is no subset list to use. You can calculate the MI within network(s) for only the evidence given in a substrate_evidence_list (must matche substrate_column of network passed in)
    
    Returns
    --------
    heatmap : pandas dataframe
        Number of substrates that overlap between kinases
    normalized : pandas dataframe
        Normalized mutual information into Jaccard Index. 
        size of intersection of two kinase networks / size of union of two kinase networks.
    heatlist : pandas dataframe
        intersection of kinase networks
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


