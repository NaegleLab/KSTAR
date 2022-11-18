import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt

from kstar import config

def numUniqueSubstrates(network = None, acc_col = 'KSTAR_ACCESSION', site_col = 'KSTAR_SITE'):
    substrates = list(network['KSTAR_ACCESSION'] + '_' + network['KSTAR_SITE'])
    return len(np.unique(substrates))
    
def averageUniqueSubstrates_KSTAR(mod_types = ['Y','ST']):
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
    
    
def experimentCoverage(experiment, network = None, all_data_columns = True, exp_cols = ['KSTAR_ACCESSION', 'KSTAR_SITE'], net_cols = ['KSTAR_ACCESSION', 'KSTAR_SITE']):
    
    net_substrates = list(network['KSTAR_ACCESSION']+'_'+network['KSTAR_SITE'])
    if !all_data_columns:
        exp_substrates = list(experiment['KSTAR_ACCESSION']+'_'+experiment['KSTAR_SITE'])
        for mod in ['Y','ST']:
            if mod == 'ST':
                mod_sub = [sub for sub in exp_substrates if '_S' in sub or '_T' in sub]
                total_substrates = len(mod_sub)
            
            else:
                mod_sub = [sub for sub in exp_substrates if '_Y' in sub]
                total_substrates = len(mod_sub)
            if total_substrates > 10:       
                for net in overlap_all[mod].keys():
                    if 'KSTAR' in net:
                        for kstar_net in kstar_names:
                            #for individual networks
                            overlap = len(set(mod_sub).intersection(set(kstar_substrates[mod][kstar_net])))
                            overlap_all[mod]['KSTAR_Individual_Networks'].append(overlap/total_substrates)
                        #across all networks
                        overlap = len(set(mod_sub).intersection(set(kstar_all[mod])))
                        overlap_all[mod]['KSTAR_Across_Networks'].append(overlap/total_substrates)
                    else:
                        overlap = len(set(mod_sub).intersection(net_substrates[net]))
                        overlap_all[mod][net].append(overlap/total_substrates)

    