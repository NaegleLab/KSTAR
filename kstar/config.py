import pickle
import os
import pandas as pd
"""
These are the GLOBAL variables for the KSTAR project and are set for default procedures for the current release
of KinPred and KSTAR. 

Update these if you want to update:
1. Where resources are placed from KSTAR and KINPRED projects: RESOURCE_DIR
2. Reference proteome: HUMAN_REF_FASTA_FILE
3. Phosphoprotoeme: HUMAN_REF_PHOSPHO
4. Networks Directory: 
5. Networks Pickle: 




"""

## BEGIN DECLARATION OF GLOBALS
#PROJECT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..') #PROJECT Directory is one level up from this config file
#PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESOURCE_DIR = f"{PROJECT_DIR}/RESOURCE_FILES" 

# RESOURCE_DIR depedencies
HUMAN_REF_FASTA_FILE = f"{RESOURCE_DIR}/Raw/HumanProteome/humanProteome_2020-02-26.fasta"  #download from KinPred https://doi.org/10.1101/2020.08.10.244426
HUMAN_REF_PHOSPHO = f"{RESOURCE_DIR}/Human_PhosphoProteome_mapped_annotated_02_26_20.csv"
NETWORK_DIR = f"{RESOURCE_DIR}/NETWORKS/Networkin/"
NETWORK_Y_PICKLE = f"{NETWORK_DIR}/network_Y.p" # created by create_networkin_pickles()
NETWORK_ST_PICKLE = f"{NETWORK_DIR}/network_ST.p" #created by create_networkin_pickles()




## END DECLARATION OF GLOBALS

def create_network_pickles():
	"""
	Given network files declared in globals, create pickles of the kstar object that can then be quickly loaded in analysis
	Assumes that the Network structure has two folders Y and ST under the NETWORK_DIR global variable and that 
	all .csv files in those directories should be loaded into a network pickle.
	"""
	phosphoTypes = ['Y', 'ST']
	for phosphoType in phosphoTypes:
		network = {}
		directory = f"{NETWORK_DIR}/{phosphoType}/"
		#get all csv files in that directory 
		for file in os.listdir(directory):
			if file.endswith(".csv"):
				#get the value of the network number
				file_noext = file.strip(".csv").split('_')
				key_name = 'nkin'+str(file_noext[-1])
				#print("Debug: key name is %s"%(key_name))
				network[key_name] = pd.read_csv(f"{directory}{file}")
		print("Loaded %d number of networks for phosphoType %s"%(len(network), phosphoType))
		pickle.dump(network, open(f"{NETWORK_DIR}/network_{phosphoType}.p", "wb"))


# network_directory='/Volumes/naegle_lab/KinaseActivity/Data/Subgraph/Modifed NetworKIN/CompendiaLimit'
# networkin = pickle.load( open(f"{network_directory}/compendia_limit.p", "rb" ) )
