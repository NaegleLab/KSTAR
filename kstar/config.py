import pickle
import os
import pandas as pd
from kstar import helpers

"""
These are the GLOBAL variables for the KSTAR project and are set for default procedures for the current release
of KinPred and KSTAR. 

Update these if you want to update:
1. Where resources are placed from KSTAR and KINPRED projects: RESOURCE_DIR
2. Reference proteome: HUMAN_REF_FASTA_FILE
3. Phosphoprotoeme: HUMAN_REF_PHOSPHO
4. Networks Directory: Downloaded from KSTAR Figshare repository, this includes heuristic networks suggested for 
5. Networks Pickle: 

"""

## BEGIN DECLARATION OF GLOBALS

KSTAR_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


#DOWNLOAD RESOURCE_FILES from Figshare to get started using networks 
#Naegle, Kristen (2021): RESOURCE_FILES. figshare. Software. https://doi.org/10.6084/m9.figshare.14885121.v3 
#Directories for the project (where to find resource file folder), resource folder, and network directories are all set by directories.txt, which
# can be updated with update_directories().
directory_file = open(f'{KSTAR_DIR}/kstar/directories.txt','r')
directories = []
for line in directory_file:
    directories.append(line.split()[0])
directory_file.close()

RESOURCE_DIR = directories[0]
NETWORK_DIR = directories[1]


# RESOURCE_DIR dependencies
HUMAN_REF_FASTA_FILE = f"{RESOURCE_DIR}/Raw/HumanProteome/humanProteome_2020-02-26.fasta"  #download from KinPred https://doi.org/10.1101/2020.08.10.244426
try:
    HUMAN_REF_SEQUENCES = helpers.process_fasta_file(HUMAN_REF_FASTA_FILE)
except FileNotFoundError:
    HUMAN_REF_SEQUENCES = None
    print('Could not find reference proteome. Please update the resource directory using config.update_directories()')

HUMAN_REF_PHOSPHO_FILE = f"{RESOURCE_DIR}/Human_PhosphoProteome_mapped_annotated_02_26_20.csv" #download from KSTAR FIGSHARE, or use helpers folder generate to create a new one
try:
    HUMAN_REF_COMPENDIA = pd.read_csv(HUMAN_REF_PHOSPHO_FILE)
except FileNotFoundError:
    HUMAN_REF_COMPENDIA = None
    print('Could not find reference phosphoproteome. Please update the resource directory using config.update_directories()')


NETWORK_Y_PICKLE = f"{NETWORK_DIR}/network_Y.p" # created by create_networkin_pickles()
NETWORK_ST_PICKLE = f"{NETWORK_DIR}/network_ST.p" #created by create_networkin_pickles()

# COLUMN NAMES USED FOR KSTAR
KSTAR_ACCESSION = 'KSTAR_ACCESSION'
KSTAR_PEPTIDE = 'KSTAR_PEPTIDE'
KSTAR_SITE = 'KSTAR_SITE'
KSTAR_KINASE = 'KSTAR_KINASE'

# Number of cores to use for parallelization, set to 1 to avoid multiprocessing
PROCESSES = 4

## END DECLARATION OF GLOBALS

def update_directories(resource = None, network = None, update_references = True, directories = directories, KSTAR_DIR = KSTAR_DIR):
    """
    Update the global variables that indicate where to find resource and network files.
    Will only update the directories that have been inputted during function use. Will also return each of the new directories
    which is the only way to update directory globals in the same iteration.
    
    Location of files should be as such:
        Reference proteome: found in project/resource
        Reference phosphoproteome: found in project/resource
        Individual Networks: found in project/resource/network
    
    Parameters
    ----------
    resource: string
        directory where resource files are kept. If it is not being updated, set to None (default).
    network: string
        directory where the network files are located. If it is not being updated, set to None (default).
    update_references: bool
        indicates whether to update where to find reference files based on provided directories. recommended to be True.
    other_parameters:
        the remaining parameters provide function access to the directories config.py is currently pointing to. Do not edit.
    """
    
    #update resource directory
    if resource is not None:
        directories[0] = resource

        
    #update network directory
    if network is not None:
        directories[1] = network

    
    with open(f'{KSTAR_DIR}/kstar/directories.txt', 'w') as d:
        d.writelines('\n'.join(directories))
    
    if update_references:
        HUMAN_REF_FASTA_FILE = f"{directories[0]}/Raw/HumanProteome/humanProteome_2020-02-26.fasta"
        try:
            HUMAN_REF_SEQUENCES = helpers.process_fasta_file(HUMAN_REF_FASTA_FILE)
        except FileNotFoundError:
            HUMAN_REF_SEQUENCES = None
            print('Still could not locate reference proteome. Please verify that provided directories are correct')
            
        HUMAN_REF_PHOSPHO_FILE = f"{directories[0]}/Human_PhosphoProteome_mapped_annotated_02_26_20.csv" #download from KSTAR FIGSHARE, or use helpers folder generate to create a new one
        try:
            HUMAN_REF_COMPENDIA = pd.read_csv(HUMAN_REF_PHOSPHO_FILE)
        except FileNotFoundError:
            HUMAN_REF_COMPENDIA = None
            print('Still could not locate reference phosphoproteome. Please verify that provided directories are correct')
        
        return directories, HUMAN_REF_SEQUENCES, HUMAN_REF_COMPENDIA
    else:
        return directories

def create_network_pickles():
	"""
	Given network files declared in globals, create pickles of the kstar object that can then be quickly loaded in analysis
	Assumes that the Network structure has two folders Y and ST under the NETWORK_DIR global variable and that 
	all .csv files in those directories should be loaded into a network pickle.
	"""
	phosphoTypes = ['Y', 'ST']
	for phosphoType in phosphoTypes:
		network = {}
		directory = f"{NETWORK_DIR}/{phosphoType}/INDIVIDUAL_NETWORKS/"
		#get all csv files in that directory 
		for file in os.listdir(directory):
			if file.endswith(".tsv"):
				#get the value of the network number
				file_noext = file.strip(".tsv").split('_')
				key_name = 'nkin'+str(file_noext[1])
				#print("Debug: key name is %s"%(key_name))
				network[key_name] = pd.read_csv(f"{directory}{file}", sep='\t')
		print("Loaded %d number of networks for phosphoType %s"%(len(network), phosphoType))
		pickle.dump(network, open(f"{NETWORK_DIR}/network_{phosphoType}.p", "wb"))
		print(f"Saved pickle file at {NETWORK_DIR}/network_{phosphoType}.p")


# network_directory='/Volumes/naegle_lab/KinaseActivity/Data/Subgraph/Modifed NetworKIN/CompendiaLimit'
# networkin = pickle.load( open(f"{network_directory}/compendia_limit.p", "rb" ) )
