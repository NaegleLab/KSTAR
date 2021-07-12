import pickle
import os
import pandas as pd
from kstar import helpers
import requests
import tarfile

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
RESOURCE_DIR = KSTAR_DIR+"/RESOURCE_FILES"

#DOWNLOAD RESOURCE_FILES from Figshare to get started using networks 
#Directories for the project (where to find resource file folder), resource folder, and network directories are all set by directories.txt, which
# can be updated with update_directories().
directory_file = open(f'{KSTAR_DIR}/kstar/directories.txt','r')
directories = []
for line in directory_file:
    directories.append(line.split()[0])
directory_file.close()


NETWORK_DIR = directories[0]



RESOURCE_URL = 'https://ndownloader.figshare.com/files/28762653' #location of corresponding release of data

# RESOURCE_DIR dependencies
HUMAN_REF_FASTA_FILE = f"{RESOURCE_DIR}/humanProteome.fasta"  #download from KinPred https://doi.org/10.1101/2020.08.10.244426
try:
    HUMAN_REF_SEQUENCES = helpers.process_fasta_file(HUMAN_REF_FASTA_FILE)
except FileNotFoundError:
    HUMAN_REF_SEQUENCES = None
    print('Could not find reference proteome. Please install the resource directory using config.install_resource_files()')

HUMAN_REF_PHOSPHO_FILE = f"{RESOURCE_DIR}/HumanPhosphoProteome.csv" #download from KSTAR FIGSHARE, or use helpers folder generate to create a new one
try:
    HUMAN_REF_COMPENDIA = pd.read_csv(HUMAN_REF_PHOSPHO_FILE)
except FileNotFoundError:
    HUMAN_REF_COMPENDIA = None
    print('Could not find reference proteome. Please install the resource directory using config.install_resource_files()')


NETWORK_Y_PICKLE = f"{NETWORK_DIR}/network_Y.p" # created by create_networkin_pickles()
NETWORK_ST_PICKLE = f"{NETWORK_DIR}/network_ST.p" #created by create_networkin_pickles()

#check to see if networks can be found
if not os.path.isdir(NETWORK_DIR):
    print('Could not find network directory. Please update using config.update_network_directory().')
elif not os.path.isfile(NETWORK_Y_PICKLE) and not os.path.isfile(NETWORK_ST_PICKLE):
    print('Could not find any network pickles. Please create network pickles using config.create_network_pickles().')

# COLUMN NAMES USED FOR KSTAR
KSTAR_ACCESSION = 'KSTAR_ACCESSION'
KSTAR_PEPTIDE = 'KSTAR_PEPTIDE'
KSTAR_SITE = 'KSTAR_SITE'
KSTAR_KINASE = 'KSTAR_KINASE'

# Number of cores to use for parallelization, set to 1 to avoid multiprocessing
PROCESSES = 4


## END DECLARATION OF GLOBALS


def install_resource_files():
    """Retrieves RESOURCE_FILES that are the companion for this version release from FigShare, unzips them to the correct directory for resource files."""

    print("Requesting file")
    r = requests.get(RESOURCE_URL, allow_redirects=True)
    outputFile = KSTAR_DIR+"/RESOURCE_FILES.tar.gz"

    open(outputFile, 'wb').write(r.content)

    
    t = tarfile.open(outputFile, 'r')


    print("Extracting %s"% outputFile)
    t.extractall(KSTAR_DIR)

#def update_network_directory(directory, *kwargs):
def update_network_directory(directory, KSTAR_DIR = KSTAR_DIR, NETWORK_DIR = NETWORK_DIR):
    """
    Update the location of network the network files, and verify that all necessary files are located in directory
    
    Parameters
    ----------
    directory: string
        path to where network files are located
    
    """
    
    #check that directory exists
    if not os.path.isdir(directory):
        print('Directory not found, so directories.txt was not updated. Please verify that directory is correct')
        return NETWORK_DIR
    
    #update network directory in directories.txt (permanent change)
    with open(f'{KSTAR_DIR}/kstar/directories.txt', 'w') as d:
        d.write(directory)
    
    #check that INDIVIDUAL_NETWORKS/ are present and they are .tsv or .csv files. Then check that network pickles exist.
    if not os.path.isdir(f'{directory}/Y/INDIVIDUAL_NETWORKS') and not os.path.isdir(f'{directory}/ST/INDIVIDUAL_NETWORKS'):
        print('Found network directory, but not the individual networks. Make sure networks are in NETWORK_DIR/{Y or ST}/INDIVIDUAL_NETWORKS')
    elif not os.path.isfile(NETWORK_Y_PICKLE) and not os.path.isfile(NETWORK_ST_PICKLE):
        print('Found individual networks, but no network pickles. Please create network pickles using config.create_network_pickles().')
    
    return directory

def create_network_pickles(phosphoTypes = ['Y','ST']):
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
        
def check_configuration():
    """
    Verify that all necessary files are downloadable and findable
    """
    ready_Y = True
    ready_ST = True
    if not os.path.isfile(HUMAN_REF_FASTA_FILE) or not os.path.isfile(HUMAN_REF_PHOSPHO_FILE):
        print('Could not find reference proteome and/or phosphoproteome. Please install the resource directory using config.install_resource_files()')
        ready_Y, ready_ST = False, False
    if not os.path.isdir(NETWORK_DIR):
        print('Could not find network directory. Please update using config.update_network_directory().')
        ready_Y, ready_ST = False, False
    else:
        if not os.path.isfile(NETWORK_Y_PICKLE) and not os.path.isfile(NETWORK_ST_PICKLE):
            print('Could not find any network pickles. Please create network pickles using config.create_network_pickles() or change network directory to where pickles are located with config.update_network_directory().')
            ready_Y, ready_ST = False, False
        elif not os.path.isfile(NETWORK_Y_PICKLE):
            ready_Y = False
        elif not os.path.isfile(NETWORK_ST_PICKLE):
            ready_ST = False
    
    if ready_Y and ready_ST:
        print('You are ready to begin! You can generate predictions for both Y and ST networks.')
    elif ready_Y:
        print('You are ready to generate predictions for Y networks, but not ST networks. If you want to generate ST predictions, create the ST network pickle with config.create_network_pickles(phosphoType = ["ST"])')
    elif ready_ST:
        print('You are ready to generate predictions for ST networks, but not Y networks. If you want to generate Y predictions, create the Y network pickle with config.create_network_pickles(phosphoType = ["Y"])')
    
        
    


# network_directory='/Volumes/naegle_lab/KinaseActivity/Data/Subgraph/Modifed NetworKIN/CompendiaLimit'
# networkin = pickle.load( open(f"{network_directory}/compendia_limit.p", "rb" ) )
