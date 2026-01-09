import os
import pandas as pd
import json
from kstar import helpers
from importlib.resources import files

"""
These are the GLOBAL variables for the KSTAR project and are set for default procedures for the current release
of KinPred and KSTAR. This is a trim version for galaxy, which will come packaged with resource files and networks.

"""

## BEGIN DECLARATION OF GLOBALS

#grab package directory (include specific folder)
DOCKER_RESOURCE_DIR = os.environ.get("KSTAR_RESOURCE_DIR", "/opt/kstar_resources/")
RESOURCE_DIR = DOCKER_RESOURCE_DIR + "/RESOURCE_FILES/"
#read information about reference files
with open(f"{RESOURCE_DIR}/reference_info.json", "r") as f:
    REFERENCE_INFO = json.load(f)

NETWORK_DIR = DOCKER_RESOURCE_DIR + "/NETWORKS/"
DEFAULT_RANDOM_ACTIVITIES_DIR = NETWORK_DIR

#specific default networks
NETWORK_NAME = {
    'Y': 'Default',
    'ST': 'Default'
}
NETWORK_SUBDIR = {'Y': f"{NETWORK_DIR}/Y/{NETWORK_NAME['Y']}/",
                  'ST': f"{NETWORK_DIR}/ST/{NETWORK_NAME['ST']}/"}


#CUSTOM_RANDOM_ACTIVITIES_DIR = directories[1]
CUSTOM_RANDOM_ACTIVITIES_DIR = None


# RESOURCE_DIR dependencies
HUMAN_REF_FASTA_FILE = RESOURCE_DIR + "humanProteome.fasta" # download from KinPred https://doi.org/10.1101/2020.08.10.244426
HUMAN_REF_SEQUENCES = helpers.process_fasta_file(HUMAN_REF_FASTA_FILE)


HUMAN_REF_PHOSPHO_FILE = RESOURCE_DIR + "HumanPhosphoProteome.csv"  # download from KSTAR FIGSHARE, or use helpers folder generate to create a new one
HUMAN_REF_COMPENDIA = pd.read_csv(HUMAN_REF_PHOSPHO_FILE)


# COLUMN NAMES USED FOR KSTAR
KSTAR_ACCESSION = 'KSTAR_ACCESSION'
KSTAR_PEPTIDE = 'KSTAR_PEPTIDE'
KSTAR_SITE = 'KSTAR_SITE'
KSTAR_KINASE = 'KSTAR_KINASE'

# PARAMETERS USED FOR RUNNING KSTAR USING PREGENERATED RANDOM ACTIVITIES
USE_PREGENERATED_RANDOM_ACTIVITIES = True
SAVE_NEW_RANDOM_ACTIVITIES = False
SAVE_RANDOM_EXPERIMENTS = False
DEFAULT_RANDOM_ACTIVITIES_DIR = NETWORK_DIR

# read in network parameters

NETWORK_INFO = {}
NETWORK_HASH = {}
#check that network directory exists
if not os.path.isdir(f"{NETWORK_DIR}"):
    print('Warning: Could not find network directory as specified in configuration file. If you have not downloaded networks, please do so using config.install_network_files(). If using your own networks, please update the configuration file using config.update_configuration() to point to correct directory and network name')
else:
    #check to make sure Y network folder exists, if it does, make sure its compatible with reference
    if not os.path.isdir(NETWORK_SUBDIR['Y']):
        print(f'Warning: Could not find tyrosine network that goes by the name {NETWORK_NAME["Y"]}. Please update the name of the network you would like to use in the configuration file using config.update_configuration()')
    else:
        #load network information, make sure its compatible with reference    
        NETWORK_INFO['Y'] = helpers.parse_network_information(NETWORK_SUBDIR['Y'])
        NETWORK_HASH['Y'] = NETWORK_INFO['Y']['unique_network_id']
        if NETWORK_INFO['Y']['unique_reference_id'] != REFERENCE_INFO['unique_reference_id']:
            print(f'Warning: The tyrosine network you have selected does not match the reference proteome used in this KSTAR installation. Please verify that the KSTAR networks were created using the same reference phosphoproteome.')

    #check to make sure ST network folder exists,if it does, make sure its compatible with reference
    if not os.path.isdir(NETWORK_SUBDIR['ST']):
        print(f'Warning: Could not find serine/threonine network that goes by the name {NETWORK_NAME["ST"]}. Please update the name of the network you would like to use in the configuration file using config.update_configuration()')
    else:
        ST_NETWORK_INFO = helpers.parse_network_information(NETWORK_SUBDIR['ST'])
        NETWORK_HASH['ST'] = ST_NETWORK_INFO['unique_network_id']
        if ST_NETWORK_INFO['unique_reference_id'] != REFERENCE_INFO['unique_reference_id']:
            print(f'Warning: The serine/threonine network you have selected does not match the reference proteome used in this KSTAR installation. Please verify that the KSTAR networks were created using the same reference phosphoproteome.')





