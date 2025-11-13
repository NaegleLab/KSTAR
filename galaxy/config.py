import os
import pandas as pd
from kstar import helpers
from importlib.resources import files

"""
These are the GLOBAL variables for the KSTAR project and are set for default procedures for the current release
of KinPred and KSTAR. This is a trim version for galaxy, which will come packaged with resource files and networks.

"""

## BEGIN DECLARATION OF GLOBALS

#grab package directory (include specific folder)
#KSTAR_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "/kstar"
#RESOURCE_DIR = KSTAR_DIR + "/RESOURCE_FILES"
#NETWORK_DIR = KSTAR_DIR + "/NETWORKS/NetworKIN/"
DOCKER_RESOURCE_DIR = os.environ.get("KSTAR_RESOURCE_DIR", "/opt/kstar_resources/")
#KSTAR_DIR = str(files("kstar"))
RESOURCE_DIR = DOCKER_RESOURCE_DIR + "/RESOURCE_FILES"
NETWORK_DIR = DOCKER_RESOURCE_DIR + "/NETWORKS/NetworKIN/"



#CUSTOM_RANDOM_ACTIVITIES_DIR = directories[1]
CUSTOM_RANDOM_ACTIVITIES_DIR = NETWORK_DIR + "/CUSTOM_PREGENERATED_RANDOM_ACTIVITIES"


# RESOURCE_DIR dependencies
HUMAN_REF_FASTA_FILE = RESOURCE_DIR + "RESOURCE_FILES/humanProteome.fasta" # download from KinPred https://doi.org/10.1101/2020.08.10.244426
HUMAN_REF_SEQUENCES = helpers.process_fasta_file(HUMAN_REF_FASTA_FILE)


HUMAN_REF_PHOSPHO_FILE = RESOURCE_DIR + "/HumanPhosphoProteome.csv"  # download from KSTAR FIGSHARE, or use helpers folder generate to create a new one
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
NETWORK_HASH_Y = 'd73e971979661777c05182efad8b314b1aa78896e242810cce8ecd51b23c85c6'
NETWORK_HASH_ST = '68328e443909aecbb6a1f6b0572da171f6b1213c95ecf1928b621ef39826c953'





