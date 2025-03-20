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
RESOURCE_DIR = KSTAR_DIR + "/RESOURCE_FILES"

# DOWNLOAD RESOURCE_FILES from Figshare to get started using networks
# Directories for the project (where to find resource file folder), resource folder, and network directories are all set by directories.txt, which
# can be updated with update_directories().
directory_file = open(f'{KSTAR_DIR}/kstar/directories.txt', 'r')
directories = []
for line in directory_file:
    directories.append(line.split()[0])
directory_file.close()

NETWORK_DIR = directories[0]

RESOURCE_URL = 'https://ndownloader.figshare.com/files/28762653'  # location of corresponding release of data
NETWORK_URL = 'https://ndownloader.figshare.com/files/52178324'

# RESOURCE_DIR dependencies
HUMAN_REF_FASTA_FILE = f"{RESOURCE_DIR}/humanProteome.fasta"  # download from KinPred https://doi.org/10.1101/2020.08.10.244426
try:
    HUMAN_REF_SEQUENCES = helpers.process_fasta_file(HUMAN_REF_FASTA_FILE)
except FileNotFoundError:
    HUMAN_REF_SEQUENCES = None
    print(
        'Could not find reference proteome. Please install the resource directory using config.install_resource_files()')

HUMAN_REF_PHOSPHO_FILE = f"{RESOURCE_DIR}/HumanPhosphoProteome.csv"  # download from KSTAR FIGSHARE, or use helpers folder generate to create a new one
try:
    HUMAN_REF_COMPENDIA = pd.read_csv(HUMAN_REF_PHOSPHO_FILE)
except FileNotFoundError:
    HUMAN_REF_COMPENDIA = None
    print(
        'Could not find reference proteome. Please install the resource directory using config.install_resource_files()')

NETWORK_Y_PICKLE = f"{NETWORK_DIR}/network_Y.p"  # created by create_networkin_pickles()
NETWORK_ST_PICKLE = f"{NETWORK_DIR}/network_ST.p"  # created by create_networkin_pickles()

# check to see if networks can be found
if not os.path.isdir(NETWORK_DIR):
    print('Could not find network directory. Please update using config.update_network_directory().')
elif not os.path.isfile(NETWORK_Y_PICKLE) and not os.path.isfile(NETWORK_ST_PICKLE):
    print('Could not find any network pickles. Please create network pickles using config.create_network_pickles().')

# COLUMN NAMES USED FOR KSTAR
KSTAR_ACCESSION = 'KSTAR_ACCESSION'
KSTAR_PEPTIDE = 'KSTAR_PEPTIDE'
KSTAR_SITE = 'KSTAR_SITE'
KSTAR_KINASE = 'KSTAR_KINASE'

# PARAMETERS USED FOR RUNNING KSTAR USING PREGENERATED EXPERIMENTS
USE_PREGEN_DATA = True
SAVE_NEW_PRECOMPUTE = False
PREGENERATED_EXPERIMENTS_PATH = NETWORK_DIR
NETWORK_HASH_Y = 'd73e971979661777c05182efad8b314b1aa78896e242810cce8ecd51b23c85c6'
NETWORK_HASH_ST = '68328e443909aecbb6a1f6b0572da171f6b1213c95ecf1928b621ef39826c953'
DIRECTORY_FOR_SAVE_PRECOMPUTE = None


## END DECLARATION OF GLOBALS


def install_resource_files():
    """Retrieves RESOURCE_FILES that are the companion for this version release from FigShare, unzips them to the correct directory for resource files."""

    print("Requesting resource file")
    r = requests.get(RESOURCE_URL, allow_redirects=True)
    outputFile = KSTAR_DIR + "/RESOURCE_FILES.tar.gz"

    open(outputFile, 'wb').write(r.content)

    t = tarfile.open(outputFile, 'r')

    print("Extracting %s" % outputFile)
    t.extractall(KSTAR_DIR)


# def update_network_directory(directory, *kwargs):

def install_network_files(target_dir=None):
    """Retrieves Network files that are the companion for this version release from FigShare, unzips them to the specified directory."""
    print("Requesting network file")
    r = requests.get(NETWORK_URL, allow_redirects=True)

    # Use target_dir if provided, otherwise default to KSTAR_DIR
    install_dir = target_dir if target_dir else KSTAR_DIR
    outputFile = os.path.join(install_dir, "NETWORKS.tar.gz")

    # Create target directory if it doesn't exist
    os.makedirs(install_dir, exist_ok=True)

    open(outputFile, 'wb').write(r.content)

    t = tarfile.open(outputFile, 'r')

    print(f"Extracting {outputFile}")
    t.extractall(install_dir)

def update_network_directory(directory, create_pickles=True, KSTAR_DIR=KSTAR_DIR, NETWORK_DIR=NETWORK_DIR):
    """
    Update the location of network the network files, and verify that all necessary files are located in directory

    Parameters
    ----------
    directory: string
        path to where network files are located


    """

    # check that directory exists
    if not os.path.isdir(directory):
        print(
            'Directory not found, so directories.txt was not updated and no pickles were generated. Please verify that directory is correct')
        return

    # update network directory in config (temporary change)
    NETWORK_DIR = directory
    # update network directory in directories.txt (permanent change)
    with open(f'{KSTAR_DIR}/kstar/directories.txt', 'w') as d:
        d.write(directory)
    print('Network directory updated.')

    # If create_pickles is true, look for individual networks and create pickles
    if create_pickles:
        if not os.path.isfile(f"{directory}/network_Y.p") and not os.path.isfile(f"{directory}/network_ST.p"):
            if os.path.isdir(f'{directory}/Y/INDIVIDUAL_NETWORKS') and os.path.isdir(
                    f'{directory}/ST/INDIVIDUAL_NETWORKS'):
                print('Found individual networks for both Y and ST, generating pickles for both')
                create_network_pickles(network_directory=directory)
            elif os.path.isdir(f'{directory}/Y/INDIVIDUAL_NETWORKS'):
                print('Individual networks found for Y, but not for ST. Generating Y pickle only.')
                create_network_pickles(['Y'], network_directory=directory)
            elif os.path.isdir(f'{directory}/ST/INDIVIDUAL_NETWORKS'):
                print('Individual networks found for ST, but not for Y. Generating ST pickle only.')
                create_network_pickles(['ST'], network_directory=directory)
        elif not os.path.isfile(f"{directory}/network_Y.p"):
            print('ST pickle already created. Generating Y pickle')
            if os.path.isdir(f'{directory}/Y/INDIVIDUAL_NETWORKS'):
                create_network_pickles(['Y'], network_directory=directory)
            else:
                print(
                    'Could not find individual networks for Y. Make sure networks are deposited in "{NETWORK_DIR}/Y/INDIVIDUAL_NETWORKS/", then rerun create_network_pickles()')
        elif not os.path.isfile(f'{directory}/network_ST.p'):
            print('Y pickle already created. Generating ST pickle')
            if os.path.isdir(f'{directory}/ST/INDIVIDUAL_NETWORKS'):
                create_network_pickles(['ST'], network_directory=directory)
            else:
                print(
                    'Could not find individual networks for ST. Make sure networks are deposited in "{NETWORK_DIR}/ST/INDIVIDUAL_NETWORKS/", , then rerun create_network_pickles()')
        else:
            print('Network pickles already generated')

    # update location of network pickles
    NETWORK_Y_PICKLE = f"{NETWORK_DIR}/network_Y.p"  # created by create_networkin_pickles()
    NETWORK_ST_PICKLE = f"{NETWORK_DIR}/network_ST.p"  # created by create_networkin_pickles()
    return NETWORK_DIR, NETWORK_Y_PICKLE, NETWORK_ST_PICKLE


# def create_network_pickles(phosphoType = ['Y','ST'], *kwargs)
def create_network_pickles(phosphoTypes=['Y', 'ST'], network_directory=NETWORK_DIR):
    """
    Given network files declared in globals, create pickles of the kstar object that can then be quickly loaded in analysis
    Assumes that the Network structure has two folders Y and ST under the NETWORK_DIR global variable and that
    all .csv files in those directories should be loaded into a network pickle.
    """

    for phosphoType in phosphoTypes:
        network = {}
        if not os.path.isfile(f"{network_directory}/network_{phosphoType}.p"):
            directory = f"{network_directory}/{phosphoType}/INDIVIDUAL_NETWORKS/"
            # get all csv files in that directory
            for file in os.listdir(directory):
                if file.endswith(".tsv"):
                    # get the value of the network number
                    file_noext = file.strip(".tsv").split('_')
                    key_name = 'nkin' + str(file_noext[1])
                    # print("Debug: key name is %s"%(key_name))
                    network[key_name] = pd.read_csv(f"{directory}{file}", sep='\t')
            print("Loaded %d number of networks for phosphoType %s" % (len(network), phosphoType))
            pickle.dump(network, open(f"{network_directory}/network_{phosphoType}.p", "wb"))
            print(f"Saved pickle file at {network_directory}/network_{phosphoType}.p")
        else:
            print(f"{phosphoType} network pickle already generated")


def check_configuration():
    """
    Verify that all necessary files are downloadable and findable
    """
    ready_Y = True
    ready_ST = True
    if not os.path.isfile(HUMAN_REF_FASTA_FILE) or not os.path.isfile(HUMAN_REF_PHOSPHO_FILE):
        print(
            'Could not find reference proteome and/or phosphoproteome. Please install the resource directory using config.install_resource_files()')
        ready_Y, ready_ST = False, False
    if not os.path.isdir(NETWORK_DIR):
        print('Could not find network directory. Please update using config.update_network_directory().')
        ready_Y, ready_ST = False, False
    else:
        if not os.path.isfile(NETWORK_Y_PICKLE) and not os.path.isfile(NETWORK_ST_PICKLE):
            print(
                'Could not find any network pickles. Please create network pickles using config.create_network_pickles() or change network directory to where pickles are located with config.update_network_directory().')
            ready_Y, ready_ST = False, False
        elif not os.path.isfile(NETWORK_Y_PICKLE):
            ready_Y = False
        elif not os.path.isfile(NETWORK_ST_PICKLE):
            ready_ST = False

    if ready_Y and ready_ST:
        print('You are ready to begin! You can generate predictions for both Y and ST networks.')
    elif ready_Y:
        print(
            'You are ready to generate predictions for Y networks, but not ST networks. If you want to generate ST predictions, create the ST network pickle with config.create_network_pickles(phosphoType = ["ST"])')
    elif ready_ST:
        print(
            'You are ready to generate predictions for ST networks, but not Y networks. If you want to generate Y predictions, create the Y network pickle with config.create_network_pickles(phosphoType = ["Y"])')

