import pickle
import os, json, re
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

#grab package directory (include specific folder)
KSTAR_DIR = os.path.dirname(os.path.abspath(__file__))
RESOURCE_DIR = KSTAR_DIR + "/RESOURCE_FILES"


###### Customizable configuration ####################
configuration_file = f'{KSTAR_DIR}/configuration.json'
with open(configuration_file, 'r') as cf:
    config_data = json.load(cf)

#extract network directory from configuration file, make sure it can be found. By default this will be located in KSTAR_DIR/NETWORKS/NetworKIN
NETWORK_DIR = config_data["network_directory"]
CUSTOM_RANDOM_ACTIVITIES_DIR = config_data["custom_pregenerated_experiments_dir"]
# PARAMETERS USED FOR RUNNING KSTAR USING PREGENERATED RANDOM ACTIVITIES
USE_PREGENERATED_RANDOM_ACTIVITIES = config_data['use_pregenerated_random_activities']
SAVE_NEW_RANDOM_ACTIVITIES = config_data['save_new_random_activities']
SAVE_RANDOM_EXPERIMENTS = config_data['save_random_experiments']
NETWORK_HASH_Y = config_data['Y_network_hash']
NETWORK_HASH_ST = config_data['ST_network_hash']


######### static configuration (user cannot change) ################
RESOURCE_URL = 'https://ndownloader.figshare.com/files/28762653'  # location of corresponding release of data
NETWORK_URL = 'https://ndownloader.figshare.com/files/52178324'

# COLUMN NAMES USED FOR KSTAR
KSTAR_ACCESSION = 'KSTAR_ACCESSION'
KSTAR_PEPTIDE = 'KSTAR_PEPTIDE'
KSTAR_SITE = 'KSTAR_SITE'
KSTAR_KINASE = 'KSTAR_KINASE'



############### Configuration setup ####################

#load resource files info
#REFERENCE_INFO_FILE = f"{RESOURCE_DIR}/reference_info.json"
#with open(REFERENCE_INFO_FILE, "r") as f:
#    REFERENCE_INFO = json.load(f)

# load reference proteome sequences
HUMAN_REF_FASTA_FILE = f"{RESOURCE_DIR}/humanProteome.fasta"  # download from KinPred https://doi.org/10.1101/2020.08.10.244426
try:
    HUMAN_REF_SEQUENCES = helpers.process_fasta_file(HUMAN_REF_FASTA_FILE)
except FileNotFoundError:
    HUMAN_REF_SEQUENCES = None
    print(
        'Could not find reference proteome. Please install the resource directory using config.install_resource_files()')

# load reference phosphoproteome
HUMAN_REF_PHOSPHO_FILE = f"{RESOURCE_DIR}/HumanPhosphoProteome.csv"  # download from KSTAR FIGSHARE, or use helpers folder generate to create a new one
try:
    HUMAN_REF_COMPENDIA = pd.read_csv(HUMAN_REF_PHOSPHO_FILE)
except FileNotFoundError:
    HUMAN_REF_COMPENDIA = None
    print(
        'Could not find reference proteome. Please install the resource directory using config.install_resource_files()')

#check network directory: if not found, set to default
if NETWORK_DIR is None or NETWORK_DIR == "":
    NETWORK_DIR = f"{KSTAR_DIR}/NETWORKS/NetworKIN"

# check to see if networks can be found
if not os.path.isdir(NETWORK_DIR):
    print('Could not find network directory. Please update using config.update_network_directory().')

DEFAULT_RANDOM_ACTIVITIES_DIR = NETWORK_DIR
Y_NETWORK_DIR = f"{NETWORK_DIR}/Y/{NETWORK_HASH_Y}/"
ST_NETWORK_DIR = f"{NETWORK_DIR}/ST/{NETWORK_HASH_ST}/"


## END DECLARATION OF GLOBALS
def update_configuration(network_dir=None, y_network_hash=None, st_network_hash=None, save_random_experiments = None, use_pregenerated_random_activities=None, save_new_random_activities=None, custom_pregenerated_experiments_dir=None):
    """
    Update configuration parameters in current iteration and save to configuration file.

    Parameters
    ----------
    use_pregenerated_random_activities : bool, optional
        Whether to use pregenerated random activities when possible, by default None
    save_new_random_activities : bool, optional
        Whether to save new random activities when they are generated, by default False
    custom_pregenerated_experiments_dir : str, optional
        Directory to save newly generated random activities for future use, by default None
    network_dir : str, optional
        Directory containing the kinase-substrate networks, by default None (which assumes it is located in kstar directory)
    y_network_hash : str, optional
        Unique identifier of the tyrosine network to use by default. 
    st_network_hash : str, optional
        Unique identifier of the serine/threonine network to use by default.
    """
    if use_pregenerated_random_activities is not None:
        global USE_PREGENERATED_RANDOM_ACTIVITIES
        USE_PREGENERATED_RANDOM_ACTIVITIES = use_pregenerated_random_activities
    if save_new_random_activities is not None:
        global SAVE_NEW_RANDOM_ACTIVITIES
        SAVE_NEW_RANDOM_ACTIVITIES = save_new_random_activities
    if custom_pregenerated_experiments_dir is not None:
        if not os.path.isdir(custom_pregenerated_experiments_dir):
            raise ValueError(f"Provided pregenerated experiments directory {custom_pregenerated_experiments_dir} does not exist. Please provide a valid directory.")
        global CUSTOM_RANDOM_ACTIVITIES_DIR
        CUSTOM_RANDOM_ACTIVITIES_DIR = custom_pregenerated_experiments_dir
    if save_random_experiments is not None:
        global SAVE_RANDOM_EXPERIMENTS
        SAVE_RANDOM_EXPERIMENTS = save_random_experiments

    if network_dir is not None:
        #check to make sure directory exists
        if not os.path.isdir(network_dir):
            raise ValueError(f"Provided network directory {network_dir} does not exist. Please provide a valid directory.")
        global NETWORK_DIR
        global Y_NETWORK_DIR
        global ST_NETWORK_DIR
        NETWORK_DIR = network_dir
        #update specifc network directories if needed
        if y_network_hash is not None:
            global NETWORK_HASH_Y
            NETWORK_HASH_Y = y_network_hash
            #Y_NETWORK_DIR = f"{NETWORK_DIR}/Y/{NETWORK_HASH_Y}/"
            Y_NETWORK_DIR = f"{NETWORK_DIR}/Y/"
        else:
            #Y_NETWORK_DIR = f"{NETWORK_DIR}/Y/{NETWORK_HASH_Y}/"
            Y_NETWORK_DIR = f"{NETWORK_DIR}/Y/"

        if st_network_hash is not None:
            global NETWORK_HASH_ST
            NETWORK_HASH_ST = st_network_hash
            #ST_NETWORK_DIR = f"{NETWORK_DIR}/ST/{NETWORK_HASH_ST}/"
            ST_NETWORK_DIR = f"{NETWORK_DIR}/ST/"
        else:
            #ST_NETWORK_DIR = f"{NETWORK_DIR}/ST/{NETWORK_HASH_ST}/"
            ST_NETWORK_DIR = f"{NETWORK_DIR}/ST/"

    #save changes to configuration file
    config_data = {
        "network_directory": NETWORK_DIR,
        'save_random_experiments': SAVE_RANDOM_EXPERIMENTS,
        "use_pregenerated_random_activities": USE_PREGENERATED_RANDOM_ACTIVITIES,
        "save_new_random_activities": SAVE_NEW_RANDOM_ACTIVITIES,
        "custom_pregenerated_experiments_dir": CUSTOM_RANDOM_ACTIVITIES_DIR,
        "Y_network_hash": NETWORK_HASH_Y,
        "ST_network_hash": NETWORK_HASH_ST
    }

    configuration_file = f'{KSTAR_DIR}/configuration.json'
    with open(configuration_file, 'w') as cf:
        json.dump(config_data, cf, indent=4)



def install_resource_files():
    """Retrieves RESOURCE_FILES that are the companion for this version release from FigShare, unzips them to the correct directory for resource files."""

    print("Requesting resource file")
    r = requests.get(RESOURCE_URL, allow_redirects=True)
    outputFile = KSTAR_DIR + "/RESOURCE_FILES.tar.gz"

    open(outputFile, 'wb').write(r.content)

    t = tarfile.open(outputFile, 'r')

    print("Extracting %s" % outputFile)
    t.extractall(KSTAR_DIR)
    t.close()
    #remove tar file
    os.remove(outputFile)

    #automatically update globals to point to new resource files and load them
    global HUMAN_REF_FASTA_FILE
    global HUMAN_REF_SEQUENCES
    global HUMAN_REF_COMPENDIA
    HUMAN_REF_FASTA_FILE = f"{RESOURCE_DIR}/humanProteome.fasta"
    HUMAN_REF_SEQUENCES = helpers.process_fasta_file(HUMAN_REF_FASTA_FILE)
    HUMAN_REF_COMPENDIA = pd.read_csv(f"{RESOURCE_DIR}/HumanPhosphoProteome.csv")


# def update_network_directory(directory, *kwargs):

def install_network_files(target_dir=None, create_pickles = False):
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
    t.close()
        # remove tar file
    os.remove(outputFile)

    #update network directory
    update_configuration(network_dir=f"{install_dir}/NETWORKS/NetworKIN")
    #update_network_directory(f"{install_dir}/NETWORKS/NetworKIN", create_pickles=False)
    #return NETWORK_DIR, NETWORK_Y_PICKLE, NETWORK_ST_PICKLE

#def update_directory_file():
#    """
#    Update the directories.txt file with new directories for network files and pregenerated #experiments

#    Parameters
#    ----------
#    network_dir: string
#        path to where network files are located
#    pregenerated_experiments_dir: string
#        path to where pregenerated experiments are located
#    custom_pregenerated_experiments_dir: string
#        path to where custom pregenerated experiments are located

#    """
#    lines = [NETWORK_DIR + '\n', CUSTOM_PREGENERATED_EXPERIMENTS_DIR + '\n']
#    # read current directories
#    with open(f'{KSTAR_DIR}/kstar/directories.txt', 'w') as d:
#        d.writelines(lines)

def find_available_network_hashes(phospho_type):
    """
    Find available network hashes in the current network directory, and return dictionary with information about them

    Returns
    -------
    available_networks : dict
        dictionary containing all available networks, in the format -> Network hash : network information dictionary
    """
    available_networks = {}

    for item in os.listdir(NETWORK_DIR + f'/{phospho_type}/'):
        item_path = os.path.join(NETWORK_DIR, item)

        #read network information
        network_info = parse_network_information_json(item_path)
        available_networks[item] = network_info
    return available_networks

def update_network_directory(directory):
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
            'Directory not found, so configuration was not updated and no pickles were generated. Please verify that directory is correct')
        return

    # update network directory in directories.txt (first line of the file)
    update_configuration(network_dir = directory)
    print('Network directory updated.')


# def create_network_pickles(phosphoType = ['Y','ST'], *kwargs)
#def create_network_pickles(phosphoTypes=['Y', 'ST'], network_directory=NETWORK_DIR):
#    """
#    Given network files declared in globals, create pickles of the kstar object that can then be quickly loaded in analysis
#    Assumes that the Network structure has two folders Y and ST under the NETWORK_DIR global variable and that
#    all .csv files in those directories should be loaded into a network pickle.
#    """

#    for phosphoType in phosphoTypes:
#        network = {}
#        if not os.path.isfile(f"{network_directory}/network_{phosphoType}.p"):
#            directory = f"{network_directory}/{phosphoType}/INDIVIDUAL_NETWORKS/"
#            # get all csv files in that directory
#            for file in os.listdir(directory):
#                if file.endswith(".tsv"):
#                    # get the value of the network number
#                    file_noext = file.strip(".tsv").split('_')
#                    key_name = 'nkin' + str(file_noext[1])
#                    # print("Debug: key name is %s"%(key_name))
#                    network[key_name] = pd.read_csv(f"{directory}{file}", sep='\t')
#            print("Loaded %d number of networks for phosphoType %s" % (len(network), phosphoType))
#            pickle.dump(network, open(f"{network_directory}/network_{phosphoType}.p", "wb"))
#            print(f"Saved pickle file at {network_directory}/network_{phosphoType}.p")
#        else:
#            print(f"{phosphoType} network pickle already generated")

#def update_pregenerated_experiments_dir(new_directory):
#    """
#    Update the directory for pregenerated experiments.
#    """
#    global CUSTOM_PREGENERATED_EXPERIMENTS_DIR
#    CUSTOM_PREGENERATED_EXPERIMENTS_DIR = new_directory
#    update_directory_file()

def parse_network_information(network_directory):
    """
    Parse the RUN_INFORMATION.txt file from network pruning run and extract its data.

    Args:
        file_path (str): Path to the RUN_INFORMATION.txt file.

    Returns:
        dict: A dictionary containing the parsed data.
    """
    file_path = os.path.join(network_directory, "RUN_INFORMATION.txt")

    try:
        with open(file_path, 'r') as file:
            content = file.read()

        return {
            "unique_id": re.search(r"Unique ID:\s+([a-fA-F0-9]+)", content).group(1),
            "date_run": re.search(r"Date Run\s+([\d-]+\s[\d:.]+)", content).group(1),
            "network_used": re.search(r"Network Used\s+([^\n]+)", content).group(1).strip(),
            "phospho_type": re.search(r"Phospho Type\s+(\w+)", content).group(1),
            "kinase_size": int(re.search(r"Kinase Size\s+(\d+)", content).group(1)),
            "site_limit": int(re.search(r"Site Limit\s+(\d+)", content).group(1)),
            "num_networks": int(re.search(r"# of Networks\s+(\d+)", content).group(1)),
            "use_compendia": re.search(r"Use Compendia\s+(\w+)", content).group(1).lower() == "yes",
            "compendia_counts": list(map(int, re.findall(r"Compendia \d+\s+(\d+)", content))),
        }
    except FileNotFoundError:
        raise FileNotFoundError(f"RUN_INFORMATION.txt file not found at: {file_path}")
    except AttributeError as e:
        raise ValueError(f"Error parsing the RUN_INFORMATION.txt file: {e}")
    
def parse_network_information_json(network_directory):
    """
    Parse the RUN_INFORMATION.json file from network pruning run and extract its data.

    Args:
        network_directory (str): Path to the network directory containing RUN_INFORMATION.json.

    Returns:
        dict: A dictionary containing the parsed data.
    """
    file_path = os.path.join(network_directory, "RUN_INFORMATION.json")

    try:
        with open(file_path, 'r') as file:
            data = json.load(file)

        return data
    except FileNotFoundError:
        raise FileNotFoundError(f"RUN_INFORMATION.json file not found at: {file_path}")
    except json.JSONDecodeError as e:
        raise ValueError(f"Error parsing the RUN_INFORMATION.json file: {e}")

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
    
    #check if network directories exist
    if not os.path.isdir(NETWORK_DIR):
        print('Could not find network directory. Please update using config.update_network_directory().')
        ready_Y, ready_ST = False, False
    else:
        if not os.path.isdir(Y_NETWORK_DIR) and not os.path.isdir(ST_NETWORK_DIR):
            print(
                'Could not find network files. Please create network pickles using config.create_network_pickles() or change network directory to where .tsv network files are located with config.update_network_directory().')
            ready_Y, ready_ST = False, False
        elif not os.path.isdir(Y_NETWORK_DIR):
            print('Could not find tyrosine network files (no Y folder in network directory)')
            ready_Y = False
        elif not os.path.isdir(ST_NETWORK_DIR):
            print('Could not find serine/threonine network files (no ST folder in network directory)')
            ready_ST = False

    #check to make sure networks match reference proteome hash
    ####### insert code here #######

    #check pregenerated experiments
    if USE_PREGENERATED_RANDOM_ACTIVITIES:
        #check tyrosine networks for pregenerated experiments

        if not os.path.isdir(f"{Y_NETWORK_DIR}/RANDOM_ACTIVITIES/"):
            print('Warning: Could not find pregenerated experiments directory for tyrosine network. This should be located within the default network directory. Without this, KSTAR will generate random experiments on the fly, which will take longer and may not be desired.')
        if not os.path.isdir(f"{ST_NETWORK_DIR}/RANDOM_ACTIVITIES/"):
            print('Warning: Could not find pregenerated experiments directory for serine/threonine network. This should be located within the default network directory. Without this, KSTAR will generate random experiments on the fly, which will take longer and may not be desired.')


    if ready_Y and ready_ST:
        print('You are ready to begin! You can generate predictions for both Y and ST networks.')
    elif ready_Y:
        print(
            'You are ready to generate predictions for Y networks, but not ST networks. If you want to generate ST predictions, create the ST network pickle with config.create_network_pickles(phosphoType = ["ST"])')
    elif ready_ST:
        print(
            'You are ready to generate predictions for ST networks, but not Y networks. If you want to generate Y predictions, create the Y network pickle with config.create_network_pickles(phosphoType = ["Y"])')

