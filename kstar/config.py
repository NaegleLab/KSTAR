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

#read information about reference files
with open(f"{RESOURCE_DIR}/reference_info.json", "r") as f:
    REFERENCE_INFO = json.load(f)


###### Customizable configuration ####################
configuration_file = f'{KSTAR_DIR}/configuration.json'
with open(configuration_file, 'r') as cf:
    config_data = json.load(cf)

#extract network directory from configuration file, make sure it can be found. By default this will be located in KSTAR_DIR/NETWORKS/
NETWORK_DIR = config_data["network_directory"]
if NETWORK_DIR is None or NETWORK_DIR == "":
    NETWORK_DIR = f"{KSTAR_DIR}/NETWORKS/"
#look for tyrosine and serine/threonine networks by name, then load network information
NETWORK_NAME = {}
NETWORK_NAME['Y'] = config_data.get('tyrosine_network_name', 'Default')
NETWORK_NAME['ST'] = config_data.get('serine_threonine_network_name', 'Default')


DEFAULT_RANDOM_ACTIVITIES_DIR = NETWORK_DIR
NETWORK_SUBDIR = {'Y': f"{NETWORK_DIR}/Y/{NETWORK_NAME['Y']}/",
                  'ST': f"{NETWORK_DIR}/ST/{NETWORK_NAME['ST']}/"}


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




######## Pregeneration parameters ################
CUSTOM_RANDOM_ACTIVITIES_DIR = config_data["custom_pregenerated_experiments_dir"]
USE_PREGENERATED_RANDOM_ACTIVITIES = config_data['use_pregenerated_random_activities']
SAVE_NEW_RANDOM_ACTIVITIES = config_data['save_new_random_activities']
SAVE_RANDOM_EXPERIMENTS = config_data['save_random_experiments']

######### static configuration (user cannot change) ################
RESOURCE_URL = 'https://figshare.com/ndownloader/files/60883930'  # location of corresponding release of data
NETWORK_URL = 'https://figshare.com/ndownloader/files/60883384'



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
    NETWORK_DIR = f"{KSTAR_DIR}/NETWORKS/"


DEFAULT_RANDOM_ACTIVITIES_DIR = {
    'Y': f"{NETWORK_DIR}/Y/{NETWORK_NAME['Y']}/RANDOM_ACTIVITIES/",
    'ST': f"{NETWORK_DIR}/ST/{NETWORK_NAME['ST']}/RANDOM_ACTIVITIES/"
}



## END DECLARATION OF GLOBALS


## USER FUNCTIONS TO UPDATE CONFIGURATION GLOBALS
def update_network_directory(network_dir = None, y_network_name = None, st_network_name = None):
    """
    Update the location of network the network files, and verify that all necessary files are located in directory

    Parameters
    ----------
    network_dir: string
        path to where network files are located
    y_network_name: string
        name of the tyrosine network to use
    st_network_name: string
        name of the serine/threonine network to use

    """
    global NETWORK_DIR
    global NETWORK_SUBDIR
    global NETWORK_NAME
    global NETWORK_INFO
    if network_dir is not None:
        NETWORK_DIR = network_dir
        #check to make sure directory exists
        if not os.path.isdir(network_dir):
            raise ValueError(f"Provided network directory {network_dir} does not exist. Please provide a valid directory.")

        if y_network_name is None:
            #check to make sure subdirectories exist
            NETWORK_SUBDIR['Y'] = f"{NETWORK_DIR}/Y/{NETWORK_NAME['Y']}/"
            if not os.path.exists(f"{NETWORK_SUBDIR['Y']}/RUN_INFORMATION.txt"):
                raise ValueError(f"Tyrosine network directory {NETWORK_SUBDIR['Y']} does not exist. Please provide a valid network name for the tyrosine network you would like to use and can be found in provided network directory.")
            NETWORK_SUBDIR['ST'] = f"{NETWORK_DIR}/ST/{NETWORK_NAME['ST']}/"

        if st_network_name is None:
            NETWORK_SUBDIR['ST'] = f"{NETWORK_DIR}/ST/{NETWORK_NAME['ST']}/"
            if not os.path.exists(f"{NETWORK_SUBDIR['ST']}/RUN_INFORMATION.txt"):
                raise ValueError(f"Serine/threonine network directory {NETWORK_SUBDIR['ST']} does not exist. Please provide a valid network name for the serine/threonine network you would like to use and can be found in provided network directory.")

    #update specifc network directories if needed
    if y_network_name is not None:
        NETWORK_NAME['Y'] = y_network_name
        #Y_NETWORK_DIR = f"{NETWORK_DIR}/Y/{NETWORK_HASH_Y}/"
        NETWORK_SUBDIR['Y'] = f"{NETWORK_DIR}/Y/{NETWORK_NAME['Y']}/"
        if not os.path.exists(f"{NETWORK_SUBDIR['Y']}/RUN_INFORMATION.txt"):
            raise ValueError(f"Tyrosine network directory {NETWORK_SUBDIR['Y']} does not exist. Please provide a valid network name for the tyrosine network you would like to use and can be found in provided network directory.")



    if st_network_name is not None:
        NETWORK_NAME['ST'] = st_network_name
        #ST_NETWORK_DIR = f"{NETWORK_DIR}/ST/{NETWORK_HASH_ST}/"
        NETWORK_SUBDIR['ST'] = f"{NETWORK_DIR}/ST/{NETWORK_NAME['ST']}/"
        if not os.path.exists(f"{NETWORK_SUBDIR['ST']}/RUN_INFORMATION.txt"):
            raise ValueError(f"Serine/threonine network directory {NETWORK_SUBDIR['ST']} does not exist. Please provide a valid network name for the serine/threonine network you would like to use and can be found in provided network directory.")
    
    #read in network info
    NETWORK_INFO['ST'] = helpers.parse_network_information(f"{NETWORK_DIR}/ST/{NETWORK_NAME['ST']}/")
    NETWORK_INFO['Y'] = helpers.parse_network_information(f"{NETWORK_DIR}/Y/{NETWORK_NAME['Y']}/")


def update_configuration(network_dir=None, y_network_name=None, st_network_name=None, save_random_experiments = None, use_pregenerated_random_activities=None, save_new_random_activities=None, custom_pregenerated_activities_dir=None):
    """
    Update configuration parameters in current iteration and save to configuration file.

    Parameters
    ----------
    use_pregenerated_random_activities : bool, optional
        Whether to use pregenerated random activities when possible, by default None
    save_new_random_activities : bool, optional
        Whether to save new random activities when they are generated, by default False
    custom_pregenerated_activities_dir : str, optional
        Directory to save newly generated random activities for future use, by default None
    network_dir : str, optional
        Directory containing the kinase-substrate networks, by default None (which assumes it is located in kstar directory)
    y_network_hash : str, optional
        Unique identifier of the tyrosine network to use by default. 
    st_network_hash : str, optional
        Unique identifier of the serine/threonine network to use by default.
    """
    global SAVE_RANDOM_EXPERIMENTS
    global CUSTOM_RANDOM_ACTIVITIES_DIR
    global SAVE_NEW_RANDOM_ACTIVITIES
    global USE_PREGENERATED_RANDOM_ACTIVITIES
    if use_pregenerated_random_activities is not None:
        USE_PREGENERATED_RANDOM_ACTIVITIES = use_pregenerated_random_activities
    if save_new_random_activities is not None:
        SAVE_NEW_RANDOM_ACTIVITIES = save_new_random_activities
    if custom_pregenerated_activities_dir is not None:
        if not os.path.isdir(custom_pregenerated_activities_dir):
            raise ValueError(f"Provided pregenerated activities directory {custom_pregenerated_activities_dir} does not exist. Please provide a valid directory.")
        CUSTOM_RANDOM_ACTIVITIES_DIR = custom_pregenerated_activities_dir
    if save_random_experiments is not None:
        SAVE_RANDOM_EXPERIMENTS = save_random_experiments

    #update network directory and names if provided
    update_network_directory(network_dir=network_dir, y_network_name=y_network_name, st_network_name=st_network_name)


    #save changes to configuration file
    config_data = {
        "network_directory": NETWORK_DIR,
        "Y_network_name": NETWORK_NAME['Y'],
        "ST_network_name": NETWORK_NAME['ST'],
        'save_random_experiments': SAVE_RANDOM_EXPERIMENTS,
        "use_pregenerated_random_activities": USE_PREGENERATED_RANDOM_ACTIVITIES,
        "save_new_random_activities": SAVE_NEW_RANDOM_ACTIVITIES,
        "custom_pregenerated_experiments_dir": CUSTOM_RANDOM_ACTIVITIES_DIR,
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




def install_network_files(target_dir=None):
    """Retrieves Network files that are the companion for this version release from FigShare, unzips them to the specified directory.
    
    Parameters
    ----------
    target_dir : str, optional
        Directory to install network files to. If None, defaults to within package location ({KSTAR_DIR}/NETWORKS/)

    """
    print("Requesting network file")
    r = requests.get(NETWORK_URL, allow_redirects=True)
    if r.status_code != 200:
        raise RuntimeError(f"Failed to download network file: HTTP {r.status_code}")

    # Use target_dir if provided, otherwise default to KSTAR_DIR
    install_dir = target_dir if target_dir else KSTAR_DIR
    outputFile = os.path.join(install_dir, "NETWORKS.tar.gz")

    # Create target directory if it doesn't exist
    os.makedirs(install_dir, exist_ok=True)

    with open(outputFile, 'wb') as f:
        f.write(r.content)

    # extract tarball
    with tarfile.open(outputFile, 'r:gz') as t:
        print(f"Extracting {outputFile}")
        t.extractall(install_dir)

    # remove tar file
    os.remove(outputFile)

    #update network directory
    update_configuration(network_dir=f"{install_dir}/NETWORKS/", y_network_name='Default', st_network_name='Default')



def find_available_networks(phospho_type):
    """
    Find available network hashes in the current network directory, and return dictionary with information about them

    Returns
    -------
    available_networks : dict
        dictionary containing all available networks, in the format -> Network hash : network information dictionary
    """
    available_networks = {}
    found_but_invalid = {}
    for item in os.listdir(NETWORK_DIR + f'/{phospho_type}/'):
        item_path = os.path.join(NETWORK_DIR, phospho_type, item)
        #check if RUN_INFORMATION.txt file exists
        if not os.path.isfile(os.path.join(item_path, "RUN_INFORMATION.txt")):
            found_but_invalid[item] = "Could not find RUN_INFORMATION.txt file"
        else:
            #read network information
            network_info = helpers.parse_network_information(item_path)
            #make sure network matches reference proteome
            if network_info['unique_reference_id'] != REFERENCE_INFO['unique_reference_id']:
                found_but_invalid[item] = "Network does not match reference proteome used in this KSTAR installation"
            elif network_info['phospho_type'] != phospho_type:
                found_but_invalid[item] = f"Network phospho type {network_info['phospho_type']} does not match requested phospho type {phospho_type}"
            elif not os.path.exists(os.path.join(item_path, f"INDIVIDUAL_NETWORKS/")):
                found_but_invalid[item] = "Could not find INDIVIDUAL_NETWORKS directory within network folder"
            elif not os.path.exists(os.path.join(item_path, f"RANDOM_ACTIVITIES/")) and USE_PREGENERATED_RANDOM_ACTIVITIES:
                found_but_invalid[item] = "Could not find RANDOM_ACTIVITIES directory within network folder"
            else:
                available_networks[item] = network_info
    return available_networks, found_but_invalid


def get_package_memory():
    #compute size of resource files
    reference_size = 0
    for root, dirs, files in os.walk(RESOURCE_DIR):
        for file in files:
            path = os.path.join(root, file)
            if os.path.isfile(path):
                reference_size += os.path.getsize(path)

    #network and pregenerated activities size    
    total_size = {'Y': 0, 'ST': 0}
    network_size = {'Y': 0, 'ST': 0}
    rand_activity_size = {'Y': 0, 'ST': 0}
    custom_activity_size = {'Y': 0, 'ST': 0}
    for phospho_type in ['Y', 'ST']:
        for root, dirs, files in os.walk(NETWORK_SUBDIR[phospho_type]):
            for file in files:
                path = os.path.join(root, file)
                if os.path.isfile(path):
                    total_size[phospho_type] += os.path.getsize(path)
                    
                    if 'INDIVIDUAL_NETWORKS' in root and file.endswith('.tsv'):
                        network_size[phospho_type] += os.path.getsize(path)
                    elif 'RANDOM_ACTIVITIES' in root and file.endswith('.tsv'):
                        rand_activity_size[phospho_type] += os.path.getsize(path)
            if CUSTOM_RANDOM_ACTIVITIES_DIR is not None:
                custom_dir = os.path.join(CUSTOM_RANDOM_ACTIVITIES_DIR, phospho_type, NETWORK_NAME[phospho_type])
                if os.path.isdir(custom_dir):
                    for root, dirs, files in os.walk(custom_dir):
                        for file in files:
                            path = os.path.join(root, file)
                            if os.path.isfile(path) and file.endswith('.tsv'):
                                custom_activity_size[phospho_type] += os.path.getsize(path)

    

    #report memory in GB
    print(f"Total reference proteome size: {reference_size / (1024**3):.2f} GB\n")
    print('Tyrosine Network Directory:')
    print(f"Total network directory size: {total_size['Y'] / (1024**3):.2f} GB")
    print(f"Total .tsv network size: {network_size['Y'] / (1024**3):.2f} GB")
    print(f"Total default random activity size: {rand_activity_size['Y'] / (1024**3):.2f} GB")
    print(f"Total custom random activity size: {custom_activity_size['Y'] / (1024**3):.2f} GB")

    print('\nSerine/Threonine Network Directory:')
    print(f"Total network directory size: {total_size['ST'] / (1024 **3):.2f} GB")
    print(f"Total .tsv network size: {network_size['ST'] / (1024 **3):.2f} GB")
    print(f"Total default random activity size: {rand_activity_size['ST'] / (1024 **3):.2f} GB")
    print(f"Total custom random activity size: {custom_activity_size['ST'] / (1024 **3):.2f} GB")






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
        if not os.path.isdir(NETWORK_SUBDIR['Y']) and not os.path.isdir(NETWORK_SUBDIR['ST']):
            print(
                'Could not find network files. Please create network pickles using config.create_network_pickles() or change network directory to where .tsv network files are located with config.update_network_directory().')
            ready_Y, ready_ST = False, False
        elif not os.path.isdir(NETWORK_SUBDIR['Y']):
            print('Could not find tyrosine network files (no Y folder in network directory)')
            ready_Y = False

            #report how many networks are available
        elif not os.path.isdir(NETWORK_SUBDIR['ST']):
            print('Could not find serine/threonine network files (no ST folder in network directory)')
            ready_ST = False

    #check to make sure networks match reference proteome hash
    ####### insert code here #######

    #check pregenerated experiments
    if USE_PREGENERATED_RANDOM_ACTIVITIES:
        #check tyrosine networks for pregenerated experiments

        if not os.path.isdir(f"{NETWORK_SUBDIR['Y']}/RANDOM_ACTIVITIES/"):
            print('Warning: Could not find pregenerated experiments directory for tyrosine network. This should be located within the default network directory. Without this, KSTAR will generate random experiments on the fly, which will take longer and may not be desired.')
        if not os.path.isdir(f"{NETWORK_SUBDIR['ST']}/RANDOM_ACTIVITIES/"):
            print('Warning: Could not find pregenerated experiments directory for serine/threonine network. This should be located within the default network directory. Without this, KSTAR will generate random experiments on the fly, which will take longer and may not be desired.')


    if ready_Y and ready_ST:
        print('You are ready to begin! You can generate predictions for both Y and ST networks.')
    elif ready_Y:
        print(
            'You are ready to generate predictions for Y networks, but not ST networks. If you want to generate ST predictions, create the ST network pickle with config.create_network_pickles(phosphoType = ["ST"])')
    elif ready_ST:
        print(
            'You are ready to generate predictions for ST networks, but not Y networks. If you want to generate Y predictions, create the Y network pickle with config.create_network_pickles(phosphoType = ["Y"])')

