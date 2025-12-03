from functools import partial
from multiprocessing.pool import Pool
import concurrent.futures
import os, json
from kstar import config, calculate, helpers
from kstar.random_experiments import generate_random_experiments
import pandas as pd
import numpy as np
import re
from tqdm import tqdm

def check_resources(phospho_type, network_dir = None, network_name = None):
    """
    Check that the reference files and networks are matched (i.e. network is built from the same reference proteome that will be used for pregenerating random activities).
    """
    #load network hashes
    if network_dir is None and network_name is None:
        #use default network directory
        network_dir = config.NETWORK_DIR
        network_info = config.NETWORK_INFO[phospho_type]
    else:
        network_subdir = os.path.join(network_dir, phospho_type, network_name)
        if not os.path.exists(network_subdir):
            raise ValueError(f"Network directory for phosphotype={phospho_type} and network_name={network_name} does not exist ({network_subdir}). Please verify that it exists and contains the necessary network files.")
        
        network_info = helpers.parse_network_information(network_subdir, file_type = 'txt')

    #compare
    if config.REFERENCE_INFO['unique_reference_id'] != network_info['unique_reference_id']:
        raise ValueError(f"{phospho_type} network does not match reference proteome. Please rebuild networks with the current reference proteome or reinstall resource files.")



def generate_random_activities_from_scratch(exp_name, exp_size, compendia_sizes, phospho_type, networks = None, network_sizes = None, filtered_compendia = None, selection_type = 'KSTAR_NUM_COMPENDIA_CLASS'):
    """
    Given an experiment size and compendia sizes, generate a random experiment and calculate KSTAR activities on it.

    Parameters
    ----------
    exp_size : int
        Size of the random experiment to generate.
    compendia_sizes : dict or list
        Sizes of each compendia class to sample from [0 (low bias), 1 (medium bias), 2 (high bias)]. If a list is provided, it is assumed to be a list of sizes for each compendium indexed by number. Will automatically adjust if fractions are provided that sum to 1 (or percents that sum to 100).
    phospho_type : str
        'Y' or 'ST' for phosphotyrosine or phosphoserine/threonine experiments.
    exp_name : str
        Name to assign to the random experiment.
    networks : dict of pd.DataFrame, optional
        Preloaded kinase-substrate networks. If None, loads from default config.
    network_sizes : dict, optional
        Number of unique sites in each kinase-substrate network. If None, calculates from networks.
    filtered_compendia : dict of pd.DataFrame, optional
        Reference phosphoproteome split into distinct dataframes based on the number of compendia each site appears in. If None, loads from default config.
    selection_type : str, optional
        Type of selection to use for filtered compendia. Default is 'KSTAR_NUM_COMPENDIA_CLASS'. Only needed if filtered_compendia is not provided
    
    
    """
    #load necessary resources (networks and compendia)
    if networks is None:
        check_resources(phospho_type)
        networks = calculate.load_networks(phospho_type)

    if network_sizes is None:
        network_sizes = {}
        for name, net in networks.items():
            network_sizes[name] = net.drop_duplicates(subset = [config.KSTAR_ACCESSION, config.KSTAR_SITE]).shape[0]

    if filtered_compendia is None:
        filtered_compendia = calculate.getFilteredCompendia(phospho_type = phospho_type, selection_type = selection_type)
    #elif selection_type not in filtered_compendia:
    #    raise ValueError(f"selection_type {selection_type} not found in provided filtered_compendia.")

    ##convert list to dictionary for provided compendia sizes if needed
    if isinstance(compendia_sizes, list):
        compendia_sizes = {i: size for i, size in enumerate(compendia_sizes)}

    #make sure there is the right number of compendia sizes
    if len(compendia_sizes) != len(filtered_compendia.keys()):
        raise ValueError(f"Number of compendia sizes provided ({len(compendia_sizes)}) does not match number of compendia classes ({selection_type} has {filtered_compendia[selection_type].nunique()} classes).")

    #make sure entries sum to exp_size
    if sum(compendia_sizes.values()) == exp_size:
        compendia_sizes_processed = compendia_sizes
    elif sum(compendia_sizes.values()) == 1:
        compendia_sizes_processed = {i: int(exp_size * size) for i, size in compendia_sizes.items()}
    else:
        compendia_sizes_processed = {i: int(exp_size * size / sum(compendia_sizes.values())) for i, size in compendia_sizes.items()}

    #generate random experiment
    rand_experiment = generate_random_experiments.build_single_filtered_experiment(compendia_sizes_processed, filtered_compendia, exp_name, selection_type = 'KSTAR_NUM_COMPENDIA_CLASS')
    act = calculate.calculate_hypergeometric_activities(rand_experiment, networks, network_sizes, exp_name)

    return act

def generate_fpr_stats_from_scratch(rand_activities, num_random_experiments = 150):
    """
    Given a dataframe of hypergeometric activities for many random experiments, calculate FPR statistics for each kinase. For each calculation, one random experiment is held out as the "real" experiment and compared to the remaining random experiments with a Mann-Whitney U test.

    Parameters
    ----------
    rand_activities : pd.DataFrame
        Dataframe of random activities with columns ['network', 'KSTAR_KINASE', 'data', 'kinase_activity'] where 'data' indicates the random experiment name. Generated by generate_default_random_activities() and saved/loaded with save_random_activities()/load_random_activities().
    num_random_experiments : int
        Number of random experiments that should be present in rand_activities. Only used to verify input data.

    Returns
    -------
    rand_stats : dict
        Dictionary of FPR statistics for each kinase.

    """
    #check that random activities dataframe has expected number of experiments
    if rand_activities['data'].nunique() != num_random_experiments:
        raise ValueError(f"Random activities dataframe has {rand_activities['data'].nunique()} unique random experiments, but expected {num_random_experiments}. Please verify the input file.")
    
    #group by kinase and random experiment
    kinases = rand_activities['KSTAR_KINASE'].unique()
    rand_activities['sample'] = rand_activities['data'].apply(lambda x: ':'.join(x.split(':')[:-1])) #grab sample associated with random expeirment
    rand_activities['rand_exp_num'] = rand_activities['data'].apply(lambda x: int(x.split(':')[-1])) #grab random experiment number
    random_grouped = rand_activities.groupby(['sample', 'KSTAR_KINASE', 'rand_exp_num'])['kinase_activity'].agg(list)

    rand_stats = {}
    for kin in kinases:
        random_kinase_activity_array = np.vstack(random_grouped.loc['data:experiment', kin])
        rand_stats[kin] = calculate.calculate_fpr_Mann_Whitney(random_kinase_activity_array)
    return rand_stats

def save_fpr_stats(rand_stats, odir):
    """
    Convert FPR stats dictionary to dataframe and save to file.

    Parameters
    ----------
    rand_stats : dict
        Dictionary of FPR statistics for each kinase. Created by generate_fpr_stats_from_scratch().
    odir : str
        Output directory to save FPR stats file. In most cases, this is the same directory as the random activities file.
    """
    #reformat to dataframe
    rand_stats_df = pd.DataFrame.from_dict(rand_stats)
    rand_stats_df.to_csv(f"{odir}/fpr_stats.tsv", sep='\t', index = False)

def load_fpr_stats(odir):
    """
    Load FPR stats from file and convert to dictionary.

    Parameters
    ----------
    odir : str
        Output directory where FPR stats file is located. In most cases, this is the same directory as the random activities file.

    Returns
    -------
    rand_stats : dict
        Dictionary of FPR statistics for each kinase.
    """
    fpr_stats_path = f"{odir}/fpr_stats.tsv"
    if not os.path.exists(fpr_stats_path):
        raise FileNotFoundError(f"FPR stats file not found at: {fpr_stats_path}")
    rand_stats_df = pd.read_csv(fpr_stats_path, sep='\t')
    rand_stats = rand_stats_df.to_dict(orient='list', index = True)
    return rand_stats

def get_exp_sizes(min, max, frac_difference = 0.2):
    """
    Given a min and maximum dataset size that should have pregenerated random experiments, return a list of experiment sizes that cover all possible sizes with approximately the given fractional difference.

    Parameters
    ----------
    min : int
        Minimum experiment size.
    max : int
        Maximum experiment size.
    frac_difference : float
        Maximum allowed fractional difference between real and pregenerated random experiment sizes.

    Returns
    -------
    exp_sizes : list
        List of experiment sizes to pregenerate.
    """
    prev_value = min

    exp_sizes = []
    while prev_value < max:
        prev_value = prev_value + int(prev_value * frac_difference)
        exp_sizes.append(prev_value)
        prev_value = prev_value + int(prev_value * frac_difference)
    return exp_sizes

def process_exp(size, compendia_sizes, phospho_type, networks, network_sizes, filtered_compendia, num_random_experiments, save_dir):
    tmp_save_dir = f"{save_dir}/{size}/"
    if not os.path.exists(tmp_save_dir):
        os.makedirs(tmp_save_dir)
    #create partial function with fixed parameters
    func = partial(generate_random_activities_from_scratch, exp_size = size, compendia_sizes = compendia_sizes, phospho_type = phospho_type, networks = networks, network_sizes = network_sizes, filtered_compendia = filtered_compendia)
    #repeat for number of random experiments (each with a different exp_name)
    exp_names = [f"data:experiment:{i}" for i in range(num_random_experiments)]
    activities_list = []
    for activities in map(func, exp_names):
        activities_list.append(activities)
    activities = pd.concat(activities_list)
    save_random_activities(activities, tmp_save_dir)

def generate_default_random_activities_v2(phospho_type, compendia_sizes, network_dir = None, network_name = None, min_size = 50, max_size = 1000, frac_difference = 0.25, networks = None, network_sizes = None, regenerate = False, filtered_compendia = None, num_random_experiments = 150, selection_type = 'KSTAR_NUM_COMPENDIA_CLASS', PROCESSES = 1):
    """
    Generate default pregenerated random activities for a given phospho_type, compendia distribution, and range of experiment sizes. Will automatically determine the experiment sizes so that all sizes between min_size and max_size are covered with approximately the given fractional difference.

    Parameters
    ----------
    phospho_type : str
        'Y' or 'ST' for phosphotyrosine or phosphoserine/threonine experiments.
    compendia_sizes : list
        Sizes of each compendia class to sample from [0 (low bias), 1 (medium bias), 2 (high bias)]. Should be provided as a list of fractions (sum to 1) or percents (sum to 100). The length of the list should match the number of compendia classes in the filtered compendia.
    network_dir : str, optional
        Directory containing the kinase-substrate networks. If None, uses default from config.
    network_name : str, optional
        Name of the kinase-substrate network to use. If None, uses default from config.
    min_size : int
        Minimum experiment size you would like to have pregenerated random activities for.
    max_size : int
        Maximum experiment size you would like to have pregenerated random activities for.
    frac_difference : float
        Maximum allowed fractional difference between real and pregenerated random experiment sizes.
    networks : dict of pd.DataFrame, optional
        Preloaded kinase-substrate networks. If None, loads from default config. will override network_dir and network_name if provided.
    network_sizes : dict, optional
        Number of unique sites in each kinase-substrate network. If None, calculates from networks.
    regenerate : bool, optional
        If True, will regenerate random activities even if they already exist on disk. Default is False.
    filtered_compendia : dict of pd.DataFrame, optional
        Reference phosphoproteome split into distinct dataframes based on the number of compendia each site appears in. If None, loads from default config.
    num_random_experiments : int, optional
        Number of random experiments to generate for each experiment size. Default is 150.
    selection_type : str, optional
        Type of selection to use for filtered compendia. Default is 'KSTAR_NUM_COMPENDIA_CLASS'. Only needed if filtered_compendia is not provided
    PROCESSES : int, optional
        Number of processes to use for multiprocessing. Default is 1 (no multiprocessing).
    """
    if network_dir is None:
        network_dir = config.NETWORK_DIR
    if network_name is None:
        network_name = config.NETWORK_NAME[phospho_type]

    full_network_path = os.path.join(network_dir, phospho_type, network_name)
    if not os.path.exists(full_network_path):
        raise ValueError(f"Network directory for phosphotype={phospho_type} and network_name={network_name} does not exist ({full_network_path}). Please create it before running generate_default_random_activities.")
    
    #check if compendia are fractions or percents
    if sum(compendia_sizes) == 1:
        #convert to percents for directory naming
        compendia_sizes = [int(size * 100) for size in compendia_sizes]
    elif sum(compendia_sizes) != 100:
        raise ValueError(f"Compendia sizes must sum to 1 (fractions) or 100 (percents). Provided sizes sum to {sum(compendia_sizes)}.")

    #create directory to save random experiments
    save_dir = f"{full_network_path}/RANDOM_ACTIVITIES/compendia={'_'.join(map(str, compendia_sizes))}/"
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    #load necessary resources (networks and compendia)
    if networks is None:
        networks = calculate.load_networks(phospho_type, network_dir=network_dir, network_name=network_name)

    if network_sizes is None:
        network_sizes = {}
        for name, net in networks.items():
            network_sizes[name] = net.drop_duplicates(subset = [config.KSTAR_ACCESSION, config.KSTAR_SITE]).shape[0]

    if filtered_compendia is None:
        filtered_compendia = calculate.getFilteredCompendia(phospho_type = phospho_type, selection_type = selection_type)

    #determine experiment sizes to pregenerate
    exp_sizes = get_exp_sizes(min_size, max_size, frac_difference)
    #iterate through and remove sizes that already exist (if regenerate is False)
    if not regenerate:
        exp_sizes = [size for size in exp_sizes if not os.path.exists(f"{save_dir}/{size}/random_enrichment.tsv")]
    


    #generate random activities for each experiment size with the given compendia distribution
    activities_dict = {}
    #for size in tqdm(exp_sizes, desc = f'Generating random activities for pregenerated experiment sizes ({[str(e) for e in exp_sizes]})'):
    exp_args = [(size, compendia_sizes, phospho_type, networks, network_sizes, filtered_compendia, num_random_experiments, save_dir) for size in exp_sizes]
    with Pool(PROCESSES) as pool:
        results = tqdm(pool.starmap(process_exp, exp_args), total = len(exp_sizes))
    return results


def generate_default_random_activities(phospho_type, compendia_sizes, network_dir = None, network_name = None, min_size = 50, max_size = 1000, frac_difference = 0.25, networks = None, network_sizes = None, regenerate = False, filtered_compendia = None, num_random_experiments = 150, selection_type = 'KSTAR_NUM_COMPENDIA_CLASS', PROCESSES = 1):
    """
    Generate default pregenerated random activities for a given phospho_type, compendia distribution, and range of experiment sizes. Will automatically determine the experiment sizes so that all sizes between min_size and max_size are covered with approximately the given fractional difference.

    Parameters
    ----------
    phospho_type : str
        'Y' or 'ST' for phosphotyrosine or phosphoserine/threonine experiments.
    compendia_sizes : list
        Sizes of each compendia class to sample from [0 (low bias), 1 (medium bias), 2 (high bias)]. Should be provided as a list of fractions (sum to 1) or percents (sum to 100). The length of the list should match the number of compendia classes in the filtered compendia.
    network_dir : str, optional
        Directory containing the kinase-substrate networks. If None, uses default from config.
    network_name : str, optional
        Name of the kinase-substrate network to use. If None, uses default from config.
    min_size : int
        Minimum experiment size you would like to have pregenerated random activities for.
    max_size : int
        Maximum experiment size you would like to have pregenerated random activities for.
    frac_difference : float
        Maximum allowed fractional difference between real and pregenerated random experiment sizes.
    networks : dict of pd.DataFrame, optional
        Preloaded kinase-substrate networks. If None, loads from default config. will override network_dir and network_name if provided.
    network_sizes : dict, optional
        Number of unique sites in each kinase-substrate network. If None, calculates from networks.
    regenerate : bool, optional
        If True, will regenerate random activities even if they already exist on disk. Default is False.
    filtered_compendia : dict of pd.DataFrame, optional
        Reference phosphoproteome split into distinct dataframes based on the number of compendia each site appears in. If None, loads from default config.
    num_random_experiments : int, optional
        Number of random experiments to generate for each experiment size. Default is 150.
    selection_type : str, optional
        Type of selection to use for filtered compendia. Default is 'KSTAR_NUM_COMPENDIA_CLASS'. Only needed if filtered_compendia is not provided
    PROCESSES : int, optional
        Number of processes to use for multiprocessing. Default is 1 (no multiprocessing).
    """
    if network_dir is None:
        network_dir = config.NETWORK_DIR
    if network_name is None:
        network_name = config.NETWORK_NAME[phospho_type]

    full_network_path = os.path.join(network_dir, phospho_type, network_name)
    if not os.path.exists(full_network_path):
        raise ValueError(f"Network directory for phosphotype={phospho_type} and network_name={network_name} does not exist ({full_network_path}). Please create it before running generate_default_random_activities.")
    
    #check if compendia are fractions or percents
    if sum(compendia_sizes) == 1:
        #convert to percents for directory naming
        compendia_sizes = [int(size * 100) for size in compendia_sizes]
    elif sum(compendia_sizes) != 100:
        raise ValueError(f"Compendia sizes must sum to 1 (fractions) or 100 (percents). Provided sizes sum to {sum(compendia_sizes)}.")

    #create directory to save random experiments
    save_dir = f"{full_network_path}/RANDOM_ACTIVITIES/compendia={'_'.join(map(str, compendia_sizes))}/"
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    #load necessary resources (networks and compendia)
    if networks is None:
        networks = calculate.load_networks(phospho_type, network_dir=network_dir, network_name=network_name)

    if network_sizes is None:
        network_sizes = {}
        for name, net in networks.items():
            network_sizes[name] = net.drop_duplicates(subset = [config.KSTAR_ACCESSION, config.KSTAR_SITE]).shape[0]

    if filtered_compendia is None:
        filtered_compendia = calculate.getFilteredCompendia(phospho_type = phospho_type, selection_type = selection_type)

    #determine experiment sizes to pregenerate
    exp_sizes = get_exp_sizes(min_size, max_size, frac_difference)
    
    #generate random activities for each experiment size with the given compendia distribution
    activities_dict = {}
    #for size in tqdm(exp_sizes, desc = f'Generating random activities for pregenerated experiment sizes ({[str(e) for e in exp_sizes]})'):
    for size in exp_sizes:
        activities_list = []
        #create new directory for this file size
        tmp_save_dir = f"{save_dir}/{size}/"
        if not os.path.exists(tmp_save_dir):
            os.makedirs(tmp_save_dir)
        #if file doesn't already exist or regenerate is True, generate random activities
        if regenerate or not os.path.exists(f"{tmp_save_dir}/random_enrichment.tsv"):
            if PROCESSES > 1:
                with Pool(processes=PROCESSES) as pool:
                    #create partial function with fixed parameters
                    func = partial(generate_random_activities_from_scratch, exp_size = size, compendia_sizes = compendia_sizes, phospho_type = phospho_type, networks = networks, network_sizes = network_sizes, filtered_compendia = filtered_compendia)
                    #repeat for number of random experiments (each with a different exp_name)
                    exp_names = [f"data:experiment:{i}" for i in range(num_random_experiments)]
                    activities_list = pool.map(func, exp_names)
                    
                    #collect results
                    #for act, _ in results:
                    #    activities_list.append(act)
            else:
                for i in tqdm(range(num_random_experiments)):
                    exp_name = f"data:experiment:{i}"
                    act = generate_random_activities_from_scratch(exp_name = exp_name, exp_size = size, compendia_sizes = compendia_sizes, phospho_type = phospho_type, networks = networks, network_sizes = network_sizes, filtered_compendia = filtered_compendia)
                    activities_list.append(act)

            activities = pd.concat(activities_list)
            save_random_activities(activities, tmp_save_dir)
        #else:
        #    print(f"Random activities for size {size} already exist and regenerate is set to False. Skipping generation.")

    #return activities_dict
def generate_default_fpr_stats(phospho_type, network_dir = None, network_name = None, regenerate = False, PROCESSES = 1):
    """
    Generate FPR statistics for all pregenerated random activities for a given phospho_type and network. Will look for existing pregenerated random activities in the default directory structure and generate FPR stats for each experiment size if they do not already exist (or if regenerate is set to True).

    Parameters
    ----------
    phospho_type : str
        'Y' or 'ST' for phosphotyrosine or phosphoserine/threonine experiments.
    network_dir : str, optional
        Directory containing the kinase-substrate networks. If None, uses default from config.
    network_name : str, optional
        Name of the KSTAR networks to use. If None, uses default from config.
    regenerate : bool, optional
        If True, will regenerate FPR stats even if they already exist on disk. Default is False.
    PROCESSES : int, optional
        Number of processes to use for multiprocessing. Default is 1 (no multiprocessing). Not currently used in this function.
    

    """
    if network_dir is None:
        network_dir = config.NETWORK_DIR
    if network_name is None:
        network_name = config.NETWORK_NAME[phospho_type]

    full_network_path = os.path.join(network_dir, phospho_type, network_name)
    if not os.path.exists(full_network_path):
        raise ValueError(f"Network directory for phosphotype={phospho_type} and network_name={network_name} does not exist ({full_network_path}). Please create it before running generate_default_random_activities.")
    if not os.path.exists(f"{full_network_path}/RANDOM_ACTIVITIES/"):
        raise ValueError(f"Random activities directory does not exist at {full_network_path}/RANDOM_ACTIVITIES/. Please generate random activities before generating FPR stats.")
    
    #iterate through existing random experiments
    compendia_folders = os.listdir(f"{full_network_path}/RANDOM_ACTIVITIES/")
    print(f"Generating FPR stats for phospho_type={phospho_type}, network_name={network_name} for compendia distributions: {compendia_folders}")
    for comp_folder in compendia_folders:
        size_folders = os.listdir(f"{full_network_path}/RANDOM_ACTIVITIES/{comp_folder}/")


        for size_folder in tqdm(size_folders, desc=f'Generating FPR stats for pregenerated experiment sizes ({comp_folder})'):
            size_folder_path = os.path.join(full_network_path, "RANDOM_ACTIVITIES",comp_folder, size_folder)
            if not os.path.exists(f"{size_folder_path}/random_enrichment.tsv"):
                print(f"Random activities for experiment size {size_folder} do not exist at {size_folder_path}/random_enrichment.tsv. Skipping FPR stats generation.")
            elif os.path.exists(f"{size_folder_path}/fpr_stats.tsv") and not regenerate:
                print(f"FPR stats for experiment size {size_folder} already exist and regenerate is set to False. Skipping generation.")
            else:
                rand_activities = load_random_activities(f"{size_folder_path}/random_enrichment.tsv")
                rand_stats = generate_fpr_stats_from_scratch(rand_activities)
                save_fpr_stats(rand_stats, size_folder_path)




def save_random_activities(activities, odir):
    """
    Given the long-form dataframe of random activities, reformat to reduce memory usage and save to file.

    Parameters
    ----------
    activities : pd.DataFrame
        Dataframe of random activities with columns ['network', 'KSTAR_KINASE', 'data', 'kinase_activity'] where 'data' indicates the random experiment name. Generated by generate_default_random_activities() and saved/loaded with save_random_activities()/load_random_activities().
    odir : str
        Output directory to save random activities file.


    """
    #reformat to reduce memory usage
    tmp_pivot = activities.pivot(columns='data', index=['network', 'KSTAR_KINASE'], values='kinase_activity').reset_index()
    #save to file
    tmp_pivot.to_csv(f"{odir}/random_enrichment.tsv", sep='\t', index = False)

def load_random_activities(file):
    """
    Load random activities from file and convert to long-form dataframe that is expected by KSTAR
    """
    rand_dataset_activities = pd.read_csv(file, delimiter='\t')
    rand_dataset_activities = rand_dataset_activities.melt(id_vars=['network', 'KSTAR_KINASE'], var_name='data', value_name='kinase_activity')
    return rand_dataset_activities


def save_random_activities_dict(activities_dict, odir):
    """
    Given a dictionary of long-form dataframes of random activities with sizes as keys, reformat to reduce memory usage and save to file.

    Parameters
    ----------
    activities_dict : dict
        Dictionary of dataframes of random activities with keys as experiment sizes and values as dataframes with columns ['network', 'KSTAR_KINASE', 'data', 'kinase_activity'] where 'data' indicates the random experiment name. Generated by generate_default_random_activities() and saved/loaded with save_random_activities()/load_random_activities().
    odir : str
        Output directory to save random activities files. Should correspond to the directory for the given phospho_type and compendia distribution.
    """
    for size, activities in activities_dict.items():
        #reformat to reduce memory usage
        tmp_pivot = activities.pivot(columns='data', index=['network', 'KSTAR_KINASE'], values='kinase_activity').reset_index()

        #make output directory if it doesn't exist
        if not os.path.exists(f"{odir}/{size}"):
            os.makedirs(f"{odir}/{size}")
        #save to file
        tmp_pivot.to_csv(f"{odir}/{size}/random_enrichment.tsv", sep='\t', index = False)


def main():
    #parse command-line arguments
    import argparse
    parser = argparse.ArgumentParser(description='Precalculate random KSTAR activities for default experiment sizes and compendia distributions.')
    parser.add_argument('--network_dir', type=str, default=None, help='Directory containing the kinase-substrate networks. If not provided, uses default from config.')
    parser.add_argument('--network_name', type=str, default=None, help='Name of the kinase-substrate network to use. If not provided, uses default from config.')
    parser.add_argument('--regenerate', action='store_true', help='If set, will regenerate pregenerated random activities and FPR stats even if they already exist on disk.')
    parser.add_argument('-p', type=int, default=1, help='Number of processes to use for multiprocessing. Default is 1 (no multiprocessing).')
    args = parser.parse_args()

    #construct default random activities for tyrosine
    min_size = 50
    max_size = 1000
    frac_difference = 0.25
    compendia_sizes = [[0,50,50], [0, 30, 70]]  # 0: low bias, 1: medium bias, 2: high bias
    for comp_size in compendia_sizes:
        activities_dict = generate_default_random_activities(phospho_type = 'Y', compendia_sizes = comp_size, min_size = min_size, max_size = max_size, frac_difference = frac_difference, num_random_experiments = 150, network_dir = args.network_dir, network_name=args.network_name, PROCESSES = args.p)



    #construct default random activities for serine/threonine
    min_size = 500
    max_size = 20000
    frac_difference = 0.25
    compendia_sizes = [[0,30,70]]
    for comp_size in compendia_sizes:
        activities_dict = generate_default_random_activities(phospho_type = 'ST', compendia_sizes = comp_size, min_size = min_size, max_size = max_size, frac_difference = frac_difference, num_random_experiments = 150, network_dir = args.network_dir, network_name=args.network_name, PROCESSES = args.p)

    print('Finished generating default pregenerated random activities. Now generating FPR statistics for all pregenerated random activities...')

    #generate FPR stats for tyrosine
    generate_default_fpr_stats(phospho_type = 'Y', network_dir = args.network_dir, network_name=args.network_name, PROCESSES = args.p)
    #generate FPR stats for serine/threonine
    generate_default_fpr_stats(phospho_type = 'ST', network_dir = args.network_dir, network_name=args.network_name, PROCESSES = args.p)
