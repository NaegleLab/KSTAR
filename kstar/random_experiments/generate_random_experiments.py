import pandas as pd
import multiprocessing
import itertools
import argparse
import os
from kstar import config

#%%
#%%
def build_filtered_experiment(experiment, compendia, filtered_compendia, num_random_experiments, name ,selection_type='KSTAR_NUM_COMPENDIA_CLASS'):
    rand_experiments = compendia[[config.KSTAR_ACCESSION, config.KSTAR_SITE]]
    if len(experiment) == 0:
        empty_columns = [f"{name}:{i}" for i in range(num_random_experiments)]
        rand_experiments = pd.concat([rand_experiments,pd.DataFrame(columns=empty_columns)])
        return rand_experiments


    sizes = experiment.groupby(selection_type).size()
    for i in range(num_random_experiments):
        rand_experiment_list = []
        for num, size in sizes.items():
            filtered = filtered_compendia[int(num)]
            filtered_random = filtered.sample(size)
            filtered_random[f"{name}:{i}"] = 1
            rand_experiment_list.append(filtered_random)
        rand_experiment = pd.concat(rand_experiment_list)
        rand_experiments = pd.merge(rand_experiments, rand_experiment, how = 'left', on = [config.KSTAR_ACCESSION, config.KSTAR_SITE])
    #rand_experiments.to_csv( f"{name}_random_experiments.tsv", index=False, sep='\t')
    return rand_experiments


def build_random_experiments(binary_evidence, compendia, num_random_experiments, phosphorylation_event, data_columns, selection_type='KSTAR_NUM_COMPENDIA_CLASS', PROCESSES = 1):
    """
    Given an experimental dataframe and the human phospho compendia, build random experiments such that each random experiment takes on the same
    distribution with respect to the study bias defined as either NUM_COMPENDIA (total number of compendia a site is annotated in) or 
    NUM_COMPENDIA_CLASS (whether it has low < 1, medium (1-3), or high study bias(>3)).
    
    Parameters
    ----------
    binary_evidence: df
        KSTAR mapped experimental dataframe that has been binarized by kstar_activity generation
    greater: Boolean
        whether to keep sites that have a numerical value >=threshold (TRUE, default) or <=threshold (FALSE)
    num_random_experiments: int
        Number of random experiments to generate for each data_column
    phosphorylation_event: {'Y', 'ST'}
        Which substrate/kinaset-type to generate random experiments for
    data_columns : list
        columns that represent experimental result
    selection_type: {'KSTAR_NUM_COMPENDIA', 'KSTAR_NUM_COMPENDIA_CLASS'}
        Whether to sample according to the absolute number of compendia or low, medium, or high study bias groupings

    Returns
    -------
    rand_experiments: pandas.DataFrame
        Dataframe of random experiments with NaN where phosphorylation sites are not selected, and 1 if they are for that experiment
    """
    #check parameters
    if selection_type != 'KSTAR_NUM_COMPENDIA':
        if selection_type != 'KSTAR_NUM_COMPENDIA_CLASS':
            raise ValueError('selection_type must be either KSTAR_NUM_COMPENDIA or KSTAR_NUM_COMPENDIA_CLASS')

    experiment = binary_evidence
    if phosphorylation_event == 'ST':
        compendia = compendia[(compendia.KSTAR_SITE.str.contains('S')) | (compendia.KSTAR_SITE.str.contains('T'))]
        experiment = experiment[(experiment.KSTAR_SITE.str.contains('S')) | (experiment.KSTAR_SITE.str.contains('T'))]
            
    elif phosphorylation_event == 'Y':
        compendia = compendia[(compendia.KSTAR_SITE.str.contains('Y'))]
        experiment = experiment[(experiment.KSTAR_SITE.str.contains('Y'))]
    else:
        raise ValueError('phosphorylation_event must be Y or ST')


    compendia = compendia[[config.KSTAR_ACCESSION, config.KSTAR_SITE, selection_type]]
    compendia = compendia.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE]).max().reset_index() #uniquify the compendia by KSTAR_ACCESSION and KSTAR_SITE
    

    #experiment = experiment.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE]).agg(agg).reset_index()
    
    sizes = compendia[selection_type].unique()
    filtered_compendia = {}
    for s in sizes:
        filtered_compendia[s] = compendia[compendia[selection_type] == s][[config.KSTAR_ACCESSION, config.KSTAR_SITE]]
    
    
    filtered_experiments  = [experiment[experiment[col] ==1 ] for col in data_columns] 


    # ************ PARALELLIZATION ************
    if PROCESSES > 1:
        pool = multiprocessing.Pool(processes = PROCESSES)
        iterable = zip(
                filtered_experiments, 
                itertools.repeat(compendia), 
                itertools.repeat(filtered_compendia), 
                itertools.repeat(num_random_experiments), 
                [col for col in data_columns] ) 
        rand_experiments_list =  pool.starmap(build_filtered_experiment, iterable)

    # ********** NO PARALLELIZATION ***********
    else:
        rand_experiments_list = []
        for experiment, data_column in zip(filtered_experiments, data_columns):
            rand_exp = build_filtered_experiment(experiment, compendia, filtered_compendia, num_random_experiments, data_column)
            rand_experiments_list.append(rand_exp)

    rand_experiments = compendia[[config.KSTAR_ACCESSION, config.KSTAR_SITE]]
    for r in rand_experiments_list:
        rand_experiments = pd.merge(rand_experiments, r, how = 'left', on = [config.KSTAR_ACCESSION, config.KSTAR_SITE])
        
    #remove sites that aren't used in any random experiments
    random_data_cols = [col for col in rand_experiments if col not in [config.KSTAR_SITE,config.KSTAR_ACCESSION]]
    rand_experiments = rand_experiments.dropna(how = 'all', subset = random_data_cols)
    return rand_experiments
        
##### Faster single experiment tests ####
def build_single_filtered_experiment(compendia_sizes, filtered_compendia, name ,selection_type='KSTAR_NUM_COMPENDIA_CLASS'):
    #rand_experiments = compendia[[config.KSTAR_ACCESSION, config.KSTAR_SITE]]
    #if len(experiment) == 0:
    #    empty_columns = [f"{name}:{i}" for i in range(num_random_experiments)]
    #    rand_experiments = pd.concat([rand_experiments,pd.DataFrame(columns=empty_columns)])
    #    return rand_experiments

    #sizes = experiment.groupby(selection_type).size()
    rand_experiment_list = []
    for num, size in compendia_sizes.items():
        filtered = filtered_compendia[int(num)]
        filtered_random = filtered.sample(size)
        filtered_random["Experiment"] = name
        rand_experiment_list.append(filtered_random)
    rand_experiment = pd.concat(rand_experiment_list)
    return rand_experiment
        
def parse_args():
    parser = argparse.ArgumentParser(description='Parse Mapping Inference Arguments')
    parser.add_argument('-e', '--exp_file', '--experiment_file', action='store', dest= 'exp_file', help='Experiment file location. csv or tsv file', required=True)
    parser.add_argument('-r','--rdir','--resource_directory', action='store',dest='rdir', help = 'resource file directory', required=True)
    parser.add_argument('-o','--odir','--output_directory', action = 'store', dest='odir', help = 'output directory for results',)
    parser.add_argument('-p','--pevent','--phosphorylation_events', action = 'store', dest='pevent', help ='phosphorylation event type', choices=['Y','S','T','ST', 'STY'], default='STY')
    parser.add_argument('-n', '--name', action = 'store', dest='name', help = 'experiment name', default='Experiment')
    parser.add_argument('-d', '--data', '--data_columns', action='store', dest='data_columns', help = 'data_columns to use', nargs='*')
    parser.add_argument('-a', '--arg', '--activity_agg', action='store', dest='agg', help = 'activity agg to use', default='count', choices =['count','mean'])
    parser.add_argument('-t', '--threshold',  action='store', dest='threshold', help = 'threshold to use for analysis', type = float, default=0.0)
    parser.add_argument('-num', '--num_random_experiments',  action='store', dest='num', help = 'Number of random experiments to generate', type = int, default=150)
    parser.add_argument('-s', '--selection_type', action='store', dest='selection_type', help='KSTAR_NUM_COMPENDIA or KSTAR_NUM_COMPENDIA_CLASS', type=str, default='KSTAR_NUM_COMPENDIA_CLASS')
    parser.add_argument('-c', '--processes', action='store', dest='processes', help='Number of processes to run in parallel', type = int, default = 1)
    results = parser.parse_args()
    return results

def process_args(results):
    # get logger
    if results.odir is None or not (path.exists(results.odir) and path.isdir(results.odir)):
        log = logger.get_logger(results.name, f"{results.name}_mapping.log")
    else:
        log = logger.get_logger(results.name, f"{results.odir}/{results.name}_mapping.log")
    
    #check if resource directory exists
    if not (path.exists(results.rdir) and path.isdir(results.rdir)):
        log.error("Please provide a valid resource directory")
        exit()
    #check if output directory exists
    if not (path.exists(results.odir) and path.isdir(results.odir)):
        log.error("Please provide a valid output directory")
        exit()
    # check if experiment file exists and is either csv or tsv file. 
    # Load experiment if valid
    if path.exists(results.exp_file) and path.isfile(results.exp_file):
        filetype = results.exp_file.split('.')[-1]
        if filetype == 'csv':
            experiment = pd.read_csv(results.exp_file)
        elif filetype == 'tsv':
            experiment = pd.read_csv(results.exp_file, sep = '\t')
        else:
            log.error("Unrecognized experiment filetype. Please use a csv or tsv file")
            exit()
    else:
        log.error("Please provide a valid experiment file")
        exit()

    # Load ProteomeScout if found
    resource_files = os.listdir(f"{results.rdir}/HUMAN_PROTEOME")
    proteomescout = None
    for f in resource_files:
        if f.startswith("ProteomeScout_KSTAR"):
            proteomescout = pd.read_csv(f"{results.rdir}/COMPENDIA_MAP/{f}")
    if proteomescout is None:
        log.error("ProteomeScout file not found. Please provide a compendia map file to use")
        exit()
    # Check all data columns provided to make sure they exist. 
    # If a column does not exist in the experiment it is removed
    data_columns = None
    if results.data_columns is not None:
        data_columns = []
        for col in results.data_columns:
            if col in columns:
                data_columns.append(col)
            else:
                log.warning(f"{col} not found in experiment columns")
        if len(data_columns) == 0:
            log.warning("No valid columns were found. Reverting to checking if 'data:' is in column name")
    if data_columns is not None and len(data_columns) == 0:
        data_columns = None
    
    return experiment, proteomescout, log, data_columns

def main():
    pool = multiprocessing.Pool(processes = config.PROCESSES)
    results = parse_args()
    experiment, proteomescout, log, data_columns  = process_args(results)
    random_experiments = build_random_experiments(experiment, proteomescout, results.agg, results.threshold, results.num, results.pevent, data_columns = None, pool=pool, selection_type='KSTAR_NUM_COMPENDIA_CLASS', PROCESSES = results.processes)
    random_experiments.to_csv(f"{results.odir}/{results.name}_random_experiments_{results.pevent}.tsv", sep = '\t', index=False)


if __name__ == "__main__":
    main()