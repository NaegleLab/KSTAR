import pandas as pd
import multiprocessing
import itertools
import argparse


def build_filtered_experiment(experiment, compendia, filtered_compendia, num_random_experiments, name):
    rand_experiments = compendia[['KSTAR_ACCESSION', 'KSTAR_SITE']].copy()
    if len(experiment) == 0:
        empty_columns = [f"{name}:{i}" for i in range(num_random_experiments)]
        rand_experiments = pd.concat([rand_experiments,pd.DataFrame(columns=empty_columns)])
        return rand_experiments


    sizes = experiment.groupby('KSTAR_NUM_COMPENDIA').size()
    for i in range(num_random_experiments):
        rand_experiment_list = []
        for num, size in sizes.iteritems():
            filtered = filtered_compendia[num]
            filtered_random = filtered.sample(size)
            filtered_random[f"{name}:{i}"] = 1
            rand_experiment_list.append(filtered_random)
        rand_experiment = pd.concat(rand_experiment_list)
        rand_experiments = pd.merge(rand_experiments, rand_experiment, how = 'left', on = ['KSTAR_ACCESSION', 'KSTAR_SITE'])
    return rand_experiments


def build_random_experiments(experiment, compendia, agg, threshold, greater, num_random_experiments, phosphorylation_event, data_columns = None):
    phosphorylation_event = tuple(phosphorylation_event)
    compendia = compendia[compendia['KSTAR_SITE'].str.startswith(phosphorylation_event)]
    experiment = experiment[experiment['KSTAR_SITE'].str.startswith(phosphorylation_event)]
    experiment = experiment.groupby(['KSTAR_ACCESSION', 'KSTAR_SITE']).agg(agg).reset_index()
    if data_columns is None:
        data_columns =[c for c in experiment.columns if c.startswith('data:')]
    sizes = compendia['KSTAR_NUM_COMPENDIA'].unique()
    filtered_compendia = {}
    for s in sizes:
        filtered_compendia[s] = compendia[compendia['KSTAR_NUM_COMPENDIA'] == s][['KSTAR_ACCESSION', 'KSTAR_SITE']]

    
    
    pool = multiprocessing.Pool()
    if greater:
        filtered_experiments = [experiment[experiment[col] >= threshold] for col in data_columns]
    else:
        filtered_experiments = [experiment[experiment[col] >= threshold] for col in data_columns]

    iterable = zip(
            filtered_experiments, 
            itertools.repeat(compendia), 
            itertools.repeat(filtered_compendia), 
            itertools.repeat(num_random_experiments), 
            [col for col in data_columns])
        
    rand_experiments_list =  pool.starmap(build_filtered_experiment, iterable)

    rand_experiments = compendia[['KSTAR_ACCESSION', 'KSTAR_SITE']].copy()
    for r in rand_experiments_list:
        rand_experiments = pd.merge(rand_experiments, r, how = 'left', on = ['KSTAR_ACCESSION', 'KSTAR_SITE'])
    return rand_experiments
        

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
        log.eror("ProteomeScout file not found. Please provide a compendia map file to use")
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
    results = parse_args()
    experiment, proteomescout, log, data_columns = process_args(results)
    random_experiments = build_random_experiments(experiment, proteomescout, results.agg, results.threshold, results.num, results.pevent, data_columns = None)
    random_experiments.to_csv(f"{results.odir}/{results.name}_random_experiments_{results.pevent}.tsv", sep = '\t', index=False)


if __name__ == "__main__":
    main()