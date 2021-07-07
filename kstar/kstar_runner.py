import logging
import pandas as pd
import argparse
from collections import defaultdict
from os import path
import os


from kstar import helpers, config, mapping, calculate

import pickle


def process_args(results):
    # get logger
    if results.odir is None or not (path.exists(results.odir) and path.isdir(results.odir)):
        log = logging.getLogger(results.name)
    else:
        log = helpers.get_logger(results.name, f"{results.odir}/{results.name}_activity.log")
    
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
    
    # Check all data columns provided to make sure they exist. 
    # If a column does not exist in the experiment it is removed
    columns = list(experiment.columns)
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

    # Map accession, peptide, site column if valid
    map_columns = defaultdict()
    if results.accession in columns:
        map_columns['accession_id'] = results.accession
    else:
        log.error(f"{results.accession} not found in experiment columns. Please provide a valid accession column")
        exit()
    if results.peptide in columns:
        map_columns['peptide'] = results.peptide
    else:
        log.error(f"{results.peptide} not found in experiment columns. Please provide a valid peptide column")
        exit()  
    if results.site is not None and results.site in columns:
        map_columns['site'] = results.site
    elif results.site is not None and results.site not in columns:
        log.error(f"{results.site} not found in experiment columns. Please provide a valid site column")
    
    greater = helpers.string_to_boolean(results.greater)
    normalize = helpers.string_to_boolean(results.normalize)
    
    return experiment, log, data_columns, map_columns, greater, normalize

def parse_args():
    parser = argparse.ArgumentParser(description='Parse KSTAR Arguments')
    parser.add_argument('--exp_file', '--experiment_file', action='store', dest= 'exp_file', help='Experiment file location. csv or tsv file', required=True)
    parser.add_argument('--odir','--output_directory', action = 'store', dest='odir', help = 'output directory for results', required=True)
    parser.add_argument('--phospho_types', action = 'store', dest='phospho_types', help ='phosphorylation event type', choices=['Y','ST'], default=['Y','ST'], nargs='*')
    parser.add_argument('--name', action = 'store', dest='name', help = 'experiment name', default='Experiment')
    parser.add_argument('--accession', action = 'store', dest='accession', help = 'Protein Accession column in experiment file', default='accession')
    parser.add_argument('--site',action = 'store', dest='site',  help='Site column in experiment file', default=None)
    parser.add_argument('--peptide', action = 'store', dest='peptide', help = 'Peptide column in experiment file', default=None)
    parser.add_argument('--window', action='store', dest='window', help = 'peptide window', type = int, default=7)
    parser.add_argument('--data_columns', action='store', dest='data_columns', help = 'data_columns to use', nargs='*', default=None)
    parser.add_argument('--activity_agg', action='store', dest='activity_agg', help = 'activity agg to use', choices =['count','mean'], default='count')
    parser.add_argument('--threshold',  action='store', dest='threshold', help = 'threshold to use for analysis', type = float, default=0.0)
    parser.add_argument('--greater', action='store', dest='greater', default='yes')
    parser.add_argument('--normalize',  action='store', dest='normalize', help = 'Whether to normalize experiments',default='yes')
    parser.add_argument('--num_random_experiments',  action='store', dest='num_random_experiments', help = 'Number of random experiments to generate', type = int, default=150)
    parser.add_argument('--target_alpha', action='store', dest='target_alpha', type=float, default=0.05)
    results = parser.parse_args()
    return results


def run_kstar_analysis(run_log, odir, name, experiment, data_columns, map_columns, window, phospho_types, activity_agg, threshold, greater, normalize,num_random_experiments, target_alpha):
    run_log.info("********** RUNNING KSTAR ANALYSIS PIPELINE **************")
    # ************ MAP EXPERIMENT *******************
    run_log.info("MAPPING EXPERIMENT")
    experiment = mapping.run_mapping(experiment, odir, name, map_columns, window, data_columns)

    data_columns = [f"data:{col}" for col in data_columns]

    # ************* RUN KINASE ACTIVITY ANALYSIS *********************
    run_log.info("RUNNING KINASE ACTIVITY ANALYSIS")
    if not os.path.exists(f"{odir}/RESULTS"): 
        os.mkdir(f"{odir}/RESULTS")  
    activity_log = helpers.get_logger(name, f'{odir}/RESULTS/{name}_kstar_activity.log')
    networks = {}
    if 'Y' in phospho_types:
        networks['Y'] = pickle.load(open(config.NETWORK_Y_PICKLE, "rb" ) )
    if 'ST' in phospho_types:
        networks['ST'] = pickle.load(open(config.NETWORK_ST_PICKLE, "rb" ) )

    kinact_dict = calculate.run_kstar_analysis(experiment, activity_log, networks, phospho_types, data_columns, activity_agg, threshold, greater)

    # ************ RUN NORMALIZATION *****************
    if normalize:  
        run_log.info("RUNNING NORMALIZATION")
        if not os.path.exists(f"{odir}/RANDOM_ANALYSIS"): 
            os.mkdir(f"{odir}/RANDOM_ANALYSIS")  
        random_log = helpers.get_logger("RANDOM", f"{odir}/RANDOM_ANALYSIS/{name}_random.log")
        calculate.normalize_analysis(kinact_dict, random_log, num_random_experiments, target_alpha)
    
    # ******************* SAVE DATA *************************
    calculate.save_kstar_slim(kinact_dict, name, odir)
    

def main():
    results = parse_args()
    experiment, run_log, data_columns, map_columns, greater, normalize,  = process_args(results)
    run_kstar_analysis(
        run_log, results.odir, results.name, 
        experiment, data_columns, map_columns, results.window, 
        results.phospho_types, results.activity_agg, results.threshold, greater, 
        normalize, results.num_random_experiments, results.target_alpha )
    
    

if __name__ == "__main__":
    main()
