import pandas as pd
from kinase_activity.src import logger
import argparse
from os import path
import os
import kstar
import pickle


def parse_args():
    parser = argparse.ArgumentParser(description='Parse KSTAR Activity Arguments')
    parser.add_argument('--exp_file', '--experiment_file', action='store', dest= 'exp_file', help='Experiment file location. csv or tsv file', required=True)
    parser.add_argument('--rdir','--resource_directory', action='store',dest='rdir', help = 'resource file directory', required=True)
    parser.add_argument('--odir','--output_directory', action = 'store', dest='odir', help = 'output directory for results',)
    parser.add_argument('--pevent','--phosphorylation_events', action = 'store', dest='pevent', help ='phosphorylation event type', choices=['Y','ST'], default=['ST','Y'], nargs='*')
    parser.add_argument('--name', action = 'store', dest='name', help = 'experiment name', default='Experiment')
    parser.add_argument('--cols', '--data_columns', action='store', dest='data_columns', help = 'data_columns to use', nargs='*')
    parser.add_argument('--agg', '--activity_agg', action='store', dest='activity_agg', help = 'activity agg to use', default='count', choices =['count','mean'])
    parser.add_argument('--thresh', '--threshold',  action='store', dest='threshold', help = 'threshold to use for analysis', type = float, default=0.0)

    results = parser.parse_args()
    return results

def process_args(results):
    # get logger
    if results.odir is None or not (path.exists(results.odir) and path.isdir(results.odir)):
        log = logger.get_logger(results.name, f"{results.name}_activity.log")
    else:
        log = logger.get_logger(results.name, f"{results.odir}/{results.name}_activity.log")
    
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
    
    # Get sequence dict from resource directory fasta file
    resource_files = os.listdir(results.rdir)

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

    
    return experiment, log, data_columns, results.pevent, results.activity_agg, results.threshold

def run_kstar_analysis(experiment, log, networks, data_columns, activity_agg, threshold):
    kinact = kstar.KinaseActivity(experiment, log)
    kinact.add_networks_batch(networks)
    kinact.calculate_kinase_activities(data_columns, agg=activity_agg, threshold=threshold)
    kinact.aggregate_activities()
    kinact.summarize_activites()
    return kinact


def save_kstar(kinact, phospho_event, name,odir):
    name = f"{name}_{phospho_event}"
    kinact.activities.to_csv(f"{odir}/{name}_activities.tsv", sep = '\t', index = False)
    kinact.agg_activities.to_csv(f"{odir}/{name}_aggregated_activities.tsv", sep = '\t', index = False)
    kinact.activity_summary.to_csv(f"{odir}/{name}_summarized_activities.tsv", sep = '\t', index = False)
    pickle.dump( kinact, open( f"{odir}/{name}_kinact.p", "wb" ) )

def main():
    results = parse_args()
    experiment, log, data_columns, phospho_events, activity_arg, threshold = process_args(results)
    

    for phospho_event in phospho_events:
        if phospho_event == 'Y':
            log.info("Running Tyrosine Kinase Activity Analysis")
            print("Running Tyrosine Kinase Activity Analysis")
            networks  = pickle.load( open( f"{results.rdir}/pTyr_networks.p", "rb" ) )
            kinact = run_kstar_analysis(experiment, log, networks, data_columns, activity_arg, threshold)
            save_kstar(kinact, phospho_event, results.name, results.odir)
            

        if phospho_event == 'ST':
            log.log("Running Serine/Threonine Kinase Activity Analysis")

    # exp_mapper.experiment.to_csv(f"{results.odir}/{results.name}_mapped.tsv", sep = '\t', index = False)
    
if __name__ == "__main__":
    main()

