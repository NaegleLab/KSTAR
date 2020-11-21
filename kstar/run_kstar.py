import logging
import pandas as pd
import argparse
from collections import defaultdict
from os import path
import os
from activity import kstar
from log.logger import get_logger
from mapper import experiment_mapper
from normalize import generate_random_experiments
import pickle
from Bio import SeqIO

def process_fasta_file(fasta_file):
    seqs = SeqIO.parse(open(fasta_file), 'fasta')

    sequences = defaultdict()
    for entry in seqs:
            seq = str(entry.seq)
            acc = entry.id.split('|')[1].strip()
            sequences[acc] = seq
    return sequences

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def process_args(results):
    # get logger
    if results.odir is None or not (path.exists(results.odir) and path.isdir(results.odir)):
        log = get_logger(results.name, f"{results.name}_activity.log")
    else:
        log = get_logger(results.name, f"{results.odir}/{results.name}_activity.log")
    
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
    
    # Get sequence dict from resource directory fasta file
    resource_files = os.listdir(f"{results.rdir}/HUMAN_PROTEOME")
    sequences = None
    for f in resource_files:
        if f.endswith(".fasta"):
            sequences = process_fasta_file(f"{results.rdir}/HUMAN_PROTEOME/{f}")
            break

    if sequences is None:
        log.eror("Fasta file not found. Please provide a UniProt fasta file to use")
        exit()

    # Load Compendia if found
    resource_files = os.listdir(f"{results.rdir}/HUMAN_PROTEOME")
    compendia = None
    for f in resource_files:
        if f.startswith("Human_PhosphoProteome_mapped_annotated"):
            compendia = pd.read_csv(f"{results.rdir}/HUMAN_PROTEOME/{f}")
    if compendia is None:
        log.eror("Compendia mp file not found. Please provide a compendia map file to use")
        exit()
    
    normalize = str2bool(results.normalize)
    
    return experiment, log, data_columns, sequences, compendia, map_columns, normalize

def parse_args():
    parser = argparse.ArgumentParser(description='Parse KSTAR Arguments')
    parser.add_argument('--exp_file', '--experiment_file', action='store', dest= 'exp_file', help='Experiment file location. csv or tsv file', required=True)
    parser.add_argument('--rdir','--resource_directory', action='store',dest='rdir', help = 'resource file directory', required=True)
    parser.add_argument('--odir','--output_directory', action = 'store', dest='odir', help = 'output directory for results',)
    parser.add_argument('--pevent','--phosphorylation_events', action = 'store', dest='pevent', help ='phosphorylation event type', choices=['Y','ST'], default=['Y','ST'], nargs='*')
    parser.add_argument('--name', action = 'store', dest='name', help = 'experiment name', default='Experiment')
    parser.add_argument('--acc', '--accession', action = 'store', dest='accession', help = 'Protein Accession column in experiment file', required=True,)
    parser.add_argument('--site',action = 'store', dest='site',  help='Site column in experiment file')
    parser.add_argument('-pep', '--peptide', action = 'store', dest='peptide', help = 'Peptide column in experiment file', required=True,)
    parser.add_argument('--window', action='store', dest='window', help = 'peptide window', type = int, default=7)
    parser.add_argument('--cols', '--data_columns', action='store', dest='data_columns', help = 'data_columns to use', nargs='*')
    parser.add_argument('--agg', '--activity_agg', action='store', dest='activity_agg', help = 'activity agg to use', default='count', choices =['count','mean'])
    parser.add_argument('--thresh', '--threshold',  action='store', dest='threshold', help = 'threshold to use for analysis', type = float, default=0.0)
    parser.add_argument('--norm', '--normalize',  action='store', dest='normalize', help = 'Whether to normalize experiments',)
    parser.add_argument('--num', '--num_random_experiments',  action='store', dest='num', help = 'Number of random experiments to generate', type = int, default=150)
    results = parser.parse_args()
    return results

def run_kstar_analysis(experiment, log, networks, data_columns, activity_agg, threshold):
    kinact = kstar.KinaseActivity(experiment, log)
    kinact.add_networks_batch(networks)
    kinact.calculate_kinase_activities(data_columns, agg=activity_agg, threshold=threshold)
    kinact.aggregate_activities()
    kinact.summarize_activites()
    return kinact


def save_kstar(kinact, phospho_event, name,odir):
    name = f"{name}:{phospho_event}"
    kinact.activities.to_csv(f"{odir}/{name}_activities.tsv", sep = '\t', index = False)
    kinact.agg_activities.to_csv(f"{odir}/{name}_aggregated_activities.tsv", sep = '\t', index = False)
    kinact.activity_summary.to_csv(f"{odir}/{name}_summarized_activities.tsv", sep = '\t', index = False)
    pickle.dump( kinact, open( f"{odir}/{name}_kinact.p", "wb" ) )

def main():
    results = parse_args()
    experiment, log, data_columns, sequences, compendia, map_columns, normalize = process_args(results)
    # ************ MAP EXPERIMENT *******************
    print("MAPPING EXPERIMENT")
    if not os.path.exists(f"{results.odir}/MAPPED_DATA"): 
        os.mkdir(f"{results.odir}/MAPPED_DATA")   
    mapping_log = get_logger(f"mapping_{results.name}", f"{results.odir}/MAPPED_DATA/mapping_{results.name}.log")
    exp_mapper = experiment_mapper.ExperimentMapper(
        experiment = experiment,
        sequences = sequences, 
        columns = map_columns, 
        logger = mapping_log, 
        compendia = compendia, 
        window = results.window, 
        data_columns = data_columns)
     
    exp_mapper.experiment.to_csv(f"{results.odir}/MAPPED_DATA/{results.name}_mapped.tsv", sep = '\t', index = False)
    experiment = exp_mapper.experiment

    if normalize:
        # ************ GENERATE RANDOM DATASETS *****************
        print("GENERATING RANDOM DATASETS")
        if not os.path.exists(f"{results.odir}/RANDOM_EXPERIMENTS"): 
            os.mkdir(f"{results.odir}/RANDOM_EXPERIMENTS")  
        random_log = get_logger("RANDOM", f"{results.odir}/RANDOM_EXPERIMENTS/{results.name}_random.log")
        random_experiments = {}
        for phospho_event in results.pevent:
            random_experiments[phospho_event] = generate_random_experiments.build_random_experiments(experiment, compendia, results.activity_agg, results.threshold, results.num, phospho_event, data_columns = None)
            random_experiments[phospho_event].to_csv(f"{results.odir}/RANDOM_EXPERIMENTS/{results.name}_random_experiments_{phospho_event}.tsv", sep = '\t', index=False)

        # ********** RUN RANDOM KINASE ACTIVITY ANALYSIS ***************
        print("RUNNING RANDOM KINASE ACTIVITY ANALYSIS")
        for phospho_event in results.pevent:
            if phospho_event == 'Y':
                log.info("Running RANDOM Tyrosine Kinase Activity Analysis")
                networks  = pickle.load( open( f"{results.rdir}/NETWORKS/TYROSINE/pTyr_networks.p", "rb" ) )
                kinact = run_kstar_analysis(random_experiments[phospho_event], random_log, networks, None, results.activity_agg, results.threshold)
                save_kstar(kinact, phospho_event, f"RANDOM_{results.name}", f"{results.odir}/RANDOM_EXPERIMENTS")
                
            if phospho_event == 'ST':
                log.log("Running Serine/Threonine Kinase Activity Analysis")

    # ************* RUN KINASE ACTIVITY ANALYSIS *********************
    print("RUNNING KINASE ACTIVITY ANALYSIS")
    if not os.path.exists(f"{results.odir}/RESULTS"): 
        os.mkdir(f"{results.odir}/RESULTS") 
    for phospho_event in results.pevent:
        if phospho_event == 'Y':
            log.info("Running Tyrosine Kinase Activity Analysis")
            print("Running Tyrosine Kinase Activity Analysis")
            networks  = pickle.load( open( f"{results.rdir}/NETWORKS/TYROSINE/pTyr_networks.p", "rb" ) )
            kinact = run_kstar_analysis(experiment, log, networks, None, results.activity_agg, results.threshold)
            save_kstar(kinact, phospho_event, results.name, f"{results.odir}/RESULTS")
            

        if phospho_event == 'ST':
            log.log("Running Serine/Threonine Kinase Activity Analysis")

    # exp_mapper.experiment.to_csv(f"{results.odir}/{results.name}_mapped.tsv", sep = '\t', index = False)
    
if __name__ == "__main__":
    main()