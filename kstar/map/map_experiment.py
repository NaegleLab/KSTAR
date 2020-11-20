
from Bio import SeqIO
import pandas as pd
from collections import defaultdict
from kinase_activity.src import logger
import argparse
from os import path
import os
import experiment_mapper


def process_fasta_file(fasta_file):
    seqs = SeqIO.parse(open(fasta_file), 'fasta')

    sequences = defaultdict()
    for entry in seqs:
            seq = str(entry.seq)
            acc = entry.id.split('|')[1].strip()
            sequences[acc] = seq
    return sequences

def parse_args():
    parser = argparse.ArgumentParser(description='Parse Mapping Inference Arguments')
    parser.add_argument('-e', '--exp_file', '--experiment_file', action='store', dest= 'exp_file', help='Experiment file location. csv or tsv file', required=True)
    parser.add_argument('-r','--rdir','--resource_directory', action='store',dest='rdir', help = 'resource file directory', required=True)
    parser.add_argument('-o','--odir','--output_directory', action = 'store', dest='odir', help = 'output directory for results',)
    # parser.add_argument('-p','--pevent','--phosphorylation_events', action = 'store', dest='pevent', help ='phosphorylation event type', choices=['Y','S','T','ST', 'STY'], default='STY')
    parser.add_argument('-n', '--name', action = 'store', dest='name', help = 'experiment name', default='Experiment')
    parser.add_argument('-a', '--accession', action = 'store', dest='accession', help = 'Protein Accession column in experiment file', required=True,)
    parser.add_argument('-s', '--site',action = 'store', dest='site',  help='Site column in experiment file')
    parser.add_argument('-p', '--peptide', action = 'store', dest='peptide', help = 'Peptide column in experiment file', required=True,)
    parser.add_argument('-w', '--window', action='store', dest='window', help = 'peptide window', type = int, default=7)
    parser.add_argument('-d', '--data', '--data_columns', action='store', dest='data_columns', help = 'data_columns to use', nargs='*')
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
    
    # Map accession, peptide, site column if valid
    columns = list(experiment.columns)
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
        exit()
    
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
    resource_files = os.listdir(f"{results.rdir}/COMPENDIA_MAP")
    compendia = None
    for f in resource_files:
        if f.startswith("Human_PhosphoProteome_mapped_annotated"):
            compendia = pd.read_csv(f"{results.rdir}/COMPENDIA_MAP/{f}")
    if compendia is None:
        log.eror("Compendia mp file not found. Please provide a compendia map file to use")
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
    
    return experiment, sequences, compendia, log, map_columns, data_columns

def main():
    results = parse_args()
    experiment, sequences, compendia, log, map_columns, data_columns = process_args(results)
    exp_mapper = experiment_mapper.ExperimentMapper(
        experiment = experiment,
        sequences = sequences, 
        columns = map_columns, 
        logger = log, 
        compendia = compendia, 
        window = 7, 
        data_columns = data_columns)
    exp_mapper.experiment.to_csv(f"{results.odir}/{results.name}_mapped.tsv", sep = '\t', index = False)
    
if __name__ == "__main__":
    main()

