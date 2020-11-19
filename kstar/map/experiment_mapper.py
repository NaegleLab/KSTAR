from Bio import SeqIO
import pandas as pd
import logging
from collections import defaultdict
import argparse
from os import path

ACCESSION_ID = 'KSTAR_ACCESSION'
SITE_ID = 'KSTAR_SITE'
PEPTIDE_ID = 'KSTAR_PEPTIDE'

class ExperimentMapper:

    def __init__(self, experiment, sequences, columns, logger, window = 7, data_columns = None):      
        """
        Given an experiment object and reference sequences, map the phosphorylation sites to the common reference.
        Inputs

        Parameters
        ----------
        name : str
            Name of experiment. Used for logging
        experiment: pandas dataframe
            Pandas dataframe of an experiment that has a reference accession, a peptide column and/or a site column. 
            The peptide column should be upper case, with lower case indicating the site of phosphorylation - this is preferred
            The site column should be in the format S/T/Y<pos>, e.g. Y15 or S345
        sequences: dict
            Dictionary of sequences. Key : accession. Value : protein sequence
        columns: dict
            Dictionary with mappings of the experiment dataframe column names for the required names 'accession_id', 'peptide', or 'site'. 
            One of 'peptide' or 'site' is required. 
        logger: Looger object
            used for logging when peptides cannot be matched and when a site location changes
        window : int
            The length of amino acids to the N- and C-terminal sides of the central phosphoprotein to map a site to.


        Returns
        -------
        experiment_mapper: ExperimentMapper object that has been mapped

        """ 
        self.experiment = experiment
        self.sequences = sequences

        if 'accession_id' not in columns.keys():
            raise ValueError('ExperimentMapper requires accession_id as a dictionary key')
        else:
            self.experiment[ACCESSION_ID] = self.experiment[columns['accession_id']] 

        if 'peptide' not in columns.keys() and 'site' not in columns.keys():
            raise ValueError('ExperimentMapper requires either site or peptide as keys in dictionary')

        self.experiment[PEPTIDE_ID] = self.experiment[columns['peptide']] if 'peptide' in columns.keys() else None
        self.experiment[SITE_ID] = self.experiment[columns['site']] if 'site' in columns.keys() else None

        self.logger = logger
        self.set_data_columns(data_columns)
        self.align_sites(window)

        
        
    def set_data_columns(self, data_columns):
        self.data_columns = []
        if data_columns is None:
            for col in self.experiment.columns:
                if col.startswith('data:'):
                    self.data_columns.append(col)
        else:
            data_rename = defaultdict()
            for col in data_columns:
                if not col.startswith('data:'):
                    data_rename[col] = f"data:{col}"
                else: 
                    data_rename[col] = col
            self.experiment.rename(columns = data_rename, inplace = True)
            self.data_columns = list(data_rename.values())

    def get_experiment(self):
        return self.experiment


    def get_sequence(self, accession):
        """
        Gets the sequence that matches the given accession
        """
        if accession in self.sequences.keys():
            return None
        sequence = self.sequences[accession]
        return sequence

    
    def align_sites(self, window = 7):
        """
        Map the peptide/sites to the common sequence reference and remove and report errors for sites that do not align as expected.
        expMapper.align_sites(window=7). Operates on the experiment dataframe of class.

        Attributes
        ----------
        window: int
            The length of amino acids to the N- and C-terminal sides of the central phosphoprotein to map a site to.

        """

        self.experiment = expand_peptide(self.experiment, PEPTIDE_ID)

        for index, row in self.experiment.iterrows():
            sequence = self.get_sequence(row[ACCESSION_ID])
            if sequence is not None:
                # If peptide provided then find site
                if row[PEPTIDE_ID] is not None:
                    
                    site =  peptide_site_number(peptide = row[PEPTIDE_ID], 
                                                site = row[SITE_ID],
                                                sequence = sequence)
                    if site is None:
                        self.logger.warning(f"SITE NOT FOUND : {row[ACCESSION_ID]}\t{row[PEPTIDE_ID]}")
                        self.experiment.loc[index, SITE_ID] = None
                        self.experiment.loc[index, PEPTIDE_ID] = None
                        continue
                    elif site != row[SITE_ID]:
                        self.logger.info(f"SITE CHANGED : {row[ACCESSION_ID]} {row[SITE_ID]} -> {site}")
                        self.experiment.loc[index, SITE_ID] = site
                    
                    peptide = get_aligned_peptide(site = site,
                                                  sequence = sequence, 
                                                  window = window)
                    if peptide is not None:
                        self.experiment.loc[index, PEPTIDE_ID] = peptide
                        
                # Peptide not provided but site is provided - build peptide
                elif row[SITE_ID] is not None:
                    peptide = get_aligned_peptide(site = row[SITE_ID],
                                                  sequence = sequence, 
                                                  window = window)
                    if peptide is not None:
                        self.experiment.loc[index, PEPTIDE_ID] = peptide
                    else:
                        self.logger.warning(f"PEPTIDE MISMATCH : {row[ACCESSION_ID]} {row[SITE_ID]}")
                        self.experiment.loc[index, SITE_ID] = None
                        self.experiment.loc[index, PEPTIDE_ID] = None
                else:
                    self.logger.warning(f"SITE NOT FOUND : {self.experiment[ACCESSION_ID]}\t{self.experiment[SITE_ID]}")
            
            # Sequence not found in compendia
            else:
                self.logger.warning(f"SEQUENCE NOT FOUND : {row[ACCESSION_ID]}")
                self.experiment.loc[index, SITE_ID] = None
                self.experiment.loc[index, PEPTIDE_ID] = None
        self.experiment.dropna(axis = 'rows', subset = [ACCESSION_ID, SITE_ID, PEPTIDE_ID], inplace = True)
        self.experiment.drop_duplicates(inplace=True)


def process_fasta_file(fasta_file):
    seqs = SeqIO.parse(open(fasta_file), 'fasta')

    sequences = defaultdict()
    for entry in seqs:
            seq = str(entry.seq)
            acc = entry.id.split('|')[1].strip()
            sequences[acc] = seq

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
    parse.add_argument('-d', '--data', '--data_columns', action='store', dest='data_columns', help = 'data_columns to use', nargs='*')
    results = parser.parse_args()
    return results

def process_args(results):
    # get logger
    logger = get_logger(results.name, filename)
    if results.odir is None or not (path.exists(results.odir) and path.isdir(results.odir)):
        logger = logging.get_logger(results.name, f"{results.name}_mapping.log")
    else:
        logger = get_logger(results.name, f"{results.odir}/{results.name}_mapping.log")
    
    #check if resource directory exists
    if not (path.exists(results.rdir) and path.isdir(results.rdir)):
        logger.error("Please provide a valid resource directory")
        exit()
    #check if output directory exists
    if not (path.exists(results.odir) and path.isdir(results.odir)):
        logger.error("Please provide a valid output directory")
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
            logger.error("Unrecognized experiment filetype. Please use a csv or tsv file")
            exit()
    else:
        logger.error("Please provide a valid experiment file")
        exit()
    
    # Map accession, peptide, site column if valid
    columns = list(experiment.columns)
    map_columns = defaultdict()
    if results.accession in columns:
        map_columns['accession_id'] = results.accession
    else:
        logger.error(f"{results.accession} not found in experiment columns. Please provide a valid accession column")
        exit()
    if results.peptide in columns:
        map_columns['peptide'] = results.peptide
    else:
        logger.error(f"{results.peptide} not found in experiment columns. Please provide a valid peptide column")
        exit()  
    if results.site is not None and results.site in columns:
        map_columns['site'] = results.site
    elif results.site is not None and results.site not in columns:
        logger.error(f"{results.site} not found in experiment columns. Please provide a valid site column")
        exit()
    
    # Get sequence dict from resource directory fasta file
    resource_files = os.listdir(results.rdir)
    sequences = None
    for f in resource_files:
        if f.endswith(".fasta"):
            sequences = process_fasta_file(f)
    if sequences is None:
        logger.eror("Fasta file not found. Please provide a UniProt fasta file to use")
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
                logger.warning(f"{col} not found in experiment columns")
        if len(data_columns) == 0:
            logger.warning("No valid columns were found. Reverting to checking if 'data:' is in column name")
    if data_columns is not None and len(data_columns) == 0:
        data_columns = None
    
    return experiment, sequences, logger, map_columns, data_columns

def main():
    results = parse_args()
    experiment, sequences, logger, map_columns, data_columns = process_args(results)
    exp_mapper = ExperimentMapper(experiment, sequences, map_columns, logger, results.window, data_columns)

    exp_mapper.experiment.to_csv(f"{results.odir}/{results.name}_mapped.tsv", sep = '\t', index = False)
    
if __name__ == "__main__":
    main()

