#!/usr/bin/env python3
from Bio import SeqIO
import logging
import argparse
import urllib.parse
import urllib.request


def process_fasta_file(fasta_file):
    """
    For configuration, to convert the global fasta sequence file into a sequence dictionary that can be used in mapping

    Parameters
    ---------
    fasta_file : str
        file location of fasta file
    
    Returns
    ---------
    sequences : dict
        {acc : sequence} dictionary generated from fasta file
    """
    seqs = SeqIO.parse(open(fasta_file), 'fasta')

    sequences = {}
    for entry in seqs:
        seq = str(entry.seq)
        acc = entry.id.split('|')[1].strip()
        sequences[acc] = seq
    return sequences


def get_logger(name, filename):
    """
    Finds and returns logger if it exists. Creates new logger if log file does not exist

    Parameters
    ----------
    name : str
    log name
    filename : str
    location to store log file
    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(filename)
    log_format = logging.Formatter('%(asctime)s\t%(name)s\t%(levelname)s:\t%(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    handler.setFormatter(log_format)
    #     if (logger.hasHandlers()):
    #         logger.handlers.clear()
    logger.addHandler(handler)
    return logger

def string_to_boolean(string):
    """
    Converts string to boolean

    Parameters
    ----------
    string :str
        input string

    Returns
    ----------
    result : bool
        output boolean
    """
    if isinstance(string, bool):
        return str
    if string.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif string.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def check_data_columns(evidence, data_columns, logger=None):
        """
        Checks data columns to make sure column is in evidence and that evidence filtered on that data column 
        has at least one point of evidence. Removes all columns that do not meet criteria
        """
        
        if logger is None:
            logger = get_logger("check_data", "check_data.log")

        if logger is None:
            lof
        new_data_columns = []
        for col in data_columns:
            if col in evidence.columns:
                if len(evidence[evidence[col] >= 1]) > 0:
                    new_data_columns.append(col)
                else:
                    logger.warning(f"{col} does not have any evidence")

            else:
                logger.warning(f"{col} not in evidence")
        data_columns = new_data_columns
        return data_columns

def set_data_columns(evidence, data_columns = None):
    """
    Sets the data columns to use in the kinase activity calculation
    If data_columns is None or an empty list then set data_columns to 
    be all columns that start with data:

    Checks all set columns to make sure columns are vaild after filtering evidence
    """
    if data_columns is None or data_columns == []:
        data_columns = []
        for col in evidence.columns:
            if col.startswith('data:'):
                data_columns.append(col)  
    else:
        data_columns = data_columns
    
    data_columns = check_data_columns(evidence, data_columns)
    return data_columns


def convert_acc_to_uniprot(df, acc_col_name, acc_col_type, acc_uni_name):
    """
    Given an experimental dataframe (df) with an accession column (acc_col_name) that is not uniprot, use uniprot 
    to append an accession column of uniprot IDS

    Parameters
    ----------
    df: pandas.DataFrame
        Dataframe with at least a column of accession of interest
    acc_col_name: string
        name of column to convert FROM
    acc_col_type: string
        Uniprot string designation of the accession type to convert FROM, see https://www.uniprot.org/help/api_idmapping
    acc_uni_name: 
        name of new column

    Returns
    -------
    appended_df: pandas.DataFrame
        Input dataframe with an appended acc column of uniprot IDs


    """
    # Convert refseq to Uniprot identifiers from the header of cptac dataset


    #get the unique identifiers 
    accVals = list(set(df[acc_col_name].values))

    #create a query string for uniprot query
    queryString = ''
    for acc in accVals:
        #remove the isoform, which maps strangely

        queryString = "%s %s"%(queryString, acc)

    url = 'https://www.uniprot.org/uploadlists/'

    params = {
    'from': acc_col_type,
    'to': 'ACC',
    'format': 'tab',
    'query': queryString#list(set(refseqList))
    }

    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
       response = f.read()
    #print(response.decode('utf-8'))

    ref_to_uni ={}

    for line in response.decode('utf-8').split('\n'):
        if line:
            l_arr = line.split('\t')
            fromVal = l_arr[0]
            if acc_col_type=='P_GI':
                fromVal = 'gi|'+str(l_arr[0]) 
            ref_to_uni[fromVal] = l_arr[1]

    #now walk through each row, create a unique and add accession
    uniprot_arr = []
    for index, row in df.iterrows():
        acc = row[acc_col_name]
        if acc in ref_to_uni:
            uniprot_arr.append(ref_to_uni[acc])
        else:
            uniprot_arr.append('NotFound')
    df[acc_uni_name]=uniprot_arr
    return df
