#from Bio import SeqIO
import os, re, json, sys, io
import logging
import argparse
import pandas as pd
import numpy as np
import urllib.parse
import urllib.request
import inspect

def suppress_print(func, *args, **kwargs):
    saved_stdout = sys.stdout
    sys.stdout = io.StringIO()
    try:
        result = func(*args, **kwargs)
    finally:
        sys.stdout = saved_stdout
    return result

def process_fasta_file(fasta_file):
	"""
	For configuration, to convert the global fasta sequence file into a sequence dictionary that can be used in mapping

	Parameters
	----------
	fasta_file : str
		file location of fasta file
	
	Returns
	---------
	sequences : dict
		{acc : sequence} dictionary generated from fasta file
	"""
	sequences = {}
	with open(fasta_file, "r") as f:
		lines = f.readlines()
		#find header lines by looking for lines that start with '>'

		for i in range(len(lines)):
			if lines[i].startswith('>'):
				accession = lines[i].split('|')[1] #get accession from header
				
				#grab all sequences
				seq = ''
				j = i + 1
				#collect sequence lines until the next header or end of file
				while j < len(lines) and not lines[j].startswith('>'):
					seq += lines[j].strip()
					j += 1
				#extract accession from header
				sequences[accession] = seq
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
	#check to make sure the logger does not already have handlers
	if logger.hasHandlers():
		logger.handlers.clear()

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

	# get the unique identifiers
	accVals = list(set(df[acc_col_name].values))

	# create a query string for uniprot query
	queryString = ''
	for acc in accVals:
		# remove the isoform, which maps strangely
		queryString = "%s %s" % (queryString, acc)

	url = 'https://www.uniprot.org/uploadlists/'

	params = {
		'from': acc_col_type,
		'to': 'ACC',
		'format': 'tab',
		'query': queryString  # list(set(refseqList))
	}

	data = urllib.parse.urlencode(params)
	data = data.encode('utf-8')
	req = urllib.request.Request(url, data)
	with urllib.request.urlopen(req) as f:
		response = f.read()
	# print(response.decode('utf-8'))

	ref_to_uni = {}

	for line in response.decode('utf-8').split('\n'):
		if line:
			l_arr = line.split('\t')
			fromVal = l_arr[0]
			if acc_col_type == 'P_GI':
				fromVal = 'gi|' + str(l_arr[0])
			ref_to_uni[fromVal] = l_arr[1]

	# now walk through each row, create a unique and add accession
	uniprot_arr = []
	for index, row in df.iterrows():
		acc = row[acc_col_name]
		if acc in ref_to_uni:
			uniprot_arr.append(ref_to_uni[acc])
		else:
			uniprot_arr.append('NotFound')
	df[acc_uni_name] = uniprot_arr
	return df

#get network hash associated with Y and ST networks
def parse_network_information(network_directory, file_type = 'txt'):
	"""
	Parse the RUN_INFORMATION.txt file from network pruning run and extract its data.

	Args:
		file_path (str): Path to the RUN_INFORMATION.txt file.

	Returns:
		dict: A dictionary containing the parsed data.
	"""
	if file_type == 'txt':
		file_path = os.path.join(network_directory, "RUN_INFORMATION.txt")
		try:
			with open(file_path, 'r') as file:
				content = file.read()
			return {
				"unique_network_id": re.search(r"Unique Network ID:\s+([a-fA-F0-9]+)", content).group(1),
				"unique_reference_id": re.search(r"Unique Reference ID:\s+([a-fA-F0-9]+)", content).group(1),
				"date_run": re.search(r"Date Run\s+([\d-]+\s[\d:.]+)", content).group(1),
				"network_used": re.search(r"Network Used\s+([^\n]+)", content).group(1).strip(),
				"phospho_type": re.search(r"Phospho Type\s+(\w+)", content).group(1),
				"kinase_size": int(re.search(r"Kinase Size\s+(\d+)", content).group(1)),
				"site_limit": int(re.search(r"Site Limit\s+(\d+)", content).group(1)),
				"num_networks": int(re.search(r"# of Networks\s+(\d+)", content).group(1)),
				"use_compendia": re.search(r"Use Compendia\s+(\w+)", content).group(1).lower() == "yes",
				"compendia_counts": list(map(int, re.findall(r"Compendia \d+\s+(\d+)", content))),
			}
		except FileNotFoundError:
			raise FileNotFoundError(f"RUN_INFORMATION.txt file not found at: {file_path}")
		except AttributeError as e:
			raise ValueError(f"Error parsing the RUN_INFORMATION.txt file: {e}")
	elif file_type == 'json':
		file_path = os.path.join(network_directory, "RUN_INFORMATION.json")

		try:
			with open(file_path, 'r') as file:
				data = json.load(file)

			return data
		except FileNotFoundError:
			raise FileNotFoundError(f"RUN_INFORMATION.json file not found at: {file_path}")
		except json.JSONDecodeError as e:
			raise ValueError(f"Error parsing the RUN_INFORMATION.json file: {e}")
	else:
		raise ValueError("file_type must be either 'txt' or 'json'")
	

def extract_relevant_kwargs(func, **kwargs):
    func_args = inspect.getfullargspec(func).args
    return {k: kwargs[k] for k in func_args if k in kwargs}

def extract_kwonlyargs(func, **kwargs):
    func_kwonlyargs = inspect.getfullargspec(func).kwonlyargs
    return {k: kwargs[k] for k in func_kwonlyargs if k in kwargs}

def calculate_jaccard_by_sets(set1, set2):
	"""
	Compares two sets and calculates the Jaccard index between them
	"""
	intersection = len(set1.intersection(set2))
	union = len(set1.union(set2))
	if union == 0:
		return 0.0
	return intersection / union

def calculate_jaccard_by_binary(set1, set2):
	"""
	Compares two binary arrays and calculates the Jaccard index between them (based on number of matches)
	"""
	#makes sure vectors are the same length
	if len(set1) != len(set2):
		raise ValueError("Input binary arrays must be of the same length.")
	
	intersection = np.logical_and(set1, set2).sum()
	union = np.logical_or(set1, set2).sum()
	if union == 0:
		return 0.0
	return intersection / union

def jaci_matrix_between_samples(evidence, samples=None):
	"""
	This function creates a looks at the similarity of evidence between samples based on Jaccard index of phosphopeptide identities

	Parameters
	----------
	evidence: pd.DataFrame
		evidence dataframe, preferably one that has been binarized
	samples: a list of sample columns 
	
	Returns
	-------
	jaccard_matrix: pd.DataFrame
		a dataframe showing the similarity of phosphopeptide identities between samples
	"""
	#from sklearn.metrics.pairwise import pairwise_distances
	if samples is None:
		samples = [col for col in evidence.columns if col.startswith("data:")]
		if len(samples) == 0:
			raise ValueError("No sample columns found in evidence dataframe. Please provide a list of sample columns.")
	else:
		samples = [s for s in samples if s in evidence.columns]

	# check if binary matrix, convert to boolean. Raise warning if not binary, but convert by removing zeros and np.nans
	if ((evidence[samples].values == 0) | (evidence[samples].values == 1)).all():
		evidence[samples] = evidence[samples].astype(bool)
	else:
		print("Warning: evidence matrix is not binary. Converting to binary by treating all non-zero, non-nan values as True.")
		evidence[samples] = evidence[samples].fillna(0)
		evidence[samples] = evidence[samples].astype(bool)


      
	#Calculate similarity between each column and construct a table to house the values
	jaccard_matrix = pd.DataFrame(np.nan, columns = samples, index = samples)
	for i in range(len(samples)-1):
		for j in range(i+1, len(samples)):
			set1 = evidence[samples[i]].values
			set2 = evidence[samples[j]].values
			similarity = calculate_jaccard_by_binary(set1,set2)
			
			#Add similarity metric to both sides of the table, so that it is symmetric
			jaccard_matrix.loc[samples[i],samples[j]] = similarity
			jaccard_matrix.loc[samples[j],samples[i]] = similarity
	
	#fill diagonal with 1s
	for s in samples:
		jaccard_matrix.loc[s,s] = 1.0

	#remove "data:" from sample names for easier readability
	jaccard_matrix.columns = [s.replace("data:", "") for s in samples]
	jaccard_matrix.index = [s.replace("data:", "") for s in samples]
	return jaccard_matrix


