from Bio import SeqIO
import logging
import argparse

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
