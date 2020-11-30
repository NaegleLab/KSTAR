from Bio import SeqIO
import logging


def process_fasta_file(fasta_file):
	'''
	For configuration, to convert the global fasta sequence file into a sequence dictionary that can be used in mapping
	'''
	seqs = SeqIO.parse(open(fasta_file), 'fasta')

	sequences = {}
	for entry in seqs:
		seq = str(entry.seq)
		acc = entry.id.split('|')[1].strip()
		sequences[acc] = seq
	return sequences


def get_logger(name, filename):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(filename)
    log_format = logging.Formatter('%(asctime)s\t%(name)s\t%(levelname)s:\t%(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    handler.setFormatter(log_format)
#     if (logger.hasHandlers()):
#         logger.handlers.clear()
    logger.addHandler(handler)
    return logger
