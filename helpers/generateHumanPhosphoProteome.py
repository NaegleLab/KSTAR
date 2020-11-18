from proteomeScoutAPI import ProteomeScoutAPI
import pandas as pd
from kinase_activity.src import experiment, logger #upddate for kstar



"""
Download the KinPred Uniprot reference from KinPred project on Figshare
Supporting Reference Data Files for KinPred (Raw) to RESOURCE_FILES directory (or update reference directory for this script)
"""
resourceDir = '../RESOURCE_FILES/'
humanProteomeFile = resourceDir+'uniprot_sequence.csv'
proteomeScoutFile = resourceDir+'Raw/ProteomeScout2020-02-07/data.tsv'
compendiaCitationsFile = resourceDir+'citations_compendia.txt'

outputFileName = resourceDir+'Human_PhosphoProteome_mapped_annotated.csv'



sequences = pd.read_csv(humanProteomeFile)


PTM_API = ProteomeScoutAPI(proteomeScoutFile)

compendia_df = pd.read_csv(compendiaCitationsFile, sep='\t')
compendiaList = list(compendia_df['Experiment ID'].values)
compendiaSet = set(compendiaList)

#to do -- copy the list files here of expiriments/compendia to count (RESOURCE_FILES) and use these to append evidence types to original dataframe, before mapping to create
# an output file (put that outputFile in RESOURCE_FILES)

#get all the header lines and parse the uniprot IDS in the reference human proteome
uniprotList = list(sequences['Entry'].values)
print("Number of human proteins: %d"%(len(uniprotList)))


#Next, get all phosphorylation sites known in those records and create a dataframe output of this information, include the number of compendia, number of experiments
# and the compendiaClass columns. 
numEvidences = 1 
array = []
count_tyr = 0
count_ser = 0
count_thr = 0


for acc in uniprotList: 
	PTMs = PTM_API.get_PTMs_withEvidence(acc)
	#mods comes back as list of lists, each with information about a site, e.g. ('37', 'Y', 'Phosphotyrosine')
	seq = PTM_API.get_sequence(acc)
	if isinstance(PTMs, list):
		for d in PTMs:
			PTM = d['mod']
			evidences = d['evidence']
			if PTM[2] in ['Phosphoserine', 'Phosphothreonine', 'Phosphotyrosine']:
				site = int(PTM[0])
				aaSite = "%s%d"%(PTM[1], site)
				pep = experiment.get_aligned_peptide(aaSite, seq, 7) #e.g. set to 7 for a 15-mer
				numCompendia = len(set(evidences).intersection(compendiaSet))

				#class is defined as low, medium or high 
				if numCompendia < 1: 
					classType = 0
				elif numCompendia <=3:
					classType = 1
				else:
					classType = 2

				array.append([acc, site, pep, PTM[1], aaSite, numCompendia, classType])

				if PTM[1]=='S':
					count_ser += 1
				elif PTM[1] =='T':
					count_thr += 1
				elif PTM[1] == 'Y':
					count_tyr += 1

print("pSer: %d\t pThr: %d\t pTyr: %d"%(count_ser, count_thr, count_tyr))
df = pd.DataFrame(array, columns=['acc', 'site', 'pep', 'type', 'typeSite', 'KSTAR_NUM_COMPENDIA', 'KSTAR_NUM_COMPENDIA_CLASS'])  
print(df.describe())
#now map that dataframe using the mapper code. 
df 

save_directory = './'
log = logger.get_logger('PROTEOMESCOUT_MAPPED', 'ProteomeScout_mapped.log')


#finally, annotate using the annotation code
colMapper = {}
colMapper['accession_id'] = 'acc'
colMapper['peptide'] = 'pep'
colMapper['site'] = 'typeSite'

e = experiment.ExperimentMapper(df, sequences, colMapper, logger=log)

   
e.experiment.to_csv(outputFileName)
