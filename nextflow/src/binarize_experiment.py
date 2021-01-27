#!/usr/bin/env python3
#%%
import pandas as pd 
import argparse
import config
import helpers

def create_binary_evidence(evidence, data_columns, agg = 'count', threshold = 1.0,  greater = True):
        """
        Returns a binary evidence data frame according to the parameters passed in for method for aggregating
        duplicates and considering whether a site is included as evidence or not

        Parameters
        ----------
        data_columns : list
            columns that represent experimental result, if None, takes the columns that start with `data:'' in experiment. 
            Pass this value in as a list, if seeking to calculate on fewer than all available data columns
        threshold : float
            threshold value used to filter rows 
        agg : {'count', 'mean'}
            method to use when aggregating duplicate substrate-sites. 
            'count' combines multiple representations and adds if values are non-NaN
            'mean' uses the mean value of numerical data from multiple representations of the same peptide.
                NA values are droped from consideration.
        greater: Boolean
            whether to keep sites that have a numerical value >=threshold (TRUE, default) or <=threshold (FALSE)
        
        Returns
        -------
        evidence_binary : pd.DataFrame
            Matches the evidence dataframe of the kinact object, but with 0 or 1 if a site is included or not.
            This is uniquified and rows that are never used are removed.
        
        """
        evidence = evidence.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE]).agg(agg).reset_index()
        
        #set the binary evidence for whether a site is included
        evidence_binary = evidence.copy()
        for col in data_columns:
            if greater:
                evidence_binary[col].mask(evidence[col] >= threshold, 1, inplace=True)
                evidence_binary[col].mask(evidence[col] < threshold, 0, inplace=True)
            else:
                evidence_binary[col].mask(evidence[col] <= threshold, 1, inplace=True)
                evidence_binary[col].mask(evidence[col] > threshold, 0, inplace=True)

        #remove phosphorylation sites that were not selected in any experiment (useful for very large experiments where removing the need to copy data reduces time)
        evidence_binary.drop(evidence_binary[evidence_binary[data_columns].sum(axis=1) == 0].index, inplace = True) 
        return evidence_binary


def parse_args():
    parser = argparse.ArgumentParser(description='Parse Normalization Arguments')
    parser.add_argument( '--evidence', action='store', dest= 'evidence', help='Experiment hypergeometric activity', required=True)
    parser.add_argument( '--data_columns', action='store', dest= 'data_columns', help='Data columns of experients', required=True, nargs='*')
    parser.add_argument('--name', action = 'store', dest='name', help = 'experiment name', default='Experiment')
    parser.add_argument('--agg', '--activity_agg', action='store', dest='activity_agg', help = 'activity agg to use', default='count', choices =['count','mean'])
    parser.add_argument('--thresh', '--threshold',  action='store', dest='threshold', help = 'threshold to use for analysis', type = float, default=1.0)
    parser.add_argument('--greater', action='store', dest='greater', help = 'whether evidence present if greater than threshold', default='yes')
    
    results = parser.parse_args()
    return results


def main():
    results = parse_args()

    evidence = pd.read_table(results.evidence)
    kstar_columns = [col for col in evidence.columns if col.startswith("KSTAR_")]
    keep_columns = kstar_columns + results.data_columns 
    evidence = evidence[keep_columns]
    greater = helpers.string_to_boolean(results.greater)
    binary_evidence = create_binary_evidence(evidence, results.data_columns, results.activity_agg, results.threshold, greater)
    binary_evidence.to_csv(f"{results.name}_binarized_experiment.tsv", sep = "\t", index=False)

    
#%%
if __name__ == "__main__":
    main()

# #%%
# evidence = pd.read_table("/Users/bj8th/Documents/GitHub/KSTAR/nextflow/test_data/BCR-ABL_mapped.tsv")
# data_columns = ["data:EOE", "data:PRE"]
# keep_columns = [config.KSTAR_ACCESSION, config.KSTAR_SITE] + data_columns 
# threshold = 0.5
# greater = True
# agg = 'mean'
# evidence = evidence[keep_columns]

# binary_evidence = create_binary_evidence(evidence, data_columns, agg, threshold, greater)
# # %%
