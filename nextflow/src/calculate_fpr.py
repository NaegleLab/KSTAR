#!/usr/bin/env python3
import pandas as  pd
import argparse

def single_experiment_kinase_fpr(pvalue_array, target_alpha):
    """
    For a single experiment and kinase calculate the FPR at the target_alpha

    Parameters
    ----------
    pvalue_array : array
        list of pvalues 
    
    Returns
    ---------

    """
    pvalue_array.sort()
    target_number = int(target_alpha * len(pvalue_array))
    pvalue = pvalue_array[target_number]
    return pvalue

def calculate_experiment_fpr(group, target_alpha, real_p_value_series):
    """
    Calcualtes the pvalue that controls FPR at target alpha and the false positive rate for the real_p_value, given the group of median pvalues
    across all experiments. 

    Parameters
    ----------
    """
    experiment_fpr = {}
    experiment_pvalue = {}
    for kinase in group.columns:
        pvalue_array = list(group[kinase])
        experiment_pvalue[kinase] = single_experiment_kinase_fpr(pvalue_array, target_alpha)
        experiment_fpr[kinase] = single_pvalue_fpr(pvalue_array, real_p_value_series[kinase])
    experiment_fpr = pd.Series(experiment_fpr)
    experiment_pvalue = pd.Series(experiment_pvalue)


    return experiment_pvalue, experiment_fpr

def generate_fpr_values(activity_summary, random_activities, target_alpha, data_column):
    """
    Given random_activities where index is Kinase and columns are random experiments kinase activity
    calculate the fpr values

    Paramters
    ---------
    activity_summary: pandas DataFrame
        index: Kinase
        columns: the experiment names
    random_activities : pands DataFrame
        index : Kinase
        columns : random experiments where each iteration is name:#
    target_alpha: float
        target alpha to use
    
    Returns
    --------
    fpr : pandas DataFrame
        index : Kinase
        column : experiment name
        value: the pvalue at which the false positive rate is controlled at target_alpha value
    """
    if target_alpha > 1 or target_alpha < 0:
        print("ERROR: Using default target_alpha of 0.05, FPR must be between 0 and 1")
        target_alpha = 0.05
    random_activities = random_activities.T



    values = calculate_experiment_fpr(random_activities, target_alpha, activity_summary[data_column])
    #assemble the dataframes
    p_value_df = pd.DataFrame(index=activity_summary.index, columns=[data_column])
    fpr_df = pd.DataFrame(index=activity_summary.index, columns=[data_column])

    p_value_df[data_column] = values[0]
    fpr_df[data_column] = values[1]

    return p_value_df, fpr_df

def single_pvalue_fpr(random_stats, p_value):
    """
    Given a list of pvalues generated from random trials and a p_value from a real test, report the positivity rate (as a fraction between 0 and 1)

    Parameters
    ----------
    randomStats: vector of floats
        list of pvalues from random tests
    p_value: float
        p_value of a real test, sorted from lowest to highest

    Returns
    -------
    fpr: float
        false positive rate of p_value as seen in random_stats

    """
    random_stats.sort()
    #find the position in random_stats that is > than that p_value and the fpr is that index/len(random_stats)
    res = len(random_stats)
    for x, val in enumerate(random_stats):
        if val > p_value:
            res = x
            break
    fpr = res/len(random_stats)
    return fpr

def parse_args():
    parser = argparse.ArgumentParser(description='Parse Normalization Arguments')
    parser.add_argument( '--experiment_activity', action='store', dest= 'experiment_activity', help='Experiment hypergeometric activity', required=True)
    parser.add_argument( '--random_activity', action='store', dest= 'random_activity', help='random hypergeometric activity', required=True)
    parser.add_argument( '--data_column', action='store', dest= 'data_column', help='Data column of random experients', required=True)
    parser.add_argument( '--fpr_alpha', action='store', dest= 'fpr_alpha', help='FPR target alpha', default=0.05, type=float)
    results = parser.parse_args()
    return results


def main():
    results = parse_args()

    experiment_activity = pd.read_table(results.experiment_activity, index_col=0)
    random_activity = pd.read_table(results.random_activity, index_col=0)

    p_value_df, fpr_df = generate_fpr_values(experiment_activity, random_activity, results.fpr_alpha, results.data_column)

    p_value_df.to_csv(f"{results.data_column}_normalization.tsv", sep = "\t")
    fpr_df.to_csv(f"{results.data_column}_fpr.tsv", sep = "\t")

if __name__ == "__main__":
    main()