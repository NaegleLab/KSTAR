#!/usr/bin/env python3
import pandas as pd 
import argparse
import config
import numpy as np

def summarize_activities(activities, method = 'median_activity'):
    """
    Builds a single combined dataframe from the provided activities such that 
    each piece of evidence is given a single column. Values are based on the method selected.
    The method must be a column in the activities

    Parameters
    ----------
    activities : dict
        hypergeometric activities that have previously been summarized by network.
        key : experiment name
        value : hypergeometric activity
    method : str
        The column in the hypergeometric activity to use for summarizing data
    
    Returns
    ---------
    activity_summary : pandas DataFrame

    """

    available_methods = list(activities.columns)
    available_methods.remove('data')
    if method not in available_methods:
        raise ValueError(f"the method '{method}' is not in the availble methods. \nAvailable methods include : {', '.join(available_methods)}")


    activity_summary = activities.pivot(index = config.KSTAR_KINASE, columns ='data', values = method).reset_index().rename_axis(None, axis=1).set_index(config.KSTAR_KINASE)
    # activity_summary = activity_summary[self.data_columns]
    return activity_summary



def parse_args():
    parser = argparse.ArgumentParser(description='Parse Normalization Arguments')
    parser.add_argument( '--activities', action='store', dest= 'activities', help='Experiment hypergeometric activity', required=True)
    parser.add_argument( '--method', action='store', dest= 'method', help='method for summarizing', default='median_activity')
    parser.add_argument( '--name', action='store', dest= 'name', help='name of activity', default='experiment')

    results = parser.parse_args()
    return results


def main():
    results = parse_args()

    activities = pd.read_table(results.activities)
    
    activity_summary = summarize_activities(activities, results.method)
    
    activity_summary.to_csv(f"{results.name}.tsv", sep = '\t', index = True)

#%%
if __name__ == "__main__":
    main()