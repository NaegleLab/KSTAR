import pandas as  pd

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

def calculate_experiment_fpr(group, target_alpha):
    """
    Calcualtes the FPR for an expeiment for all kinases
    """
    experiment_fpr = {}
    for kinase in group.columns:
        pvalue_array = list(group[kinase])
        experiment_fpr[kinase] = single_experiment_kinase_fpr(pvalue_array, target_alpha)
    experiment_fpr = pd.Series(experiment_fpr)
    return experiment_fpr

def generate_fpr_values(random_activities, target_alpha):
    """
    Given random_activites where index is Kinase and columns are random experiments kianse activity
    calculate the fpr values

    Paramters
    ---------
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
    """
    if target_alpha > 1 or target_alpha < 0:
        print("ERROR: Using default target_alpha of 0.05, FPR must be between 0 and 1")
        target_alpha = 0.05
    df_rename = {col:':'.join(col.split(':')[:-1]) for col in random_activities.columns}
    random_activities = random_activities.T
    random_activities.rename(index  = df_rename, inplace = True)

    fpr = random_activities.groupby(random_activities.index).apply(lambda group: calculate_experiment_fpr(group, target_alpha))
    fpr = fpr.T
    return fpr

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



