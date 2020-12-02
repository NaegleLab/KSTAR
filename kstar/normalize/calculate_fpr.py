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
    df_rename = {col:':'.join(col.split(':')[:-1]) for col in random_activities.columns}
    random_activities = random_activities.T
    random_activities.rename(index  = df_rename, inplace = True)

    fpr = random_activities.groupby(random_activities.index).apply(lambda group: calculate_experiment_fpr(group, 0.05))
    fpr = fpr.T
    return fpr


