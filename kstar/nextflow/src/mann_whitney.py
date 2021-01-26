import logging
import pandas as pd 
import numpy as np 

from itertools import repeat
import concurrent.futures

# def calculate_Mann_Whitney_activities_sig(self, log, number_sig_trials = 100):
#         """
#         For a kinact_dict, where random generation and activity has already been run for the phospho_types of interest, 
#         this will calculate the Mann-Whitney U test for comparing the array of p-values for real data 
#         to those of random data, across the number of networks used.
#         It will also calculate the false positive rate for a pvalue, given observations of a random bootstrapping analysis
        
#         Parameters
#         ----------
#         kinact_dict: dictionary
#             A dictionary of kinact objects, with keys 'Y' and/or 'ST'
#         log: logger 
#             Logger for logging activity messages
#         phospho_types: {['Y', 'ST'], ['Y'], ['ST']}
#             Which substrate/kinaset-type to run activity for: Both ['Y, 'ST'] (default), Tyrosine ['Y'], or Serine/Threonine ['ST']
#         number_sig_trials: int
#             Maximum number of significant trials to run

            
        
#         Returns
#         -------

#         """
#         #First, check that objects are correct and values can be found
#         if not isinstance(self.random_kinact.activities_list, pd.DataFrame):
#             raise ValueError("Random activities do not exist, please run kstar_activity.normalize_analysis")

#         if number_sig_trials > self.num_random_experiments:
#             log.info("Warning: number of trials for Mann Whitney exceeds number available, using %d instead of %d"%(self.num_random_experiments, number_sig_trials))
#             number_sig_trials = self.num_random_experiments

#         self.activities_mann_whitney = pd.DataFrame(index=self.activities_normalized.index, columns=self.activities_normalized.columns)
#         self.fpr_mann_whitney = pd.DataFrame(index=self.activities_normalized.index, columns=self.activities_normalized.columns)
#         #for every kinase and every dataset, calculate and assemble dataframes of activities and significance values

#         for exp in self.data_columns:
#             log.info("MW Working on %s: "%(exp))

#             #Get a subset of the random and real activites for this experiment
#             activities_sub = self.activities_list[self.activities_list['data']==exp]
#             rand_activities_sub = self.random_kinact.activities_list[self.random_kinact.activities_list['data'].str.startswith(exp)]

#             pval_arr = []
#             fpr_arr = []
#             with concurrent.futures.ProcessPoolExecutor(max_workers=config.PROCESSES) as executor:
#                 for pval, fpr in executor.map(calculate_MannWhitney_one_experiment_one_kinase, repeat(activities_sub), repeat(rand_activities_sub), repeat(self.num_networks), self.activities_normalized.index, repeat(exp), repeat(number_sig_trials)):
#                     pval_arr.append(pval)
#                     fpr_arr.append(fpr)
#                 #print(pval_arr)
#             self.activities_mann_whitney[exp] = pval_arr
#             self.fpr_mann_whitney[exp] = fpr_arr
#             #for kinase in self.normalized_summary.index:
#             #    log.info("\tKinase: %s"%(kinase))
#             #    self.activities_mann_whitney.at[kinase, exp], self.significance_mann_whitney.at[kinase, exp] = calculate_MannWhitney_one_experiment_one_kinase(
#             #        activities_sub, rand_activities_sub, len(self.networks), kinase, exp, number_sig_trials, target_alpha=target_alpha)      


def run_Mann_Whitney_pipeline(activity_list, random_activity_list, normalized_activities, num_networks, num_sig_trials, num_random_experiments, experiment_name, max_cpus):
    if number_sig_trials > num_random_experiments:
        logging.info("Warning: number of trials for Mann Whitney exceeds number available, using %d instead of %d"%(num_random_experiments, number_sig_trials))
        number_sig_trials = num_random_experiments
    activities_mann_whitney = pd.DataFrame(index=activities_normalized.index, columns=activities_normalized.columns)
    fpr_mann_whitney = pd.DataFrame(index=activities_normalized.index, columns=activities_normalized.columns)

    pval_arr = []
    fpr_arr = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_cpus) as executor:
        for pval, fpr in executor.map(calculate_MannWhitney_one_experiment_one_kinase, repeat(activity_list), repeat(random_activity_list), repeat(num_networks), activities_normalized.index, repeat(experiment_name), repeat(number_sig_trials)):
            pval_arr.append(pval)
            fpr_arr.append(fpr)
        #print(pval_arr)
    self.activities_mann_whitney[exp] = pval_arr
    self.fpr_mann_whitney[exp] = fpr_arr

def calculate_fpr_Mann_Whitney(random_kinase_activity_array, number_sig_trials):
    """
    Given an mxn array of kinase activities from m random experiments across n networks
    use bootstrapping to calculate an empirical p-value at which the false positive rate is controlled. 
    This function takes one of m random experiments and calculates the Mann Whitney U pvalue 
    then finds the pvalue at which the target_alpha is achieved
    Parameters
    ----------
    random_kinase_activity_array: np.array
        See calculate_MannWhitney_one_experiment_one_kinase for unwrapping all activities for a kinase and experiment
    number_sig_trials:
        Number of random trials to perform, where one random set is tested number_sig_trials against the full background
    
    Returns
    -------
    random_stats: np.array 
        A vector of Mann Whitney p-values that is m long, representing pvalues from m bootstrap tests
    
    """
    #calculate the significance by taking each experiment 
    [m, n] = random_kinase_activity_array.shape
    if number_sig_trials > m: 
        print("Warning, using %d, maximum number for significance"%(m))
        number_sig_trials = m
    random_stats = np.empty([number_sig_trials])
    for i in range(0, number_sig_trials):
        #take out one vector as real
        sample = random_kinase_activity_array[i,:]
        bgnd = np.delete(random_kinase_activity_array, i, 0) #remove the sample before testing
        [stat, random_stats[i]] = stats.mannwhitneyu(-np.log10(sample), -np.log10(bgnd.reshape(bgnd.size)), alternative='greater')
    return random_stats

def calculate_MannWhitney_one_experiment_one_kinase(kinact_activities, rand_activities, number_networks, kinase, experiment, number_sig_trials):
    """
    For a given kinact object, where random generation and activity has already been run, this will calculate
    the Mann-Whitney U test between the p-values across all networks for the given experiment name 
    and from the random networks. It will also calculate the significance value for the given test 
    based on the target_alpha value by using each random set as a real set to bootstrap. 
    
    Parameters
    ----------
    kinact: kinact object
        A kinact object (not a dictionary)
    kinase: str
        Kinase name to measure significance for
    experiment: str
        Experiment name to measure significance for

    Returns
    -------
    p-value: float
        p-value that results from Mann Whitney U test
    fpr_value: float
        the false positive rate where the p_value for the real experiment lies, given the random experiments
    """
    
    
    kinase_activity_list = kinact_activities[(kinact_activities[config.KSTAR_KINASE]==kinase) & (kinact_activities['data']==experiment)].kinase_activity.values
    
    random_kinase_activity_array = np.empty([number_sig_trials, number_networks])

    for i in range(0, number_sig_trials):
        # get the kinase activity values for all networks for a random set
        experiment_name = experiment+':'+str(i)
        random_kinase_activity_array[i,:]=rand_activities[(rand_activities[config.KSTAR_KINASE]==kinase) & (rand_activities['data']==experiment_name)].kinase_activity.values
       
    [stat, p_value] = stats.mannwhitneyu(-np.log10(kinase_activity_list), -np.log10(random_kinase_activity_array.reshape(random_kinase_activity_array.size)), alternative='greater')
    
    randomStats = calculate_fpr_Mann_Whitney(random_kinase_activity_array, number_sig_trials)
    
    #sig_value = calculate_fpr.single_experiment_kinase_fpr(randomStats, target_alpha)
    fpr_value = calculate_fpr.single_pvalue_fpr(randomStats, p_value)
    #sig_value = 1
    
    return p_value, fpr_value



def parse_args():
    parser = argparse.ArgumentParser(description='Parse Normalization Arguments')
    parser.add_argument( '--activity_list', action='store', dest= 'activity_list', required=True)

    parser.add_argument( '--random_activity_list', action='store', dest= 'random_activity_list', required=True)
    parser.add_argument( '--normalized_activities', action='store', dest= 'normalized_activities', required=True)
    parser.add_argument( '--num_networks', action='store', dest= 'num_networks', required=True,type=int)
    parser.add_argument( '--num_sig_trials', action='store', dest= 'num_sig_trials', required=True,type=int)
    parser.add_argument( '--num_random_experiments', action='store', dest= 'num_random_experiments', required=True, type=int)
    parser.add_argument( '--experiment_name', action='store', dest= 'experiment_name', required=True)
    parser.add_argument( '--max_cpus', action='store', dest= 'max_cpus', required=True, type=int, default=1)


    results = parser.parse_args()
    return results


def main():
    results = parse_args()

    activity_list = pd.read_table(results.activity_list, index_col=0)
    random_activity_list = pd.read_table(results.random_activity_list, index_col=0)
    normalized_activities = pd.read_table(results.normalized_activities, index_col=0)

    
    

    normalized_activity_list.to_csv(f"{data_column}_normalized_activity_list.tsv", sep = "\t")
    normalized_agg_activities.to_csv(f"{data_column}_normalized_aggregate_activity.tsv", sep = "\t")
    activities_normalized.to_csv(f"{data_column}_normalized_activities.tsv", sep = "\t")

if __name__ == "__main__":
    main()