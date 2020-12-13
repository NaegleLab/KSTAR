import pandas as pd
import numpy as np
import scipy.stats as stats
from statsmodels.stats import multitest
import multiprocessing
from collections import defaultdict
import pickle
import os
from datetime import datetime
import itertools
from itertools import repeat
import concurrent.futures

from kstar import config
from kstar.normalize import generate_random_experiments, calculate_fpr

class KinaseActivity:
    def __init__(self, evidence, logger, phospho_type = 'Y'):
        """
        Kinase Activity calculates the estimated activity of kinases given an experiment using hypergeometric distribution.
        Hypergeometric distribution examines the number of protein sites found to be active in evidence compared to the 
        number of protein sites attributed to a kinase on a provided network.

        Parameters
        ----------
        evidence : pandas df
            a dataframe that contains (at minimum, but can have more) data columms as evidence to use in analysis and KSTAR_ACCESSION and KSTAR_SITE
        evidence_columns : dict
            columns corresponding to substrate and site
            required keys : substrate, site
        logger : Logger object
        """
        self.phospho_type = phospho_type
        self.set_evidence(evidence)
        

        self.networks = defaultdict()
        self.network_sizes = defaultdict()
        self.logger = logger
        self.normalizers = defaultdict()

        self.activities = None
        self.agg_activities = None
        self.activity_summary = None

        self.normalized = False
        self.random_experiments = None
        self.random_kinact = None
        self.normalized_activities = None
        self.normalized_agg_activities = None
        self.normalized_summary = None

        self.aggregate = 'count'
        self.threshold = 1.0
        self.greater = True

        self.run_date = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    

    def check_data_columns(self):
        """
        Checks data columns to make sure column is in evidence and that evidence filtered on that data column 
        has at least one point of evidence. Removes all columns that do not meet criteria
        """
        new_data_columns = []
        evidence = self.evidence.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE]).agg(self.aggregate).reset_index()
        for col in self.data_columns:
            if col in self.evidence.columns:
                if self.greater:
                    if len(evidence[evidence[col] >= self.threshold]) > 0:
                        new_data_columns.append(col)
                    else:
                        self.logger.warning(f"{col} does not have any evidence")
                else:
                    if len(evidence[evidence[col] <= self.threshold]) > 0:
                        new_data_columns.append(col)
                    else:
                        self.logger.warning(f"{col} does not have any evidence")
            else:
                self.logger.warning(f"{col} not in evidence")
        self.data_columns = new_data_columns

    def set_data_columns(self, data_columns = None):
        """
        Sets the data columns to use in the kinase activity calculation
        If data_columns is None or an empty list then set data_columns to 
        be all columns that start with data:

        Checks all set columns to make sure columns are vaild after filtering evidence
        """
        if data_columns is None or data_columns == []:
            self.data_columns = []
            for col in self.evidence.columns:
                if col.startswith('data:'):
                    self.data_columns.append(col)  
        else:
            self.data_columns = data_columns
        self.check_data_columns()      
    
    def run_normalization(self, logger, num_random_experiments=150, target_alpha=0.05):
        """
        Run entire normaliation pipeline 
        """
        self.logger.info("Running Normalization Pipeline")
        self.normalized = True
        self.num_random_experiments = num_random_experiments
        self.logger.info("Generating random experiments")
        self.random_experiments = generate_random_experiments.build_random_experiments(
            self.evidence_binary, 
            config.HUMAN_REF_COMPENDIA, 
            num_random_experiments,
            self.phospho_type, 
            self.data_columns )

        
        self.logger.info("Calculating random kinase activities")
        self.random_kinact = KinaseActivity(self.random_experiments, logger, phospho_type=self.phospho_type)
        self.random_kinact.add_networks_batch(self.networks)
        self.random_kinact.calculate_kinase_activities( agg='count', threshold=1.0, greater=True )
        self.random_kinact.aggregate_activities()
        self.random_kinact.activity_summary = self.random_kinact.summarize_activities()

        self.logger.info("Normalizing Activities")
        self.normalizers = calculate_fpr.generate_fpr_values(self.random_kinact.activity_summary, target_alpha)

        self.normalize_activities()
        self.normalized_summary = self.summarize_activities(self.normalized_agg_activities,'median_normalized_activity',normalized=True)

    def add_networks_batch(self,networks):
        for nid, network in networks.items():
            self.add_network(nid, network)

    def add_network(self, network_id, network, network_size = None):
        """
        Add network to be analyzed

        Parameters
        ----------
        network_id : str
            name of the network
        network : pandas DataFrame
            network with columns substrate_id, site, kinase_id
        """
        self.networks[network_id] = network
        if network_size is not None and isinstance(network_size, int):
            self.network_sizes[network_id] = network_size
        else:
            self.network_sizes[network_id] = len(network.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE]).size())
            self.logger.info(f'ADD NETWORK : Number of Accession Sites : {self.network_sizes[network_id]}')
    
    def get_run_date(self):
        return self.run_date

    def set_evidence(self, evidence):
        """
        Evidence to use in analysis

        Parameters
        ----------
        evidence : pandas DataFrame
            substrate sites with activity seen. 
            columns : dict for column mapping
                substrate : Uniprot ID (P12345)
                site : phosphorylation site (Y123)
        """
        columns = list(evidence.columns)
        if config.KSTAR_ACCESSION in columns and config.KSTAR_SITE in columns:
            phospho_type = tuple(self.phospho_type)
            self.evidence = evidence[evidence[config.KSTAR_SITE].str.startswith(phospho_type)].copy()
        else:
            self.logger.warning(f"Evidence not set. Evidence columns must include '{config.KSTAR_ACCESSION}' and '{config.KSTAR_SITE}' keys")
    
    def set_normalizers(self, normalizers):
        """
        Sets the Kinase normalization values to be used in calculating the updated p-values
        Updated p-values normalized via original p-value / normalization-factor * normalization_factor

        Parameters
        ----------
        normalizers : dict or pandas df
            Kinase Activity Normalizers
            index : kinase
            columns : data_column normalization values
        """
        self.normalizers = normalizers

    # def calculate_hypergeometric_single_network(self, evidence, network_id):
    #     """
    #     Hypergeometric Cumulative Distribution Function calculated for each kinase given evidence
    #         k : number of times kinase seen in evidence
    #         M : number of unique sites in network
    #         n : number of times kinase seen in network
    #         N : size of evidence
        
    #     Parameters
    #     ----------
    #     evidence : pandas df
    #         subset of kstar evidence that has been filtered to only include evidence associated with experimetn
    #     network_id : str
    #         network to use for analysis
        
    #     Returns
    #     -------
    #     results : pandas df
    #         Hypergeometric results of evidence for given network
    #         index : kinase_id
    #         columns
    #             frequency : number of times kinase seen in network
    #             kinase_activity : activity derived from hypergometric cdf 
    #     """
        
    #     network = self.networks[network_id]
    #     intersect = pd.merge(network, evidence, how='inner',
    #         on=[config.KSTAR_ACCESSION, config.KSTAR_SITE])
    #     counts = intersect.groupby(config.KSTAR_KINASE).size()
    #     N = len(intersect.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE]).size())
    #     results = pd.DataFrame(counts, columns = ['frequency'])
    #     results['kinase_activity'] = 1.0

    #     K = network.groupby(config.KSTAR_KINASE).size()
        
    #     kinases = counts.index
    #     for kin in kinases:
    #         k = 0
    #         if counts.loc[kin] > 0:
    #             k = counts.loc[kin] - 1
    #         prb = stats.hypergeom.sf(
    #             k = int(k), 
    #             M = int(self.network_sizes[network_id]), 
    #             n = int(K.loc[kin]), 
    #             N = int(N)
    #             )
    #         results.at[kin, 'kinase_activity'] = prb
    #     results['network'] = network_id
    #     return results
        
    # def calculate_hypergeometric_activities(self, evidence, name):
    #     """
    #     Perform hypergeometric kinase activity analysis given evidence on all networks
        
    #     Parameters
    #     ----------
    #     evidence : pandas df
    #         subset of class evidence variable where data is filtered based on experiment
    #     name : str
    #         name of experiment being performed
            
    #     Returns
    #     ---------
    #     fdr_act : pd DataFrame
    #         network : network name, from networks key
    #         frequency : number of times kinase was seen in subgraph of evidence and network
    #         kinase_activity : hypergeometric kinase activity
    #         fdr_corrected_kinase_activity : kinase activity after fdr correction
    #         significant : whether kinase activity is significant based on fdr alpha
    #     combined : pd DataFrame
    #         significant : number of networks where kinase was found to be significant
    #         fraction_significant : fraction of networks kinase was found to be significant through FDR correction
    #         avg_p_value : combined p-values of kinase using mean
    #         median_p_value : combined p-values of kinase using median
    #     """
        
    #     self.logger.info(f"Running hypergeometric analysis on {name}")
        
    #     results = []
    #     for network_id in self.networks.keys(): # calculate kinase activity within each network 
    #         result = self.calculate_hypergeometric_single_network(evidence, network_id) 
    #         results.append(result)

    #     # combine results into single dataframe
    #     hyp_act = pd.concat(results)
    #     hyp_act = hyp_act.reset_index()
    #     hyp_act['data'] = name
    #     return hyp_act

    def create_binary_evidence(self, data_columns = None, agg = 'count', threshold = 1.0,  greater = True):
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

        evidence = self.evidence.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE]).agg(agg).reset_index()
        
        #set the binary evidence for whether a site is included
        evidence_binary = evidence.copy()
        for col in self.data_columns:
            if greater:
                evidence_binary[col].mask(evidence[col] >= threshold, 1, inplace=True)
                evidence_binary[col].mask(evidence[col] < threshold, 0, inplace=True)
            else:
                evidence_binary[col].mask(evidence[col] <= threshold, 1, inplace=True)
                evidence_binary[col].mask(evidence[col] > threshold, 0, inplace=True)

        #remove phosphorylation sites that were not selected in any experiment (useful for very large experiments where removing the need to copy data reduces time)
        evidence_binary.drop(evidence_binary[evidence_binary[self.data_columns].sum(axis=1) == 0].index, inplace = True) 
        return evidence_binary


    def calculate_kinase_activities(self, data_columns = None, agg = 'count', threshold = 1.0,  greater = True):
        """
        Calculates combined activity of experiments based that uses a threshold value to determine if an experiment sees a site or not
        To use values use 'mean' as agg
            mean aggregation drops NA values from consideration
        To use count use 'count' as agg - present if not na

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
        activities : dict
            key : experiment
            value : pd DataFrame
                network : network name, from networks key
                kinase : kinase examined
                frequency : number of times kinase was seen in subgraph of evidence and network
                kinase_activity : hypergeometric kinase activity
        
        """

        #for each phosphoType, check that the 
        self.aggregate = agg
        self.threshold = threshold
        self.greater = greater
        self.set_data_columns(data_columns)

        self.evidence_binary = self.create_binary_evidence(data_columns = self.data_columns, agg = self.aggregate, threshold = self.threshold,  greater = self.greater)
        
        # if no data columns are provided use all columns that start with data:
        # data columns that filtered have no evidence are removed
        self.logger.info(f"Kinase Activity will be run on the following data columns: {','.join(self.data_columns)}")
        

        # MULTIPROCESSING
        if config.PROCESSES > 1:
            manager = multiprocessing.Manager()
            pool = multiprocessing.Pool(processes = config.PROCESSES)
            filtered_evidence_list  = [self.evidence_binary[self.evidence_binary[col] ==1 ] for col in self.data_columns] 
            networks = itertools.repeat(self.networks)  
            network_sizes = itertools.repeat(self.network_sizes)
            iterable = zip(filtered_evidence_list, networks, network_sizes, self.data_columns)
            activities_list = pool.starmap(calculate_hypergeometric_activities, iterable)
        
        # SINGLE CORE PROCESSING
        else:
            activities_list =[]
            for col in self.data_columns:
                filtered_evidence = self.evidence_binary[self.evidence_binary[col] == 1]
                act = calculate_hypergeometric_activities(filtered_evidence, self.networks, self.network_sizes, col)
                act['data'] = col
                activities_list.append(act)

        self.activities = pd.concat(activities_list)
        return self.activities

    def summarize_activities(self, activities = None, method = 'median_activity', normalized = False):
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
            The column in the hypergeometric activity to use for summaring data
        
        Returns
        ---------
        activity_summary : pandas DataFrame

        """
        if activities is None:
            activities = self.agg_activities
        available_methods = list(activities.columns)
        available_methods.remove('data')
        if method not in available_methods:
            raise ValueError(f"the method '{method}' is not in the availble methods. \nAvailable methods include : {', '.join(available_methods)}")


        activity_summary = activities.pivot(index = config.KSTAR_KINASE, columns ='data', values = method).reset_index().rename_axis(None, axis=1).set_index(config.KSTAR_KINASE)
        activity_summary = activity_summary[self.data_columns]
        return activity_summary

    def aggregate_activities(self, activities = None):
        """
        Aggregate network activity using median for all activities

        Parameters
        ---------
        activities : dict
            key : Experiment
            value : kinase activity result
        
        Returns 
        ---------
        summaries : dict
            key : experiment
            value : summarized kinase activities accross networks
        """
        if activities is None:
            activities = self.activities
        self.agg_activities = activities.groupby(['data', config.KSTAR_KINASE ]).agg(
            median_activity = ('kinase_activity', 'median'),
        ).reset_index()
        return self.agg_activities

    def normalize_activities(self, activities = None, default_normalization = 0.05, normalization_multiplier = 0.05):
        """
        Generates the normalized activities and corresponding summary statistics

        Parameters
        ----------
        activities : dict
            hypergeometric activities generated by KSTAR algorithm
            key : experiment
            value : hypergeometric result
        
        Returns
        --------
        normalized_activities : dict
            normalized Kinase Activity results
            key : experiment
            value : normalized results
        summarized_activities : dict
            summary statistics of normalized kinase activity
            key : experiment
            value : summary statistics
        """
        self.normalized = True
        
        if activities is None:
            activities = self.activities
        normalized_activities = []
        for data in self.data_columns:
            if type(self.normalizers) is dict:
                normalizers = self.normalizers
            elif type(self.normalizers) is pd.DataFrame and data in self.normalizers.columns:
                normalizers = self.normalizers[data].to_dict()
            else:
                normalizers = {}
            activity = activities[activities['data'] == data]
            normalized_activities.append(
                self.calculate_normalized_activity(
                    activity, 
                    normalizers, 
                    default_normalization, 
                    normalization_multiplier))
        self.normalized_activities = pd.concat(normalized_activities)
        self.aggregate_normalized_activities(self.normalized_activities)
        
        return self.normalized_activities, self.normalized_agg_activities

    def calculate_normalized_activity(self, kinase_activity, normalizers, default_normalization = 0.05, normalization_multiplier = 0.05):
        """
        Activity is normalized based on kinase normalization factors. 
        Added Columns : 
            Normalization Factor : normalization factor for given Kinase
            Significant : 1 if Kinase Activity <= Normalization Factor
            Normalized Activity : ( Kinase Activity ) / ( Normalization Factor ) * Normalization Multiplier

        Parameters
        ----------
        kinase_activity : pandas df
            hypergeometric activity calculated through KSTAR algorithm
        default_normalization : float
            Normalization factor to use if one not provided
        normalization_multiplier : float
            normalization multiplier to use for calculated normalized kinase activity
        """
        normalized_activity = kinase_activity.copy()
        normalized_activity['Normalization Factor'] = normalized_activity[config.KSTAR_KINASE].apply(lambda name: normalizers[name] if name in normalizers.keys() else default_normalization)
        if len(self.networks) > 0:
            normalized_activity['Significant'] = normalized_activity.apply(lambda row: (row['kinase_activity'] <= row['Normalization Factor'] / len( self.networks ) ) * 1, axis = 1)
        else:
            normalized_activity['Significant'] = normalized_activity.apply(lambda row: (row['kinase_activity'] <= row['Normalization Factor']) * 1, axis = 1)
        normalized_activity['Normalized Activity'] = normalized_activity.apply(lambda row: np.clip( row['kinase_activity'] / row['Normalization Factor'] * normalization_multiplier, 0.0, 1.0 ), axis = 1)
        return normalized_activity

    def aggregate_normalized_activities(self, normalized_activities):
        """
        Summarizes normalized kinase activity results of all networks by kinase
        Summaries provided : 
            median original activity
            average original activity
            median normalized activity
            average normalized activity
            count significant
            fraction significant
        
        Parameters
        ----------
        normalized_activity : pandas df
            Normalized kinase activty
        
        Returns
        --------
        summary : pandas df
            summarized data of all networks by kinase
        """
        self.normalized_agg_activities = normalized_activities.groupby(['data',config.KSTAR_KINASE ]).agg(
            median_original_activity = ('kinase_activity', 'median'),
            median_normalized_activity = ('Normalized Activity', 'median'),
            count_significant = ('Significant', 'sum'),
            fraction_significant = ('Significant', 'mean')
        ).reset_index()
        return self.normalized_agg_activities

    def  find_pvalue_limits(self, data_columns, agg = 'count', threshold = 1.0):
        """
        For each data column and network find the lowest p-value achievable and how many 
        seen sites are required to get to that limit.
        Assumptions
            - kinase size in network is same for all kinases

        Parameters
        ----------
        data_columns : list
            what columns in evidence to compare
        agg : str
            aggregate function - what function to use for determining if site is present
                count : use when using activity_count
                mean : use when using activity_threshold
        threshold : float
            threshold to use in determining if site present in evidence
        
        Returns
        -------
        all_limits : pandas DataFrame
            p-value limits of each column for each network
            columns:
                evidence        evidence data column
                network         network being compared
                kinase          kinase being evaluated
                evidence_size   size of evidence
                limit_size      number of sites to get non-zero p-value
                p-value          p-value generated
        limit_summary : pandas DataFrame
            summary of all_limits by taking average over by evidence
        """

        
        evidence = self.evidence.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE]).agg(agg)
        all_rows = []

        kinase_sizes = {}
        for nid, network in self.networks.items():
            # kinase_sizes[nid] = network.groupby(config.KSTAR_KINASE).size()
            kinase_sizes[nid] = int(network.groupby(config.KSTAR_KINASE).size().mean())
            

        for col in data_columns:
            print(col)
            filtered_evidence = evidence[evidence[col] > threshold]
            N = len(filtered_evidence.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE]).size())
            for nid, M in self.network_sizes.items():
                # kinases = list(kinase_sizes[nid].index)
                # for kinase in kinases:
                # n = int(kinase_sizes[nid].loc[kinase])
                n = kinase_sizes[nid]

                for k in range(N,0,-1):
                    pvalue = stats.hypergeom.sf(
                        k = k, 
                        M = M, 
                        n = n, 
                        N = N
                    )
                    # pvalue = 1 - prb
                    if pvalue > 0:
                        break
                row = {
                    'evidence' : col,
                    'network' : nid,
                    # 'kinase' : kinase,
                    'evidence_size' : N,
                    'limit_size' : k,
                    'kinase_size' : n,
                    'network_size' : M,
                    'p-value' : pvalue
                    }
                all_rows.append(row)
        all_limits = pd.DataFrame(all_rows)
        limit_summary = all_limits.groupby('evidence').mean()
        return all_limits, limit_summary

    def calculate_Mann_Whitney_activities_sig(self, log, number_sig_trials = 100, target_alpha=0.05):
        """
        For a kinact_dict, where random generation and activity has already been run for the phospho_types of interest, 
        this will calculate the Mann-Whitney U test for comparing the array of p-values for real data 
        to those of random data, across the number of networks used.
        It will also calculate the significance value for the given test 
        based on the target_alpha value by using each random set as a real set to bootstrap, using number_sig_trials.
        
        Parameters
        ----------
        kinact_dict: dictionary
            A dictionary of kinact objects, with keys 'Y' and/or 'ST'
        log: logger 
            Logger for logging activity messages
        phospho_types: {['Y', 'ST'], ['Y'], ['ST']}
            Which substrate/kinaset-type to run activity for: Both ['Y, 'ST'] (default), Tyrosine ['Y'], or Serine/Threonine ['ST']
        number_sig_trials: int
            Maximum number of significant trials to run
        target_alpha: float
            target false positive rate, value between 0 and 1, default =0.05
            
        
        Returns
        -------

        """
        #First, check that objects are correct and values can be found
        if not isinstance(self.random_kinact.activities, pd.DataFrame):
            raise ValueError("Random activities do not exist, please run kstar_activity.normalize_analysis")

        if number_sig_trials > self.num_random_experiments:
            log.info("Warning: number of trials for Mann Whitney exceeds number available, using %d instead of %d"%(self.num_random_experiments, number_sig_trials))
            number_sig_trials = self.num_random_experiments

        self.activities_mann_whitney = pd.DataFrame(index=self.normalized_summary.index, columns=self.normalized_summary.columns)
        self.significance_mann_whitney = pd.DataFrame(index=self.normalized_summary.index, columns=self.normalized_summary.columns)
        #for every kinase and every dataset, calculate and assemble dataframes of activities and significance values

        for exp in self.data_columns:
            log.info("MW Working on %s: "%(exp))

            #Get a subset of the random and real activites for this experiment
            activities_sub = self.activities[self.activities['data']==exp]
            rand_activities_sub = self.random_kinact.activities[self.random_kinact.activities['data'].str.startswith(exp)]

            pval_arr = []
            sig_arr = []
            with concurrent.futures.ProcessPoolExecutor(max_workers=config.PROCESSES) as executor:
                for pval, sig in executor.map(calculate_MannWhitney_one_experiment_one_kinase, repeat(activities_sub), repeat(rand_activities_sub), repeat(len(self.networks)), self.normalized_summary.index, repeat(exp), repeat(number_sig_trials), repeat(target_alpha)):
                    pval_arr.append(pval)
                    sig_arr.append(sig)
                #print(pval_arr)
            self.activities_mann_whitney[exp] = pval_arr
            self.significance_mann_whitney[exp] = sig_arr
            #for kinase in self.normalized_summary.index:
            #    log.info("\tKinase: %s"%(kinase))
            #    self.activities_mann_whitney.at[kinase, exp], self.significance_mann_whitney.at[kinase, exp] = calculate_MannWhitney_one_experiment_one_kinase(
            #        activities_sub, rand_activities_sub, len(self.networks), kinase, exp, number_sig_trials, target_alpha=target_alpha)      




def calculate_hypergeometric_single_network(evidence, network, network_size, network_id ):
        """
        Hypergeometric Cumulative Distribution Function calculated for each kinase given evidence
            k : number of times kinase seen in evidence
            M : number of unique sites in network
            n : number of times kinase seen in network
            N : size of evidence
        
        Parameters
        ----------
        evidence : pandas df
            subset of kstar evidence that has been filtered to only include evidence associated with experiment
        network_id : str
            network to use for analysis
        
        Returns
        -------
        results : pandas df
            Hypergeometric results of evidence for given network
            index : kinase_id
            columns
                frequency : number of times kinase seen in network
                kinase_activity : activity derived from hypergometric cdf 
        """
        
        
        intersect = pd.merge(network, evidence, how='inner',
            on=[config.KSTAR_ACCESSION, config.KSTAR_SITE])
        
        counts = intersect.groupby(config.KSTAR_KINASE).size()
        N = len(intersect.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE]).size())
        results = pd.DataFrame(counts, columns = ['frequency'])

        results['kinase_activity'] = 1.0

        K = network.groupby(config.KSTAR_KINASE).size()
        
        kinases = counts.index
        for kin in kinases:
            k = 0
            if counts.loc[kin] > 0:
                k = counts.loc[kin] - 1
            prb = stats.hypergeom.sf(
                k = int(k), 
                M = int(network_size), 
                n = int(K.loc[kin]), 
                N = int(N)
                )
            results.at[kin, 'kinase_activity'] = prb

        kinases = network[config.KSTAR_KINASE].unique()
        for kin in kinases:
            if kin not in results.index:
                results.at[kin,'frequency'] = 0
                results.at[kin,'kinase_activity'] = 1.0
        results['network'] = network_id
        return results
        
def calculate_hypergeometric_activities(evidence, networks, network_sizes, name):
        """
        Perform hypergeometric kinase activity analysis given evidence on all networks
        
        Parameters
        ----------
        evidence : pandas df
            subset of class evidence variable where data is filtered based on experiment
        name : str
            name of experiment being performed
            
        Returns
        ---------
        fdr_act : pd DataFrame
            network : network name, from networks key
            frequency : number of times kinase was seen in subgraph of evidence and network
            kinase_activity : hypergeometric kinase activity
            fdr_corrected_kinase_activity : kinase activity after fdr correction
            significant : whether kinase activity is significant based on fdr alpha
        combined : pd DataFrame
            significant : number of networks where kinase was found to be significant
            fraction_significant : fraction of networks kinase was found to be significant through FDR correction
            avg_p_value : combined p-values of kinase using mean
            median_p_value : combined p-values of kinase using median
        """
        
        # self.logger.info(f"Running hypergeometric analysis on {name}")
        
        results = []
        for network_id in networks.keys(): # calculate kinase activity within each network 
            result =calculate_hypergeometric_single_network(evidence, networks[network_id], network_sizes[network_id], network_id) 
            results.append(result)

        # combine results into single dataframe
        hyp_act = pd.concat(results)
        hyp_act = hyp_act.reset_index()
        hyp_act['data'] = name
        return hyp_act

"""
****************************************
Methods for Mann Whitney analysis
****************************************
"""

def calculate_fpr_Mann_Whitney(random_kinase_activity_array, number_sig_trials, target_alpha=0.05):
    """
    Given an mxn array of kinase activities from m random experiments across n networks
    use bootstrapping to calculate an empirical p-value at which the false positive rate is controlled. 
    This function takes one of m random experiments and calculates the Mann Whitney U pvalue 
    then finds the pvalue at which the target_alpha is achieved
    Parameters
    ----------
    random_kinase_activity_array: np.array
        See calculate_MannWhitney_one_experiment_one_kinase for unwrapping all activities for a kinase and experiment
    target_alpha: float
        false positive rate between 0 and 1
    
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

def calculate_MannWhitney_one_experiment_one_kinase(kinact_activities, rand_activities, number_networks, kinase, experiment, number_sig_trials, target_alpha=0.05):
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
    target_alpha: float
        target false positive rate, value between 0 and 1, default =0.05
        
    
    Returns
    -------
    p-value: float
        p-value that results from Mann Whitney U test
    fpr_significance: float
        the p-value where target_alpha is achieved
    """
    
    
    kinase_activity_list = kinact_activities[(kinact_activities[config.KSTAR_KINASE]==kinase) & (kinact_activities['data']==experiment)].kinase_activity.values
    
    random_kinase_activity_array = np.empty([number_sig_trials, number_networks])

    for i in range(0, number_sig_trials):
        # get the kinase activity values for all networks for a random set
        experiment_name = experiment+':'+str(i)
        random_kinase_activity_array[i,:]=rand_activities[(rand_activities[config.KSTAR_KINASE]==kinase) & (rand_activities['data']==experiment_name)].kinase_activity.values
       
    [stat, p_value] = stats.mannwhitneyu(-np.log10(kinase_activity_list), -np.log10(random_kinase_activity_array.reshape(random_kinase_activity_array.size)), alternative='greater')
    
    randomStats = calculate_fpr_Mann_Whitney(random_kinase_activity_array, number_sig_trials, target_alpha=target_alpha)
    
    sig_value = calculate_fpr.single_experiment_kinase_fpr(randomStats, target_alpha)
    #sig_value = 1
    
    return p_value, sig_value

"""
****************************************
Methods for running KSTAR pipeline
****************************************
"""
def run_kstar_analysis(experiment, log, networks, phospho_types =['Y', 'ST'], data_columns = None, agg = 'count', threshold = 1.0,  greater = True):
    """
    A super method to establish a kstar KinaseActivity object from an experiment with an activity log
    add the networks, calculate, aggregate, and summarize the activities into a final activity object

    Parameters
    ----------
    experiment: pandas df
        experiment dataframe that has been mapped, includes KSTAR_SITE, KSTAR_ACCESSION, etc.
    log: logger object
        Log to write activity log error and update to
    networks: dictionary of dictionaries
        Outer dictionary keys are 'Y' and 'ST'.
        Establish a network by loading a pickle of desired networks. See the helpers and config file for this.
        If downloaded from FigShare, then the GLOBAL network pickles in config file can be loaded
        For example: networks['Y'] = pickle.load(open(config.NETWORK_Y_PICKLE, "rb" ))
    phospho_types: {['Y', 'ST'], ['Y'], ['ST']}
        Which substrate/kinaset-type to run activity for: Both ['Y, 'ST'] (default), Tyrosine ['Y'], or Serine/Threonine ['ST']
    data_columns : list
        columns that represent experimental result, if None, takes the columns that start with `data:'' in experiment. 
        Pass this value in as a list, if seeking to calculate on fewer than all available data columns
    agg : {'count', 'mean'}
        method to use when aggregating duplicate substrate-sites. 
        'count' combines multiple representations and adds if values are non-NaN
        'mean' uses the mean value of numerical data from multiple representations of the same peptide.
            NA values are droped from consideration.
    threshold : float
        threshold value used to filter rows
    greater: Boolean
        whether to keep sites that have a numerical value >=threshold (TRUE, default) or <=threshold (FALSE)
    
    Returns
    -------
    kinactDict: dictionary of Kinase Activity Objects
        Outer keys are phosphoTypes run 'Y' and 'ST'
        Includes the activities dictionary (see calculate_kinase_activities)
        aggregation of activities across networks (see aggregate activities)
        activity summary (see summarize_activities)

    """


    kinact_dict = {}

    # For each phosphoType of interest, establish a kinase activity object on a filtered dataset and run, aggregate, and summarize activity
    for phospho_type in phospho_types:
        #first check that networks for the phosphotypes were passed in
        if phospho_type not in networks:
            print("ERROR: Please pass networks as dictionary with phosphotype key")
        #filter the experiment (log how many are of that type)
        if phospho_type == 'ST':
            experiment_sub = experiment[(experiment.KSTAR_SITE.str.contains('S')) | (experiment.KSTAR_SITE.str.contains('T'))]
            log.info("Running Serine/Threonine Kinase Activity Analysis")
        elif phospho_type == 'Y':
            experiment_sub = experiment[(experiment.KSTAR_SITE.str.contains('Y'))]
            log.info("Running Tyrosine Kinase Activity Analysis")

        else:
            print("ERROR: Did not recognize phosphoType %s, which should only include 'Y' or 'ST' "%(phospho_type))
            return
        kinact = KinaseActivity(experiment_sub, log, phospho_type=phospho_type)
        kinact.add_networks_batch(networks[phospho_type])
        kinact.calculate_kinase_activities(data_columns, agg=agg, threshold=threshold, greater=greater)
        kinact.aggregate_activities()
        kinact.activity_summary = kinact.summarize_activities()
        kinact_dict[phospho_type] = kinact
    return kinact_dict

def normalize_analysis(kinact_dict, log, num_random_experiments=150, target_alpha = 0.05):
    """
    Creates random experiments, drawn from the human phosphoproteome, according to the distribution of the number of compendia
    that each data column in the experiment has for num_random_experiments. Kinase activity calculation is then run on every random experiment
    and then assembled for each data column to calculate the pvalue at which a kinase hits the target_alpha value of false positives. 
    Finally, it takes these empirically defined pvalue corrections and normalizes the kinase activities according to these, such that
    the real experimental data also has the target_alpha value. 

    Params
    ------
    kinact_dict: KinaseActivities dictionary
        Has keys ['Y'] and/or ['ST'] and values that are KinaseActivity objects. These objects are modified to add normalization
    log: logger 
        Logger for logging activity messages
    num_random_experiments: int
        Number of random experiments, for each data column, to create and run activity from
    target_alpha: float 
        Value between 0 and 1 that is the target false positive rate. 

    Returns
    -------


    """

    if target_alpha > 1 or target_alpha < 0: 
        raise ValueError('ERROR: target_alpha must be value between 0 and 1')

    for phospho_type, kinact in kinact_dict.items():
        kinact.run_normalization(log, num_random_experiments, target_alpha)
    
def save_kstar(kinact_dict, name, odir, PICKLE=True):
        """
        Having performed kinase activities (run_kstar_analyis), save each of the important dataframes to files and the final pickle
        Saves an activities, aggregated_activities, summarized_activities tab-separated files
        Saves a pickle file of dictionary 
        
        Parameters
        ----------
        kinact_dict: dictionary of Kinase Activity Objects
            Outer keys are phosphoTypes run 'Y' and 'ST'
            Includes the activities dictionary (see calculate_kinase_activities)
            aggregation of activities across networks (see aggregate activities)
            activity summary (see summarize_activities)
        name: string 
            The name to use when saving activities
        odir:  string
            Outputdirectory to save files and pickle to
        PICKLE: boolean
            Whether to save the entire pickle file 

        Returns
        -------
        Nothing

        """

        if not os.path.exists(f"{odir}/RESULTS"): 
            os.mkdir(f"{odir}/RESULTS") 

        if not os.path.exists(f"{odir}/RANDOM_ANALYSIS"): 
            os.mkdir(f"{odir}/RANDOM_ANALYSIS") 
        
        for phospho_type in kinact_dict:
            kinact = kinact_dict[phospho_type]
            name_out = f"{name}_{phospho_type}"
            kinact.activities.to_csv(f"{odir}/RESULTS/{name_out}_activities.tsv", sep = '\t', index = False)
            kinact.agg_activities.to_csv(f"{odir}/RESULTS/{name_out}_aggregated_activities.tsv", sep = '\t', index = False)
            kinact.activity_summary.to_csv(f"{odir}/RESULTS/{name_out}_summarized_activities.tsv", sep = '\t', index = False)
            kinact.evidence_binary.to_csv(f"{odir}/RESULTS/{name_out}_binarized_experiment.tsv", sep='\t', index=False)

            if kinact.normalized:

                kinact.random_kinact.activities.to_csv(f"{odir}/RANDOM_ANALYSIS/{name_out}_random_activities.tsv", sep = '\t', index = False)
                kinact.random_kinact.agg_activities.to_csv(f"{odir}/RANDOM_ANALYSIS/{name_out}_random_aggregated_activities.tsv", sep = '\t', index = False)
                kinact.random_kinact.activity_summary.to_csv(f"{odir}/RANDOM_ANALYSIS/{name_out}_random_summarized_activities.tsv", sep = '\t', index = False)

                kinact.normalized_activities.to_csv(f"{odir}/RESULTS/{name_out}_normalized_activities.tsv", sep = '\t', index = False)
                kinact.normalized_agg_activities.to_csv(f"{odir}/RESULTS/{name_out}_normalized_aggregated_activities.tsv", sep = '\t', index = False)
                kinact.normalized_summary.to_csv(f"{odir}/RESULTS/{name_out}_normalized_summarized_activities.tsv", sep = '\t', index = False)

                kinact.random_experiments.to_csv(f"{odir}/RANDOM_ANALYSIS/{name_out}_random_experiments.tsv", sep = '\t', index = False)

            if hasattr(kinact, 'activities_mann_whitney'):
                kinact.activities_mann_whitney.to_csv(f"{odir}/RESULTS/{name_out}_mann_whitney_activities.tsv", sep='\t', index = False) 
                kinact.significance_mann_whitney.to_csv(f"{odir}/RESULTS/{name_out}_mann_whitney_significance.tsv", sep='\t', index = False) 


        pickle.dump( kinact_dict, open( f"{odir}/RESULTS/{name}_kinact.p", "wb" ) )

