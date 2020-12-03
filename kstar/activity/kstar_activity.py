import pandas as pd
import numpy as np
import scipy.stats as stats
from statsmodels.stats import multitest
import multiprocessing
from collections import defaultdict
import pickle
import os
from datetime import datetime

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
        
        self.set_data_columns()

        self.run_date = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    

    def check_data_columns(self, data_columns):
        new_data_columns = []
        evidence = self.evidence.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE]).agg(self.aggregate).reset_index()
        for col in data_columns:
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
        return new_data_columns

    def set_data_columns(self):
        self.data_columns = []
        for col in self.evidence.columns:
            if col.startswith('data:'):
                self.data_columns.append(col)        

    def set_network_columns(self, network_columns):
        """
        sets the network columns to use in analysis
        Parameters
        ----------
        network_columns : dict
            required keys : substrate, site, kinase
        """
        if 'substrate' in network_columns.keys() and 'site' in network_columns.keys() and 'kinase' in network_columns.keys():
            self.network_columns = network_columns
            self.logger.info("Network Columns Changed : " + network_columns)
        else:
            self.logger.warning("Network Columns Not Changed. Columns must include substrate, site, and kinase keys")
    
    def run_normalization(self, logger, num_random_experiments=150, target_alpha=0.05):
        """
        Run entire normaliation pipeline 
        """
        self.logger.info("Running Normalization Pipeline")
        self.normalized = True
        self.logger.info("Generation random experiments")
        self.random_experiments = generate_random_experiments.build_random_experiments(
            self.evidence, 
            config.HUMAN_REF_COMPENDIA, 
            self.aggregate, 
            self.threshold,
            self.greater, 
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
            key : str
                Kinase Name
            value : float
                Normalization Factor
        """
        self.normalizers = normalizers


    def calculate_hypergeometric_single_network(self, evidence, network_id):
        """
        Hypergeometric Cumulative Distribution Function calculated for each kinase given evidence
            k : number of times kinase seen in evidence
            M : number of unique sites in network
            n : number of times kinase seen in network
            N : size of evidence
        
        Parameters
        ----------
        evidence : pandas df
            subset of kstar evidence that has been filtered to only include evidence associated with experimetn
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
        
        network = self.networks[network_id]
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
                M = int(self.network_sizes[network_id]), 
                n = int(K.loc[kin]), 
                N = int(N)
                )
            results.at[kin, 'kinase_activity'] = prb
        results['network'] = network_id
        return results
        

    def calculate_hypergeometric_activities(self, evidence, name):
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
        
        self.logger.info(f"Running hypergeometric analysis on {name}")
        
        results = []
        for network_id in self.networks.keys(): # calculate kinase activity within each network 
            result = self.calculate_hypergeometric_single_network(evidence, network_id) 
            results.append(result)

        # combine results into single dataframe
        hyp_act = pd.concat(results)
        hyp_act = hyp_act.reset_index()
        hyp_act['data'] = name
        return hyp_act


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

        if data_columns is None:
            data_columns = self.data_columns
        
        self.data_columns = self.check_data_columns(self.data_columns)

        #log here what datacolumns will be operated on, report number of experiments found

        

        evidence = self.evidence.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE]).agg(agg).reset_index()

        # MULTIPROCESSING
        num_cores_available = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes = num_cores_available - 2)
        if greater:
            filtered_evidence_list  = [evidence[evidence[col] >= threshold] for col in self.data_columns]
            iterable = zip(filtered_evidence_list, self.data_columns)
        else:
            filtered_evidence_list  = [evidence[evidence[col] <= threshold] for col in self.data_columns]
            iterable = zip(filtered_evidence_list, self.data_columns)
        activities_list = pool.starmap(self.calculate_hypergeometric_activities, iterable)
        
        # SINGLE CORE PROCESSING
        # for col in data_columns:
        #     if greater:
        #         filtered_evidence = evidence[evidence[col] > threshold]
        #     else:
        #         filtered_evidence = evidence[evidence[col] < threshold]
        #     act = self.calculate_hypergeometric_activities(filtered_evidence)
        #     act['data'] = col
        #     activities_list.append(act)
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
    

    def calculate_fdr(self, activities, default_alpha = 0.05):
        """
        Given raw kinase activities and the normalizaiton factors, FDR correction is 
        performed on activities for all networks using normalizaiton factor
        as FDR alpha.
        """
        exp_activities = pd.concat(activities, names = ['Data','__drop__']).reset_index().drop(columns='__drop__')
        kinases = exp_activities[config.KSTAR_KINASE].unique()
        
        norm_activities = []

        for dataset in self.data_columns:
            activity = self.activities[self.activities['data'] == dataset]
            for kinase in kinases:
                if type(self.normalizers) is pd.DataFrame:
                    if kinase in self.normalizers.index and dataset in self.normalizers.columns:
                        alpha = self.normalizers.loc[kinase][dataset]
                    else:
                        alpha = default_alpha
                        self.logger.warning(f"FDR Correction : Normalization not found {dataset}:{kinase}. Default alpha of {default_alpha} used")
                else:
                    alpha = self.normalizers[kinase] if kinase in self.normalizers.keys() else default_alpha
                
                norm_activity = activity[activity[config.KSTAR_KINASE] == kinase]
                rej, cor, _, _ = multitest.multipletests(norm_activity['kinase_activity'], method = 'fdr_tsbky', alpha=alpha)
                norm_activity['fdr_activity'] = cor
                norm_activity['fdr_significant'] = rej * 1
                norm_activities.append(norm_activity)
        
        norm_activities = pd.concat(norm_activities, ignore_index=True)
        return norm_activities
        

    def aggregate_fdr(self, activities):
        agg_activities =defaultdict()
        for key, activity in activities.items():
            agg_activities[key] = activity.groupby(['data',config.KSTAR_KINASE]).agg(
            median_activity = ('kinase_activity', 'median'),
            median_fdr_activity = ('fdr_activity','median')).reset_index()
        return agg_activities

    def complete_fdr_response(self, response, alpha):
        """
        FDR correction performed on response dataframe by first melting all experimental columns before 
        performing FDR correction through Two-stage Benjamini, Krieger, & Yekutieli FDR procedure

        Parameters
        ----------
        response : pandas df
        alpha : float
            fdr alpha to use for significance
        
        Returns
        ----------
        fdr : pandas df
            fdr corrected values
        significant : pandas df
            binary significant/not based on FDR 
        """
        reshaped = pd.melt(response.reset_index(), id_vars='index')

        rej, cor, _, _ = multitest.multipletests(reshaped['value'], method = 'fdr_tsbky', alpha=alpha)
        reshaped['significant'] = rej * 1
        reshaped['fdr'] = cor

        pivot = reshaped.pivot(index = 'index', columns = 'variable')

        fdr = pivot['fdr']
        fdr.columns.name = None
        
        significant = pivot['significant']
        significant.columns.name = None

        return fdr, significant


    def experiment_fdr_response(self, response, alpha):
        fdr = pd.DataFrame(index = response.index)
        significant = pd.DataFrame(index = response.index)
        for experiment in response.columns:
            rej, cor, _, _ = multitest.multipletests(response[experiment], method = 'fdr_tsbky', alpha=alpha)
            fdr[experiment] = cor
            significant[experiment] = rej*1
        return fdr, significant


    def fdr_corrected_activity(self, kinase_activity, alpha):
        """
        Performs False Discovery Rate correction using Two-stage Benjamini, Krieger, & Yekutieli FDR procedure
        Benjamini, Yoav, Abba M. Krieger, and Daniel Yekutieli. 2006. 
        “Adaptive Linear Step-up Procedures That Control the False Discovery Rate.” 
        Biometrika 93 (3) (September 1): 491–507. doi:10.1093/biomet/93.3.491.
    
        Parameters
        ----------
        kinase_activity : pandas df
            kinase activity p-values
        alpha : flaot
            FWER, family-wise error rate
        
        Returns 
        -------
        results : pandas df
            kinase_activity dataframe with added columns 
                significant : 1 if p-value is statistically significant given alpha
                fdr_corrected_kinase_activity : fdr corrected activity by tsbky method
        """
        results = kinase_activity.copy()
        rej, cor, _, _ = multitest.multipletests(results['kinase_activity'], method = 'fdr_tsbky', alpha=alpha)
        results['fdr_significant'] = rej * 1
        results['fdr_activity'] = cor
        return results

def run_kstar_analysis(experiment, log, networks, phospho_types =['Y', 'ST'], data_columns = None, agg = 'count', threshold = 1.0,  greater = True):
    """
    A super method to establish a kstar KinaseActivity object from an experiment with an activity log
    add the networks, calculate, aggregate, and summarize the activities into a final activity object

    Parameters
    ----------
    experiment: pandas df
        experiment dataframe that has been mapped, includes KSTAR_SITE, KSTAR_ACCESSION, etc.
    log: string
        File name to write activity log error and update to
    networks: dictionary of dictionaries
        Outer dictionary keys are 'Y' and 'ST'.
        Establish a network by loading a pickle of desired networks. See the helpers and config file for this.
        If downloaded from FigShare, then the GLOBAL network pickles in config file can be loaded
        For example: networks['Y'] = pickle.load(open(config.NETWORK_Y_PICKLE, "rb" ))
    phosphoTypes: {['Y', 'ST'], ['Y'], ['ST']}
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
            log.log("Running Serine/Threonine Kinase Activity Analysis")
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
    


def save_kstar(kinact_dict, name, odir):
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

            if kinact.normalized:

                kinact.random_kinact.activities.to_csv(f"{odir}/RANDOM_ANALYSIS/{name_out}_random_activities.tsv", sep = '\t', index = False)
                kinact.random_kinact.agg_activities.to_csv(f"{odir}/RANDOM_ANALYSIS/{name_out}_random_aggregated_activities.tsv", sep = '\t', index = False)
                kinact.random_kinact.activity_summary.to_csv(f"{odir}/RANDOM_ANALYSIS/{name_out}_random_summarized_activities.tsv", sep = '\t', index = False)

                kinact.normalized_activities.to_csv(f"{odir}/RESULTS/{name_out}_normalized_activities.tsv", sep = '\t', index = False)
                kinact.normalized_agg_activities.to_csv(f"{odir}/RESULTS/{name_out}_normalized_aggregated_activities.tsv", sep = '\t', index = False)
                kinact.normalized_summary.to_csv(f"{odir}/RESULTS/{name_out}_normalized_summarized_activities.tsv", sep = '\t', index = False)

                kinact.random_experiments.to_csv(f"{odir}/RANDOM_ANALYSIS/{name_out}_random_experiments.tsv", sep = '\t', index = False)

        pickle.dump( kinact_dict, open( f"{odir}/RESULTS/{name}_kinact.p", "wb" ) )

