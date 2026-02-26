import json
import os
import re
import pickle
import itertools
import numbers
import shutil

import pandas as pd
import numpy as np
import scipy.stats as stats
import multiprocessing
import numbers

import seaborn as sns
import matplotlib.pyplot as plt
import concurrent.futures

from tqdm import tqdm
from datetime import datetime
from collections import defaultdict
from kstar import config, helpers
from kstar.random_experiments import generate_random_experiments, calculate_fpr
from kstar.plot import DotPlot, KSTAR_PDF, plot_jaci_between_samples


class DatasetSizeError(Exception):
    """Raised when dataset size is too small for analysis"""
    pass

class KinaseActivity:
    """
    Kinase Activity calculates the estimated activity of kinases given an experiment using hypergeometric distribution.
    Hypergeometric distribution examines the number of protein sites found to be active in evidence compared to the
    number of protein sites attributed to a kinase on a provided network.

    Parameters
    ----------
    evidence : pandas df
        a dataframe that contains (at minimum, but can have more) data columms as evidence to use in analysis and KSTAR_ACCESSION and KSTAR_SITE
    odir : string
        output directory where results will be saved
    name : string
        name of the experiment, used to label output files. Default is 'experiment'
    kinases : list or None
        list of kinases to predict activity for. If None, will use all kinases found in the provided networks
    network_dir : string or None
        directory where pruned KSTAR networks are located. If None, will use config.NETWORK_DIR. If network files were downloaded with config.install_network_files(), this directory should already be set and does not need to be provided.
    network_name : string or None
        name of the network to use. If None, will use the default network name from config based on phospho_type
    data_columns: list
        list of the columns containing the abundance values, which will be used to determine which sites will be used as evidence for activity prediction in each sample
    phospho_type: string, either 'Y' or 'ST'
        indicates the phospho modification of interest
    logger : Logger object or None
        keeps track of kstar analysis, including any errors that occur. If None, a new logger will be created automatically
    min_dataset_size_for_pregenerated: int
        minimum dataset size required to use pregenerated random activities (by number of sites used as evidence). Default is 150
    max_diff_from_pregenerated: float
        maximum percent difference between dataset size and pregenerated random activity size to use pregenerated data. Default is 0.20 (i.e. 20%)
    seed : int or None
        random seed to use for random number generation. If None, seed will be set to current time

    Attributes
    ----------
    -------------------
    Upon Initialization
    -------------------
    evidence: pandas dataframe
        inputted dataset used for kinase activity calculation
    networks: dict
        dictionary of pruned kinase substrate networks, with keys as network ids and values as pandas dataframes
    data_columns: list
        list of columns containing abundance values, which will be used to determine which sites will be used as evidence. If inputted data_columns parameter was None, this lists includes in column in evidence prefixed by 'data:'
    logger : Logger object
        keeps track of kstar analysis, including any errors that occur
    aggregate: string
        the type of aggregation to use when determining binary evidence, either 'count' or 'mean'. Default is 'count'.
    run_date: string
        indicates the date that kinase activity object was initialized
    random_seed: int
        random seed used for activity calculation. Only relevant if not using pregenerated random activities
    network_info: dict
        metadata about the loaded networks
    network_hash: string
        unique identifier for the loaded networks
    kinases: list
        list of kinases to predict activity for

    ---------------------------------
    After Hypergeometric Calculations
    ---------------------------------
    real_enrichment: pandas dataframe
        p-values obtained for all pruned networks indicating statistical enrichment of a kinase's substrates for each network, based on hypergeometric
        tests
    activities: pandas dataframe
        median p-values obtained from the real_enrichment object for each experiment/kinase
    agg_activities: pandas dataframe

    -----------------------------------
    After Random Enrichment Calculation
    -----------------------------------
    random_experiments: pandas dataframe
        contains information about the sites randomly sampled for each random experiment. Will only be saved if save_random_experiments=True.
    random_enrichment: KinaseActivity object
        KinaseActivity object containing random activities predicted from each of the random experiments
    data_columns_from_scratch: list
        list of data columns which generated random activities from scratch
    data_columns_with_pregenerated: list
        list of data columns which generated random activities from pregenerated random activities

    ---------------------------
    After Mann Whitney Analysis
    ---------------------------
    activities_mann_whitney: pandas dataframe
        p-values obtained from comparing the real distribution of p-values to the distribution of p-values from random datasets, based
        the Mann Whitney U-test
    fpr_mann_whitney: pandas dataframe
        false positive rates for predicted kinase activities

    """

    def __init__(self, evidence, odir, name = 'experiment', data_columns=None, phospho_type='Y', kinases = None, network_dir = None, logger = None, network_name = None, seed = None):
        self.odir = odir
        self.name = name
        self.phospho_type = phospho_type
        self.set_evidence(evidence)
        self.networks = defaultdict()
        self.network_sizes = defaultdict()

        #if no seed is set, set seed to current time
        if seed is None:
            seed = int(datetime.now().timestamp())
        np.random.seed(seed)
        self.random_seed = seed

        #set up logger
        #if directory doesn't exist yet, create it
        if not os.path.exists(f"{odir}/RESULTS/{phospho_type}"):
            if not os.path.exists(f"{odir}/RESULTS"):
                os.mkdir(f"{odir}/RESULTS")
            os.mkdir(f"{odir}/RESULTS/{phospho_type}")

        if logger is not None:
            self.logger = logger
        else:
            self.logger = helpers.get_logger(f"mapping_{self.name}", f"{odir}/RESULTS/{phospho_type}/activity_{self.name}.log")
        # self.normalizers = defaultdict()


        
        if network_dir is not None:
            #set network directory based on user input
            self.network_name = network_name if network_name is not None else config.NETWORK_NAME[self.phospho_type]
            self.network_directory = network_dir + f"/{self.phospho_type}/{self.network_name}/"
        else:
            #use default network directory from config
            self.network_name = network_name if network_name is not None else config.NETWORK_NAME[self.phospho_type]
            self.network_directory = config.NETWORK_SUBDIR[self.phospho_type]

        #check to make sure network directory exists
        if not os.path.exists(self.network_directory):
            raise FileNotFoundError(f"Network directory not found at: {self.network_directory}. Please download networks using config.install_network_files() or indicate where networks are found with config.update_network_directory().")
        

        #load network and meta information
        self.add_networks_from_directory()
        self.network_info = helpers.parse_network_information(self.network_directory)
        self.network_hash = self.network_info['unique_network_id']
        #make sure network matches reference phosphoproteome
        if self.network_info['unique_reference_id'] != config.REFERENCE_INFO['unique_reference_id']:
            raise TypeError("Network was not built with the same reference phosphoproteome as the current configuration. Please either update the reference phosphoproteome or use a network built with the current reference phosphoproteome.")

        #filter networks to only include specified kinases, if provided, otherwise get all kinases in network
        if kinases is not None:
            self.kinases = kinases
            #filter networks to only include specified kinases
            warned_about_missing_kinases = False
            for nid, net in self.networks.items():
                self.networks[nid] = net[net['KSTAR_KINASE'].isin(kinases)]

                #check for missing kinases, warn user only once if any are missing
                missing_kinases = [k for k in kinases if k not in net['KSTAR_KINASE'].unique()]
                if len(missing_kinases) > 0 and not warned_about_missing_kinases:
                    print(f"Warning: The following kinases were not found in network: {', '.join(missing_kinases)}. You can use the get_available_kinases function to see which kinases are available in the networks.")
                    warned_about_missing_kinases = True

        else:
            #get all kinases in network
            nkeys = list(self.networks.keys())
            self.kinases = self.networks[nkeys[0]]['KSTAR_KINASE'].unique().tolist()


        self.num_networks = None

        self.real_enrichment = None
        # self.agg_activities = None
        # self.activities = None

        # Location of random experiments
        self.random_experiments = None

        # added fields for pregenerated_random
        self.compendia_distribution = None
        self.data_columns_from_scratch = None
        self.pregenerated_experiments_path = self.network_directory + "/RANDOM_ACTIVITIES/" 



        self.compendia_paths = {
            ('Y', True): 'compendia=0_30_70',
            ('Y', False): 'compendia=0_50_50',
            ('ST', None): 'compendia=0_30_70'
        }
        
        # end of new fields for pregenerated_random

        self.aggregate = 'mean'
        self.threshold = None
        self.greater = True

        self.run_date = datetime.now().strftime("%m/%d/%Y %H:%M:%S")

        self.dropped_columns = []
        # if data columns is None, set data columns to be columns with data: in front
        self.set_data_columns(data_columns=data_columns)
        #calculate the compendia distriubtion for each column
        #self.get_compendia_distribution(selection_type = "KSTAR_NUM_COMPENDIA_CLASS")
        #get available pregenerated sizes
        #self.pregenerated_sizes = self.check_file_sizes_for_pregenerated()

    def _report_warning(self, message):
        """
        Reports a warning message to both the logger and standard output

        Parameters
        ----------
        message : str
            warning message to report
        """
        print("WARNING:", message)
        self.logger.warning(message)

    def _check_single_threshold(self, threshold=None, greater = True, agg ='mean', min_evidence_size=0, allow_column_loss = True):
        if not isinstance(threshold, numbers.Number):
            return False
        
        evidence = self.evidence.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE])[self.data_columns].agg(agg)
        binary_evidence = evidence >= threshold if greater else evidence <= threshold
        num_sites = binary_evidence.sum()
        if allow_column_loss and (num_sites > min_evidence_size).any():
            return True
        elif not allow_column_loss and (num_sites > min_evidence_size).all():
            return True
        else:
            return False

    def check_valid_threshold(self, threshold=None, greater = True, agg ='mean', min_evidence_size=0, allow_column_loss = True):
        """
        Given a threshold value or a list of threshold values, make sure that at least one threshold(s) result in at least one data column having evidence (or all data columns if allow_column_loss = False)
        
        Parameters
        ----------
        threshold: float or list
            cutoff for inclusion as evidence for prediction. If greater = True, sites with quantification greater than the threshold are used as evidence.
        greater: bool
            whether to use sites greater (True) or less (False) than the threshold
        agg: str
            how to combine sites with multiple instances in experiment
        min_evidence_size: int
            minimum number of sites required for a data column to be considered for activity calculation
        allow_column_loss: bool
            whether to allow some data columns to have no evidence (True) or require all data columns to have evidence (False)
        """
        if isinstance(threshold, list):
            for t in threshold:
                is_valid = self._check_single_threshold(t, greater, agg, min_evidence_size, allow_column_loss)
                #remove invalid thresholds from list
                if not is_valid:
                    threshold.remove(t)

            if len(threshold) == 0:
                return False
            else:
                return True
        else:
            is_valid =self._check_single_threshold(threshold, greater, agg, min_evidence_size, allow_column_loss)
            return is_valid
        


    def check_data_columns(self, min_evidence_size = 0):
        """
        Checks data columns to make sure column is in evidence and that evidence filtered on that data column
        has at least one point of evidence (or minimum set by min_evidence_size). Removes all columns that do not meet criteria

        Parameters
        ----------
        min_evidence_size : int
            minimum number of sites required for a data column to be considered for activity calculation
        """
        new_data_columns = []
        evidence = self.evidence.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE])[self.data_columns].agg(
            self.aggregate).reset_index()
        for col in self.data_columns:
            if col in self.evidence.columns:
                if self.threshold is not None:
                    if isinstance(self.threshold, list):
                        #use smallest threshold if multiple provided (largest if self.greater = False)
                        threshold = min(self.threshold) if self.greater else max(self.threshold)
                    else:
                        threshold = self.threshold

                    if self.greater:
                        if len(evidence[evidence[col] >= threshold]) > min_evidence_size:
                            new_data_columns.append(col)
                        else:
                            self._report_warning(f"{col} does not have sufficient evidence, and will not be used.")
                            self.dropped_columns.append(col)
                    else:
                        if len(evidence[evidence[col] <= threshold]) > min_evidence_size:
                            new_data_columns.append(col)
                        else:
                            self._report_warning(f"{col} does not have sufficient evidence, and will not be used.")
                            self.dropped_columns.append(col)
                else:
                    if ~evidence[col].isna().all():
                        new_data_columns.append(col)
                    else:
                        self._report_warning(f"{col} does not have any evidence, and will not be used.")
                        self.dropped_columns.append(col)
            else:
                self._report_warning(f"{col} not in evidence, and will not be used.")
                self.dropped_columns.append(col)

        #make sure there is at least one data column
        if len(new_data_columns) == 0:
            raise DatasetSizeError(f"No valid data columns found after filtering evidence, as no columns have more than the requested minimum evidence size after thresholding ({min_evidence_size} sites). Please check that data columns are in evidence and that threshold is appropriate.")
        self.data_columns = new_data_columns

    def set_data_columns(self, data_columns=None):
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
            #make sure inputted data columns are in evidence (check if any have 'data:') in front, if not add it
            final_data_columns = []
            for col in data_columns:
                if col in self.evidence.columns:
                    final_data_columns.append(col)
                elif 'data:' + col in self.evidence.columns:
                    final_data_columns.append('data:' + col)
                else:
                    self._report_warning(f"{col} not found in evidence, and will not be used as a data column.")

            if len(final_data_columns) == 0:
                raise KeyError("Could not find any of the inputted data columns in the evidence dataframe. Please check that the data columns are in the evidence dataframe and try again.")

            self.data_columns = final_data_columns
            #print(self.data_columns)
        self.check_data_columns()

    def test_threshold(self, threshold, agg='mean', greater=True, plot=False, return_evidence_sizes=False, min_evidence_size = 0):
        """
        Given a threshold value, calculate the distribution of evidence sizes (i.e. number of sites used in prediction for each sample in the experiment).

        Parameters
        ----------
        threshold: float
            cutoff for inclusion as evidence for prediction. If greater = True, sites with quantification greater than the threshold are used as evidence.
        agg: str
            how to combine sites with multiple instances in experiment
        greater: bool
            whether to use sites greater (True) or less (False) than the threshold
        plot: bool
            whether to plot a histogram of the evidence sizes used and heatmap of Jaccard similarity between samples
        return_evidence_sizes: bool
            indicates whether to return the evidence sizes for all samples or not
        min_evidence_size: int
            minimum number of sites required for a data column to be considered for activity calculation

        Returns
        -------
        Outputs the minimum, maximum, and median evidence sizes across all samples. May return evidence sizes of all samples as pandas series
        """
        evidence_binary = self.create_binary_evidence(agg=agg, threshold=threshold, greater=greater, drop_empty_columns=False, min_evidence_size = min_evidence_size)
        num_sites = evidence_binary[self.data_columns].sum()
        print('Number of Sites That Will Be Used As Evidence For Each Sample:')
        print('Maximum =', num_sites.max())
        print('Minimum =', num_sites.min())
        print('Median =', num_sites.median())

        #get similarity between samples
        print('\n')
        jaci = helpers.jaci_matrix_between_samples(evidence_binary, self.data_columns)
        print('Average Jaccard Similarity Between Samples:', round(jaci[jaci != 1].mean().mean(), 3))
        print('Maximum Jaccard Similarity Between Samples:', round(jaci[jaci != 1].max().max(), 3))

        #number of columns without enough evidence
        print('\n')
        num_columns_without_enough_evidence = (num_sites <= min_evidence_size).sum()
        print(f"Number of columns without enough evidence (min = {min_evidence_size}): {num_columns_without_enough_evidence}")



        if plot:
            fig, ax = plt.subplots(figsize = (12,4), ncols = 2, width_ratios = [1,0.7])
            fig.subplots_adjust(wspace=0.8)
            plt_data = num_sites.reset_index()
            sns.histplot(data=plt_data, x=0, bins = 20, ax = ax[0])
            ax[0].set_xlabel('Number of Sites Used As Evidence')
            ax[0].set_ylabel('Frequency')
            ax[0].set_title('Distribution of Evidence\nSizes Across Samples')

            #plot jaccard similarity heatmap
            plot_jaci_between_samples(evidence_binary, self.data_columns, ax = ax[1], title = 'Jaccard Similarity\nBetween Samples', cmap = 'coolwarm', annot = False)

        if return_evidence_sizes:
            return num_sites
        
    def test_threshold_range(self, min_threshold, max_threshold, step=0.1, agg = 'mean', greater = True, min_evidence_size = 0, desired_evidence_size = None, desired_similarity = None, show_recommended = False):
        """
        Given a range of threshold values, calculate the distribution of evidence sizes (i.e. number of sites used in prediction for each sample in the experiment) and Jaccard similarity between samples at each threshold

        Parameters
        ----------
        min_threshold: float
            minimum cutoff for inclusion as evidence for prediction. If greater = True, sites with quantification greater than the threshold are used as evidence.
        max_threshold: float
            maximum cutoff for inclusion as evidence for prediction. If greater = True, sites with quantification greater than the threshold are used as evidence.
        step: float
            step size to use when iterating through threshold range
        agg: str
            how to combine sites with multiple instances in experiment
        greater: bool
            whether to use sites greater (True) or less (False) than the threshold
        min_evidence_size: int
            minimum number of sites required for a data column to be considered for activity calculation
        desired_evidence_size: int or None
            target evidence size to use for plotting. If None, will use 150 for phospho_type 'Y' and 1500 for phospho_type 'ST'
        desired_similarity: float or None
            target similarity to use for plotting. If None, will use 0.6 for phospho_type 'Y' and 0.7 for phospho_type 'ST'
        show_recommended: bool
            whether to show recommended evidence size and similarity lines on the plots
        """
        #size information
        median_sizes = []
        min_sizes = []
        max_sizes = []
        average_jac = []
        max_jac = []
        num_valid_columns = []

        #set range of thresholds
        if greater:
            range_values = np.arange(min_threshold, max_threshold + step, step)
        else:
            range_values = np.arange(max_threshold, min_threshold - step, -step)

        #iterate through threholds and calculate evidence sizes and jaccard similarity
        plt_range = []
        for threshold in range_values:
            #calculate binary evidence, checking if no evidence is found
            try:
                evidence_binary = helpers.suppress_print(self.create_binary_evidence, agg=agg, threshold=threshold, greater=greater, drop_empty_columns = False, min_evidence_size = min_evidence_size)
            except Exception as e:
                print('Error encountered when calculating evidence at threshold', threshold)
                print('Error message:', str(e))
                break


            #calculate evidence sizes
            num_sites = evidence_binary[self.data_columns].sum()
            if num_sites.min() == 0:
                direction = 'greater' if greater else 'less'
                print(f"Warning: Some data columns have no evidence at thresholds {direction} than {threshold}. Stopping further calculations.")
                break

            plt_range.append(threshold)
            median_sizes.append(num_sites.median())
            min_sizes.append(num_sites.min())
            max_sizes.append(num_sites.max())

            #get number of valid columns
            num_valid_columns.append((num_sites > min_evidence_size).sum())

            if len(self.data_columns) > 1:
                #calculate jaccard similarity between samples based on phosphosite identity
                jaci = helpers.jaci_matrix_between_samples(evidence_binary, self.data_columns)
                average_jac.append(jaci[jaci != 1].mean().mean())
                max_jac.append(jaci[jaci != 1].max().max())
        
        #set desired evidence size if not provided (for plotting)
        if desired_evidence_size is None:
            desired_evidence_size = {'Y':150, 'ST':1500}[self.phospho_type]
        if desired_similarity is None:
            desired_similarity = 0.6
        
        #setup figure
        if len(self.data_columns) <= 1:
            fig, size_ax = plt.subplots(figsize = (4,4), nrows=1, ncols=1)
        else:
            fig, axes = plt.subplots(figsize = (8,4), nrows=1, ncols=2)
            size_ax = axes[0]
            sim_ax = axes[1]
            fig.subplots_adjust(wspace=0.6)
        #make evidence size plot
        
        size_ax.plot(plt_range, median_sizes, linestyle='-', c = 'black', label = 'Median Across Samples')
        size_ax.plot(plt_range, min_sizes, linestyle = ':', c = 'black', label = 'Min/Max Across Samples')
        size_ax.plot(plt_range, max_sizes, linestyle = ':', c = 'black', label = 'Min/Max Across Samples')
        size_ax.set_xlabel('Threshold')
        size_ax.set_ylabel('Evidence Size')
        if show_recommended:
            size_ax.axhline(y=desired_evidence_size, color='red', linestyle='-', alpha = 0.5, label = 'Target Sample Size')
        size_ax.legend(loc=(0, 1.05))

        #make jaccard similarity plot
        if len(self.data_columns) > 1:
            sim_ax.plot(plt_range, average_jac, linestyle = '-', c = 'black', label = 'Average Jaccard Similarity')
            sim_ax.plot(plt_range, max_jac, linestyle = ':', c = 'black', label = 'Maximum Jaccard Similarity')
            if show_recommended:
                sim_ax.axhline(y=desired_similarity, color='red', linestyle='-', alpha = 0.5, label = 'Target Maximum Similarity')
            sim_ax.set_xlabel('Threshold')
            sim_ax.set_ylabel('Jaccard Similarity')
            sim_ax.legend(loc = (0, 1.05))

        #if any data columns were removed due to lack of evidence, print a warning
        if any([n < len(self.data_columns) for n in num_valid_columns]):
            #find first threshold where number of valid columns is less than total data columns
            for i, n in enumerate(num_valid_columns):
                if n < len(self.data_columns):
                    warning_threshold = range_values[i]
                    break
            direction = 'greater' if greater else 'less'
            print(f"Warning: Some data columns were removed due to lack of evidence at thresholds {direction} than {warning_threshold}.")


        
    def recommend_threshold_based_on_similarity(self, max_similarity = 0.7, min_threshold = 0, max_threshold = np.inf, step = 0.1, greater = True, agg = 'mean', pick_best_by = 'max'):
        """
        Given a minimum threshold (and max if provided), iteratively test thresholds to find the threshold that results in average Jaccard similarity between samples below the desired similarity
        
        Parameters
        ----------
        max_similarity: float
            maximum average Jaccard similarity allowed between samples
        min_threshold: float
            minimum cutoff for inclusion as evidence for prediction. If greater = True, sites with quantification greater than the threshold are used as evidence.
        max_threshold: float
            maximum cutoff for inclusion as evidence for prediction. If greater = True, sites with quantification
            greater than the threshold are used as evidence.
        step: float
            step size to use when iterating through threshold range
        greater: bool
            whether to use sites greater (True) or less (False) than the threshold
        agg: str
            how to combine sites with multiple instances in experiment
        pick_best_by: str
            method to use when aggregating Jaccard similarity values across samples, either 'max' or 'median'

        Returns
        -------
        threshold: float
            recommended threshold that results in average Jaccard similarity between samples below max_similarity
        ave_jaci: float
            average Jaccard similarity between samples at the recommended threshold
        """
        #make sure there are multiple data columns to compare similarity between
        if len(self.data_columns) <= 1:
            raise ValueError("Cannot recommend threshold based on similarity with only one data column, since there are no other samples to compare to")
        if greater:
            threshold = min_threshold - step
            ave_jaci = 1
            while ave_jaci >= max_similarity and threshold <= max_threshold:
                threshold += step
                evidence_binary = self.create_binary_evidence(agg=agg, threshold=threshold, greater=greater, drop_empty_columns=False)
                jaci = helpers.jaci_matrix_between_samples(evidence_binary, self.data_columns)
                #calculate average jaccard similarity excluding self-comparisons
                ave_jaci = helpers.agg_jaccard(jaci, pick_best_by)
        else:
            threshold = max_threshold + step
            ave_jaci = 1
            while ave_jaci >= max_similarity and threshold >= min_threshold:
                threshold -= step
                evidence_binary = self.create_binary_evidence(agg=agg, threshold=threshold, greater=greater, drop_empty_columns=False)
                jaci = helpers.jaci_matrix_between_samples(evidence_binary, self.data_columns)
                ave_jaci = helpers.agg_jaccard(jaci, pick_best_by)

        return threshold, ave_jaci

    def recommend_threshold_by_size(self, desired_evidence_size = None, min_threshold = 1, max_threshold = np.inf, step = 0.1, pick_best_by = 'median', greater = True, agg = 'mean'):
        """
        Given a minimum threshold (and max if provided), iteratively test thresholds to find the threshold that results in desired evidence size (by default based on median)

        Parameters
        ----------
        desired_evidence_size: int
            target evidence size to use when recommending threshold
        min_threshold: float
            minimum cutoff for inclusion as evidence for prediction. If greater = True, sites with quantification greater than the threshold are used as evidence.
        max_threshold: float
            maximum cutoff for inclusion as evidence for prediction. If greater = True, sites with quantification greater than the threshold are used as evidence.
        step: float
            step size to use when iterating through threshold range
        pick_best_by: str
            method to use when aggregating evidence size values across samples, either 'max' or 'median'
        greater: bool
            whether to use sites greater (True) or less (False) than the threshold
        agg: str
            how to combine sites with multiple instances in experiment
        
        Returns
        -------
        threshold: float
            recommended threshold that results in desired evidence size
        evidence_size: int
            evidence size at the recommended threshold
        """
        #if not indicated, use default desired evidence size
        if desired_evidence_size is None:
            desired_evidence_size = 1500 if self.phospho_type == 'ST' else 150 

        # find evidence with minimum threshold
        if greater:
            threshold = min_threshold - step
            evidence_size = np.inf
            # iterate through thresholds until evidence size is less than desired evidence size or threshold is greater than max_threshold
            while evidence_size > desired_evidence_size and (threshold + step) <= max_threshold:
                threshold += step
                evidence_binary = self.create_binary_evidence(agg=agg, threshold=threshold, greater=greater, drop_empty_columns=False)
                num_sites = evidence_binary[self.data_columns].sum()
                evidence_size = num_sites.agg(pick_best_by)
        else:
            threshold = max_threshold + step
            evidence_size = np.inf
            # iterate through thresholds until evidence size is less than desired evidence size or threshold is less than min_threshold
            while evidence_size > desired_evidence_size and (threshold - step) >= min_threshold:
                threshold -= step
                evidence_binary = self.create_binary_evidence(agg=agg, threshold=threshold, greater=greater, drop_empty_columns=False)
                num_sites = evidence_binary[self.data_columns].sum()
                # aggregate counts
                evidence_size = num_sites.agg(pick_best_by)

        return threshold, evidence_size
    
    def get_allowable_threshold(self, greater = True, agg='mean', min_evidence_size = 20, allow_column_loss = False):
        """
        Determine the minimum/maximum threshold that still results in all data columns having evidence

        Parameters
        ----------
        greater: bool
            whether to use sites greater (True) or less (False) than the threshold
        agg: str
            how to combine sites with multiple instances in experiment
        min_evidence_size: int
            minimum number of sites required for a data column to be considered for activity calculation

        Returns
        -------
        allowable threshold: float
            maximum or minimum threshold that still results in all data columns having evidence (or at least one if min_evidence_size = None)
        """
        min_evidence_size = min_evidence_size if min_evidence_size is not None else 0
        #aggregate evidence by site (same as is done in create_binary_evidence)
        agg_evidence = self.evidence.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE])[self.data_columns].agg(agg).reset_index()
        if greater:
            if min_evidence_size == 0 and allow_column_loss:
                #set allowable threshold to maximum value in data columns (value that would cause no evidence to be found in any column)
                allowable_threshold = agg_evidence[self.data_columns].max().max()
            elif min_evidence_size == 0 and not allow_column_loss:
                #if no minimum evidence size, set allowable threshold to largest value that still results in at least one evidence point per data column
                allowable_threshold = agg_evidence[self.data_columns].max().min()
            elif min_evidence_size > 0:
                #identify the nth largest value in each data column, where n = min_evidence_size
                nth_largest_values = []
                dropped_columns = []
                for col in self.data_columns:
                    sorted_values = agg_evidence[col].dropna().sort_values(ascending=False).reset_index(drop=True)
                    if len(sorted_values) >= min_evidence_size:
                        nth_largest = sorted_values.iloc[min_evidence_size - 1]
                        nth_largest_values.append(nth_largest)
                    else:
                        dropped_columns.append(col)

                #make sure there is at least one column with enough evidence
                if len(nth_largest_values) == 0:
                    raise DatasetSizeError(f"No columns met the minimum evidence size requirement of {min_evidence_size} sites.")
                if not allow_column_loss and len(dropped_columns) > 0:
                    raise DatasetSizeError(f"The following data columns do not have enough evidence to meet the minimum evidence size requirement of {min_evidence_size}: {', '.join(dropped_columns)}. Please lower the min_evidence_size parameter or set allow_column_loss = True to allow for some data columns to be removed from analysis.")
            
                elif len(dropped_columns) > 0:
                    self._report_warning(f"Warning: The following data columns do not have enough evidence to meet the minimum evidence size requirement of {min_evidence_size} and will be ignored when calculating allowable threshold: {', '.join(dropped_columns)}.")

                #set smallest of these values as the allowable threshold
                allowable_threshold = min(nth_largest_values)
            elif min_evidence_size is None:
                #set allowable threshold to maximum value in data columns (value that would cause no evidence to be found in any column)
                allowable_threshold = agg_evidence[self.data_columns].max().max()
        else:
            if min_evidence_size == 0:
                allowable_threshold = agg_evidence[self.data_columns].min().max()
            elif min_evidence_size > 0:
                #identify the nth smallest value in each data column, where n = min_evidence_size
                nth_smallest_values = []
                dropped_columns = []
                for col in self.data_columns:
                    sorted_values = agg_evidence[col].dropna().sort_values(ascending=True).reset_index(drop=True)
                    if len(sorted_values) >= min_evidence_size:
                        nth_smallest = sorted_values.iloc[min_evidence_size - 1]
                        nth_smallest_values.append(nth_smallest)
                    else:
                        dropped_columns.append(col)


                #make sure there is at least one column with enough evidence
                if len(nth_smallest_values) == 0:
                    raise DatasetSizeError(f"No columns met the minimum evidence size requirement of {min_evidence_size} sites.")
                elif not allow_column_loss and len(dropped_columns) > 0:
                    raise DatasetSizeError(f"The following data columns do not have enough evidence to meet the minimum evidence size requirement of {min_evidence_size} sites: {', '.join(dropped_columns)}. Please lower the min_evidence_size parameter or set allow_column_loss = True to allow for some data columns to be removed from analysis.")
                elif len(dropped_columns) > 0:
                    self._report_warning(f"Warning: The following data columns do not have enough evidence to meet the minimum evidence size requirement of {min_evidence_size} sites and will be ignored when calculating allowable threshold: {', '.join(dropped_columns)}.")
                allowable_threshold = max(nth_smallest_values)
            elif min_evidence_size is None:
                #set allowable threshold to minimum value in data columns
                allowable_threshold = self.evidence[self.data_columns].min().min()

        return allowable_threshold
    
    def recommend_threshold(self, desired_evidence_size = None, max_similarity = 0.7, consider_size = True, consider_similarity = True, min_threshold = -np.inf, max_threshold = np.inf, step = 0.1, pick_best_size_by = 'median', pick_best_similarity_by = 'max', greater = True, agg = 'mean', min_evidence_size = 20, allow_column_loss = False):
        """
        Recommend a threshold, one based on desired evidence size and one based on maximum average Jaccard similarity between samples. Will report the characteristics of the resulting evidences for both thresholds

        Parameters
        ----------
        desired_evidence_size: int
            target evidence size to use when recommending threshold
        max_similarity: float
            maximum average Jaccard similarity between samples to use when recommending threshold. Default is 0.7
        consider_size: bool
            whether to consider evidence size when recommending threshold
        consider_similarity: bool
            whether to consider similarity between data columns when recommending threshold
        min_threshold: float
            minimum threshold to consider when recommending threshold. Must be provided if greater = True. Default is -infinity
        max_threshold: float
            maximum threshold to consider when recommending threshold. Must be provided if greater = False. Default is infinity
        step: float
            step size to use when iterating through thresholds
        pick_best_size_by: str
            method to use when aggregating evidence size values across samples, recommended to be either 'min', 'max', or 'median'
        pick_best_similarity_by: str
            method to use when aggregating Jaccard similarity values across samples, recommended to be either 'max' or 'median'
        greater: bool
            whether to use sites greater (True) or less (False) than the threshold
        agg: str
            how to combine sites with multiple instances in experiment
        min_evidence_size: int
            minimum number of sites required for a data column to be considered for activity calculation
        allow_column_loss: bool
            whether to allow some data columns to be lost when recommending threshold based on size. If False, will raise an error if min/max thresholds provided result in loss of any data columns

        Returns
        -------
        float
            recommended threshold value
        """
        #make sure min and max thresholds are set appropriately
        if greater and min_threshold == -np.inf:
            raise ValueError("When greater = True, min_threshold must be set to a finite value.")
        if not greater and max_threshold == np.inf:
            raise ValueError("When greater = False, max_threshold must be set to a finite value.")
        if not consider_size and not consider_similarity:
            raise ValueError("At least one of consider_size or consider_similarity must be set to True.")


        #if not allowing column loss, check that provided min/max thresholds will not result in any data columns being lost
        if not allow_column_loss:
            #determine the minimum/maximum threshold that still results in all data columns having evidence, based on min_evidence_size
            allowable_threshold = self.get_allowable_threshold(greater=greater, min_evidence_size=min_evidence_size, allow_column_loss=allow_column_loss)
            #raise errors or warnings if provided min/max thresholds are outside of allowable range
            if greater:
                if allowable_threshold < min_threshold:
                    raise ValueError(f'Min threshold of {min_threshold} is too high and results in no evidence for at least on data column. Your options include:\n1) Setting min_threshold to lower than {allowable_threshold}\n2) Lowering the `min_evidence_size` parameter in recommend_threshold to allow for smaller dataset sizes\n3) Set alow_column_loss = True if you do not care if some data columns are lost).')
                elif max_threshold == np.inf:
                    max_threshold = round(allowable_threshold, 2)
                elif allowable_threshold < max_threshold:
                    print(f'Warning: Max threshold of {max_threshold} is too high, resulting in no evidence for at least one data column. Maximum allowable threshold to retain evidence for all data columns is {allowable_threshold}. Setting max_threshold to {allowable_threshold}. You can lower the `min_evidence_size` parameter in recommend_threshold to allow higher thresholds (set to negative if you do not care if some data columns are lost).\n')
                    max_threshold = round(allowable_threshold,2)

            else:
                if allowable_threshold > max_threshold:
                    raise ValueError(f'Max threshold of {max_threshold} is too high and results in no evidence for at least on data column. Your options include:\n1) Setting max_threshold to higher than {allowable_threshold}\n2) Lowering the `min_evidence_size` parameter in recommend_threshold to allow for smaller dataset sizes\n3) Set allow_column_loss = True if you do not care if some data columns are lost).')
                elif min_threshold == -np.inf:
                    min_threshold = round(allowable_threshold, 2)
                elif allowable_threshold > min_threshold:
                    print(f'Warning: Min threshold of {min_threshold} is too low, resulting in no evidence for at least one data column. Minimum allowable threshold to retain evidence for all data columns is {allowable_threshold}. Setting min_threshold to {allowable_threshold}. You can lower the `min_evidence_size` parameter in recommend_threshold to allow lower thresholds (set allow_column_loss=True if you do not care if some data columns are lost).\n')
                    min_threshold = round(allowable_threshold,2)


        #identify and test optimal threshold based on evidence size
        if consider_size:
            threshold_by_size, evidence_size = self.recommend_threshold_by_size(desired_evidence_size=desired_evidence_size, pick_best_by=pick_best_size_by, min_threshold=min_threshold, max_threshold=max_threshold, step=step, greater=greater, agg=agg)
            print(f'Recommended threshold based on {pick_best_size_by} evidence_size:', round(threshold_by_size, 2), f'({pick_best_size_by} sites = {evidence_size})')


        #if also considering similarity between data columns, check here
        if consider_similarity:
            if len(self.data_columns) <= 1:
                print('Warning: Cannot recommend threshold based on similarity with only one data column, since there are no other samples to compare to. Skipping similarity-based threshold recommendation.')
                consider_similarity = False
            else:
                #identify and test optimal threshold based on similarity
                threshold_by_similarity, similarity = self.recommend_threshold_based_on_similarity(max_similarity=max_similarity, min_threshold=min_threshold, max_threshold=max_threshold, step=step, greater=greater, agg=agg, pick_best_by=pick_best_similarity_by)
                print(f'Recommended threshold based on {pick_best_similarity_by} similarity:', round(threshold_by_similarity, 2), f'({pick_best_similarity_by} Jaccard similarity = {round(similarity, 2)})')  

        #choose the more stringent of the two thresholds
        if consider_size and consider_similarity:
            if greater:
                recommended_threshold = round(max(threshold_by_size, threshold_by_similarity), 2)
            else:
                recommended_threshold = round(min(threshold_by_size, threshold_by_similarity), 2)

            #check if recommended thresholds are 
        elif consider_size:
            recommended_threshold = threshold_by_size
        elif consider_similarity:
            recommended_threshold = threshold_by_similarity

        print(f"Final recommended threshold = {recommended_threshold})")


        return recommended_threshold



    def network_check_for_pregeneration(self, pregenerated_path):
        """
        Check if the network hash matches a pre-generated network in a pregen_experiments directory.
        and verifies RUN_INFORMATION.txt within the hash subdirectory.

        Returns:
            bool: True if the data matches, False otherwise.
        """
        # Check if this is a directory and matches the target hash
        if os.path.isdir(pregenerated_path):
            # Path to RUN_INFORMATION.txt in the hash subdirectory
            default_file_path = os.path.join(pregenerated_path, "RUN_INFORMATION.txt")

            if not os.path.exists(default_file_path):
                return False

            # Parse and compare the RUN_INFORMATION.txt files
            default_values = helpers.parse_network_information(pregenerated_path)
            if default_values['unique_network_id'] != self.network_hash:
                return False

            return True
        else:
            # If no matching hash directory is found
            return False


    def get_compendia_distribution(self, data_columns=None,
                                   selection_type='KSTAR_NUM_COMPENDIA_CLASS'):
        """
        Get the compendia distribution for each data column.

        Parameters
        ----------
        with_pregenerated_evidence : pandas DataFrame
            KSTAR mapped experimental dataframe that has been binarized by kstar_activity generation.
        data_columns : list
            Columns that represent experimental results.
        selection_type : str, optional
            The type of compendia selection, by default 'KSTAR_NUM_COMPENDIA_CLASS'.

        Returns
        -------
        dict
            Dictionary containing the compendia distribution for each data column.
        """
        if data_columns is None:
            data_columns = self.data_columns

        compendia_distribution = {}
        for col in data_columns:
            filtered_data = self.evidence_binary[self.evidence_binary[col] == 1]
            compendia_counts = filtered_data[selection_type].value_counts(normalize=True).mul(100).to_dict()
            compendia_distribution[col] = {k: round(compendia_counts.get(k, 0), 2) for k in range(3)}
        return compendia_distribution
    


    def check_file_sizes_for_pregenerated(self):
        """
        Check the sizes of pre-generated files for the given datasets.

        This function identifies the appropriate pre-generated file sizes based on the dataset size and returns a list of these sizes.

        Returns
        -------
        pregenerated_sizes_files : list
            List of sizes of pre-generated files that match the dataset size criteria.
        """
        pregenerated_sizes_files = {}
        #grab compendia distributions
        default_compendia_distributions = os.listdir(self.pregenerated_experiments_path)
        for dist in default_compendia_distributions:
            #experiment sizes will be the directories within each compendia distribution
            file_sizes = os.listdir(os.path.join(self.pregenerated_experiments_path, dist))
            #convert to int
            file_sizes = [int(size) for size in file_sizes]
            pregenerated_sizes_files[dist.split('=')[1]] = file_sizes

        #if custom pregenerated path is provided, check that as well
        if not self.default_pregen_only:
            custom_compendia_distributions = os.listdir(self.custom_pregenerated_path)
            for dist in custom_compendia_distributions:
                tmp_path = os.path.join(self.custom_pregenerated_path, dist)
                #make sure it's a directory
                if os.path.isdir(tmp_path):
                    file_sizes = os.listdir(tmp_path)
                    #convert to int
                    file_sizes = [int(size) for size in file_sizes]
                    pregenerated_sizes_files[dist.split('=')[1]] = file_sizes
                    #experiment sizes will be the directories within each compendia distribution
                    if dist.split('=')[1] not in pregenerated_sizes_files:
                        pregenerated_sizes_files[dist.split('=')[1]] = file_sizes
                    else:
                        pregenerated_sizes_files[dist.split('=')[1]].extend(file_sizes)

        return pregenerated_sizes_files
    
    def determine_pregen_compendia_dist(self, high_study_bias_perc):
        """
        Based on the phospho_type and compendia class counts, determine the appropriate pregenerated compendia distribution to use.
        
        Parameters
        ----------
        high_study_bias_perc : int
            percentage of sites in experiment that have high study bias (KSTAR_COMPENDIA_CLASS = 2)
        """
        if self.phospho_type == 'Y':
            pregen_compendia_dist = '0_30_70' if high_study_bias_perc > 60 else '0_50_50'
        elif self.phospho_type == 'ST':
            pregen_compendia_dist = '0_30_70' 
        else:
            raise ValueError(f"ERROR: Unrecognized phosphoType '{self.phospho_type}'. Must be 'Y' or 'ST'.")
        return pregen_compendia_dist

    def determine_if_pregen(self, dataset, max_diff_from_pregenerated=0.25, min_dataset_size_for_pregenerated=150):
        """
        Given a dataset size and compendia distribution, determine if pregenerated data should be used, and output the location of the file to use.

        Parameters
        ----------
        size : int
            size of actual experiment/sample (i.e number of unique phosphosites)
        high_study_bias_perc : int
            percentage of sites in experiment that have high study bias (KSTAR_COMPENDIA_CLASS = 2)
        pregenerated_sizes : dict
            contains the sizes of the available random experiments for different compendia distributions (generated by `check_file_sizes_for_pregenerated()` function)
        custom_pregenerated_path : str
            path to custom pregenerated experiments. Only needed if pregenerated_sizes is not provided
        """
        #get real experiment characteristics
        high_study_bias_perc = int(self.compendia_distribution[dataset][2])
        size = self.dataset_sizes[dataset]

        pregen_compendia_dist = self.determine_pregen_compendia_dist(high_study_bias_perc)


        #get experiment sizes associated with compendia distribution
        if not hasattr(self, 'pregenerated_sizes'):
            self.pregenerated_sizes = self.check_file_sizes_for_pregenerated()

        available_sizes = self.pregenerated_sizes[pregen_compendia_dist]
        #make sure dataset is large enough to use pregenerated data
        large_enough = size >= min_dataset_size_for_pregenerated
        #make sure pregenerated data is comparable in size to dataset
        close_enough = any(
            abs(pregen_size - size) / size < max_diff_from_pregenerated for pregen_size in
            available_sizes)
        # make sure pregenerated files match network
        use_pregen = large_enough and close_enough
        self.logger.info(f"Use pregen for {dataset}: {use_pregen} (Dataset Size: {size}, Large Enough: {large_enough}, Close Enough: {close_enough})")
        if use_pregen:
            closest_size = min(available_sizes, key=lambda x: abs(x - size))
            matched_experiment_dir = self.get_pregenerated_experiment_dir(closest_size, high_study_bias_perc)
        else:
            matched_experiment_dir = None
        return use_pregen, matched_experiment_dir
    
    def get_pregenerated_subdirectory(self, size, high_study_bias_perc):
        """
        Given a column of interest, determine the specific pregenerated path to use based on the compendia distribution and dataset size. This will output the subdirectory to append to either the default or custom pregenerated path.
        
        Parameters
        ----------
        col : str
            Column of interest in the dataset.
        """
        # Determine the compendia directory based on phospho_type and class counts
        #get real experiment characteristics
        pregen_compendia_dist = self.determine_pregen_compendia_dist(high_study_bias_perc)

        # Save activities to file
        specific_path = os.path.join(f'compendia={pregen_compendia_dist}', str(size))
        return specific_path

    def get_pregenerated_experiment_dir(self, size, high_study_bias_perc):
        """
        Given a the compendia distribution and size, get the path to the pregenerated experiment directory. Will first check default directory, then custom directory if provided.
        """
        pregenerated_subdirectory = self.get_pregenerated_subdirectory(size, high_study_bias_perc)
        default_exp_file_path = os.path.join(
            str(self.pregenerated_experiments_path), pregenerated_subdirectory
        )
        custom_exp_file_path = os.path.join(
            str(self.custom_pregenerated_path), pregenerated_subdirectory
        )

        if os.path.exists(default_exp_file_path):
            return default_exp_file_path
        elif os.path.exists(custom_exp_file_path):
            return custom_exp_file_path
        else:
            raise FileNotFoundError(f"Could not find pregenerated experiment directory located at {pregenerated_subdirectory} in either default or custom pregenerated paths.")
        
    def check_pregenerated(self, default_pregen_only=False, custom_pregen_dir=None, save_new_random_activities = False):
        """
        Make sure that either default pregenerated data or custom pregenerated data exists and matches the network hash.
        
        Parameters
        ----------
        default_pregen_only : bool
            Whether to only use the default pregenerated data, by default False.
        custom_pregen_dir : str, optional
            Path to the custom pregenerated data directory, by default None.
        """
        #default path should be located in network directory, check to make sure it and any provided custom directory exists
        custom_pregen_dir = custom_pregen_dir if custom_pregen_dir is not None else config.CUSTOM_RANDOM_ACTIVITIES_DIR

        self.default_pregen_only = default_pregen_only
        if not default_pregen_only and custom_pregen_dir is None:
            #self._report_warning("Warning: custom_pregen_dir is not provided, so will only using default pregenerated data.")
            self.default_pregen_only = True
            self.custom_pregenerated_path = None
            if not os.path.exists(self.pregenerated_experiments_path):
                raise ValueError(f'Could not find pregenerated random activities in expected network directory ({self.pregenerated_experiments_path}). Please pregenerate activities with the pregenerate module or set use_pregen = False.')
            else:
                self.network_check = True

        elif default_pregen_only and custom_pregen_dir is not None:
            self._report_warning("Warning: custom_pregen_dir is provided but default_pregen_only is set to True. Ignoring custom pregenerated directory and only using default pregenerated data.")
            self.custom_pregenerated_path = None
            #make sure default pregenerated path exists
            if not os.path.exists(self.pregenerated_experiments_path):
                raise ValueError(f'Could not find pregenerated random activities in expected network directory ({self.pregenerated_experiments_path}). Please pregenerate activities with the pregenerate module or set use_pregen = False.')
            else:
                self.network_check = True
        elif not default_pregen_only and custom_pregen_dir is not None:
            #set and check custom pregenerated path
            if not os.path.exists(custom_pregen_dir):
                raise ValueError(f"Provided custom pregenerated directory '{custom_pregen_dir}' does not exist. Please provide a valid directory or set use_default_pregen_only to True.")
            
            custom_pregenerated_path = os.path.join(custom_pregen_dir, self.phospho_type, self.network_name)

            #check if custom pregenerated directory exists
            if not os.path.exists(custom_pregenerated_path):
                if not save_new_random_activities:
                    self._report_warning(f'Found provided pregenerated directory ({custom_pregen_dir}), but it does not contain a subdirectory for this network. Will only use default pregenerated data. To save new pregenerated data to this directory, set save_new_random_activities = True.')
                else:
                    self._report_warning(f'Provided pregenerated directory ({custom_pregen_dir}) does not contain a subdirectory for this network. Creating new pregenerated directory at {custom_pregenerated_path} to save new random activities.')
                    #add run information to custom pregenerated directory
                    os.makedirs(custom_pregenerated_path)
                    #copy run information file from network directory
                    run_info_file = os.path.join(self.network_directory, "RUN_INFORMATION.txt")
                    shutil.copyfile(run_info_file, os.path.join(custom_pregenerated_path, "RUN_INFORMATION.txt"))

                self.custom_pregenerated_path = custom_pregenerated_path
                self.network_check = True
            else:
                #check to make sure network hash matches
                match = self.network_check_for_pregeneration(custom_pregenerated_path)
                if not match:
                    print("Warning: Provided custom pregenerated directory either could not be find or does not match the network hash (i.e. was created using a different network). Please either set use_default_pregen_only to True or fix directory")
                    self.custom_pregenerated_path = None
                    self.default_pregen_only = True
                if (not match and not os.path.exists(self.pregenerated_experiments_path)):
                    raise ValueError('Could not find pregenerated random activities (either default or in the provided custom directory. Please pregenerate activities with the pregenerate module or set use_pregen = False.)')
                else:
                    self.network_check = True
                    self.custom_pregenerated_path = custom_pregenerated_path
        else:
            self.custom_pregenerated_path = None
            if  not os.path.exists(self.pregenerated_experiments_path):
                raise ValueError(f'Could not find pregenerated random activities in expected network directory ({self.pregenerated_experiments_path}). Please pregenerate activities with the pregenerate module or set use_pregen = False.')
            else:
                self.network_check = True
        
    def setup_pregenerated(self, max_diff_from_pregenerated=0.25, min_dataset_size_for_pregenerated=150):
        """
        For each dataset, determine if pregenerated data can be used based on dataset size and compendia distribution, and identify the specific pregenerated experiment directory to use for each experiment

        Parameters
        ----------
        max_diff_from_pregenerated : float
            Maximum allowed difference in size between the dataset and pregenerated data to use pregenerated data, by default 0.25.
        min_dataset_size_for_pregenerated : int
            Minimum dataset size required to use pregenerated data, by default 150.
        """
        # Classify datasets
        self.check_file_sizes_for_pregenerated()
        #network_check = self.network_check_for_pregeneration()
        if not self.network_check:
            self.logger.warning("Network used does not match any pre-generated networks. All datasets will be calculated from scratch.")
            print("Network used does not match any pre-generated networks. All datasets will be calculated from scratch.")
            self.data_columns_with_pregenerated = []
            self.data_columns_from_scratch = self.data_columns
            self.selected_pregenerated_paths = {}
        else:
            self.data_columns_with_pregenerated = []
            self.data_columns_from_scratch = []
            self.selected_pregenerated_paths = {}

            #make sure parameters are valid and store if so
            if isinstance(max_diff_from_pregenerated, float) and 0 <= max_diff_from_pregenerated <= 1:
                self.max_diff_from_pregenerated = max_diff_from_pregenerated
            else:
                raise ValueError("max_diff_from_pregenerated must be a float between 0 and 1.")
            
            if isinstance(min_dataset_size_for_pregenerated, int) and min_dataset_size_for_pregenerated >= 0:
                self.min_dataset_size_for_pregenerated = min_dataset_size_for_pregenerated
            else:
                raise ValueError("min_dataset_size_for_pregenerated must be a non-negative integer.")

            for dataset in self.data_columns:
                #determine if pregenerated data can be used for this dataset, and get matched experiment directory if so
                use_pregen, matched_experiment_dir  = self.determine_if_pregen(dataset, max_diff_from_pregenerated=max_diff_from_pregenerated, min_dataset_size_for_pregenerated=min_dataset_size_for_pregenerated)
                self.selected_pregenerated_paths[dataset] = matched_experiment_dir
                if use_pregen:
                    self.data_columns_with_pregenerated.append(dataset)
                else:
                    self.data_columns_from_scratch.append(dataset)


            

    def get_random_activities(self, num_random_experiments=150, use_pregenerated_random_activities=None, default_pregen_only = False, save_new_random_activities=None,  custom_pregenerated_path=None, save_random_experiments=None, require_pregenerated = False, max_diff_from_pregenerated=0.25, min_dataset_size_for_pregenerated=150, show_taskbar = True, PROCESSES=1):
        """
        Generate random experiments and calculate kinase activities.Either uses pre-generated activity lists or
        generates new random experiments based on the provided parameters.

        Parameters
        ----------
        num_random_experiments : int, optional
            Number of random experiments to generate, by default 150.
        use_pregenerated_random_activities : bool, optional
            Whether to use pre-generated data, by default None and will use configuration value.
        default_pregen_only : bool, optional
            Whether to only use the default pregenerated data found in the network directory folder, by default False.
        save_new_random_activities : bool, optional
            Whether to save new pregenerated data, by default None and will use configuration value
        custom_pregenerated_path : str, optional
            Directory to save new precomputed data, by default None and will use configuration value.
        save_random_experiments : bool, optional
            Whether to save the generated random experiments, by default None and will use configuration value.
        require_pregenerated : bool, optional
            Whether to require using pre-generated data for all datasets, by default False. This is will ensure fast run times, but may result in some datasets not being processed if they do not have matching pre-generated data (most commonly due to smaller samples).
        max_diff_from_pregenerated : float, optional
            Maximum allowed difference in size between the dataset and pregenerated data to use pregenerated data, by default 0.25.
        min_dataset_size_for_pregenerated : int, optional
            Minimum dataset size required to use pregenerated data, by default 150.
        PROCESSES : int, optional
            Number of processes to use for parallel computation, by default 1.

        
        New Attributes
        --------------
        self.random_enrichment : pandas.DataFrame
            DataFrame containing the kinase activities for each random experiment.
        self.random_experiments : pandas.DataFrame
            DataFrame containing the generated random experiments. Only saved if save_random_experiments is set to True and not using pregenerated data.

        """
        self.logger.info("Running Randomization Pipeline")
        self.use_pregen_data = use_pregenerated_random_activities if use_pregenerated_random_activities is not None else config.USE_PREGENERATED_RANDOM_ACTIVITIES
        #if use_pregen_data is True, force num_random_experiments to 150
        if self.use_pregen_data and num_random_experiments != 150:
            print("Note: When using pregenerated random activities, num_random_experiments must be set to 150 to match pregenerated data.")
            num_random_experiments = 150
        
        custom_pregenerated_path = custom_pregenerated_path if custom_pregenerated_path is not None else config.CUSTOM_RANDOM_ACTIVITIES_DIR
        if save_new_random_activities and custom_pregenerated_path is None:
            raise ValueError("If save_new_random_activities is set to True, must provide a `custom_pregenerated_path` to save new pregenerated data.")
        

        self.require_pregenerated = require_pregenerated
        if self.require_pregenerated and not self.use_pregen_data:
            self._report_warning("require_pregenerated is set to True but use_pregen_data is set to False. Setting use_pregen_data to True.")
            self.use_pregen_data = True

        if self.use_pregen_data:
            #make sure pregenerated data exists and matches network
            self.check_pregenerated(default_pregen_only=default_pregen_only, custom_pregen_dir=custom_pregenerated_path)

            #initialize other parameters
            self.save_new_pregenerated = save_new_random_activities if save_new_random_activities is not None else config.SAVE_NEW_RANDOM_ACTIVITIES
            


            # Prepare for pregenerated data by determining which datasets can use pregenerated data
            self.setup_pregenerated(max_diff_from_pregenerated=max_diff_from_pregenerated, min_dataset_size_for_pregenerated=min_dataset_size_for_pregenerated)


            # Process pre-generated data, one dataset at a time
            print('Loading pre-generated random activities for datasets where applicable...')
            #report on which datasets will use pregenerated data
            if len(self.data_columns_from_scratch) > 0 and not self.require_pregenerated:
                print(f"{len(self.data_columns_from_scratch)} out of {len(self.data_columns)} columns did not have appropriate pregenerated random enrichment and will calculate from scratch: {', '.join(self.data_columns_from_scratch)}")

                self.logger.info(f"Generating random experiments for: {', '.join(self.data_columns_from_scratch)}")
                self.calculate_random_enrichment(num_random_experiments, selection_type='KSTAR_NUM_COMPENDIA_CLASS',
                    save_random_experiments = save_random_experiments, save_new_random_activities = save_new_random_activities, show_taskbar=show_taskbar, PROCESSES=PROCESSES
                )
            elif len(self.data_columns_with_pregenerated) == 0 and self.require_pregenerated:
                raise ValueError("No datasets had appropriate pregenerated random enrichment to use. Please set require_pregenerated to False to calculate random activities from scratch for datasets without matched pregenerated data.")
            elif len(self.data_columns_from_scratch) > 0 and self.require_pregenerated:
                print(f"The following columns did not have a matched pregenerated dataset and will not be processed {', '.join(self.data_columns_from_scratch)}. \nTo calculate random activities for these datasets automatically, please set require_pregenerated to False.")
                self.logger.info(f"Skipping activity calculation for the following columns that did not have matched pregenerated dataset: {', '.join(self.data_columns_from_scratch)}")
                #remove columns without pregenerated data from data_columns
                self.data_columns = self.data_columns_with_pregenerated
                self.data_columns_from_scratch = []

            #load pregenerated data for datasets that can use it
            if len(self.data_columns_with_pregenerated) > 0:
                # Process pre-generated data, one dataset at a time
                self.load_pregenerated_random_enrichment()

            # Add pre-generated data to random enrichment
            self.add_pregenerated_to_random_enrichment()

        # Process all from-scratch data if use_pregen_data is False
        else:
            self.data_columns_from_scratch = self.data_columns
            self.data_columns_with_pregenerated = []
            self.logger.info(f"Generating random experiments for: {', '.join(self.data_columns)}")
            self.calculate_random_enrichment(num_random_experiments, selection_type='KSTAR_NUM_COMPENDIA_CLASS',
                save_random_experiments = save_random_experiments, show_taskbar=show_taskbar, PROCESSES=PROCESSES
            )


        # Add pre-generated data to random enrichment
        #self.add_pregenerated_to_random_enrichment()


    def calculate_random_enrichment(self, num_random_experiments, selection_type='KSTAR_NUM_COMPENDIA_CLASS', save_random_experiments=False, save_new_random_activities=False, show_taskbar = True, PROCESSES=1):
        """
        Calculate the kinase activities for each random experiment and discard the resulting random experiment unless save_random_experiments is set to True.

        Parameters
        ----------
        num_random_experiments : int
            Number of random experiments to generate.
        selection_type : str, optional
            The type of compendia selection, by default 'KSTAR_NUM_COMPENDIA_CLASS'.
        save_random_experiments : bool, optional
            Whether to save the generated random experiments, by default False.
        show_taskbar : bool, optional
            Whether to show a progress bar, by default True.
        PROCESSES : int, optional
            Number of processes to use for parallel computation, by default 1.

        Updated Attributes
        ------------------
        self.random_enrichment : pandas.DataFrame
            DataFrame containing the kinase activities for each random experiment.
        """
        if self.use_pregen_data:
            num_random_experiments = 150
        else:
            num_random_experiments = num_random_experiments
        """
        Change to calculate_random_enrichment to immediately calculate the kinase activities for each random experiment and discard the resulting random experiment.
        """
        #if save_new_random_activities is True, make sure that pregenerated paths are setup
        if save_new_random_activities and not hasattr(self, 'custom_pregenerated_path'):
            raise ValueError("If save_new_random_activities is set to True, custom_pregenerated_path must be known. Run self.check_pregenerated() first to setup pregenerated paths.")
        

        self.logger.info("Running Randomization Pipeline")
        self.num_random_experiments = num_random_experiments
        #grab sites filtered by compendia bias
        filtered_compendia = getFilteredCompendia(phospho_type=self.phospho_type, selection_type=selection_type)

        if PROCESSES > 1:
            pool = multiprocessing.Pool(processes=PROCESSES)
            # prepare iterable for all 150 experiments to be generated
            filtered_compendia = itertools.repeat(filtered_compendia)
            networks = itertools.repeat(self.networks)
            network_sizes = itertools.repeat(self.network_sizes)
            rand_exp_numbers = list(range(num_random_experiments))
            save_experiments = itertools.repeat(save_random_experiments)
            rand_experiments = []

            combined_activities_list = []
            for col in tqdm(self.data_columns_from_scratch, desc="Calculating activities from random experiments for each column not using pregenerated random activities"):
                activities_list = []
                compendia_sizes = self.evidence_binary[self.evidence_binary[col] == 1].groupby(
                            selection_type).size()
                compendia_sizes = itertools.repeat(compendia_sizes)
                col_name = itertools.repeat(col)
                iterable = zip(compendia_sizes, filtered_compendia, networks, network_sizes, col_name,
                                rand_exp_numbers, save_experiments)
                result = pool.starmap(calculate_random_activity_singleExperiment, iterable)
                # extract results
                if save_random_experiments:
                    activities, experiments = zip(*result)
                    rand_experiments.extend(experiments)
                    activities_list.extend(activities)
                else:
                    activities_list = activities_list + result


                if save_new_random_activities and self.dataset_sizes[col] >= self.min_dataset_size_for_pregenerated:
                    self.save_new_pregenerated = save_new_random_activities
                    activities_list_df = pd.concat(activities_list).reset_index(drop=True)
                    self.save_new_pregenerated_random_activities(activities_list_df, col)
                combined_activities_list.extend(activities_list)

            # save data
            self.random_enrichment = pd.concat(combined_activities_list).reset_index(drop=True)
            if save_random_experiments:
                rand_experiments = pd.concat(rand_experiments).reset_index(drop=True)
                rand_experiments['weight'] = 1
                self.random_experiments = rand_experiments.pivot(index=[config.KSTAR_ACCESSION, config.KSTAR_SITE],columns='Experiment',
                    values='weight').reset_index()

        else:
            all_rand_experiments = []
            combined_activities_list = []
            rand_exp_numbers = list(range(num_random_experiments))
            for col in tqdm(self.data_columns_from_scratch, desc="Calculating activities from random experiments for each dataset not using pregenerated random activities", disable=not show_taskbar):
                activities_list = []
                #group evidence by compendia sizes (number of compendia each site is found in)
                compendia_sizes = self.evidence_binary[self.evidence_binary[col] == 1].groupby(
                        selection_type).size()
                
                # generate random experiments and calculate activities
                for i in range(num_random_experiments):
                    results = calculate_random_activity_singleExperiment(compendia_sizes, filtered_compendia, self.networks, self.network_sizes, col, i, save_random_experiments)
                    # extract results (either activity alone or experiment)
                    if save_random_experiments:
                        act, exp = results
                        all_rand_experiments.append(exp)
                    else:
                        act = results
                    activities_list.append(act)

                # if want to save new precomputed data for each experiment
                if self.save_new_pregenerated and self.dataset_sizes[col] >= self.min_dataset_size_for_pregenerated:
                    activities_list_df = pd.concat(activities_list).reset_index(drop=True)
                    self.save_new_pregenerated_random_enrichment(activities_list_df, col)
                combined_activities_list.extend(activities_list)

                # combine all random activities
            self.random_enrichment = pd.concat(combined_activities_list).reset_index(drop=True)

            # reformat random experiments into single matrix, if wanting to save
            if save_random_experiments:
                all_rand_experiments = pd.concat(all_rand_experiments)
                all_rand_experiments['weight'] = 1
                self.random_experiments = all_rand_experiments.pivot(index=[config.KSTAR_ACCESSION, config.KSTAR_SITE],columns='Experiment',
                    values='weight').reset_index()



    def load_pregenerated_random_enrichment(self):
        """
        Load pre-generated random enrichment for the datasets that have been identified for use with pregenerated data.

        This function processes datasets that have associated pre-generated enrichment from random experiments. It identifies the appropriate pre-generated file based on the size of the dataset and appends the activities to the provided list.

        Updated Attributes
        ------------------
        self.pregenerated_random_activities : pandas.DataFrame
            DataFrame containing the loaded pre-generated random enrichment activities.
        

        """
                    #if the specific pregenerated experiments to use have not been determined, get them
        if not hasattr(self, 'selected_pregenerated_paths'):
            self.setup_pregenerated()

        if len(self.data_columns_with_pregenerated) > 0:
            #initialize lists
            self.num_random_experiments = 150
            pregen_activities_list = []

            for dataset in self.data_columns_with_pregenerated:
                #grab random enrichment file path, then load random activities
                matched_file_path = os.path.join(self.selected_pregenerated_paths[dataset], 'random_enrichment.tsv')
                rand_dataset_activities = load_random_enrichment(matched_file_path)

                #replace the generic prefix (most likely 'data:experiment') with the dataset of interest name
                pattern = re.compile(r':\d+$')
                if not rand_dataset_activities['data'].str.match(f"^{dataset}:\d+$").all():
                    rand_dataset_activities['data'] = rand_dataset_activities['data'].apply(
                        lambda x: f"{dataset}{pattern.search(x).group()}" if pattern.search(x) else x
                    )
                #add random activities to list
                pregen_activities_list.append(rand_dataset_activities)

            pregen_activities_list = pd.concat(pregen_activities_list)
            self.pregenerated_random_activities = pregen_activities_list
        else:
            self.pregenerated_random_activities = pd.DataFrame()

    def add_pregenerated_to_random_enrichment(self):
        """
        Combine pre-generated random activities with random enrichment, sort based on the "data" column,
        and reorganize the combined DataFrame based on the original column order in self.data_columns.

        If use_pregen_data is True and data_columns_from_scratch is None, uses only pre-generated activities.
        If use_pregen_data is True and data_columns_from_scratch exists, combines both pre-generated and
        newly calculated random activities.
        If use_pregen_data is False, uses only newly calculated random activities.

        Returns
        -------
        None
            Updates self.random_enrichment with the combined and sorted activities
        """

        if self.use_pregen_data:
            if len(self.data_columns_from_scratch) == 0:
                self.random_enrichment = self.pregenerated_random_activities
            elif len(self.data_columns_with_pregenerated) == 0:
                self.random_enrichment = self.random_enrichment
            else:
                # Combine the two DataFrames
                combined_df = pd.concat([self.random_enrichment, self.pregenerated_random_activities],
                                        ignore_index=True)
                combined_df['suffix'] = combined_df['data'].apply(lambda x: int(x.split(':')[-1]))
                combined_df['data_prefix'] = combined_df['data'].apply(lambda x: ':'.join(x.split(':')[:-1]))
                combined_df = combined_df.sort_values(by=["data_prefix", "suffix"])
                combined_df = combined_df.drop(columns=['suffix', 'data_prefix'])
                self.random_enrichment = combined_df.reset_index(drop=True)
        else:
            self.random_enrichment = self.random_enrichment


    def save_new_pregenerated_random_activities(self, activities_list_df, col):
        """
        Save the new precomputed random enrichment activities to a file.

        This function saves the provided DataFrame of random enrichment activities to a file, using the specified column name.

        Parameters
        ----------
        activities_list_df : pandas.DataFrame
            DataFrame containing the random enrichment activities to be saved.
        col : str
            Column name to be used for saving the activities.

        Returns
        -------
        None
        """
        # Determine the compendia directory based on phospho_type and class counts
        #get real experiment characteristics
        high_study_bias_perc = int(self.compendia_distribution[col][2])
        size = self.dataset_sizes[col]


        pregen_subdirectory = self.get_pregenerated_subdirectory(size, high_study_bias_perc)
        custom_exp_path = os.path.join(self.custom_pregenerated_path, pregen_subdirectory)


        #if path does not exist, create it
        if not os.path.exists(custom_exp_path):
            os.makedirs(custom_exp_path)


        save_random_enrichment(activities_list_df, custom_exp_path)
        file_path = os.path.join(custom_exp_path, 'random_enrichment.tsv')
        self.logger.info(f"New precomputed random enrichment data for {col} saved to {file_path}")

        # Save run information content, if not already present
        run_info_save_path = os.path.join(self.custom_pregenerated_path, "RUN_INFORMATION.txt")
        if not os.path.exists(run_info_save_path):
            run_info_file = os.path.join(self.network_directory, "RUN_INFORMATION.txt")
            shutil.copyfile(run_info_file, run_info_save_path)

            self.logger.info(f"RUN_INFORMATION.txt copied to {run_info_save_path}")

    def get_run_information_content(self):
        """
        Retrieve network information from RUN_INFORMATION.txt based on phospho_type.

        Reads the RUN_INFORMATION.txt file from the appropriate network directory based on
        the phospho_type ('Y' or 'ST'). The file contains network configuration details
        including unique ID, date, network specifications, and compendia counts.

        Returns
        -------
        Contents of RUN_INFORMATION.txt if found.
        'RUN_INFORMATION.txt file not found.' if the file doesn't exist.

        """
        base_path = self.network_directory
        if self.phospho_type == 'Y':
            run_info_path = os.path.join(base_path, 'Y', "RUN_INFORMATION.txt")
        elif self.phospho_type == 'ST':
            run_info_path = os.path.join(base_path, 'ST', "RUN_INFORMATION.txt")
        else:
            raise ValueError("Invalid phospho_type")

        try:
            with open(run_info_path, 'r') as file:
                return file.read()
        except FileNotFoundError:
            self.logger.warning(f"RUN_INFORMATION.txt file not found at: {run_info_path}")
            return "RUN_INFORMATION.txt file not found."


    def add_networks_batch(self, networks):
        for nid, network in networks.items():
            self.add_network(nid, network)
        self.num_networks = len(networks)

    def add_networks_from_directory(self):
        """
        Add all networks from the specified directory.

        Parameters
        ----------
        network_directory : str
            Path to the directory containing network files.
        """
        individual_network_dir = os.path.join(self.network_directory, 'INDIVIDUAL_NETWORKS')
        for filename in os.listdir(individual_network_dir):
            net_num = filename.split('_')[1]
            if filename.endswith('.tsv'):
                file_path = os.path.join(individual_network_dir, filename)
                network = pd.read_csv(file_path, sep='\t')
                #restrict to specified kinases if provided
                self.add_network(net_num, network)


    def add_network(self, network_id, network, network_size=None):
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
            self.network_sizes[network_id] = network.drop_duplicates(subset=[config.KSTAR_ACCESSION, config.KSTAR_SITE]).shape[0]
            self.logger.info(f'ADD NETWORK : Number of Accession Sites : {self.network_sizes[network_id]}')

    def get_run_date(self):
        """
        return date that kinase activities were run
        """
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
            raise ValueError("Invalid evidence DataFrame. Make sure the dataset has been mapped to the reference phosphoproteome (see KSTAR tutorial) prior to activity calculation. Evidence dataframe must contain 'KSTAR_ACCESSION' and 'KSTAR_SITE' columns.")
            #self.logger.warning(
            #    f"Evidence not set. Evidence columns must include '{config.KSTAR_ACCESSION}' and '{config.KSTAR_SITE}' keys")

    def create_binary_evidence(self, agg='mean', threshold=1.0, evidence_size=None, greater=True, min_evidence_size = 0, drop_empty_columns = True):
        """
        Returns a binary evidence data frame according to the parameters passed in for method for aggregating
        duplicates and considering whether a site is included as evidence or not

        Parameters
        ----------
        threshold : float, list
            threshold value used to filter rows
        evidence_size: None or int
            the number of sites to use for prediction for each sample. If a value is provided, this will override the threshold, and will instead obtain the N sites with the greatest abundance within each sample.
        agg : {'count', 'mean'}
            method to use when aggregating duplicate substrate-sites.
            'count' combines multiple representations and adds if values are non-NaN
            'mean' uses the mean value of numerical data from multiple representations of the same peptide.
                NA values are droped from consideration.
        greater: Boolean
            whether to keep sites that have a numerical value >=threshold (TRUE, default) or <=threshold (FALSE)
        min_evidence_size : int
            minimum number of sites required for a data column to be considered for activity calculation
        drop_empty_columns : bool
            whether to drop data columns with fewer than min_evidence_size sites

        Returns
        -------
        evidence_binary : pd.DataFrame
            Matches the evidence dataframe of the kinact object, but with 0 or 1 if a site is included or not.
            This is uniquified and rows that are never used are removed.


        """
        #make sure parameters are valid
        if evidence_size is not None and isinstance(evidence_size, int):
            if evidence_size < 0:
                raise ValueError("evidence_size must be a positive integer.")
            
            # if evidence size is given, ignore threshold (report this to user)
            if threshold is not None:
                self.logger.warning("Both evidence_size and threshold were provided. Ignoring threshold and using evidence_size.")
            #set threshold to None to indicate evidence size is being used
            threshold = None
        elif evidence_size is not None and not isinstance(evidence_size, int):
            raise ValueError("evidence_size must be an integer or None.")
        elif isinstance(threshold, list):
            if not all(isinstance(t, numbers.Number) for t in threshold):
                raise ValueError("All elements in threshold list must be integers or floats.")
        elif not isinstance(threshold, numbers.Number):
            raise ValueError("Threshold must be an integer, float, or list of integers or floats.")

        if agg not in ['mean', 'min','max', 'median']:
            raise ValueError("Aggregation method must be one of 'mean', 'min', 'max', or 'median'.")
        

        self.threshold = threshold
        self.evidence_size = evidence_size
        self.greater = greater

        #check data columns are valid and have evidence
        if drop_empty_columns:
            self.check_data_columns(min_evidence_size = min_evidence_size)

        # collapse sites into single row based on agg parameter
        evidence = self.evidence.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE])[self.data_columns].agg(
            agg).reset_index()

        # set the binary evidence for whether a site is included
        evidence_binary = evidence.copy()
        if evidence_size is None:
            if isinstance(threshold, list):
                #if multiple thresholds are provided, create new columns for each threshold
                new_data_columns = []
                for col in self.data_columns:
                    for thresh in threshold:
                        new_col = f"{col}_threshold={thresh}"
                        if greater:
                            evidence_binary[new_col] = (evidence_binary[col] >= thresh).astype(int)
                        else:
                            evidence_binary[new_col] = (evidence_binary[col] <= thresh).astype(int)
                        new_data_columns.append(new_col)

                #remove original data columns
                evidence_binary.drop(columns=self.data_columns, inplace=True)
                if drop_empty_columns:
                    #remove any columns that have less than min_evidence_size sites
                    for col in new_data_columns:
                        if evidence_binary[col].sum() < min_evidence_size:
                            evidence_binary.drop(columns=col, inplace=True)
                            new_data_columns.remove(col)
                
                self.data_columns = new_data_columns

            else:
                for col in self.data_columns:
                    if greater:
                        evidence_binary[col] = (evidence_binary[col] >= threshold).astype(int)
                    else:
                        evidence_binary[col] = (evidence_binary[col] <= threshold).astype(int)
        else:
            for col in self.data_columns:
                # check how many non_nan sites there (if less than N, set n to be equal to number of sites available)
                num_sites_available = evidence_binary[col].count()
                if num_sites_available >= evidence_size:
                    n = evidence_size
                else:
                    n = num_sites_available

                if greater:
                    max_indices = np.argsort(-evidence_binary[col].values)[0:n]
                    evidence_binary[col] = 0
                    col_loc = np.where(evidence_binary.columns == col)[0][0]
                    evidence_binary.iloc[max_indices, col_loc] = 1
                else:
                    min_indices = np.argsort(evidence_binary[col].values)[0:n]
                    evidence_binary[col] = 0
                    col_loc = np.where(evidence_binary.columns == col)[0][0]
                    evidence_binary.iloc[min_indices, col_loc] = 1

        # remove phosphorylation sites that were not selected in any experiment (useful for very large experiments where removing the need to copy data reduces time)
        evidence_binary.drop(evidence_binary[evidence_binary[self.data_columns].sum(axis=1) == 0].index, inplace=True)

        # add back compendia/study bias information to binary evidence
        compendia = config.HUMAN_REF_COMPENDIA[
            ['KSTAR_ACCESSION', 'KSTAR_SITE', 'KSTAR_NUM_COMPENDIA', 'KSTAR_NUM_COMPENDIA_CLASS']]
        evidence_binary = evidence_binary.merge(compendia, on=[config.KSTAR_ACCESSION, config.KSTAR_SITE], how='left')
        return evidence_binary

    def calculate_kinase_activities(self, agg='mean', threshold=1.0, evidence_size=None, greater=True, min_evidence_size = 0, PROCESSES=1):
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
        evidence_size : int or None
            the number of sites to use for prediction for each sample. If a value is provided, this will override the threshold, and will instead obtain the N sites with the greatest abundance within each sample (or lowest if greater=False).
        agg : {'count', 'mean'}
            method to use when aggregating duplicate substrate-sites.
            'count' combines multiple representations and adds if values are non-NaN
            'mean' uses the mean value of numerical data from multiple representations of the same peptide.
                NA values are droped from consideration.
        greater: Boolean
            whether to keep sites that have a numerical value >=threshold (TRUE, default) or <=threshold (FALSE)
        min_evidence_size : int
            minimum number of sites required for a data column to be considered for activity calculation
        PROCESSES : int
            number of processes to use for multiprocessing

        New attributes
        --------------
        activities : pd.DataFrame
            All hypergeometric activities for each sample, network, and kinase combination

        """

        # for each phosphoType, check that the
        self.aggregate = agg
        # if evidence size is given, ignore threshold, otherwise record threshold
        if evidence_size is None:
            self.threshold = threshold
            self.evidence_size = None
        else:
            self.threshold = None
            self.evidence_size = evidence_size
        self.greater = greater
        # self.set_data_columns(data_columns)

        self.evidence_binary = self.create_binary_evidence(agg=self.aggregate, threshold=self.threshold,
                                                           evidence_size=self.evidence_size, greater=self.greater, min_evidence_size = min_evidence_size)
        #get dataset sizes and compendia distribution of binary evidence
        self.dataset_sizes = {col: self.evidence_binary[col].sum() for col in self.data_columns}
        self.compendia_distribution = self.get_compendia_distribution()
        
        # if no data columns are provided use all columns that start with data:
        # data columns that filtered have no evidence are removed
        self.logger.info(f"Kinase Activity will be run on the following data columns: {','.join(self.data_columns)}")

        # MULTIPROCESSING
        if PROCESSES > 1:
            pool = multiprocessing.Pool(processes=PROCESSES)

            filtered_evidence_list = [self.evidence_binary[self.evidence_binary[col] == 1] for col in self.data_columns]
            networks = itertools.repeat(self.networks)
            network_sizes = itertools.repeat(self.network_sizes)
            iterable = zip(filtered_evidence_list, networks, network_sizes, self.data_columns)
            real_enrichment = pool.starmap(calculate_hypergeometric_activities, iterable)

        # SINGLE CORE PROCESSING
        else:
            real_enrichment = []
            for col in self.data_columns:
                filtered_evidence = self.evidence_binary[self.evidence_binary[col] == 1]
                act = calculate_hypergeometric_activities(filtered_evidence, self.networks, self.network_sizes, col)
                act['data'] = col
                real_enrichment.append(act)
        self.num_networks = len(self.network_sizes)

        self.real_enrichment = pd.concat(real_enrichment)
        return self.real_enrichment

    def summarize_activities(self, activities=None, method='median_activity', normalized=False):
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
        if activities is None:
            activities = self.agg_activities
        available_methods = list(activities.columns)
        available_methods.remove('data')
        if method not in available_methods:
            raise ValueError(
                f"the method '{method}' is not in the availble methods. \nAvailable methods include : {', '.join(available_methods)}")

        activity_summary = activities.pivot(index=config.KSTAR_KINASE, columns='data',
                                            values=method).reset_index().rename_axis(None, axis=1).set_index(
            config.KSTAR_KINASE)
        activity_summary = activity_summary[self.data_columns]
        return activity_summary

    def aggregate_activities(self, activities=None):
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
            activities = self.real_enrichment
        self.agg_activities = activities.groupby(['data', config.KSTAR_KINASE]).agg(
            median_activity=('kinase_activity', 'median'),
        ).reset_index()
        return self.agg_activities

    def find_pvalue_limits(self, data_columns, agg='count', threshold=1.0):
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
            filtered_evidence = evidence[evidence[col] > threshold]
            N = len(filtered_evidence.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE]).size())
            for nid, M in self.network_sizes.items():
                # kinases = list(kinase_sizes[nid].index)
                # for kinase in kinases:
                # n = int(kinase_sizes[nid].loc[kinase])
                n = kinase_sizes[nid]

                for k in range(N, 0, -1):
                    pvalue = stats.hypergeom.sf(
                        k=k,
                        M=M,
                        n=n,
                        N=N
                    )
                    # pvalue = 1 - prb
                    if pvalue > 0:
                        break
                row = {
                    'evidence': col,
                    'network': nid,
                    # 'kinase' : kinase,
                    'evidence_size': N,
                    'limit_size': k,
                    'kinase_size': n,
                    'network_size': M,
                    'p-value': pvalue
                }
                all_rows.append(row)
        all_limits = pd.DataFrame(all_rows)
        limit_summary = all_limits.groupby('evidence').mean()
        return all_limits, limit_summary

    def calculate_Mann_Whitney_activities_sig(self, show_taskbar=True, PROCESSES=1):
        """
        If random enrichment has already occurred, calculate the Mann-Whitney U test for comparing the array of p-values for real data to those of random data, across the number of networks used. It will also calculate the false positive rate for a pvalue by treating one of the randome experiments as 'real' and calculating how many random experiments have a more significant pvalue than the 'real' experiment.

        Parameters
        ----------
        show_taskbar : bool
            whether to show the tqdm taskbar for progress through datasets
        PROCESSES : int
            Number of processes to use for multiprocessing
        

        Returns
        -------

        """
        if not isinstance(self.random_enrichment, pd.DataFrame):
            raise ValueError("Random activities do not exist, please run kstar_activity.randomized_analysis")
        elif self.random_enrichment.empty:
            raise ValueError("Random activities not properly generated, please run KinaseActivity.get_random_activities()")


        #initialize output dataframes
        self.activities_mann_whitney = pd.DataFrame(index=self.kinases,columns=self.data_columns)
        self.activities_mann_whitney.index.name = 'KSTAR_KINASE'
        self.fpr_mann_whitney = pd.DataFrame(index=self.kinases, columns=self.data_columns)
        self.fpr_mann_whitney.index.name = 'KSTAR_KINASE'



        #group activities by sample and kinase
        real_grouped = self.real_enrichment.groupby(['data', 'KSTAR_KINASE'])['kinase_activity'].agg(list)
        self.random_enrichment['sample'] = self.random_enrichment['data'].apply(lambda x: ':'.join(x.split(':')[:-1])) #grab sample associated with random expeirment
        self.random_enrichment['rand_exp_num'] = self.random_enrichment['data'].apply(lambda x: int(x.split(':')[-1])) #grab random experiment number
        random_grouped = self.random_enrichment.groupby(['sample', 'KSTAR_KINASE', 'rand_exp_num'])['kinase_activity'].agg(list)


        # for every kinase and every dataset, calculate and assemble dataframes of activities and significance values
        for exp in tqdm(self.data_columns, desc='Calculating final activities with the mann whitney U test', disable = not show_taskbar):
            self.logger.info("MW Working on %s: " % (exp))
            pval_arr = []
            fpr_arr = []
            rand_stats = []
            #if using pregenerated data, load pregenerated fpr stats for FPR calculation
            if exp in self.data_columns_with_pregenerated:
                pregenerated_fpr_stats_path = os.path.join(self.selected_pregenerated_paths[exp], 'fpr_stats.tsv')
                if os.path.exists(pregenerated_fpr_stats_path):
                    pregenerated_fpr_stats = pd.read_csv(pregenerated_fpr_stats_path, sep='\t').to_dict(orient='list', index = True)
                else:
                    self.logger.warning(f"Pregenerated FPR stats file not found for {exp} at {pregenerated_fpr_stats_path}. FPR will be calculated without pregenerated data.")
                    pregenerated_fpr_stats = None
            else:
                pregenerated_fpr_stats = None

            if PROCESSES > 1 or pregenerated_fpr_stats is not None:  #note: its slower to multiprocess without pregenerated data due to overhead, so just use single processor
                #split dataframe ahead of time for each kinase
                real_enrichment_exp_list = [real_grouped.loc[(exp, kinase)] for kinase in self.kinases]
                random_enrichment_exp_list = [random_grouped.loc[(exp, kinase)] for kinase in self.kinases]
                pregenerated_fpr_stats_list = [pregenerated_fpr_stats[kinase] for kinase in self.kinases] if pregenerated_fpr_stats is not None else [None for _ in self.kinases]
                #iterate through kinases in parallel
                with concurrent.futures.ProcessPoolExecutor(max_workers=PROCESSES) as executor:
                    for pval, fpr, tmp_stats in executor.map(calculate_MannWhitney_one_experiment_one_kinase, real_enrichment_exp_list, random_enrichment_exp_list, pregenerated_fpr_stats_list):
                        pval_arr.append(pval)
                        fpr_arr.append(fpr)
                        rand_stats.append(tmp_stats)
                            
            else:
                for kinase in self.kinases:
                    tmp_real_grouped = real_grouped.loc[(exp, kinase)]
                    tmp_random_grouped = random_grouped.loc[(exp, kinase)]
                    pregenerated_fpr_stats_kinase = pregenerated_fpr_stats.get(kinase, None) if pregenerated_fpr_stats is not None else None
                    
                    pval, fpr, tmp_stats = calculate_MannWhitney_one_experiment_one_kinase(tmp_real_grouped, tmp_random_grouped, pregenerated_fpr_stats=pregenerated_fpr_stats_kinase)
                    pval_arr.append(pval)
                    fpr_arr.append(fpr)
                    rand_stats.append(tmp_stats)
                

            self.activities_mann_whitney[exp] = pval_arr
            self.fpr_mann_whitney[exp] = fpr_arr
            if self.save_new_pregenerated and self.dataset_sizes[exp] >= self.min_dataset_size_for_pregenerated:

                #get save directory
                pregen_subdirectory = self.get_pregenerated_subdirectory(
                    self.dataset_sizes[exp], 
                    int(self.compendia_distribution[exp][2])
                )
                custom_exp_path = os.path.join(self.custom_pregenerated_path, pregen_subdirectory)

                #convert to dataframe
                rand_stats = {kinase: stats for kinase, stats in zip(self.kinases, rand_stats)}
                rand_stats = pd.DataFrame.from_dict(rand_stats)
                #save
                rand_stats.to_csv(os.path.join(custom_exp_path, 'fpr_stats.tsv'), sep='\t', index=False)
                self.logger.info(f"New pregenerated FPR stats for {exp} saved to {os.path.join(custom_exp_path, 'fpr_stats.tsv')}")


    def get_param_dict(self, params_to_ignore = ['network_sizes', 'pregenerated_experiments_path', 'mann_whitney']):
        """
        Get a dictionary of important parameters needed to reinstantiate the KSTAR object
        """
        params = {}
        #iterate through attributes and save those needed to reinstantiate object
        for attr in vars(self):
            #check for basic types only
            if attr in params_to_ignore:
                pass
            elif isinstance(getattr(self, attr), (int, float, str, bool, type(None))):
                params[attr] = getattr(self, attr)

            elif isinstance(getattr(self, attr), list):
                #make sure list is of basic types
                if all(isinstance(i, (int, float, str, bool)) for i in getattr(self, attr)):
                    params[attr] = getattr(self, attr)

            elif isinstance(getattr(self, attr), dict):
                #make sure dict is of basic types
                if all(isinstance(k, str) and isinstance(v, (int, float, str, bool, type(None))) for k,v in getattr(self, attr).items()):
                    params[attr] = getattr(self, attr)
        return params

    def make_summary_pdf(self, regenerate_plots = False):
        """
        Create a summary PDF of the kinase activity results

        Parameters
        ----------
        regenerate_plots : bool
            Whether to regenerate plots even if they already exist
        """
        #construct param_dict
        param_dict = self.get_param_dict()
        summary_pdf = KSTAR_PDF(activities = self.activities_mann_whitney, fpr = self.fpr_mann_whitney, odir = self.odir + f'/RESULTS/{self.phospho_type}/', name = self.name, binarized_experiment = self.evidence_binary, param_dict = param_dict)
        summary_pdf.generate(regenerate_plots = regenerate_plots)

    def make_dotplot(self, include_evidence_sizes=True, **kwargs):
        """
        Create a dotplot of the kinase activity results

        Parameters
        ----------
        include_evidence_sizes : bool
            Whether to include evidence sizes in the dotplot
        **kwargs
            Additional keyword arguments to pass to the DotPlot initialization and make_complete_dotplot methods
        """
        #isolate kwargs unique to dotplot init
        init_kwargs = helpers.extract_relevant_kwargs(DotPlot.__init__, **kwargs)
        other_kwargs = {k:v for k,v in kwargs.items() if k not in init_kwargs}
        #initialize dotplot
        dotplot = DotPlot(values = self.activities_mann_whitney, fpr = self.fpr_mann_whitney, **init_kwargs)
        #make dotplot
        if include_evidence_sizes:
            dotplot.make_complete_dotplot(binary_evidence=self.evidence_binary, include_recommendations=True, phospho_type = self.phospho_type, **other_kwargs)
        else:
            dotplot.make_complete_dotplot(include_recommendations=True, phospho_type = self.phospho_type, **other_kwargs)





def get_available_kinases(phospho_type, network_dir=None):
    """
    Get all kinases available in the networks of the specified phospho_type

    Parameters
    ----------
    phospho_type : str
        'Y' or 'ST' for network type
    network_dir : str
        Path to the directory containing network files.

    Returns
    -------
    kinases : set
        Set of all kinases available in the networks
    """
    if network_dir is None:
        network_dir = config.NETWORK_DIR
    if phospho_type == 'Y':
        network_directory = os.path.join(network_dir, 'Y/INDIVIDUAL_NETWORKS/')
    elif phospho_type == 'ST':
        network_directory = os.path.join(network_dir, 'ST/INDIVIDUAL_NETWORKS/')

    #load in one network file to get list of kinases
    files = os.listdir(network_directory)
    for filename in files:
        if filename.endswith('.tsv'):
            network = pd.read_csv(os.path.join(network_directory, filename), sep='\t')
            kinases = np.sort(network['KSTAR_KINASE'].unique().tolist())
            break
    return kinases


def load_networks(phospho_type, network_dir = None, network_name = None, kinases = None):
    """
    Load all networks from the specified directory. into a dictionary

    Parameters
    ----------
    phospho_type : str
        'Y' or 'ST' for network type
    network_dir : str
        Path to the directory containing network files.
    kinases : list
        list of kinases to obtain edges for, if None, load all kinases
    """
    if network_dir is None:
        network_dir = config.NETWORK_DIR

    if network_name is None:
        network_name = config.NETWORK_NAME[phospho_type]

    #combine network name and phosphotype info
    network_directory = os.path.join(network_dir, phospho_type, network_name, 'INDIVIDUAL_NETWORKS')
    warned = False
    networks = {}
    for filename in os.listdir(network_directory):
        net_num = filename.split('_')[1]
        if filename.endswith('.tsv'):
            #load network file
            file_path = os.path.join(network_directory, filename)
            network = pd.read_csv(file_path, sep='\t')
            #if kinase list is provided, filter network to only include those kinases
            if kinases is not None:
                network = network[network['KSTAR_KINASE'].isin(kinases)]
                
                #if some kinases are provided that are not in the network, print warning
                missing_kinases = set(kinases) - set(network['KSTAR_KINASE'].unique())
                if len(missing_kinases) > 0 and not warned:
                    warned = True
                    print(f"Warning: The following kinases were not found in network {net_num}: {', '.join(missing_kinases)}. You can use the get_available_kinases function to see which kinases are available in the networks.")
            #add to network dictionary
            networks[net_num] = network
    return networks

def getFilteredCompendia(phospho_type, selection_type='KSTAR_NUM_COMPENDIA_CLASS'):
    """
    Get phosphorylation sites binned based on selection type
    """

    compendia = config.HUMAN_REF_COMPENDIA[
        config.HUMAN_REF_COMPENDIA[config.KSTAR_SITE].str.startswith(tuple(phospho_type))]

    compendia = compendia[[config.KSTAR_ACCESSION, config.KSTAR_SITE, selection_type]]
    compendia = compendia.groupby([config.KSTAR_ACCESSION,
                                    config.KSTAR_SITE]).max().reset_index()  # uniquify the compendia by KSTAR_ACCESSION and KSTAR_SITE
    sizes = compendia[selection_type].unique()
    filtered_compendia = {}
    for s in sizes:
        filtered_compendia[s] = compendia[compendia[selection_type] == s][
            [config.KSTAR_ACCESSION, config.KSTAR_SITE]]

    return filtered_compendia

def load_random_enrichment(file):
    """
    Given a file path, load random enrichment activities and melt into long format for use with KSTAR analysis
    """
    rand_dataset_enrichment = pd.read_csv(file, delimiter='\t')
    rand_dataset_enrichment = rand_dataset_enrichment.melt(id_vars=['network', 'KSTAR_KINASE'], var_name='data', value_name='kinase_activity')
    return rand_dataset_enrichment

def save_random_enrichment(rand_enrichment, odir):
    #reformat to reduce memory usage
    tmp_pivot = rand_enrichment.pivot(columns='data', index=['network', 'KSTAR_KINASE'], values='kinase_activity').reset_index()
    #save to file
    tmp_pivot.to_csv(f"{odir}/random_enrichment.tsv", sep='\t', index = False)


def calculate_hypergeometric_single_network(evidence, network, network_size, network_id):
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
    results = pd.DataFrame(counts, columns=['frequency'])

    results['kinase_activity'] = 1.0

    K = network.groupby(config.KSTAR_KINASE).size()

    kinases = counts.index
    for kin in kinases:
        k = 0
        if counts.loc[kin] > 0:
            k = counts.loc[kin] - 1
        prb = stats.hypergeom.sf(
            k=int(k),
            M=int(network_size),
            n=int(K.loc[kin]),
            N=int(N)
        )
        results.at[kin, 'kinase_activity'] = prb

    kinases = network[config.KSTAR_KINASE].unique()
    for kin in kinases:
        if kin not in results.index:
            results.at[kin, 'frequency'] = 0
            results.at[kin, 'kinase_activity'] = 1.0
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
    for network_id in networks.keys():  # calculate kinase activity within each network
        result = calculate_hypergeometric_single_network(evidence, networks[network_id], network_sizes[network_id],
                                                         network_id)
        results.append(result)

    # combine results into single dataframe
    hyp_act = pd.concat(results)
    hyp_act = hyp_act.reset_index()
    hyp_act['data'] = name
    return hyp_act


#def calculate_random_activity_singleExperiment(filtered_evidence, filtered_compendia, networks, network_sizes, data_col,
#                                               num_random_experiments=150, save_experiments=False):
#    if save_experiments:
#        real_enrichment = []
#        all_rand_experiments = []
#        for i in range(num_random_experiments):
#            name = f'{data_col}:{i}'
#            rand_experiment = generate_random_experiments.build_single_filtered_experiment(filtered_evidence,
#                                                                                           filtered_compendia, name,
#                                                                                           selection_type='KSTAR_NUM_COMPENDIA_CLASS')
#            all_rand_experiments.append(rand_experiment)
#            act = calculate_hypergeometric_activities(rand_experiment, networks, network_sizes, name)
#            real_enrichment.append(act)
#        return real_enrichment, all_rand_experiments
#    else:
#        real_enrichment = []
#        for i in range(num_random_experiments):
#            name = f'{data_col}:{i}'
#            rand_experiment = generate_random_experiments.build_single_filtered_experiment(filtered_evidence,
#                                                                                           filtered_compendia, name,
#                                                                                           selection_type='KSTAR_NUM_COMPENDIA_CLASS')
#            act = calculate_hypergeometric_activities(rand_experiment, networks, network_sizes, name)
#            real_enrichment.append(act)
#        return real_enrichment



def calculate_random_activity_singleExperiment(compendia_sizes, filtered_compendia, networks, network_sizes, data_col,
                                                exp_num, 
                                                save_experiments=False):
    """
    Given the properites of a given experiment, generate a single random experiment with the same properties and calculate its
    hypergeometric activities across all networks.
    """
    name = f'{data_col}:{exp_num}'
    rand_experiment = generate_random_experiments.build_single_filtered_experiment(compendia_sizes, filtered_compendia,
                                                                                   name,
                                                                                   selection_type='KSTAR_NUM_COMPENDIA_CLASS')
    act = calculate_hypergeometric_activities(rand_experiment, networks, network_sizes, name)
    if save_experiments:
        return act, rand_experiment
    else:
        return act

"""
****************************************
Methods for Mann Whitney analysis
****************************************
"""


def calculate_fpr_Mann_Whitney(random_kinase_activity_array):
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
    # calculate the significance by taking each experiment
    [number_sig_trials, n] = random_kinase_activity_array.shape
    #if number_sig_trials > m:
    #    print("Warning, using %d, maximum number for significance" % (m))
    #    number_sig_trials = m
    random_stats = np.empty([number_sig_trials])
    for i in range(0, number_sig_trials):
        # take out one vector as real
        sample = random_kinase_activity_array[i, :]
        bgnd = np.delete(random_kinase_activity_array, i, 0)  # remove the sample before testing
        [stat, random_stats[i]] = stats.mannwhitneyu(-np.log10(sample), -np.log10(bgnd.reshape(bgnd.size)),
                                                     alternative='greater')
    return random_stats


def calculate_MannWhitney_one_experiment_one_kinase(real_activities_grouped_one_kinase_one_exp, rand_activities_grouped_one_kinase_one_exp, pregenerated_fpr_stats = None):
    """
    For a given kinact object, where random generation and activity has already been run, this will calculate the Mann-Whitney U test between the p-values across all networks for the given experiment name and from the random networks. It will also calculate the significance value for the given test based on the target_alpha value by using each random set as a real set to bootstrap.

    Parameters
    ----------
    real_activities_grouped_one_kinase_one_exp : pandas.Series
        Multi-index series containing kinase activities for a given experiment grouped by kinase. Should be indexed by kinase, with each entry being a list of activities across all networks. You can obtain this by performing a groupby on the real enrichment DataFrame.
    rand_activities_grouped_one_kinase_one_exp : pandas.Series
        Multi-index series object containing kinase activities for random experiments associated with the same sample as kinact_activities_sub, grouped by kinase and random experiment number. Should be indexed by kinase and rand experiment number, with each entry being a list of activities across all networks. You can obtain this by performing a groupby on the random enrichment DataFrame.
    kinase : str
        Kinase name to measure significance for.
    experiment : str
        Experiment name to measure significance for.
    pregenerated_fpr_stats : dict or None
        Pre-generated random statistics for FPR calculation, matching the random activities provided. If None, will generate new random statistics using calculate_fpr_Mann_Whitney.

    Returns
    -------
    p_value : float
        p-value that results from Mann-Whitney U test.
    fpr_value : float
        The false positive rate where the p_value for the real experiment lies, given the random experiments.
    """
    kinase_activity_list = real_activities_grouped_one_kinase_one_exp
    #reformat random activities into array
    random_kinase_activity_array = np.vstack(rand_activities_grouped_one_kinase_one_exp)

    #remove one random experiment to make real and random Mann Whitney tests comparable (same number of random experiments used for comparison)
    i = np.random.randint(0, random_kinase_activity_array.shape[0])
    bgnd = np.delete(random_kinase_activity_array, i, 0)

    #compare real enrichment to random background (size = num_random_experiments - 1)
    [stat, p_value] = stats.mannwhitneyu(-np.log10(kinase_activity_list),
        -np.log10(np.concatenate(bgnd)),alternative='greater')
    
    # Calculate FPR using the helper function
    if pregenerated_fpr_stats is not None:
        randomStats = pregenerated_fpr_stats
    else:
        randomStats = calculate_fpr_Mann_Whitney(random_kinase_activity_array)
    fpr_value = calculate_fpr.single_pvalue_fpr(randomStats, p_value)
    return p_value, fpr_value, randomStats


"""
****************************************
Methods for running KSTAR pipeline
****************************************
"""


def enrichment_analysis(experiment, odir, name='experiment', phospho_types=['Y', 'ST'], data_columns=None, agg='mean',threshold=1.0, evidence_size=None, greater=True, min_evidence_size=0, allow_column_loss = True, kinases=None, PROCESSES=1, **kwargs):
    """
    Function to establish a kstar KinaseActivity object from an experiment with an activity log
    add the networks, calculate, aggregate, and summarize the hypergeometric enrichment into a final activity object. Should be followed by
    randomized_analyis, then Mann_Whitney_analysis.

    Parameters
    ----------
    experiment: pandas df
        experiment dataframe that has been mapped, includes KSTAR_SITE, KSTAR_ACCESSION, etc.
    odir : str
        path to where you would like logger and output saved
    name : str
        name to use for outputs
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
    threshold : float or dict
        threshold value used to filter rows. If provided as a dictionary, keys should be 'Y' and/or 'ST' with float values for each phospho_type.
    evidence_size : int or dict
        size of evidence to use for filtering. If provided as a dictionary, keys should be 'Y' and/or 'ST' with int values for each phospho_type. Will overide threshold if both provided.
    min_evidence_size : int
        minimum size of evidence to run kinase activity on. Default 0, meaning any data column with at least one site will be run on
    greater: Boolean
        whether to keep sites that have a numerical value >=threshold (TRUE, default) or <=threshold (FALSE)
    PROCESSES : int
        number of processes to use for parallel computation, by default 1.
    **kwargs
        Additional keyword arguments to pass to the KinaseActivity class

    Returns
    -------
    kinactDict: dictionary of Kinase Activity Objects
        Outer keys are phosphoTypes run 'Y' and 'ST'
        Includes the activities dictionary (see calculate_kinase_activities)
        aggregation of activities across networks (see aggregate activities)
        activity summary (see summarize_activities)

    """
    #make sure either threshold or evidence_size is provided
    if threshold is None and evidence_size is None:
        raise ValueError("Either threshold or evidence_size must be provided to run enrichment_analysis.")
    
    #check to make sure threshold and evidence_size are properly formatted
    if isinstance(threshold, dict):
        for pt in phospho_types:
            if pt not in threshold:
                raise ValueError("When providing threshold as a dictionary, it must include keys for all phospho_types being run.")
    else:
        threshold = {pt: threshold for pt in phospho_types}
    
    if isinstance(evidence_size, dict):
        for pt in phospho_types:
            if pt not in evidence_size:
                raise ValueError("When providing evidence_size as a dictionary, it must include keys for all phospho_types being run.")
    elif evidence_size is None or isinstance(evidence_size, int):
        evidence_size = {pt: evidence_size for pt in phospho_types}
    else:
        raise TypeError("evidence_size must be provided as an int or as a dictionary with keys for all phospho_types being run.")

    if isinstance(kinases, dict):
        for pt in phospho_types:
            if pt not in kinases:
                raise ValueError("When providing evidence_size as a dictionary, it must include keys for all phospho_types being run.")
    else:
        kinases = {pt: kinases for pt in phospho_types}

    if not isinstance(greater, bool):
        raise TypeError("greater parameter must be a boolean value of True or False.")

    #make sure odir exists
    if not os.path.exists(odir):
        raise ValueError(f"Output directory {odir} does not exist. Please create it before running enrichment_analysis.")
    
    kinact_dict = {}
    # For each phosphoType of interest, establish a kinase activity object on a filtered dataset and run, aggregate, and summarize activity
    if phospho_types is None:
        raise TypeError("phospho_types must be provided as a list containing 'Y' and/or 'ST' (['Y', 'ST'], ['Y'], or ['ST'])")
    elif not isinstance(phospho_types, list):
        raise TypeError("phospho_types must be provided as a list containing 'Y' and/or 'ST' (['Y', 'ST'], ['Y'], or ['ST'])")
    elif len(phospho_types) == 0:
        raise ValueError("phospho_types must contain at least one of 'Y' or 'ST'")

    
    for phospho_type in phospho_types:
        # filter the experiment (log how many are of that type)
        if phospho_type == 'ST':
            experiment_sub = experiment[
                (experiment.KSTAR_SITE.str.contains('S')) | (experiment.KSTAR_SITE.str.contains('T'))]
        #    log.info("Running Serine/Threonine Kinase Activity Analysis")
        elif phospho_type == 'Y':
            experiment_sub = experiment[(experiment.KSTAR_SITE.str.contains('Y'))]
        #    log.info("Running Tyrosine Kinase Activity Analysis")

        else:
            raise ValueError("ERROR: Did not recognize phosphoType %s, which should only include 'Y' or 'ST' " % (phospho_type))
        
        #make sure there are any sites of that type
        if experiment_sub.shape[0] == 0:
            print(f"No phosphorylation sites of type {phospho_type} found in dataset, skipping kinase activity calculations for this phospho_type.")
            continue

        #initialize KinaseActivity object

        kinact = KinaseActivity(experiment_sub, odir = odir, name = name, data_columns=data_columns, phospho_type=phospho_type, kinases=kinases[phospho_type], **kwargs)
        #make sure thresholding parameters are valid (unless evidence_size is provided)
        if evidence_size[phospho_type] is None:
            good_params = kinact.check_valid_threshold(threshold = threshold[phospho_type], greater = greater, agg = agg, min_evidence_size = min_evidence_size, allow_column_loss=allow_column_loss)
            if not good_params:
                print(f"Skipping kinase activity calculations for {phospho_type} kinases due to no data columns with sufficient evidence after thresholding. Please consider adjusting threshold.")
                continue

        #perform hypergeometric enrichment across all networks
        kinact.calculate_kinase_activities(agg=agg, threshold=threshold[phospho_type], evidence_size=evidence_size[phospho_type], greater=greater, min_evidence_size=min_evidence_size, PROCESSES=PROCESSES)
        kinact_dict[phospho_type] = kinact

        
    if len(kinact_dict) == 0:
        raise ValueError("No KinaseActivity objects were successfully created. Please check your input experiment and parameters.")

    return kinact_dict

def randomized_analysis(kinact_dict, show_taskbar = True, **kwargs):
    """
    Perform randomized analysis on kinase activity data.

    Parameters
    ----------
    kinact_dict : dict
        Dictionary containing kinase activity data.
    kwargs : keyword arguments
        Additional keyword arguments for random activity generation passed to KinaseActivity.get_random_activities method.

        These can include:
        num_random_experiments : int, optional
            Number of random experiments to generate, by default 150.
        use_pregen_data : bool, optional
            Whether to use pre-generated data, by default False.
        max_diff_from_pregenerated : float, optional
            Maximum fractional difference allowed from pre-generated data, by default 0.25.
        min_dataset_size_for_pregenerated : int, optional
            Minimum dataset size to use pre-generated data, by default 150.
        default_pregen_only : bool, optional
            Whether to only use default pre-generated data (and not any activities in custom path), by default False.
        require_pregenerated : bool, optional
            Whether to require pre-generated data, by default False. This will ensure fast performance, but may result in some data columns being dropped
        custom_pregenerated_path : str, optional
            Directory to save new precomputed data, by default None.
        save_random_experiments : bool, optional
            Whether to save the generated random experiments, by default None.
        save_new_random_activities : bool, optional
            Whether to save new precomputed data, by default None.
        PROCESSES : int, optional
            Number of processes to use for parallel computation, by default 1.

    Returns
    -------
    None
    """
    for phospho_type, kinact in kinact_dict.items():
        #if network_hash is None:
        #    if phospho_type == 'Y':
        #        network_hash = config.NETWORK_HASH_Y
        #    elif phospho_type == 'ST':
        #        network_hash = config.NETWORK_HASH_ST
        #if not re.fullmatch(r'[a-fA-F0-9]{64}', network_hash):
        #    raise ValueError("network_hash must be a valid SHA-256 hash")

        kinact.get_random_activities(show_taskbar=show_taskbar, **kwargs)
        #except Exception as e:
        # Ensure `random_experiments` is stored in `kinact_dict`
        #kinact_dict[phospho_type].random_experiments = kinact.random_experiments

def Mann_Whitney_analysis(kinact_dict, show_taskbar=True, PROCESSES=1):
    """
    For a kinact_dict, where random generation and activity has already been run for the phospho_types of interest,
    this will calculate the Mann-Whitney U test for comparing the array of p-values for real data
    to those of random data, across the number of networks used.
    It will also calculate the false positive rate for a pvalue, given observations of a random bootstrapping analysis

    Parameters
    ----------
    kinact_dict: dictionary
        A dictionary of kinact objects, with keys 'Y' and/or 'ST'
    PROCESSES: int
        number of processes to use for parallel computation, by default 1.
    """

    for phospho_type, kinact in kinact_dict.items():
        kinact.calculate_Mann_Whitney_activities_sig(show_taskbar=show_taskbar, PROCESSES=PROCESSES)

def run_kstar_analysis(experiment, odir, name='experiment', phospho_types=['Y', 'ST'], data_columns=None, threshold=1.0, evidence_size=None, greater=True, save_output=True, show_taskbar = True, PROCESSES=1, **kwargs):
    """
    Given a mapped experiment, run the KSTAR analysis pipeline.

    Parameters
    ----------
    experiment: DataFrame
        Mapped experiment data
    odir: string
        Output directory
    name: string
        Name of the experiment
    phospho_types: list
        List of phospho types to analyze
    network_dir: string
        Directory containing network data
    data_columns: list
        Columns to use from the data
    agg: string
        Aggregation method
    threshold: float
        Threshold for analysis
    evidence_size: int
        Size of evidence
    greater: bool
        Whether to use greater comparison
    PROCESSES: int
        Number of processes to use
    **kwargs
        Additional keyword arguments for enrichment_analysis, randomized_analysis, and save_kstar functions.

    """
    #find relevant kwargs
    enrichment_kwargs = helpers.extract_relevant_kwargs(enrichment_analysis, **kwargs)
    randomized_kwargs = helpers.extract_relevant_kwargs(KinaseActivity.get_random_activities, **kwargs)
    save_kwargs = helpers.extract_relevant_kwargs(save_kstar, **kwargs)
    identified_kwargs = set(enrichment_kwargs.keys()).union(set(randomized_kwargs.keys())).union(set(save_kwargs.keys()))

    #throw error if any unexpected kwargs are provided
    unexpected_kwargs = set(kwargs.keys()).difference(identified_kwargs)
    if len(unexpected_kwargs) > 0:
        raise ValueError(f"Unexpected keyword arguments provided: {', '.join(unexpected_kwargs)}. Please check the function parameters for enrichment_analysis, randomized_analysis, and save_kstar.")
    
    #start enrichment analysis
    print('Starting kinase-substrate enrichment analysis...')
    kinact_dict = enrichment_analysis(experiment, odir = odir, name = name, phospho_types = phospho_types, threshold = threshold, data_columns=data_columns, greater = greater, evidence_size=evidence_size, PROCESSES = PROCESSES, **enrichment_kwargs)
    #
    print('Starting calculation of random activities...')
    randomized_analysis(kinact_dict, show_taskbar = show_taskbar, PROCESSES=PROCESSES, **randomized_kwargs)
    #
    print('Comparing kinase-substrate enrichment from the real experiment to random experiments...')
    Mann_Whitney_analysis(kinact_dict, show_taskbar=show_taskbar, PROCESSES = PROCESSES)


    if save_output:
        print(f'Saving KSTAR analysis results in {odir}/RESULTS/')
        save_kstar(kinact_dict, name, odir, **save_kwargs)

    print('Done.')
    return kinact_dict

    
    



def save_kstar(kinact_dict, name, odir, minimal = True, ftype = 'tsv', param_format = 'json'):
    """
    Having performed kinase activities (run_kstar_analyis), save each of the important dataframes, minimizing the memory storage needed to get back
    to a rebuilt version for plotting results and analysis. For each phospho_type in the kinact_dict, at a minimum, this will save the binarized evidence, mann whitney activities and fpr dataframes, and parameters used during run. If you would like to save all files (hypergeometric and random enrichment intermediate files), set minimal = False.

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
    minimal: bool
        Whether to save only minimal files or all intermediate files
    ftype: {'tsv', 'csv'}
        Format to save dataframes in, either tsv or csv
    param_format: {'pickle', 'json'}
        Format to save parameter dictionary in, either pickle or json. Json is recommended for easier human readability
    

    Returns
    -------
    Nothing

    """

    if not os.path.exists(f"{odir}/RESULTS"):
        os.mkdir(f"{odir}/RESULTS")

    param_dict = {}

    for phospho_type in kinact_dict:

        kinact = kinact_dict[phospho_type]
        name_out = f"{name}_{phospho_type}"

        param_temp = kinact.get_param_dict()

        if ftype == 'tsv':
            kinact.evidence_binary.to_csv(f"{odir}/RESULTS/{phospho_type}/{name_out}_binarized_experiment.tsv", sep='\t', index=False)
        elif ftype == 'csv':
            kinact.evidence_binary.to_csv(f"{odir}/RESULTS/{phospho_type}/{name_out}_binarized_experiment.csv", index=False)

        if hasattr(kinact, 'activities_mann_whitney'):
            param_temp['mann_whitney'] = True
            if ftype == 'tsv':
                kinact.activities_mann_whitney.to_csv(f"{odir}/RESULTS/{phospho_type}/{name_out}_mann_whitney_activities.tsv", sep='\t',
                                                    index=True)
                kinact.fpr_mann_whitney.to_csv(f"{odir}/RESULTS/{phospho_type}/{name_out}_mann_whitney_fpr.tsv", sep='\t', index=True)
            elif ftype == 'csv':
                kinact.activities_mann_whitney.to_csv(f"{odir}/RESULTS/{phospho_type}/{name_out}_mann_whitney_activities.csv", index=True)
                kinact.fpr_mann_whitney.to_csv(f"{odir}/RESULTS/{phospho_type}/{name_out}_mann_whitney_fpr.csv", index=True)

        param_dict[phospho_type] = param_temp

        #if indicated, save additional files (this is really just the hypergeometric and random enrichment intermediate files)
        if not minimal:
            if ftype == 'tsv':
                kinact.real_enrichment.to_csv(f"{odir}/RESULTS/{phospho_type}/{name_out}_hypergeometric_enrichment.tsv", sep='\t')
                if hasattr(kinact, 'random_experiments') and kinact.random_experiments is not None:
                    kinact.random_experiments.to_csv(f"{odir}/RESULTS/{phospho_type}/{name_out}_random_experiments.tsv", sep='\t', index=False)
                if hasattr(kinact, 'random_enrichment'):
                    kinact.random_enrichment.to_csv(f"{odir}/RESULTS/{phospho_type}/{name_out}_random_enrichment.tsv", sep='\t')
            elif ftype == 'csv':
                kinact.real_enrichment.to_csv(f"{odir}/RESULTS/{phospho_type}/{name_out}_hypergeometric_enrichment.csv")
                if hasattr(kinact, 'random_experiments') and kinact.random_experiments is not None:
                    kinact.random_experiments.to_csv(f"{odir}/RESULTS/{phospho_type}/{name_out}_random_experiments.csv", index=False)
                if hasattr(kinact, 'random_enrichment'):
                    kinact.random_enrichment.to_csv(f"{odir}/RESULTS/{phospho_type}/{name_out}_random_enrichment.csv")

    # save the parameters in a pickle or json file for reinstantiating object information
    if param_format == 'pickle':
        pickle.dump(param_dict, open(f"{odir}/RESULTS/{name}_params.p", "wb"))
    elif param_format == 'json':
        with open(f"{odir}/RESULTS/{name}_params.json", 'w') as json_file:
            json.dump(param_dict, json_file, indent=4)
    else:
        raise ValueError("param_format must be either 'pickle' or 'json'")
    
def save_kstar_slim(kinact_dict, name, odir):
    raise DeprecationWarning("save_kstar_slim is deprecated, please use save_kstar with minimal=True instead.")

def load_param_dict(name, odir):
    """
    Given the name and output directory of a saved kstar analyis, load the parameters used for KSTAR analysis, which should be stored as either a pickle or json file.

    Parameters
    ----------
    name: string
        The name to used when saving activities and mapped data
    odir:  string
        Output directory of saved files and parameter pickle
    """
    # First check for the param file
    if os.path.exists(f"{odir}/RESULTS/{name}_params.json"):
        with open(f"{odir}/RESULTS/{name}_params.json", 'r') as json_file:
            param_dict = json.load(json_file)
    elif os.path.exists(f"{odir}/RESULTS/{name}_params.p"):
        param_dict = pickle.load(open(f"{odir}/RESULTS/{name}_params.p", "rb"))
    elif not os.path.exists(f"{odir}/RESULTS/") and not os.path.exists(f"{odir}/RESULTS/"):
        print(f"ERROR: Cannot find RESULTS directory in output directory: {odir}/RESULTS/")
        param_dict = None
    else:
        print(f"ERROR: Cannot find parameter dictionary file in RESULTS: {odir}/RESULTS/{name}_params.p or {odir}/RESULTS/{name}_params.json")
        param_dict = None
    return param_dict
    

def from_kstar(name, odir, ftype = 'tsv'):
    """
    Given the name and output directory of a saved kstar analyis, load the parameters and minimum dataframes needed for reinstantiating a kinact object
    This minimum list will allow you to repeat normalization or mann whitney at a different false positive rate threshold and plot results.

    Parameters
    ----------
    name: string
        The name to used when saving activities and mapped data
    odir:  string
        Output directory of saved files and parameter pickle
    """

    # First check for the param file
    param_dict = load_param_dict(name, odir)
    if param_dict is None:
        return
    
    kinact_dict = {}
    for phospho_type in param_dict.keys():
        params = param_dict[phospho_type]
        print(f"Loading data for phospho_type: {phospho_type}")
        name_out = f"{name}_{phospho_type}"


        # check that the minimum file set exists so we can use binary_evidence file as the experiment
        if ftype == 'tsv':
            evidence_binary = pd.read_csv(f"{odir}/RESULTS/{phospho_type}/{name_out}_binarized_experiment.tsv", sep='\t')
        elif ftype == 'csv':
            evidence_binary = pd.read_csv(f"{odir}/RESULTS/{phospho_type}/{name_out}_binarized_experiment.csv")
        else:
            raise ValueError("ftype must be either 'tsv' or 'csv'")
        #grab additional parameters needed to reinstate object
        #network_dir = params.get('network_directory', None)
        #network_name = params.get('network_name', None)


        #load activity logger
        kinact = KinaseActivity(evidence_binary, odir=odir, name=name, phospho_type=phospho_type)

        #kinact.real_enrichment = pd.read_csv(f"{odir}/RESULTS/{name_out}_real_enrichment.tsv", sep='\t', index_col=0)
        kinact.evidence_binary = evidence_binary

        # read mann_whitney and load
        if ftype == 'tsv':
            kinact.activities_mann_whitney = pd.read_csv(f"{odir}/RESULTS/{phospho_type}/{name_out}_mann_whitney_activities.tsv",
                                                            sep='\t', index_col=config.KSTAR_KINASE)
            kinact.fpr_mann_whitney = pd.read_csv(f"{odir}/RESULTS/{phospho_type}/{name_out}_mann_whitney_fpr.tsv", sep='\t',
                                                    index_col=config.KSTAR_KINASE)
        elif ftype == 'csv':
            kinact.activities_mann_whitney = pd.read_csv(f"{odir}/RESULTS/{phospho_type}/{name_out}_mann_whitney_activities.csv",
                                                            index_col=config.KSTAR_KINASE)
            kinact.fpr_mann_whitney = pd.read_csv(f"{odir}/RESULTS/{phospho_type}/{name_out}_mann_whitney_fpr.csv", index_col=config.KSTAR_KINASE)
        else:
            raise ValueError("ftype must be either 'tsv' or 'csv'")
        #params.pop('mann_whitney', None)


        for param_name in params:
            setattr(kinact, param_name, params[param_name])
        kinact_dict[phospho_type] = kinact
    return kinact_dict


def from_kstar_nextflow(name, odir, log=None):
    """
    Given the name and output directory of a saved kstar analyis from the nextflow pipeline, load the results into new kinact object with
    the minimum dataframes required for analysis (binary experiment, hypergeometric activities, normalized activities, mann whitney activities)

    Parameters
    ----------
    name: string
        The name to used when saving activities and mapped data
    odir:  string
        Output directory of saved files
    log: logger
        logger used when loading nextflow data into kinase activity object. If not provided, new logger will be created.
    """
    # create new logger if not provided
    if log is None:
        log = helpers.get_logger(f"activity_{name}", f"{odir}/activity_{name}.log")

    kinact_dict = {}
    for phospho_type in ['Y', 'ST']:
        # check to see if folder for phosphotype is present (have Y or ST activities been calculated)
        if os.path.isdir(f'{odir}/{phospho_type}'):
            # check that the minimum file set exists so we can use binary_evidence file as the experiment
            evidence_binary = pd.read_csv(f"{odir}/{phospho_type}/binary_experiment/{name}_binarized_experiment.tsv",
                                          sep='\t')

            kinact = KinaseActivity(evidence_binary, log, phospho_type=phospho_type)

            kinact.real_enrichment = pd.read_csv(
                f"{odir}/{phospho_type}/hypergeometric_activity/{name}_real_enrichment.tsv", sep='\t', index_col=0)
            kinact.activities = pd.read_csv(f"{odir}/{phospho_type}/hypergeometric_activity/{name}_activities.tsv",
                                            sep='\t', index_col=config.KSTAR_KINASE)
            kinact.evidence_binary = evidence_binary

            # read mann_whitney and load
            kinact.activities_mann_whitney = pd.read_csv(
                f"{odir}/{phospho_type}/mann_whitney/{name}_mann_whitney_activities.tsv", sep='\t',
                index_col=config.KSTAR_KINASE)
            kinact.fpr_mann_whitney = pd.read_csv(f"{odir}/{phospho_type}/mann_whitney/{name}_mann_whitney_fpr.tsv",
                                                  sep='\t', index_col=config.KSTAR_KINASE)

            kinact_dict[phospho_type] = kinact

    # check to see if any files were loaded
    if kinact_dict.keys() == 0:
        print('KSTAR results were not found for either Y or ST')

    return kinact_dict

