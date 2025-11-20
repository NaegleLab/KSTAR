import json
import os
import re
import pickle
import itertools
import warnings

import pandas as pd
import numpy as np
import scipy.stats as stats
import multiprocessing

import seaborn as sns
import matplotlib.pyplot as plt
import concurrent.futures
from functools import partial

from tqdm import tqdm
from datetime import datetime
from itertools import repeat
from collections import defaultdict
from kstar import config, helpers, mapping
from kstar.random_experiments import generate_random_experiments, calculate_fpr


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

    Attributes
    ----------
    -------------------
    Upon Initialization
    -------------------
    evidence: pandas dataframe
        inputted evidence column
    data_columns: list
        list of columns containing abundance values, which will be used to determine which sites will be used as evidence. If inputted data_columns parameter was None, this lists includes in column in evidence prefixed by 'data:'
    logger : Logger object
        keeps track of kstar analysis, including any errors that occur
    phospho_type: string
        indicated phosphomod of interest
    network_directory: string
        directory where kinase substrate networks can be downloaded, as indicated in config.py
    normalized: bool
        indicates whether normalization analysis has been performed
    aggregate: string
        the type of aggregation to use when determining binary evidence, either 'count' or 'mean'. Default is 'count'.
    threshold: float
        cutoff to use when determining what sites to use for each experiment
    greater: bool
        indicates whether sites with greater or lower abundances than the threshold will be used
    run_data: string
        indicates the date that kinase activity object was initialized
    min_dataset_size_for_pregenerated: int
        minimum dataset size required to use pregenerated random activities
    max_diff_from_pregenerated: float
        maximum percent difference between dataset size and pregenerated random activity size to use pregenerated data

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
        contains information about the sites randomly sampled for each random experiment

    random_kinact: KinaseActivity object
        KinaseActivity object containing random activities predicted from each of the random experiments

    ---------------------------
    After Mann Whitney Analysis
    ---------------------------
    activities_mann_whitney: pandas dataframe
        p-values obtained from comparing the real distribution of p-values to the distribution of p-values from random datasets, based
        the Mann Whitney U-test
    fpr_mann_whitney: pandas dataframe
        false positive rates for predicted kinase activities

    """

    def __init__(self, evidence, odir, name = 'experiment', data_columns=None, phospho_type='Y', kinases = None, network_dir = None, use_pregen = None, default_pregen_only = True, custom_pregen_dir = None, min_dataset_size_for_pregenerated=150, max_diff_from_pregenerated=0.20, logger = None, network_name = None):
        self.odir = odir
        self.name = name
        self.phospho_type = phospho_type
        self.set_evidence(evidence)

        self.networks = defaultdict()
        self.network_sizes = defaultdict()
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
            self.network_directory = config.NETWORK_SUBDIR[self.phospho_type]

        #check to make sure network directory exists
        if not os.path.exists(self.network_directory):
            raise FileNotFoundError(f"Network directory not found at: {self.network_directory}. Please download networks using config.install_network_files() or indicate where networks are found with config.update_network_directory().")
        

        #load network and meta information
        self.network_name = network_name if network_name is not None else config.NETWORK_NAME[self.phospho_type]
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
        self.randomized = False
        self.random_experiments = None

        # added fields for pregenerated_random
        self.min_dataset_size_for_pregenerated = min_dataset_size_for_pregenerated
        self.max_diff_from_pregenerated = max_diff_from_pregenerated
        self.random_activities_list = None
        self.compendia_distribution = None
        self.data_columns_from_scratch = None
        self.use_pregen_data = use_pregen if use_pregen is not None else config.USE_PREGENERATED_RANDOM_ACTIVITIES
        self.save_new_precompute = config.SAVE_NEW_RANDOM_ACTIVITIES
        self.pregenerated_experiments_path = self.network_directory + "/RANDOM_ACTIVITIES/" 


        #default path should be located in network directory, check to make sure it and any provided custom directory exists
        self.default_pregen_only = default_pregen_only
        if not default_pregen_only:
            if custom_pregen_dir is None:
                self.custom_pregenerated_path = os.path.join(config.CUSTOM_RANDOM_ACTIVITIES_DIR, self.network_name)
            else:
                 self.custom_pregenerated_path = os.path.join(custom_pregen_dir, self.network_name)
            match = self.network_check_for_pregeneration(self.custom_pregenerated_path)
            if not match:
                print("Warning: Provided custom pregenerated directory either could not be find or does not match the network hash (i.e. was created using a different network). Please either set use_default_pregen_only to True or fix directory")
                self.custom_pregenerated_path = None
                self.default_pregen_only = True
            
            if (not match and not os.path.exists(self.pregenerated_experiments_path)) and self.use_pregen_data:
                raise ValueError('Could not find pregenerated random activities (either default or in the provided custom directory. Please pregenerate activities with the pregenerate module or set use_pregen = False.)')
            else:
                self.network_check = True
        else:
            self.custom_pregenerated_path = None
            if  not os.path.exists(self.pregenerated_experiments_path) and self.use_pregen_data:
                raise ValueError(f'Could not find pregenerated random activities in expected network directory ({self.pregenerated_experiments_path}). Please pregenerate activities with the pregenerate module or set use_pregen = False.')
            else:
                self.network_check = True

        self.compendia_paths = {
            ('Y', True): 'compendia=0_30_70',
            ('Y', False): 'compendia=0_50_50',
            ('ST', None): 'compendia=0_30_70'
        }
        
        # end of new fields for pregenerated_random

        self.aggregate = 'mean'
        self.threshold = None
        self.greater = True

        self.run_date = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

        # if data columns is None, set data columns to be columns with data: in front
        self.set_data_columns(data_columns=data_columns)
        #calculate the compendia distriubtion for each column
        #self.get_compendia_distribution(selection_type = "KSTAR_NUM_COMPENDIA_CLASS")
        #get available pregenerated sizes
        self.pregenerated_sizes = self.check_file_sizes_for_pregenerated()

    def check_data_columns(self):
        """
        Checks data columns to make sure column is in evidence and that evidence filtered on that data column
        has at least one point of evidence. Removes all columns that do not meet criteria
        """
        new_data_columns = []
        evidence = self.evidence.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE])[self.data_columns].agg(
            self.aggregate).reset_index()
        for col in self.data_columns:
            if col in self.evidence.columns:
                if self.threshold is not None:
                    if self.greater:
                        if len(evidence[evidence[col] >= self.threshold]) > 0:
                            new_data_columns.append(col)
                        else:
                            print(f"{col} does not have any evidence, and will not be used.")
                            self.logger.warning(f"{col} does not have any evidence")
                    else:
                        if len(evidence[evidence[col] <= self.threshold]) > 0:
                            new_data_columns.append(col)
                        else:
                            print(f"{col} does not have any evidence, and will not be used.")
                            self.logger.warning(f"{col} does not have any evidence")
                else:
                    if ~evidence[col].isna().all():
                        new_data_columns.append(col)
                    else:
                        print(f"{col} does not have any evidence, and will not be used.")
                        self.logger.warning(f"{col} does not have any evidence")
            else:
                print(f"{col} not in evidence, and will not be used.")
                self.logger.warning(f"{col} not in evidence")

        #make sure there is at least one data column
        if len(new_data_columns) == 0:
            raise ValueError("ERROR: No valid data columns found after filtering evidence. Please check that data columns are in evidence and that evidence has at least one point of evidence per data column after filtering if threshold is set.")
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
            self.data_columns = data_columns
            #print(self.data_columns)
        self.check_data_columns()

    def test_threshold(self, threshold, agg='mean', greater=True, plot=False, return_evidence_sizes=False):
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
            whether to plot a histogram of the evidence sizes used
        return_site_nums: bool
            indicates whether to return the evidence sizes for all samples or not

        Returns
        -------
        Outputs the minimum, maximum, and median evidence sizes across all samples. May return evidence sizes of all samples as pandas series
        """
        evidence_binary = self.create_binary_evidence(agg=agg, threshold=threshold, greater=greater)
        num_sites = evidence_binary[self.data_columns].sum()
        print('Number of Sites That Will Be Used As Evidence For Each Sample:')
        print('Maximum =', num_sites.max())
        print('Minimum =', num_sites.min())
        print('Median =', num_sites.median())

        if plot:
            plt_data = num_sites.reset_index()
            sns.histplot(data=plt_data, x=0)
            plt.xlabel('Number of Sites Used As Evidence')

        if return_evidence_sizes:
            return num_sites



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
            default_values = helpers.parse_network_information(default_file_path)
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
    
    #def find_closest_compendia_distribution(self, high_study_fraction):


    #def check_file_sizes_for_pregenerated(self):
    #    """
    #    Check the sizes of pre-generated files for the given datasets.

    #    This function identifies the appropriate pre-generated file sizes based on the dataset size and returns a list of these sizes.

    #    Returns
    #    -------
    #    pregenerated_sizes_files : list
    #        List of sizes of pre-generated files that match the dataset size criteria.
    #    """

    #    pregenerated_sizes_files = []
    #    compendia_class_counts = self.get_compendia_distribution(self.data_columns)
    #    for dataset in self.data_columns:
    #        # Determine path key (phospho_type + class condition)
    #        key = (self.phospho_type, compendia_class_counts[dataset][2] > 60 if self.phospho_type == 'Y' else None)

    #        if key not in self.compendia_paths:
    #            raise ValueError(f"ERROR: Unrecognized phosphoType '{self.phospho_type}'. Must be 'Y' or 'ST'.")

    #        compendia_file_path = os.path.join(
    #            str(self.pregenerated_experiments_path), str(self.compendia_paths[key])
    #        )

    #        if os.path.exists(compendia_file_path):
    #            pregenerated_sizes = os.listdir(compendia_file_path)
    #        else:
    #            raise FileNotFoundError(f"Directory not found: {compendia_file_path}")
    #    return pregenerated_sizes

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
                file_sizes = os.listdir(os.path.join(self.custom_pregenerated_path, dist))
                #convert to int
                file_sizes = [int(size) for size in file_sizes]
                pregenerated_sizes_files[dist.split('=')[1]] = file_sizes
                #experiment sizes will be the directories within each compendia distribution
                if dist.split('=')[1] not in pregenerated_sizes_files:
                    pregenerated_sizes_files[dist.split('=')[1]] = file_sizes
                else:
                    pregenerated_sizes_files[dist.split('=')[1]].extend(file_sizes)

        return pregenerated_sizes_files

    
    def determine_if_pregen(self, size, high_study_bias_perc, pregenerated_sizes = None, custom_pregenerated_path = None):
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
        if pregenerated_sizes is None:
            pregenerated_sizes = self.check_file_sizes_for_pregenerated(custom_pregenerated_path=custom_pregenerated_path)

        if self.phospho_type == 'Y':
            pregen_compendia_dist = '0_30_70' if high_study_bias_perc > 60 else '0_50_50'
        else:
            pregen_compendia_dist = '0_30_70' 

        #get experiment sizes associated with compendia distribution
        available_sizes = pregenerated_sizes[pregen_compendia_dist]
        #make sure dataset is large enough to use pregenerated data
        large_enough = size >= self.min_dataset_size_for_pregenerated
        #make sure pregenerated data is comparable in size to dataset
        close_enough = any(
            abs(pregen_size - size) / size < self.max_diff_from_pregenerated for pregen_size in
            available_sizes)
        # make sure pregenerated files match network
        use_pregen = large_enough and close_enough
        return use_pregen
            

    def calculate_random_activities(self, num_random_experiments=150, use_pregenerated_random_activities=None, save_new_random_activities=None,  custom_pregenerated_path=None, save_random_experiments=None, PROCESSES=1):
        """
        Generate random experiments and calculate kinase activities.Either uses pre-generated activity lists or
        generates new random experiments based on the provided parameters.

        Parameters
        ----------
        logger : Logger object
            Logger to record the progress and any issues during the randomization pipeline.
        num_random_experiments : int, optional
            Number of random experiments to generate, by default 150.
        use_pregen_data : bool, optional
            Whether to use pre-generated data, by default None.
        save_new_precompute : bool, optional
            Whether to save new precomputed data, by default None.
        pregenerated_experiments_path : str, optional
            Path to the directory containing pre-generated experiments, by default None.
        directory_for_save_precompute : str, optional
            Directory to save new precomputed data, by default None.
        network_hash : str, optional
            Hash of the network used, by default None.
        save_random_experiments : bool, optional
            Whether to save the generated random experiments, by default None.
        PROCESSES : int, optional
            Number of processes to use for parallel computation, by default 1.

        Returns
        -------
        None
        """
        self.logger.info("Running Randomization Pipeline")
        if use_pregenerated_random_activities is not None:
            self.use_pregen_data = use_pregenerated_random_activities
        if save_new_random_activities is not None:
            self.save_new_precompute = save_new_random_activities
        if custom_pregenerated_path is not None:
            self.directory_for_save_precompute = custom_pregenerated_path


        if use_pregenerated_random_activities:
            # Classify datasets
            pregenerated_sizes = self.check_file_sizes_for_pregenerated()
            #network_check = self.network_check_for_pregeneration()
            if not self.network_check:
                self.logger.warning("Network used does not match any pre-generated networks. All datasets will be calculated from scratch.")
                print("Network used does not match any pre-generated networks. All datasets will be calculated from scratch.")
                self.data_columns_with_pregenerated = []
                self.data_columns_from_scratch = self.data_columns
            else:
                self.data_columns_with_pregenerated = []
                self.data_columns_from_scratch = []
                for dataset in self.data_columns:
                    high_study_bias_perc = int(self.compendia_distribution[dataset][2])
                    size = self.dataset_sizes[dataset]
                    use_pregen = self.determine_if_pregen(size = size, high_study_bias_perc=high_study_bias_perc, pregenerated_sizes = pregenerated_sizes)
                    if use_pregen:
                        self.data_columns_with_pregenerated.append(dataset)
                    else:
                        self.data_columns_from_scratch.append(dataset)

            # Process pre-generated data, one dataset at a time
            print('Loading pre-generated random activities for datasets where applicable...')
            self.load_pregenerated_random_activities()

            #report on which datasets will use pregenerated data
            if len(self.data_columns_from_scratch) > 0:
                print(f"{len(self.data_columns_from_scratch)} out of {len(self.data_columns)} columns did not have appropriate pregenerated random activities and will calculate from scratch: {', '.join(self.data_columns_from_scratch)}")

                self.logger.info(f"Generating random experiments for: {', '.join(self.data_columns_from_scratch)}")
                self.calculate_random_enrichment(num_random_experiments, selection_type='KSTAR_NUM_COMPENDIA_CLASS',
                    save_random_experiments = save_random_experiments, PROCESSES=PROCESSES
                )

            # Add pre-generated data to random enrichment
            self.add_pregenerated_to_random_enrichment()

        # Process all from-scratch data if use_pregen_data is False
        else:
            self.data_columns_from_scratch = self.data_columns
            self.logger.info(f"Generating random experiments for: {', '.join(self.data_columns)}")
            self.calculate_random_enrichment(num_random_experiments, selection_type='KSTAR_NUM_COMPENDIA_CLASS',
                save_random_experiments = save_random_experiments, PROCESSES=PROCESSES
            )


        # Add pre-generated data to random enrichment
        #self.add_pregenerated_to_random_enrichment()


    def calculate_random_enrichment(self, num_random_experiments, selection_type='KSTAR_NUM_COMPENDIA_CLASS', save_random_experiments=False, PROCESSES=1):
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
        PROCESSES : int, optional
            Number of processes to use for parallel computation, by default 1.

        Returns
        -------
        None
        """
        if self.use_pregen_data:
            num_random_experiments = 150
        else:
            num_random_experiments = num_random_experiments
        """
        Change to calculate_random_enrichment to immediately calculate the kinase activities for each random experiment and discard the resulting random experiment.
        """
        self.logger.info("Running Randomization Pipeline")
        self.num_random_experiments = num_random_experiments
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


                if self.save_new_precompute:
                    activities_list_df = pd.concat(activities_list).reset_index(drop=True)
                    self.save_new_precomputed_random_enrichment(activities_list_df, col)
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
            for col in tqdm(self.data_columns_from_scratch, desc="Calculating activities from random experiments for each dataset not using pregenerated random activities"):
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
                if self.save_new_precompute:
                    activities_list_df = pd.concat(activities_list).reset_index(drop=True)
                    self.save_new_precomputed_random_enrichment(activities_list_df, col)
                combined_activities_list.extend(activities_list)

                # combine all random activities
            self.random_enrichment = pd.concat(combined_activities_list).reset_index(drop=True)

            # reformat random experiments into single matrix, if wanting to save
            if save_random_experiments:
                all_rand_experiments = pd.concat(all_rand_experiments)
                all_rand_experiments['weight'] = 1
                self.random_experiments = all_rand_experiments.pivot(index=[config.KSTAR_ACCESSION, config.KSTAR_SITE],columns='Experiment',
                    values='weight').reset_index()



    def load_pregenerated_random_activities(self):
        """
        Load pre-generated random activities for the given datasets.

        This function processes datasets that have pre-generated random experiments. It identifies the appropriate
        pre-generated file based on the size of the dataset and appends the activities to the provided list.

        Parameters
        ----------
        with_pregenerated_evidence : pandas.DataFrame
            DataFrame containing the evidence for the datasets with pre-generated random experiments.
        with_pregenerated : list
            List of dataset names that have pre-generated random experiments.
        random_activities_list : list
            List to which the concatenated activities of each dataset will be appended.

        Returns
        -------
        None
        """
        if len(self.data_columns_with_pregenerated) > 0:
            #initialize lists
            self.num_random_experiments = 150
            pregen_activities_list = []

            for dataset in self.data_columns_with_pregenerated:
                with_pregenerated_evidence = self.evidence_binary[['KSTAR_ACCESSION', 'KSTAR_SITE','KSTAR_NUM_COMPENDIA_CLASS', 'KSTAR_NUM_COMPENDIA', dataset]]
                compendia_class_counts = self.compendia_distribution[dataset]
                key = (self.phospho_type, compendia_class_counts[2] > 60 if self.phospho_type == 'Y' else None)

                if key not in self.compendia_paths:
                    raise ValueError(f"ERROR: Unrecognized phosphoType '{self.phospho_type}'. Must be 'Y' or 'ST'.")

                compendia_file_path = os.path.join(
                    str(self.pregenerated_experiments_path), str(self.compendia_paths[key])
                )

                if not os.path.exists(compendia_file_path):
                    raise FileNotFoundError(f"Directory not found: {compendia_file_path}")

                pregenerated_files = os.listdir(compendia_file_path)

                size = self.dataset_sizes[dataset]
                closest_size = min(pregenerated_files, key=lambda x: abs(x - size))
                matched_file_path = os.path.join(compendia_file_path, str(closest_size), 'random_enrichment.tsv')
                rand_dataset_activities = load_random_activities(matched_file_path)

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

    def save_new_precomputed_random_enrichment(self, activities_list_df, col):
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
        run_info_content = self.get_run_information_content()
        size = self.evidence_binary[col].sum()

        # Determine the compendia directory based on phospho_type and class counts
        compendia_class_counts = self.get_compendia_distribution(self.evidence_binary, [col])
        compendia_paths = {
            ('Y', True): 'compendia_0_0_1_30_2_70',
            ('Y', False): 'compendia_0_5_1_50_2_45',
            ('ST', None): 'compendia_0_0_1_30_2_70'
        }

        key = (self.phospho_type, compendia_class_counts[col][2] > 60 if self.phospho_type == 'Y' else None)
        if key not in compendia_paths:
            raise ValueError("Invalid phospho_type or compendia class counts")

        compendia_file_path = os.path.join(str(self.directory_for_save_precompute), str(self.phospho_type), str(self.network_hash),
                                           str(compendia_paths[key]))
        os.makedirs(compendia_file_path, exist_ok=True)

        # Save activities to file
        file_name = f"{size}.tsv"
        file_path = os.path.join(compendia_file_path, file_name)

        #if path does not exist, create it
        if not os.path.exists(file_path):
            os.makedirs(os.path.dirname(file_path))

        pd.DataFrame(activities_list_df).to_csv(file_path, index=True, sep='\t')
        self.logger.info(f"New precomputed random enrichment data for {col} saved to {file_path}")

        # Save run information content
        run_info_save_path = os.path.join(self.directory_for_save_precompute, self.phospho_type, self.network_hash,
                                          "RUN_INFORMATION.txt")
        with open(run_info_save_path, 'w') as file:
            file.write(run_info_content)
        self.logger.info(f"RUN_INFORMATION.txt saved to {run_info_save_path}")

    def get_run_information_content(self):
        """
        Retrieve network information from RUN_INFORMATION.txt based on phospho_type.

        Reads the RUN_INFORMATION.txt file from the appropriate network directory based on
        the phospho_type ('Y' or 'ST'). The file contains network configuration details
        including unique ID, date, network specifications, and compendia counts.

        Returns
        -------
        str
            Contents of RUN_INFORMATION.txt if found.
            'RUN_INFORMATION.txt file not found.' if the file doesn't exist.

        Raises
        ------
        ValueError
            If phospho_type is not 'Y' or 'ST'.
        """
        base_path = config.NETWORK_DIR
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
        network_directory = os.path.join(config.NETWORK_SUBDIR[self.phospho_type], 'INDIVIDUAL_NETWORKS/')

        for filename in os.listdir(network_directory):
            net_num = filename.split('_')[1]
            if filename.endswith('.tsv'):
                file_path = os.path.join(network_directory, filename)
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

    def create_binary_evidence(self, agg='mean', threshold=1.0, evidence_size=None, greater=True):
        """
        Returns a binary evidence data frame according to the parameters passed in for method for aggregating
        duplicates and considering whether a site is included as evidence or not

        Parameters
        ----------
        threshold : float
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

        Returns
        -------
        evidence_binary : pd.DataFrame
            Matches the evidence dataframe of the kinact object, but with 0 or 1 if a site is included or not.
            This is uniquified and rows that are never used are removed.


        """
        #make sure parameters are valid
        if not isinstance(threshold, (int, float)):
            raise ValueError("Threshold must be an integer or float.")
        if evidence_size is not None and not isinstance(evidence_size, int):
            raise ValueError("evidence_size must be an integer or None.")
        elif evidence_size is not None and evidence_size <= 0:
            raise ValueError("evidence_size must be a positive integer.")

        if agg not in ['mean', 'min','max', 'median']:
            raise ValueError("Aggregation method must be one of 'mean', 'min', 'max', or 'median'.")
        

        self.threshold = threshold
        self.evidence_size = evidence_size
        self.greater = greater
        self.check_data_columns()

        # collapse sites into single row based on agg parameter
        evidence = self.evidence.groupby([config.KSTAR_ACCESSION, config.KSTAR_SITE])[self.data_columns].agg(
            agg).reset_index()

        # set the binary evidence for whether a site is included
        evidence_binary = evidence.copy()
        if evidence_size is None:
            for col in self.data_columns:
                if greater:
                    evidence_binary[col] = (evidence_binary[col] >= threshold).astype(int)
                else:
                    evidence_binary[col] = (evidence_binary[col] <= threshold).astype(int)
        else:
            for col in self.data_columns:
                # check how many non_nan sites there (if less than N, set n to be equal to number of sites available)
                num_sites_available = evidence_binary.dropna().shape[0]
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

    def calculate_kinase_activities(self, agg='mean', threshold=1.0, evidence_size=None, greater=True, PROCESSES=1):
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
                                                           evidence_size=self.evidence_size, greater=self.greater)
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

    def calculate_Mann_Whitney_activities_sig(self, PROCESSES=1):
        """
        For a kinact_dict, where random generation and activity has already been run for the phospho_types of interest,
        this will calculate the Mann-Whitney U test for comparing the array of p-values for real data
        to those of random data, across the number of networks used.
        It will also calculate the false positive rate for a pvalue, given observations of a random bootstrapping analysis

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

        Returns
        -------

        """
        if not isinstance(self.random_enrichment, pd.DataFrame):
            raise ValueError("Random activities do not exist, please run kstar_activity.randomized_analysis")

        number_sig_trials = self.num_random_experiments - 1

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
        for exp in tqdm(self.data_columns, desc='Calculating final activities with the mann whitney U test'):
            self.logger.info("MW Working on %s: " % (exp))
            pval_arr = []
            fpr_arr = []
            if PROCESSES > 1:
                #setup partial function for unchanging parameters
                partial_func = partial(calculate_MannWhitney_one_experiment_one_kinase, real_grouped, random_grouped, exp)
                #iterate through kinases in parallel
                with concurrent.futures.ProcessPoolExecutor(max_workers=PROCESSES) as executor:
                    for pval, fpr in executor.map(partial_func, self.kinases):
                        pval_arr.append(pval)
                        fpr_arr.append(fpr)
            else:
                for kinase in self.kinases:
                    pval, fpr = calculate_MannWhitney_one_experiment_one_kinase(real_grouped, random_grouped,exp, kinase)
                    pval_arr.append(pval)
                    fpr_arr.append(fpr)

            self.activities_mann_whitney[exp] = pval_arr
            self.fpr_mann_whitney[exp] = fpr_arr

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

def load_random_activities(file):
    rand_dataset_activities = pd.read_csv(file, delimiter='\t')
    rand_dataset_activities = rand_dataset_activities.melt(id_vars=['network', 'KSTAR_KINASE'], var_name='data', value_name='kinase_activity')
    return rand_dataset_activities


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


def calculate_MannWhitney_one_experiment_one_kinase(real_activities_grouped, rand_activities_grouped, experiment, kinase):
    """
    For a given kinact object, where random generation and activity has already been run, this will calculate the Mann-Whitney U test between the p-values across all networks for the given experiment name and from the random networks. It will also calculate the significance value for the given test based on the target_alpha value by using each random set as a real set to bootstrap.

    Parameters
    ----------
    real_activities_grouped : pandas.Series
        Multi-index series containing kinase activities for a given experiment grouped by kinase. Should be indexed by kinase, with each entry being a list of activities across all networks. You can obtain this by performing a groupby on the real enrichment DataFrame.
    rand_activities_grouped : pandas.Series
        Multi-index series object containing kinase activities for random experiments associated with the same sample as kinact_activities_sub, grouped by kinase and random experiment number. Should be indexed by kinase and rand experiment number, with each entry being a list of activities across all networks. You can obtain this by performing a groupby on the random enrichment DataFrame.
    kinase : str
        Kinase name to measure significance for.
    experiment : str
        Experiment name to measure significance for.

    Returns
    -------
    p_value : float
        p-value that results from Mann-Whitney U test.
    fpr_value : float
        The false positive rate where the p_value for the real experiment lies, given the random experiments.
    """
    kinase_activity_list = real_activities_grouped.loc[experiment,kinase]
    #grab 
    random_kinase_activity_array = np.vstack(rand_activities_grouped.loc[experiment, kinase])

    #remove one random experiment to make real and random Mann Whitney tests comparable (same number of random experiments used for comparison)
    i = np.random.randint(0, random_kinase_activity_array.shape[0])
    bgnd = np.delete(random_kinase_activity_array, i, 0)

    #compare real enrichment to random background (size = num_random_experiments - 1)
    [stat, p_value] = stats.mannwhitneyu(-np.log10(kinase_activity_list),
        -np.log10(np.concatenate(bgnd)),alternative='greater')
    
    # Calculate FPR using the helper function
    randomStats = calculate_fpr_Mann_Whitney(random_kinase_activity_array)
    fpr_value = calculate_fpr.single_pvalue_fpr(randomStats, p_value)
    return p_value, fpr_value


"""
****************************************
Methods for running KSTAR pipeline
****************************************
"""


def enrichment_analysis(experiment, odir, name='experiment', phospho_types=['Y', 'ST'], network_dir = None, data_columns=None, agg='mean',
                        threshold=1.0, evidence_size=None, greater=True, PROCESSES=1, logger = None, min_dataset_size_for_pregenerated=150, max_diff_from_pregenerated=0.20):
    """
    Function to establish a kstar KinaseActivity object from an experiment with an activity log
    add the networks, calculate, aggregate, and summarize the hypergeometric enrichment into a final activity object. Should be followed by
    randomized_analyis, then Mann_Whitney_analysis.

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
    #make sure odir exists
    if not os.path.exists(odir):
        raise ValueError(f"Output directory {odir} does not exist. Please create it before running enrichment_analysis.")
    kinact_dict = {}
    # For each phosphoType of interest, establish a kinase activity object on a filtered dataset and run, aggregate, and summarize activity
    if not isinstance(phospho_types, list):
        raise ValueError("phospho_types must be a list containing 'Y' and/or 'ST' (['Y', 'ST'], ['Y'], or ['ST'])")
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
            
        kinact = KinaseActivity(experiment_sub, odir = odir, name = name, data_columns=data_columns, phospho_type=phospho_type, network_dir=network_dir, logger = logger, min_dataset_size_for_pregenerated=min_dataset_size_for_pregenerated, max_diff_from_pregenerated=max_diff_from_pregenerated,)

        kinact.calculate_kinase_activities(agg=agg, threshold=threshold, evidence_size=evidence_size, greater=greater, 
                                           PROCESSES=PROCESSES)
        #kinact.aggregate_activities()
        #kinact.activities = kinact.summarize_activities()
        kinact_dict[phospho_type] = kinact
    return kinact_dict

def randomized_analysis(kinact_dict, num_random_experiments=150, use_pregen_data = None,
                        save_new_random_activities = None, pregenerated_experiments_path = None,
                        custom_pregenerated_dir = None, save_random_experiments = None, PROCESSES = 1):
    """
    Perform randomized analysis on kinase activity data.

    Parameters
    ----------
    kinact_dict : dict
        Dictionary containing kinase activity data.
    log : Logger object
        Logger to record the progress and any issues during the randomization pipeline.
    num_random_experiments : int, optional
        Number of random experiments to generate, by default 150.
    use_pregen_data : bool, optional
        Whether to use pre-generated data, by default False.
    save_new_precompute : bool, optional
        Whether to save new precomputed data, by default None.
    pregenerated_experiments_path : str, optional
        Path to the directory containing pre-generated experiments, by default None.
    directory_for_save_precompute : str, optional
        Directory to save new precomputed data, by default None.
    network_hash : str, optional
        Hash of the network used, by default None.
    save_random_experiments : bool, optional
        Whether to save the generated random experiments, by default None.
    PROCESSES : int, optional
        Number of processes to use for parallel computation, by default 1.

    Returns
    -------
    None
    """
    # Validate input parameters
    if not isinstance(use_pregen_data, bool) and use_pregen_data is not None:
        raise ValueError("use_pregen_data must be True or False")
    if not isinstance(save_new_random_activities, bool) and save_new_random_activities is not None:
        raise ValueError("save_new_random_activities must be True or False")
    if not isinstance(pregenerated_experiments_path, str) and pregenerated_experiments_path is not None:
        raise ValueError("pregenerated_experiments_path must be a string")
    if save_new_random_activities:
        #will want to chang
        if not isinstance(custom_pregenerated_dir, str) and custom_pregenerated_dir is not None:
            raise ValueError("When save_new_precompute is True, directory_for_save_precompute must be provided as a non-empty string")



    for phospho_type, kinact in kinact_dict.items():
        #if network_hash is None:
        #    if phospho_type == 'Y':
        #        network_hash = config.NETWORK_HASH_Y
        #    elif phospho_type == 'ST':
        #        network_hash = config.NETWORK_HASH_ST
        #if not re.fullmatch(r'[a-fA-F0-9]{64}', network_hash):
        #    raise ValueError("network_hash must be a valid SHA-256 hash")
    
        kinact.calculate_random_activities(num_random_experiments, use_pregenerated_random_activities=use_pregen_data,save_new_random_activities = save_new_random_activities, custom_pregenerated_path=custom_pregenerated_dir,
                                           save_random_experiments=save_random_experiments, PROCESSES=PROCESSES)
        # Ensure `random_experiments` is stored in `kinact_dict`
        #kinact_dict[phospho_type].random_experiments = kinact.random_experiments

def Mann_Whitney_analysis(kinact_dict, PROCESSES=1):
    """
    For a kinact_dict, where random generation and activity has already been run for the phospho_types of interest,
    this will calculate the Mann-Whitney U test for comparing the array of p-values for real data
    to those of random data, across the number of networks used.
    It will also calculate the false positive rate for a pvalue, given observations of a random bootstrapping analysis

    Parameters
    ----------
    kinact_dict: dictionary
        A dictionary of kinact objects, with keys 'Y' and/or 'ST'
    log: logger
        Logger for logging activity messages
    number_sig_trials: int
        Maximum number of significant trials to run
    """

    for phospho_type, kinact in kinact_dict.items():
        kinact.calculate_Mann_Whitney_activities_sig(PROCESSES=PROCESSES)

def run_kstar_analysis(experiment, odir, name='experiment', phospho_types=['Y', 'ST'], network_dir = None, data_columns=None, agg='mean',
                        threshold=1.0, evidence_size=None, greater=True, num_random_experiments = 100, use_pregen_data = None,
                        save_new_precompute = None, pregenerated_experiments_path = None,
                        custom_pregen_dir = None, save_random_experiments = None, logger = None, mapped = True, map_dict = None, save = True, PROCESSES=1):
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
    num_random_experiments: int
        Number of random experiments to run
    use_pregen_data: bool
        Whether to use pre-generated data
    save_new_precompute: bool
        Whether to save new precomputed data
    pregenerated_experiments_path: string
        Path to pregenerated experiments
    directory_for_save_precompute: string
        Directory to save precomputed data
    save_random_experiments: bool
        Whether to save random experiments
    logger: logger
        Logger for logging messages
    PROCESSES: int
        Number of processes to use

    """
    if not mapped:
        if map_dict is None:
            raise ValueError("map_dict must be provided if mapped is False. This should indicate the name of accession and site/peptide columns ({'accession': 'accession_column_name', 'site': 'site_column_name', 'peptide':'peptide_column_name'})")
        print('Mapping experiment to reference phosphoprotome...')
        exp_mapper = mapping.ExperimentMapper(experiment, map_dict, odir = odir, name = name)
        if save:
            print(f'Saving mapped experiment and relevant stats in {odir}/MAPPED_DATA/')
            exp_mapper.save_experiment()
        experiment = exp_mapper.experiment
    
    #start enrichment analysis
    print('Starting kinase-substrate enrichment analysis...')
    kinact_dict = enrichment_analysis(experiment, odir = odir, name = name, phospho_types = phospho_types, network_dir = network_dir, threshold = threshold, agg = agg, data_columns=data_columns, greater = greater, logger = logger, evidence_size=evidence_size,PROCESSES = PROCESSES)
    #
    print('Starting calculation of random activities...')
    randomized_analysis(kinact_dict, num_random_experiments=num_random_experiments, use_pregen_data=use_pregen_data, save_new_precompute=save_new_precompute, pregenerated_experiments_path=pregenerated_experiments_path,custom_pregen_dir=custom_pregen_dir, save_random_experiments=save_random_experiments, PROCESSES=PROCESSES)
    #
    print('Comparing kinase-substrate enrichment from the real experiment to random experiments...')
    Mann_Whitney_analysis(kinact_dict, PROCESSES = PROCESSES)

    if save:
        print(f'Saving KSTAR analysis results in {odir}/RESULTS/')
        save_kstar_slim(kinact_dict, name, odir)

    print('Done.')

    
    


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

    for phospho_type in kinact_dict:
        kinact = kinact_dict[phospho_type]
        name_out = f"{name}_{phospho_type}"
        kinact.real_enrichment.to_csv(f"{odir}/RESULTS/{name_out}_real_enrichment.tsv", sep='\t')
        kinact.evidence_binary.to_csv(f"{odir}/RESULTS/{name_out}_binarized_experiment.tsv", sep='\t', index=False)
        if hasattr(kinact, 'random_experiments') and kinact.random_experiments is not None:
            kinact.random_experiments.to_csv(f"{odir}/RESULTS/{name_out}_random_experiments.tsv", sep='\t', index=False)
        if hasattr(kinact, 'random_enrichment'):
            kinact.random_enrichment.to_csv(f"{odir}/RESULTS/{name_out}_random_enrichment.tsv", sep='\t')
        if hasattr(kinact, 'activities_mann_whitney'):
            kinact.activities_mann_whitney.to_csv(f"{odir}/RESULTS/{name_out}_mann_whitney_activities.tsv", sep='\t',
                                                  index=True)
            kinact.fpr_mann_whitney.to_csv(f"{odir}/RESULTS/{name_out}_mann_whitney_fpr.tsv", sep='\t', index=True)

    if PICKLE:
        pickle.dump(kinact_dict, open(f"{odir}/RESULTS/{name}_kinact.p", "wb"))


def save_kstar_slim(kinact_dict, name, odir, param_format = 'pickle'):
    """
    Having performed kinase activities (run_kstar_analyis), save each of the important dataframes, minimizing the memory storage needed to get back
    to a rebuilt version for plotting results and analysis. For each phospho_type in the kinact_dict, this will save three .tsv files for every activities
    analysis run, two additional if random analysis was run, and two more if Mann Whitney based analysis was run. It also creates a readme file of the parameter values
    used

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

    param_dict = {}

    for phospho_type in kinact_dict:

        kinact = kinact_dict[phospho_type]
        name_out = f"{name}_{phospho_type}"

        param_temp = {}
        #iterate through attributes and save those needed to reinstantiate object
        for attr in vars(kinact):
            #check for basic types only
            if isinstance(getattr(kinact, attr), (int, float, str, bool, type(None))):
                param_temp[attr] = getattr(kinact, attr)

       # param_temp['run_date'] = kinact.run_date
       # param_temp['network_directory'] = kinact.network_directory
       # param_temp['network_name'] = kinact.network_name
       # param_temp['network_hash'] = kinact.network_hash
       # param_temp['num_networks'] = kinact.num_networks
       # param_temp['data_columns'] = kinact.data_columns
       # param_temp['threshold'] = kinact.threshold
       # param_temp['aggregate'] = kinact.aggregate
       # param_temp['greater'] = kinact.greater
       # param_temp['num_random_experiments'] = kinact.num_random_experiments
       # param_temp['use_pregen_data'] = kinact.use_pregen_data
        param_temp['mann_whitney'] = False
       # param_temp['eataset_sizes'] = kinact.dataset_sizes
       # param_temp['compendia_distribution'] = kinact.compendia_distribution



        #kinact.real_enrichment.to_csv(f"{odir}/RESULTS/{name_out}_real_enrichment.tsv", sep='\t')
        # kinact.activities.to_csv(f"{odir}/RESULTS/{name_out}_activities.tsv", sep = '\t', index = True)
        kinact.evidence_binary.to_csv(f"{odir}/RESULTS/{name_out}_binarized_experiment.tsv", sep='\t', index=False)

        #if hasattr(kinact, 'random_enrichment'):
        #    param_temp['randomized'] = True
        #    param_temp['num_random_experiments'] = kinact.num_random_experiments
        #    kinact.random_enrichment.to_csv(f"{odir}/RESULTS/{name_out}_random_enrichment.tsv", sep='\t')
        #    # if hasattr(kinact, 'random_experiments'):

        if hasattr(kinact, 'activities_mann_whitney'):
            param_temp['mann_whitney'] = True
            kinact.activities_mann_whitney.to_csv(f"{odir}/RESULTS/{name_out}_mann_whitney_activities.tsv", sep='\t',
                                                  index=True)
            kinact.fpr_mann_whitney.to_csv(f"{odir}/RESULTS/{name_out}_mann_whitney_fpr.tsv", sep='\t', index=True)

        param_dict[phospho_type] = param_temp

    # save the parameters in a pickle file for reinstantiating object information
    if param_format == 'pickle':
        pickle.dump(param_dict, open(f"{odir}/RESULTS/{name}_params.p", "wb"))
    elif param_format == 'json':
        with open(f"{odir}/RESULTS/{name}_params.json", 'w') as json_file:
            json.dump(param_dict, json_file, indent=4)
    else:
        raise ValueError("param_format must be either 'pickle' or 'json'")


def from_kstar_slim(name, odir):
    """
    Given the name and output directory of a saved kstar analyis, load the parameters and minimum dataframes needed for reinstantiating a kinact object
    This minimum list will allow you to repeat normalization or mann whitney at a different false positive rate threshold and plot results.

    Parameters
    ----------
    name: string
        The name to used when saving activities and mapped data
    odir:  string
        Output directory of saved files and parameter pickle
    log: logger
        Logger for logging activity messages
    """

    # First check for the param file
    if os.path.exists(f"{odir}/RESULTS/{name}_params.json"):
        with open(f"{odir}/RESULTS/{name}_params.json", 'r') as json_file:
            param_dict = json.load(json_file)
    elif os.path.exists(f"{odir}/RESULTS/{name}_params.p"):
        param_dict = pickle.load(open(f"{odir}/RESULTS/{name}_params.p", "rb"))
    else:
        print(f"ERROR: Cannot find parameter dictionary file in RESULTS: {odir}/RESULTS/{name}_params.p")
        return
    
    kinact_dict = {}
    for phospho_type in param_dict.keys():
        params = param_dict[phospho_type]
        name_out = f"{name}_{phospho_type}"


        # check that the minimum file set exists so we can use binary_evidence file as the experiment
        evidence_binary = pd.read_csv(f"{odir}/RESULTS/{name_out}_binarized_experiment.tsv", sep='\t')
        #grab additional parameters needed to reinstate object
        network_dir = params.get('network_directory', None)
        network_name = params.get('network_name', None)


        #load activity logger
        kinact = KinaseActivity(evidence_binary, odir=odir, name=name, phospho_type=phospho_type)

        #kinact.real_enrichment = pd.read_csv(f"{odir}/RESULTS/{name_out}_real_enrichment.tsv", sep='\t', index_col=0)
        kinact.evidence_binary = evidence_binary

        if 'mann_whitney' in params.keys():
            # read mann_whitney and load
            kinact.activities_mann_whitney = pd.read_csv(f"{odir}/RESULTS/{name_out}_mann_whitney_activities.tsv",
                                                         sep='\t', index_col=config.KSTAR_KINASE)
            kinact.fpr_mann_whitney = pd.read_csv(f"{odir}/RESULTS/{name_out}_mann_whitney_fpr.tsv", sep='\t',
                                                  index_col=config.KSTAR_KINASE)
            params.pop('mann_whitney', None)


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

