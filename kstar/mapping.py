
import pandas as pd
from collections import defaultdict
import re
import os
import numpy as np
from kstar import config, helpers



class ExperimentMapper:
    """
    Given an experiment object and reference sequences, map the phosphorylation sites to the common reference.
    Inputs

    Parameters
    ----------
    name : str
        Name of experiment. Used for logging
    experiment: pandas dataframe
        Pandas dataframe of an experiment that has a reference accession, a peptide column and/or a site column. 
        The peptide column should be upper case, with lower case indicating the site of phosphorylation - this is preferred
        The site column should be in the format S/T/Y<pos>, e.g. Y15 or S345
    columns: dict
        Dictionary with mappings of the experiment dataframe column names for the required names 'accession_id', 'peptide', or 'site'. 
        One of 'peptide' or 'site' is required. 
    logger: Logger object
        used for logging when peptides cannot be matched and when a site location changes
    sequences: dict
        Dictionary of sequences. Key : accession. Value : protein sequence. 
        Default is imported from kstar.config
    compendia: pd.DataFrame
        Human phosphoproteome compendia, mapped to KinPred and annotated with number of compendia. 
        Default is imported from kstar.config
    window : int
        The length of amino acids to the N- and C-terminal sides of the central phosphoprotein to map a site to.
        Default is 7.
    data_columns: list, or empty
        The list of data columns to use. If this is empty, logger will look for anything that starts with statement data: and those values
        Default is None.

    
    Attributes
    ----------
    experiment: pandas dataframe
        mapped experiment, which for each peptide, no contains the mapped accession, site, peptide, number of compendia, compendia type
    sequences: dict
        Dictionary of sequences passed into the class
    compendia: pandas dataframe
        compendia dataframe passed into the class
    data_columns: list
        indicates which columns will be used as data
    

    """ 
    #for documentation purposes convert non required parameters to kwargs
    #def __init__(self, experiment, columns, logger, *kwargs): 
    def __init__(self, experiment, columns, logger, sequences=config.HUMAN_REF_SEQUENCES, compendia=config.HUMAN_REF_COMPENDIA, window = 7, data_columns = None): 
        self.experiment = experiment
        self.sequences = sequences
        self.compendia = compendia
        def set_accession_id(accession):
            acc = accession.split('-')
            if len(acc) > 1:
                acc = acc[:-1]
            return '-'.join(acc)

        if 'accession_id' not in columns.keys():
            raise ValueError('ExperimentMapper requires accession_id as a dictionary key')
        else:
            self.experiment[config.KSTAR_ACCESSION] = self.experiment[columns['accession_id']].apply(set_accession_id) 

        if 'peptide' not in columns.keys() and 'site' not in columns.keys():
            raise ValueError('ExperimentMapper requires either site or peptide as keys in dictionary')

        self.experiment[config.KSTAR_PEPTIDE] = self.experiment[columns['peptide']] if 'peptide' in columns.keys() else None
        self.experiment[config.KSTAR_SITE] = self.experiment[columns['site']] if 'site' in columns.keys() else None

        self.logger = logger
        self.set_data_columns(data_columns)
        self.align_sites(window)

        compendia = self.compendia[[config.KSTAR_ACCESSION, config.KSTAR_SITE, 'KSTAR_NUM_COMPENDIA', 'KSTAR_NUM_COMPENDIA_CLASS']]
        self.experiment = pd.merge(self.experiment, compendia, how = 'inner', on = [config.KSTAR_ACCESSION, config.KSTAR_SITE] )
        #after merge, NUM_COMPENDIA/CLASS as become a float, likely due to mismatches, so assume those evidences are 0
        self.experiment['KSTAR_NUM_COMPENDIA'] = self.experiment['KSTAR_NUM_COMPENDIA'].fillna(0.0).astype(int)
        self.experiment['KSTAR_NUM_COMPENDIA_CLASS'] = self.experiment['KSTAR_NUM_COMPENDIA_CLASS'].fillna(0.0).astype(int)

    def set_data_columns(self, data_columns):
        """
        Identifies which columns in the experiment should be used as data columns. If data_columns is provided,
        then 'data:' is added to the front and experiment dataframe is renamed. Otherwise, function will look for columns
        with 'data:' in front and this to the data_columns attribute.
        """
        self.data_columns = []
        if data_columns is None:
            for col in self.experiment.columns:
                if col.startswith('data:'):
                    self.data_columns.append(col)
        else:
            data_rename = defaultdict()
            for col in data_columns:
                if not col.startswith('data:'):
                    data_rename[col] = f"data:{col}"
                else: 
                    data_rename[col] = col
            self.experiment.rename(columns = data_rename, inplace = True)
            self.data_columns = list(data_rename.values())

    def get_experiment(self):
        """
        Return the mapped experiment dataframe
        """
        return self.experiment


    def get_sequence(self, accession):
        """
        Gets the sequence that matches the given accession
        """
        if accession in self.sequences.keys():
            return self.sequences[accession]
        return None

    
    def align_sites(self, window = 7):
        """
        Map the peptide/sites to the common sequence reference and remove and report errors for sites that do not align as expected.
        expMapper.align_sites(window=7). Operates on the experiment dataframe of class.

        Parameters
        ----------
        window: int
            The length of amino acids to the N- and C-terminal sides of the central phosphoprotein to map a site to.

        """

        self.experiment = expand_peptide(self.experiment, config.KSTAR_PEPTIDE)

        for index, row in self.experiment.iterrows():
            sequence = self.get_sequence(row[config.KSTAR_ACCESSION])
            if sequence is not None:
                # If peptide provided then find site
                if row[config.KSTAR_PEPTIDE] is not None:
                    
                    site =  peptide_site_number(peptide = row[config.KSTAR_PEPTIDE], 
                                                site = row[config.KSTAR_SITE],
                                                sequence = sequence)
                    if site is None:
                        self.logger.warning(f"SITE NOT FOUND : {row[config.KSTAR_ACCESSION]}\t{row[config.KSTAR_PEPTIDE]}")
                        self.experiment.loc[index, config.KSTAR_SITE] = None
                        self.experiment.loc[index, config.KSTAR_PEPTIDE] = None
                        continue
                    elif site != row[config.KSTAR_SITE]:
                        self.logger.info(f"SITE CHANGED : {row[config.KSTAR_ACCESSION]} {row[config.KSTAR_SITE]} -> {site}")
                        self.experiment.loc[index, config.KSTAR_SITE] = site
                    
                    peptide = get_aligned_peptide(site = site,
                                                  sequence = sequence, 
                                                  window = window)
                    if peptide is not None:
                        self.experiment.loc[index, config.KSTAR_PEPTIDE] = peptide
                        
                # Peptide not provided but site is provided - build peptide
                elif row[config.KSTAR_SITE] is not None:
                    peptide = get_aligned_peptide(site = row[config.KSTAR_SITE],
                                                  sequence = sequence, 
                                                  window = window)
                    if peptide is not None:
                        self.experiment.loc[index, config.KSTAR_PEPTIDE] = peptide
                    else:
                        self.logger.warning(f"PEPTIDE MISMATCH : {row[config.KSTAR_ACCESSION]} {row[config.KSTAR_SITE]}")
                        self.experiment.loc[index, config.KSTAR_SITE] = None
                        self.experiment.loc[index, config.KSTAR_PEPTIDE] = None
                else:
                    self.logger.warning(f"SITE NOT FOUND : {self.experiment[config.KSTAR_ACCESSION]}\t{self.experiment[config.KSTAR_SITE]}")
            
            # Sequence not found in compendia
            else:
                self.logger.warning(f"SEQUENCE NOT FOUND : {row[config.KSTAR_ACCESSION]}")
                self.experiment.loc[index, config.KSTAR_SITE] = None
                self.experiment.loc[index, config.KSTAR_PEPTIDE] = None
        self.experiment.dropna(axis = 'rows', subset = [config.KSTAR_ACCESSION, config.KSTAR_SITE, config.KSTAR_PEPTIDE], inplace = True)
        self.experiment.drop_duplicates(inplace=True)


def expand_peptide(df, peptide_column):
    """
    Expands a dataframe based on the peptide column such that each 
    peptide contains one modification type. 
    Example
        Accession  Peptide  Info
        P0001      aBCDeF   Good
        P0002      TUvWXy   Bad

        Becomes
        Accession  Peptide  Info
        P0001      aBCDEF   Good
        P0001      ABCDeF   Good
        P0002      TUvWXY   Bad
        P0002      TUVWXy   Bad
    
    Parameters
    --------
    df : pandas df
        Experiment pandas dataframe
    peptide_column: str
        peptide column to expand
    
    Returns
    ---------
    df : pandas DataFrame
        expanded dataframe

    """
    df_list = []
    for index, row in df.iterrows():
        if row[peptide_column] is not None:
            mods = find_modified_sites(row[peptide_column])
            for mod in mods:
                tmp = row.copy()
                tmp[peptide_column] = list(tmp[peptide_column].upper())
                tmp[peptide_column][mod[1]] = tmp[peptide_column][mod[1]].lower()
                tmp[peptide_column] = ''.join(tmp[peptide_column])
                df_list.append(tmp)
        else:
            df_list.append(row)
    df = pd.DataFrame(df_list).reset_index(drop=True)
    return df
        
def find_modified_sites(peptide, modification_types = None):
        """
        Finds all modification sites, indicated by a lower letter in the peptide column 
        and returns as a list of lists 
        Parameters
        ----------
        row : pandas row
        peptide_column : column that contains peptide with mod site
            mod site identified by being lower case

        Returns : [[modtype, relative_location]]
        1st element is the modification type
        2nd element is the relative modification location
        """

        if modification_types is not None:
            mod_types = ''.join(modification_types)
            mod_types = '[' + mod_types + ']'
            return [[m.group() , m.start()] for m in re.finditer(mod_types, peptide)]
        return [[m.group() , m.start()] for m in re.finditer('[a-z]', peptide)]

def find_peptide_locations(peptide, sequence):
    """
    Finds all start locations of a peptide in a given sequence. 
    Sequence starts at location 1.
    Example: Sequence: ABCDEFGABCDHIJK, Peptide: BCD
    would return [2, 9], as peptide is found at location 2 and 9
    """
    return [m.start() + 1 for m in re.finditer(peptide, sequence)]

def peptide_site_number(peptide, sequence, site = None, modification_types = None):
    """
    Finds the site number by finding the modification site in sequence and the peptide location
    in the sequence. A modification site is defined as where in a peptide a lower case letter appears.
    If a site is provided with alphabetical character the alphabetical character is
    the modification type checked for in the peptide.
    If a site contains numbers the closest match of the petide in the sequence to the site number is used as the site
    If no modification sites are found or peptide is not found in sequence None is returned
    
    Example: sequence is ZZZZABCDEFGZZZABCDEFGHIZZ, peptide is AbCdEfG
    If no site is provided and no modification types are provided then B6 is returned
    If no site is provided and modification_types is ['d','f'] or 'df' then D8 is returned
    If site is 30 and modification_types is ['d','f'] then F20 is returned
    If site is D23 then D18 is returned, modification_type is ignored
    If site is X50 then X50 is returned as modification_type x is not found in peptide
    Parameters
    ----------

    peptide :str
        peptide string where modification site is first lower-case letter found
    sequence : str
        sequence string where all characters are uppercase
    site : str
        site number to check against. can either have letters or just numbers
    modification_types: list
        list of modification types to check against while checking peptide
        ex: ['s', 't', 'y'] would ingore all modifications except for s,t,y
        not used if site contains alphabetical characters
    
    Returns
    ---------
    site : str
        site of peptide in sequence dataframe in format S123
        closest to original site if provided, otherwise first location found
        if no modification sites are found, the sequence cannot be found, or the peptide
            cannot be found in the sequence then the original site is returned
    """
    if pd.isnull(peptide):
        return None

    peptide = ''.join(filter(str.isalpha, peptide)) # keep only alphabetical letters
    all_mod_sites = find_modified_sites(peptide, None) # get all mod sites
        
    if isinstance(site, str):
        modification_types = site[0].lower()
    
    mod_sites = []
    if modification_types is not None:
        for mod_site in all_mod_sites:
            if mod_site[0] in modification_types:
                mod_sites.append(mod_site)
    else:
        mod_sites = all_mod_sites
    
    # if no mod sites are found return None
    if len(mod_sites) == 0:
        return None

    peptide_locations = find_peptide_locations(peptide.upper(), sequence) 


    relative_location = mod_sites[0][1] # only looking at first mod site location
    site_type = mod_sites[0][0].upper() # first mod type (STY)
    
    # get all locations where peptide appears in sequence
    peptide_locations = find_peptide_locations(peptide.upper(), sequence) 
    if len(peptide_locations)>0:
        if site is not None:
            # get closest peptide location to provided site
            site_num = int(re.sub("[^0-9]", "", str(site))) #remove all letters and make int
            closeness = np.inf
            closest_site = np.inf
            for loc in peptide_locations:
                true_location = loc + relative_location
                if abs(true_location - site_num) < closeness:
                    closeness = abs(true_location - site_num)
                    closest_site = true_location
            return site_type + str(closest_site)
        else:
            # no site info provided : return first instance
            return site_type + str(peptide_locations[0] + relative_location)  
    return None

def peptide_site_number_df(peptide, accession, df_sequence, site = None, modification_types = None):
    """
    Finds the site number by finding the modification site in sequence and the peptide location
    in the sequence from the given dataframe.
    If the sequence cannot be found then the original site is returned.
    See peptide_site_number(peptide,sequence,site,modification_types) 
    for details of how site is found once sequence is found
    
    Parameters
    ----------

    peptide :str
        peptide string where modification site is first lower-case letter found
    accession : str
        UniProt accession number 
    df_sequence : pandas DataFrame
        sequence dataframe pulled from UniProt with columns Entry, Sequence
    site : str
        site number to check against. can either have letters or just numbers
    modification_types: list
        list of modification types to check against while checking peptide
        ex: ['s', 't', 'y'] would ingore all modifications except for s,t,y
    
    Returns
    ---------
    site : str
        site of peptide in sequence dataframe in format S123
        closest to row[site_column] if provided, otherwise first location found
        if no modification sites are found, the sequence cannot be found, or the peptide
            cannot be found in the sequence then the original site is returned
    """
    sequence = df_sequence[df_sequence['Entry'] == accession]
    if len(sequence) > 0: # get first sequence for accession
        sequence = sequence.iloc[0]['Sequence']
        return peptide_site_number(peptide, sequence, site, modification_types)
    return site

def is_valid_site(site, sequence):
    """
    Checks whether the site provided is a valid site in the given sequence
    Valid sites are sites where site has same modification as sequence location
    
    Parameters
    ----------
    site: str
        site AA#, e.g. Y12
    sequence : str
        sequence to check
    
    Returns
    --------
    is_valid : bool
        True if site is valid. False otherwise
    """
    site_type  = ''.join(filter(str.isalpha, site)).upper()
    site_number = int( ''.join( filter(lambda x: x.isdigit(), site) ) )
    if site_number > 0 and site_number <= len(sequence):
        return site_type == sequence[site_number - 1]
    return False

def get_aligned_peptide(site, sequence, window):
    """
    Returns aligned peptide where modificaiton site is lowercase. 
    Aligned peptide contains +/- window AA surrounding modification site. 
    No buffer characters are added.
    If the site is not a valid location in the sequence then None is returned
    
    Parameters
    ----------
    site : str
        modification site. AA#. e.g. Y12
    sequence : str
        sequence to use for generating aligned peptide
    window : int
        number of amino acids to pull surrounding modification site
    
    Returns
    ----------
    peptide : str
        peptide of AA surrounding mod site. Mod site is lowercase. 
        None is returned if not a valid site
            e.g. ABCDeFGHI where site=E5, window=4
            e.g. AbCDEF    where site=B2, window=4 
    """
    if is_valid_site(site, sequence): 
        site_number = int( ''.join( filter(lambda x: x.isdigit(), site) ) ) - 1
        start = max(0, site_number - window)
        stop = min(len(sequence), site_number + window + 1)
        sequence = list(sequence)
        sequence[site_number] = sequence[site_number].lower()
        peptide = sequence[start:stop]
        peptide = ''.join(peptide)
        return peptide
    return None

def run_mapping(experiment, odir, name, map_columns, window=7, data_columns=None):
    if not os.path.exists(f"{odir}/MAPPED_DATA"): 
        os.mkdir(f"{odir}/MAPPED_DATA")   
    mapping_log = helpers.get_logger(f"mapping_{name}", f"{odir}/MAPPED_DATA/mapping_{name}.log")
    exp_mapper = ExperimentMapper(
        experiment = experiment,
        sequences = config.HUMAN_REF_SEQUENCES, 
        columns = map_columns, 
        logger = mapping_log, 
        compendia = config.HUMAN_REF_COMPENDIA, 
        window = window, 
        data_columns = data_columns)
     
    experiment = exp_mapper.experiment 
    experiment.to_csv(f"{odir}/MAPPED_DATA/{name}_mapped.tsv", sep = '\t', index = False)
    
    return experiment