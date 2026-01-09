import re
import numpy as np
import string

def check_for_modification(peptide):
    """
    Check if a peptide sequence contains any modified residues (lowercase s, t, or y)

    Parameters
    ----------
    peptide : str
        Peptide sequence to check

    Returns
    -------
    bool
        True if modified residues are found, False otherwise
    """
    if re.search(r'[sty]', peptide):
        return True
    else:
        return False

def change_mod_annotation(peptide, indicator = '*', after = True):
    """
    If a peptide uses a special character or series of characters to indicate phosphorylation (e.g. *), this function will find the indicator and lowercase the appropriate residue.

    Parameters
    ----------
    peptide : str
        Peptide sequence with modification indicator
    indicator : str
        String that indicates modification location
    after : bool
        Whether the indicator appears before the modified residue (False) or after (True, default)

    Returns
    -------
    str or np.nan
        Reformatted peptide sequence with lowercase modified residues, or np.nan if no modifications found
    """
    if not after:
        peptide = re.sub(re.escape(indicator) + r'([STY])', lambda x: x.group(1).lower(), peptide)
    else:
        #reformat peptide to lowercase modified residues based on indicator that comes after modified residue
        peptide = re.sub(r'([STY])' + re.escape(indicator), lambda x: x.group(1).lower(), peptide)

    #if indicator includes parentheses or brackets, remove any other content in those
    if '(' in indicator:
        #remove anything else in parentheses
        peptide = re.sub(r'\(([^)]*)\)', '', peptide)

    if '[' in indicator:
        #remove anything else in brackets
        peptide = re.sub(r'\[([^]]*)\]', '', peptide)

    #if tryptic fragments are included with periods, remove them
    if '.' in peptide:
        peptide = re.sub(r'^[A-Z]\.', '', peptide)  #remove leading amino acid and period
        peptide = re.sub(r'\.[A-Z]$', '', peptide)  #remove trailing amino acid and period
    
    #remove any other characters that are not letters
    peptide = re.sub(r'[^A-Za-z]', '', peptide)

    #check if there are any modifications left, if not, return np.nan
    #if not check_for_modification(peptide):
    #    return np.nan

    return peptide

def fix_pX_annotation(pep : str):
    """
    For cases where peptide sequences have modifications donated with lowercase letters (pY, pS, pT, or mK for other mods), this function will identify each phosphoresidue, lowercase it, and uppercase all other residues (removing any other modification annotations).
    """
    # First, remove any lowercased letters that are not p
    pep = re.sub(r'[a-oq-z]', '', pep)
    # Now replace pY, pS, and pT with lowercase y, s, and t
    pep = pep.replace('pY', 'y').replace('pS', 's').replace('pT', 't')
    
    return pep

def fix_centered_uppercase_peptide_seq(pep):
    """
    This finds the phosphoresidue that is at the center of the peptide sequence and lowercases that single peptide.
    """
    if pep != '':
        pep_length = len(pep) - 1
        mod_res = int(pep_length/2)
        pep_pieces = [pep[0:mod_res], pep[mod_res].lower(), pep[mod_res+1:]]
    else:
        pep_pieces = ['']
    return ''.join(pep_pieces)


def fix_unimod_ids(pep : str):
    """
    Finds Unimod accession IDs for modifications (e.g. Unimod_XX) and removes the placeholder. If it is Unimod_21 then lowercases the immediately previous character as that stands for phosphorylation.
    """

    phos = '(UniMod_21)'
    if phos in pep:
        
        pep_str = pep.split('(UniMod_21)')
        lim = len(pep_str) -1
        peplist = pep_str
        for i,substr in enumerate(pep_str):
            if i < lim:
                peplist[i] = peplist[i][0:-1] + peplist[i][-1].lower()
        new_pep = ''.join(peplist)
    else:
        new_pep = pep

    # Now remove all the other peptide modifications
    if '(' in new_pep:
        pep_str = re.split(r'[()]', new_pep)
        peplist = [x for x in pep_str if 'UniMod' not in x]
        new_pep = ''.join(peplist)
        
    return new_pep


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
    

def fix_prob_sequence(pep, thres = 0.75):
    ''' This cleans up sequences from MaxQuant that has the probabilities of the modification for a specific residue, by removing the probability and lowercasing the residue.'''
    
    # First break up the sequence by finding the (
    broke_pep = re.split(r'[()]', pep)

    # Now look for the numbers and then convert those to numbers
    res_probs = [float(i)  for i in broke_pep if is_number(i)]
    pep_parts = [i for i in broke_pep if not is_number(i)]
    


    for phos_indx, val in enumerate(res_probs):
        if val >= thres:
            # Now get the string part that is that index and lowercase the final character
            phos_pep = pep_parts[phos_indx]
            phos_pep_part1 = phos_pep[0:len(phos_pep)-1]
            phos_res = phos_pep[len(phos_pep)-1].lower()
            phos_pep = ''.join([phos_pep_part1, phos_res])
            pep_parts[phos_indx] = phos_pep

    
    new_pep = ''.join(pep_parts)
    return new_pep

def check_and_correct_nonphospho(pep: str):
    acceptable_mods = ['y','t','s']
    new_pep = [x.upper() if x not in acceptable_mods else x for x in pep ]
    new_pep = ''.join(new_pep)
    return new_pep

def detect_peptide_format(peptide:str, ignore_unrecognized = False):
    '''
    This takes a peptide sequence and from a predefined list of possibilities to determine how to adjust the sequence into an acceptable format for downstream analysis.

    Parameters
    ----------
    - peptide: str
        The peptide sequence of interest

    Returns
    -------
    - pep_format: str
        The format that the peptide falls under (e.g. acceptable, uppercase, annotated, maxquant_perc, unimod)
    '''

    # If there is pY, pS, or pT anywhere in the sequence, it is lowercase_p format 
    pY_format = any([x in peptide for x in ['pY','pS','pT']])
    # Check for lowercase letters
    lowercase = any([x.islower() for x in peptide])
    #check for parentheses or brackets
    parentheses = '(' in peptide
    bracket = '[' in peptide
    #check for possible indicators of phosphorylation
    other_indicator_symbols = ['*', '#']
    other_indicators = any([x in peptide for x in other_indicator_symbols])
    #check for other special characters that are not parentheses/brackets/underscores/other indicators
    others_chars = set(string.punctuation).difference(['(',')','.','[',']', '_'] + other_indicator_symbols)
    other_special = any([x in others_chars for x in peptide])

    if pY_format:
        pep_format = 'lowercase_p'
    elif parentheses or bracket:
        if 'Unimod' in peptide:
            pep_format = 'unimod'
        elif any([is_number(x) for x in peptide]): #check if there are numbers in the parentheses
            #check if numbers are alone in parentheses and between 0 and 1, indicating maxquant probability format
            paren_contents = re.findall(r'\(([^)]*)\)', peptide)
            #try to convert contents to float and check range
            if all([is_number(x) for x in paren_contents]):
                numbers = [float(x) for x in paren_contents]
                if all([0 <= x <= 1 for x in numbers]):
                    pep_format = 'maxquant_perc'
                else:
                    if not ignore_unrecognized:
                        raise PeptideSequenceError('New format detected with numbers in parentheses that are not between 0 and 1, cannot continue.')
                    else:
                        pep_format = 'unrecognized'
            else:
                if not ignore_unrecognized:
                    raise PeptideSequenceError('Numbers partially detected in parentheses but not all contents are numbers, unclear format, cannot continue.')
                else:
                    pep_format = 'unrecognized'
        elif any([x in peptide for x in ['Phospho', 'ph']]):
            pep_format = 'annotated'
        else:
            if not ignore_unrecognized:
                raise PeptideSequenceError('New format detected with information encompassed within the peptide and cannot continue.')
            else:
                pep_format = 'unrecognized'
    elif other_indicators:
        pep_format = 'annotated'
    elif other_special:
        other_chars_found = set([x for x in peptide if x in others_chars])
        if not ignore_unrecognized:
            raise PeptideSequenceError(f'New format detected with special characters ({', '.join(other_chars_found)}) and cannot continue.')
        else:
            pep_format = 'unrecognized'
    elif lowercase:
        # Check that there are y,s, or t
        if any([aa in peptide for aa in ['y','s','t']]):
            pep_format = 'acceptable'
        else:
            if not ignore_unrecognized:
                raise PeptideSequenceError('Lowercase letters detected but none correspond to phosphorylated residues (s, t, or y).')
            else:
                pep_format = 'unrecognized'
    else:
        pep_format = 'uppercase'
    return pep_format

def detect_most_common_format(peptide_list):
    """
    This function takes a list of peptide sequences and detects the most common format among them.
    
    :param peptide_list: Description
    """
    format_counts = {}
    for pep in peptide_list:
        try:
            pep_format = detect_peptide_format(pep)
            format_counts[pep_format] = format_counts.get(pep_format, 0) + 1
        except PeptideSequenceError:
            format_counts['unrecognized'] = format_counts.get('unrecognized', 0) + 1

    # if any formats were recognized, determine the most common one and report
    if format_counts:
        most_common_format = max(format_counts, key=format_counts.get)
        #raise error if majority are unrecognized
        if most_common_format == 'unrecognized':
            raise PeptideSequenceError('Majority of peptide formats are unrecognized.')
        #warn if some are unrecognized
        elif most_common_format != 'unrecognized' and 'unrecognized' in format_counts:
            print(f"Warning: Some peptides has an unrecognized format ({format_counts['unrecognized']}), although the majority are recognized.")
            del format_counts['unrecognized']

        #report all detected formats if multiple
        if len(format_counts) > 1:
            print("Multiple peptide formats detected in the provided list:")
            for key in format_counts.keys():
                print(f"Detected {format_counts[key]} peptides with format: {key}")

        #only use uppercase format if no other formats are detected (otherwise this may throw errors if there are too many sequences without modifications)
        if len(format_counts) > 1 and 'uppercase' in format_counts:
            print('Many peptides only contain uppercase letters, but other sequences detected different formats, suggesting that these uppercase-only sequences do not have modifications. Ignoring uppercase-only sequences for format detection.')
            del format_counts['uppercase']

        most_common_format = max(format_counts, key=format_counts.get)
        
        #if most common format is "uppercase", check to see if all sequences are of equal length (indicating centered modification)
        if most_common_format == 'uppercase':
            lengths = [len(pep) for pep in peptide_list]
            if len(set(lengths)) == 1:
                most_common_format = 'centered_uppercase'
            else:
                raise PeptideSequenceError('Provided sequences are all uppercased but are not centered, so it is not clear which residue is modified.')

        print(f'Most common peptide format detected: {most_common_format}')
        return most_common_format
    else:
        raise PeptideSequenceError('No valid peptide formats detected in the provided list.')
    
def format_peptide(peptide, format = 'acceptable', indicator = None, after = True):
    """
    Given a peptide and its current format, reformat it so that modified residues are lowercase s, t, or y, and any other symbols/annotations are removed.
    
    Parameters
    ---------- 
    peptide : str
        Peptide sequence to reformat
    format : str
        Current format of the peptide sequence (acceptable, uppercase, annotated, maxquant_perc, unimod, centered_uppercase)
    indicator : str
        If format is 'annotated', the string that indicates modification location
    after : bool
        If format is 'annotated', whether the indicator appears before the modified residue (False) or after (True, default)
    
    Returns
    -------
    str
        Reformatted peptide sequence with modified residues as lowercase s, t, or y, and other symbols/annotations removed.
    """
    if format == 'acceptable':
        return peptide
    elif format == 'centered_uppercase':
        new_pep = fix_centered_uppercase_peptide_seq(peptide)
    elif format == 'annotated':
        if indicator is None:
            raise ValueError('If sequence is annotated, an indicator string for the phosphorylation site must be provided.')
        new_pep = change_mod_annotation(peptide, indicator, after)
    elif format == 'maxquant_perc':
        new_pep = fix_prob_sequence(peptide)
    elif format == 'unimod':
        new_pep = fix_unimod_ids(peptide)
    elif format == 'lowercase_p':
        new_pep = fix_pX_annotation(peptide)
    else:
        raise PeptideSequenceError(f'Unrecognized peptide format: {format}')
    
    # Final check to ensure only acceptable modifications are present
    if new_pep == new_pep:
        new_pep = check_and_correct_nonphospho(new_pep)

    return new_pep

def format_peptide_list(peptide_list, format = 'annotated', indicator = '*', after = True):
    """
    Given a list of peptides and their current format, reformat them so that modified residues are lowercase s, t, or y, and any other symbols/annotations are removed.
    
    Parameters
    ---------- 
    peptide_list : list of str
        List of peptide sequences to reformat
    format : str
        Current format of the peptide sequences (acceptable, uppercase, annotated, maxquant_perc, unimod, centered_uppercase)
    indicator : str
        If format is 'annotated', the string that indicates modification location
    after : bool
        If format is 'annotated', whether the indicator appears before the modified residue (False) or after (True, default)
    
    Returns
    -------
    list of str
        List of reformatted peptide sequences with modified residues as lowercase s, t, or y, and other symbols/annotations removed.
    """
    formatted_peptides = []
    for pep in peptide_list:
        new_pep = format_peptide(pep, format=format, indicator=indicator, after=after)
        formatted_peptides.append(new_pep)
    return formatted_peptides

def detect_annotation(peptide):
    """
    Identify the indicator/annotation used in a peptide sequence to denote phosphorylation sites.
    
    :param peptide: Description
    """
    possible_indicators = ['*', '#', '(Phospho)', '(ph)', '[Phospho]', '[ph]']
    present_indicators = set([indicator for indicator in possible_indicators if indicator in peptide])
    if len(present_indicators) == 0:
        raise PeptideSequenceError('No recognized modification indicators found in the peptide sequence, you will need to manually indicate what indicator to use.')
    elif len(present_indicators) > 1:
        raise PeptideSequenceError(f'Multiple modification indicators found in the peptide sequence: {", ".join(present_indicators)}. You will need to  manually indicate what indicator to use.')
    else:
        #determine if indicator is before or after modified residue
        indicator = list(present_indicators)[0]
        indicator_index = peptide.index(indicator)
        if indicator_index == 0:
            #Indicator found at start of peptide, assuming it appears before modified residue.
            after = False
        elif indicator_index == len(peptide) - len(indicator):
            #Indicator found at end of peptide, assuming it appears after modified residue.
            after = True
        else:
            #Indicator found in middle of peptide, need to check surrounding characters
            if peptide[indicator_index - 1] in ['S', 'T', 'Y']:
                after = True
            elif peptide[indicator_index + len(indicator)] in ['S', 'T', 'Y']:
                after = False
            else:
                raise PeptideSequenceError('Unable to determine if modification indicator appears before or after modified residue. You will need to manually indicate this.')

        return indicator, after


def format_peptide_from_df(df, peptide_col, new_peptide_col = None, autodetect = True, **kwargs):
    """
    Given a phosphoproteomics dataset with a column of peptide sequences, reformat the peptides so that modified residues are lowercase s, t, or y, and any other symbols/annotations are removed. The reformatted peptides will be added as a new column in the dataframe.

    Can either autodetect the peptide format or have it specified. If autodetecting, the most common format among the peptides will be used, although a warning will be issued to check the results alongside a column indicating the detected format for each peptide.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe containing phosphoproteomics data
    peptide_col : str
        Name of the column containing peptide sequences
    new_peptide_col : str, optional
        Name of the new column to add with reformatted peptides. If None, defaults to 'Formatted Peptide'.
    autodetect : bool, optional
        Whether to autodetect the peptide format (default True). If False, the format must be specified via the 'format' keyword argument.
    **kwargs : dict
        Additional arguments to pass to the formatting functions, such as 'indicator' and 'after' if format is 'annotated'.
    """
    
    peptides = df[peptide_col].tolist()
    if autodetect:
        print('Autodetecting peptide format...')
        detected_format = detect_most_common_format(peptides)
        print('Proceeding to format peptides...')
    elif 'format' in kwargs:
        detected_format = kwargs['format']
    else:
        raise ValueError('If autodetect is False, a peptide format must be provided. Options are: acceptable, uppercase, annotated, maxquant_perc, unimod, centered_uppercase.')
    
    if detected_format == 'annotated':
        if ('indicator' in kwargs and 'after' not in kwargs) or ('after' in kwargs and 'indicator' not in kwargs):
            raise ValueError('If one of `indicator` or `after` is provided, both must be provided together.')
        if 'indicator' not in kwargs or 'after' not in kwargs:
            #make sure multiple indicators are not present and detect which one is used
            indicators_found = set()
            for pep in peptides:
                if any([x in pep for x in ['*', '#', '(Phospho)', '(ph)', '[Phospho]', '[ph]']]):
                    indicator, after = detect_annotation(pep)
                    indicators_found.add((indicator, after))

            if len(indicators_found) > 1:
                raise PeptideSequenceError('Multiple modification indicators detected in the peptide list, you will need to manually indicate what indicator to use.')


    
    formatted_peptides = format_peptide_list(peptides, format=detected_format, **kwargs)

    #add formatted peptides to dataframe
    new_peptide_col = new_peptide_col if new_peptide_col is not None else 'Formatted Peptide'
    df[new_peptide_col] = formatted_peptides

    #report warning to check peptides if autodetect was used
    if autodetect:
        format_col = f'Detected Peptide Format'
        df[format_col] = df[peptide_col].apply(lambda x: detect_peptide_format(x, ignore_unrecognized=True))
        print('\nPeptide formatting complete. Please verify that the formatted peptides are correct before proceeding.\nYou can see the detected format for each peptide in the column: ', format_col)

    return df



class PeptideSequenceError(Exception):
    pass