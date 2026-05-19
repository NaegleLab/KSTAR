import re
import time
import json
import zlib
from xml.etree import ElementTree
from urllib.parse import urlparse, parse_qs, urlencode
import requests
from requests.adapters import HTTPAdapter, Retry
import pandas as pd
import numpy as np

#fetch uniprot_fields from uniprot.org
import requests

uniprot_url = 'https://rest.uniprot.org/configure/idmapping/fields'
response = requests.get(uniprot_url)
if response.status_code == 200:
    uniprot_fields = response.json()
else:
    print(f'Warning: Unable to fetch uniprot fields from {uniprot_url}. Status code: {response.status_code}. Will have to hardcode fields if needed.')

#get all valid uniprot fields
uniprot_field_names = []
for group in uniprot_fields['groups']:
    for item in group['items']:
        uniprot_field_names.append(item['name'])


##### Code adapted from https://www.uniprot.org/help/id_mapping_prog #####
POLLING_INTERVAL = 3
API_URL = "https://rest.uniprot.org"


retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

#get ids found in kstar reference
from kstar import config
kstar_ids = config.HUMAN_REF_COMPENDIA['KSTAR_ACCESSION'].unique()


def check_response(response):
    try:
        response.raise_for_status()
    except requests.HTTPError:
        print(response.json())
        raise


def submit_id_mapping(from_db, to_db, ids, taxonID=None):
    if taxonID is None:
        request = requests.post(
            f"{API_URL}/idmapping/run",
            data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
        )
    else:
        request = requests.post(
            f"{API_URL}/idmapping/run",
            data={"from": from_db, "to": to_db, "ids": ",".join(ids), "taxId": str(taxonID)},
        )
    check_response(request)
    return request.json()["jobId"]


def get_next_link(headers):
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def check_id_mapping_results_ready(job_id):
    while True:
        request = session.get(f"{API_URL}/idmapping/status/{job_id}")
        check_response(request)
        j = request.json()
        if "jobStatus" in j:
            if j["jobStatus"] == "RUNNING":
                print(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(request["jobStatus"])
        else:
            return bool(j["results"] or j["failedIds"])


def get_batch(batch_response, file_format, compressed):
    batch_url = get_next_link(batch_response.headers)
    while batch_url:
        batch_response = session.get(batch_url)
        batch_response.raise_for_status()
        yield decode_results(batch_response, file_format, compressed)
        batch_url = get_next_link(batch_response.headers)


def combine_batches(all_results, batch_results, file_format):
    if file_format == "json":
        for key in ("results", "failedIds"):
            if key in batch_results and batch_results[key]:
                all_results[key] += batch_results[key]
    elif file_format == "tsv":
        return all_results + batch_results[1:]
    else:
        return all_results + batch_results
    return all_results


def get_id_mapping_results_link(job_id):
    url = f"{API_URL}/idmapping/details/{job_id}"
    request = session.get(url)
    check_response(request)
    return request.json()["redirectURL"]


def decode_results(response, file_format, compressed):
    if compressed:
        decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
        if file_format == "json":
            j = json.loads(decompressed.decode("utf-8"))
            return j
        elif file_format == "tsv":
            return [line for line in decompressed.decode("utf-8").split("\n") if line]
        elif file_format == "xlsx":
            return [decompressed]
        elif file_format == "xml":
            return [decompressed.decode("utf-8")]
        else:
            return decompressed.decode("utf-8")
    elif file_format == "json":
        return response.json()
    elif file_format == "tsv":
        return [line for line in response.text.split("\n") if line]
    elif file_format == "xlsx":
        return [response.content]
    elif file_format == "xml":
        return [response.text]
    return response.text


def get_xml_namespace(element):
    m = re.match(r"\{(.*)\}", element.tag)
    return m.groups()[0] if m else ""


def merge_xml_results(xml_results):
    merged_root = ElementTree.fromstring(xml_results[0])
    for result in xml_results[1:]:
        root = ElementTree.fromstring(result)
        for child in root.findall("{http://uniprot.org/uniprot}entry"):
            merged_root.insert(-1, child)
    ElementTree.register_namespace("", get_xml_namespace(merged_root[0]))
    return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)


def print_progress_batches(batch_index, size, total):
    n_fetched = min((batch_index + 1) * size, total)
    print(f"Fetched: {n_fetched} / {total}")


def get_id_mapping_results_search(url):
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    if "size" in query:
        size = int(query["size"][0])
    else:
        size = 500
        query["size"] = size
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    parsed = parsed._replace(query=urlencode(query, doseq=True))
    url = parsed.geturl()
    request = session.get(url)
    check_response(request)
    results = decode_results(request, file_format, compressed)
    total = int(request.headers["x-total-results"])
    print_progress_batches(0, size, total)
    for i, batch in enumerate(get_batch(request, file_format, compressed), 1):
        results = combine_batches(results, batch, file_format)
        print_progress_batches(i, size, total)
    if file_format == "xml":
        return merge_xml_results(results)
    return results


def get_id_mapping_results_stream(url):
    if "/stream/" not in url:
        url = url.replace("/results/", "/results/stream/")
    request = session.get(url)
    check_response(request)
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    return decode_results(request, file_format, compressed)

def construct_id_map(accessions, from_id = 'UniProtKB_AC-ID', to_id = 'UniProtKB', taxonID = 9606):
    """
    Construct a mapping from one type of ID to another using UniProt mapping services

    Parameters
    ----------
    accessions : list
        List of accessions to map
    from_id : str
        Name of the database in UniProt rest api to convert from. Options include:
        'UniProtKB_AC-ID' : UniProt Accession or Entry name
        'Gene_Name' : Gene name
    to_id : str
        Name of the database in UniProt rest api to convert to. Options include:
        'UniProtKB' : UniProt Accession or Entry name
        'Gene_Name' : Gene name
    taxonID : int
        Taxonomic ID for the species of interest

    Returns
    -------
    dict
        Dictionary mapping from_id to to_id
    """
    job_id = submit_id_mapping(
        from_db=from_id, to_db=to_id, ids=accessions, taxonID=taxonID
    )

    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        results = get_id_mapping_results_search(link)

    # create mapping dictionary
    id_to_acc = {}
    for entry in results['results']:
        from_id_entry = entry['from']
        if to_id == 'Gene_Name':
            to_id_entry = entry['to']
        else:
            to_id_entry = entry['to']['primaryAccession']
        id_to_acc[from_id_entry] = to_id_entry

    return id_to_acc


def check_if_in_reference(accession):
    return accession.split('-')[0] in kstar_ids

def get_all_ids_in_df(df, accession_col, id_sep = None):
    """
    Get a list of all unique accessions in a dataframe given the column containing the accessions and the separator (if multiple accessions are in the same cell)

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe containing the dataset with accession column
    accession_col : str
        Name of the column containing the protein accessions
    id_sep : str, optional
        Separator for multiple accessions in the same cell (default is None, meaning there is only one accession per cell)

    Returns
    -------
    list
        List of all unique accessions in the dataframe
    """
    if id_sep is not None:
        #check if separator exists in the column
        if df[accession_col].str.contains(id_sep).any():
            ids = df[accession_col].str.split(id_sep)
            ids = ids.explode().unique().tolist()
        else:
            ids = df[accession_col].unique().tolist()
    else:
        ids = df[accession_col].unique().tolist()

    return ids

def check_and_convert_UniProtKB_accessions(accessions):
    """
    Given a list of UniProtKB accessions, check if any are non-reviewed accessions and if so, convert them to their associated reviewed accession using uniprot mapping services. This is important since non-reviewed accessions may not be found in the KSTAR reference but their associated reviewed accessions may be.

    Parameters
    ----------
    accessions : list
        List of accessions to check and convert

    Returns
    -------
    id_to_acc : dict
        Dictionary mapping from original accession to reviewed accession (if conversion was needed)
    """
    #
    found_accessions = [i for i in accessions if check_if_in_reference(i)]
    not_found_accessions = [i for i in accessions if not check_if_in_reference(i)]

    #for entries not found in reference, check if there is a more updated accession associated with the same gene
    id_map_to_gene = construct_id_map(not_found_accessions, from_id = 'UniProtKB_AC-ID', to_id = 'Gene_Name')
    #construct map to convert from gene name back to swissprot accession
    genes = list(id_map_to_gene.values())
    id_map_to_acc = construct_id_map(genes, from_id = 'Gene_Name', to_id = 'UniProtKB-Swiss-Prot', taxonID = 9606)
    #combine maps into one
    id_to_acc = {k:id_map_to_acc[v] for k, v in id_map_to_gene.items() if v in id_map_to_acc}
    #add known accessions found in reference to mapping dictionary
    for acc in found_accessions:
        id_to_acc[acc] = acc
    return id_to_acc

def convert_row_with_separator(row, accession_col, id_sep, id_to_acc):
    """
    Given a row with multiple accessions in the same cell separated by a separator, convert each accession to its associated reviewed accession (if needed) and return a string with the converted accessions separated by the same separator

    Parameters
    ----------
    row : pandas.Series
        Row of the dataframe containing the accession column
    accession_col : str
        Name of the column containing the protein accessions
    id_sep : str
        Separator for multiple accessions in the same cell
    id_to_acc : dict
        Dictionary mapping from original accession to reviewed accession (if conversion was needed)

    Returns
    -------
    str
        String with converted accessions separated by the same separator
    """
    if id_sep not in row[accession_col]:
        return id_to_acc.get(row[accession_col], np.nan)
    else:
        accessions = row[accession_col].split(id_sep)
        converted_accessions = np.unique([id_to_acc[acc] for acc in accessions if acc in id_to_acc])
        if len(converted_accessions) == 0:
            return np.nan
        else:
            return id_sep.join(converted_accessions)

def convert_ids_from_df(df, accession_col, from_id = 'UniProtKB_AC-ID', to_id = 'UniProtKB-Swiss-Prot', id_sep = ';', taxonID = 9606, remove_unmapped = True, report_missed = True):
    """
    Given a dataset and the location of the accession column, use uniprot mapping services to convert from one type of id to another

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe containing the dataset with accession column
    accession_col : str
        Name of the column containing the protein accessions
    from_id : str
        Name of the database in UniProt rest api to convert from. Options include:
        'UniProtKB_AC-ID' : UniProt Accession or Entry name
        'Gene_Name' : Gene name
    to_id : str
        Name of the database in UniProt rest api to convert to (default is 'UniProtKB') 
    remove_unmapped : bool
        Whether to remove rows with unmapped accessions (default is True). Otherwise these rows will be returned with NaN accessions and can be accessed through the returned missing_rows dataframe.
    """
    #check from_id 
    if from_id not in ['UniProtKB_AC-ID', 'Gene_Name', 'GI_number', 'RefSeq_Protein', 'RefSeq_Nucleotide', 'Ensembl', 'Ensembl_Protein', 'Ensembl_Transcript', 'CCDS', 'UniParc', 'UniRef50', 'UniRef90', 'UniRef100', 'UniProtKB']:
        raise ValueError(f'from_id {from_id} not supported. Supported options at this time include: UniProtKB_AC-ID, Gene_Name')
    
    if (from_id == 'Gene_Name') and to_id == 'UniProtKB':
        print(f'Warning: Mapping from {from_id} may result in many-to-one mappings. Switching to "UniProtKB-Swiss-Prot" as `to_id`.')
        to_id = 'UniProtKB-Swiss-Prot'
    elif (from_id == 'UniProtKB_AC-ID') and to_id == 'UniProtKB-Swiss-Prot':
        #won't return anything if try to map to swiss prot
        to_id = 'UniProtKB'
    

    #check to make sure to_id is valid with from_id
    valid_to_ids = get_valid_to_ids(from_id)
    if from_id == 'UniProtKB_AC-ID':
        valid_to_ids.append('Gene_Name')
    if to_id not in valid_to_ids:
        raise ValueError(f'to_id {to_id} not valid for from_id {from_id}. Valid to_id options for this from_id include: {valid_to_ids}')
    
    #check if taxonID is needed
    if from_id == 'Gene_Name' and taxonID is None:
        raise ValueError('taxonID is required when from_id is Gene_Name since gene names are not unique across species. Please provide a taxonID to ensure accurate mapping.')
    elif from_id != 'Gene_Name' and taxonID is not None:
        taxonID = None
    
    #remove rows without accession ids
    df = df.dropna(subset = [accession_col]).copy()
    ids = get_all_ids_in_df(df, accession_col, id_sep)

    if from_id == 'UniProtKB':
        #for uniprot accessions, check if in kstar reference, if not attempt to convert to reviewed accession from swissprot
        id_to_acc = check_and_convert_UniProtKB_accessions(ids)
    else:
        
        #check if separator exists
        #if id_sep is not None:
        #    if df[accession_col].str.contains(id_sep).any():
        ids = df[accession_col].unique().tolist()
        if from_id == 'GI_number':
            #remove 'gi|' prefix from gi numbers for uniprot mapping
            ids = [id.strip('gi|') for id in ids]
        id_to_acc = construct_id_map(ids, from_id = from_id, to_id = to_id, taxonID = taxonID)

        #if from_id is GI_number, need to add 'gi|' prefix back to unmapped ids since uniprot mapping requires removing the prefix but our dataset may still have it for unmapped ids
        if from_id == 'GI_number':
            id_to_acc = {f'gi|{key}': value for key, value in id_to_acc.items()}

    if id_sep is None:
        df['Accession'] = df[accession_col].map(id_to_acc)
    else:
        df['Accession'] = df.apply(lambda row: convert_row_with_separator(row, accession_col, id_sep, id_to_acc), axis=1)
    num_rows_without_acc = df['Accession'].isna().sum()
    if num_rows_without_acc > 0:
        if report_missed:
            print(f'Warning: {num_rows_without_acc} rows could not be mapped to UniProtKB accessions and will be removed from dataset')
        missing_rows = df[df['Accession'].isna()]
        if remove_unmapped:
            df = df.dropna(subset = ['Accession'])
    else:
        missing_rows = pd.DataFrame()

    return df, missing_rows

def automatic_id_conversion(df, accession_col, taxonID = 9606, keep_isoform_info = False,remove_unmapped = True, id_sep = None):
    """
    Given a dataframe with an accession column, automatically detect the type of accession and convert it to UniProtKB accessions using uniprot mapping services. This is useful for datasets that may have a mix of different types of accessions or when the type of accession is not specified.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe containing the dataset with accession column
    accession_col : str
        Name of the column containing the protein accessions
    taxonID : int
        Taxonomic ID for the species of interest (default is 9606 for human). Only needed if converting from gene names since gene names are not unique across species.
    remove_unmapped : bool
        Whether to remove rows with unmapped accessions (default is True). Otherwise these rows will be returned with NaN accessions and can be accessed through the returned missing_rows dataframe.
    """
    #remove rows without accession ids
    df = df.dropna(subset = [accession_col]).copy()
    #if no separator is provided, check if there might be and report
    if id_sep is None:
        if df[accession_col].apply(lambda x: ';' in x or ',' in x if isinstance(x, str) else False).any():
            print(f'Warning: It looks like there may be multiple accessions in the {accession_col} column. If so, please provide the separator using the `id_sep` argument to ensure accurate mapping. Attempting to automatically detect separator...')
            #count instances of semicolons and commas to determine separator if there are multiple accessions in the same cell
            semicolon_count = df[accession_col].apply(lambda x: x.count(';') if isinstance(x, str) else 0).sum()
            comma_count = df[accession_col].apply(lambda x: x.count(',') if isinstance(x, str) else 0).sum()

            #based on counts, determine likely separator
            if semicolon_count > comma_count and semicolon_count > df.shape[0] * 0.1: #if more than 10% of rows have semicolons, assume that is the separator
                id_sep = ';'
                print(f'Automatically detected ";" as separator for multiple accessions in the same cell.')
            elif comma_count > semicolon_count and comma_count > df.shape[0] * 0.1: #if more than 10% of rows have commas, assume that is the separator
                id_sep = ','
                print(f'Automatically detected "," as separator for multiple accessions in the same cell.')
            else:
                print(f'Could not automatically detect separator for multiple accessions in the same cell. Proceeding with assumption that there is only one accession per cell. If this is not the case, please provide the separator using the `id_sep` argument to ensure accurate mapping.')




    #get accession types and whether they are in reference
    df['Original Accession Type'] = df.apply(lambda x: identify_accession_type_in_row(x, accession_col, id_sep = id_sep), axis=1)
    accession_types = df['Original Accession Type'].unique().tolist()


    #for accessions not found in reference, split by accession type and convert each type separately
    split_dfs = {acc_type: df[df['Original Accession Type'] == acc_type] for acc_type in df['Original Accession Type'].unique()}

    missed = {}
    for acc_type, split_df in split_dfs.items():
        if acc_type in ['UniProtKB_ACC-ID']:
            split_dfs[acc_type], missed[acc_type] = convert_ids_from_df(split_df, accession_col, from_id = acc_type, to_id = 'UniProtKB', taxonID = taxonID, remove_unmapped = remove_unmapped, report_missed=False)
        elif acc_type == 'genename/unrecognized':
            #assume these are gene names and try to convert from gene names to swissprot accessions
            print("Couldn't identify accession type for some accessions. Attempting to convert these using gene names...")
            split_dfs[acc_type], missed[acc_type] = convert_ids_from_df(split_df, accession_col, from_id = 'Gene_Name', to_id = 'UniProtKB-Swiss-Prot', taxonID = taxonID, remove_unmapped = remove_unmapped, report_missed=False)
        else:
            split_dfs[acc_type], missed[acc_type] = convert_ids_from_df(split_df, accession_col, from_id = acc_type, to_id = 'UniProtKB-Swiss-Prot', taxonID = taxonID, remove_unmapped = remove_unmapped, report_missed=False)

    #combine split dfs and missing df
    df = pd.concat(list(split_dfs.values()), ignore_index = True)
    missing_rows = pd.concat(list(missed.values()), ignore_index = True)

    #remove isoform information if keep_isoform_information is false (remove '-1', '-2', etc. from accessions)
    if not keep_isoform_info:
        df['Accession'] = df['Accession'].apply(lambda x: get_unique_uniprot_ids(x, id_sep = id_sep) if isinstance(x, str) else x)


    num_mapped = df['Accession'].notna().sum()
    num_missed = missing_rows.shape[0]
    print(f'After automatic id conversion, {num_missed} rows could not be mapped to UniProtKB accessions. Overall, {num_mapped} ({num_mapped / (num_mapped + num_missed) * 100:.1f}%) rows have valid accessions after conversion.')
    return df, missing_rows



def get_valid_from_ids():
    """
    Return a list of valid from_id options for uniprot mapping services

    Returns
    -------
    list
        List of valid from_id options
    """
    from_ids = []
    for group in uniprot_fields['groups']:
        for item in group['items']:
            if item['from']:
                from_ids.append(item['name'])
    return from_ids

def get_valid_to_ids(from_id = None):
    """
    Return a list of valid to_id options for uniprot mapping services

    Parameters
    ----------
    from_id : str, optional
        Name of the database in UniProt rest api to convert from. If provided, only return to_id options that are compatible with the given from_id.

    Returns
    -------
    list
        List of valid to_id options
    """
    if from_id is None:
        to_ids = []
        for group in uniprot_fields['groups']:
            for item in group['items']:
                if item['to']:
                    to_ids.append(item['name'])
    else:
        #find from id information
        field_info = get_field_information(from_id)
        ruleId = field_info['ruleId']
        #if no rule, grab all to ids
        if ruleId is None:
            to_ids = []
            for group in uniprot_fields['groups']:
                for item in group['items']:
                    if item['to']:
                        to_ids.append(item['name'])
        else:
            #if rule exists, grab to ids for that rule
            to_ids = uniprot_fields['rules'][ruleId - 1]['tos']

    return to_ids

def get_field_information(field_name):
    """
    Given a uniprot field name, return the field information from uniprot.org

    Parameters
    ----------
    field_name : str
        Name of the uniprot field

    Returns
    -------
    dict
        Field information dictionary
    """
    for group in uniprot_fields['groups']:
        for item in group['items']:
            if item['name'] == 'UniParc':
                return item
    
    raise ValueError(f'Field name {field_name} not found in uniprot fields')

def check_if_id_supported(from_id):
    """
    Check if the given from_id is supported for uniprot mapping services

    Parameters
    ----------
    from_id : str
        Name of the database in UniProt rest api to convert from. Options include:

    """
    field_info = get_field_information(from_id)
    return field_info['from'], field_info['to']

def is_swissprot_accession(s):
    """
    Returns True if the string matches a Swiss-Prot (UniProtKB/Swiss-Prot) accession format. This is defined as either a 6-character or 10-character alphanumeric string with specific patterns.
    """
    #indicate which characters are currently found at the start of swissprot accessions based on uniprot field information
    start_charcters = r'[ABCDEFGHIOPQ]'
    # Swiss-Prot (UniProtKB/Swiss-Prot) accession regex
    SWISSPROT_REGEX = re.compile(r'^(?:' + start_charcters + r'[0-9][A-Z0-9]{3}[0-9]|A0A[0-9][A-Z0-9]{6})$')
    return bool(SWISSPROT_REGEX.match(s))



def identify_accession_type(acc):
    """
    Given a protein accession, identify the type of accession it is
    
    Parameters
    ----------
    acc : str
        Protein accession

    Returns
    -------
    str
        Type of accession. Options include:
        'UniProtKB_AC-ID' : UniProt Accession or Entry name
        'Ensembl Gene' : Ensembl Gene ID
        'Ensembl Protein' : Ensembl Protein ID
        'Ensembl Transcript' : Ensembl Transcript ID
        'RefSeq_Nucleotide' : RefSeq Nucleotide ID
        'RefSeq_Protein' : RefSeq Protein ID
        'CCDS' : CCDS ID
    """
    #determine if accessions are uniprot
    if acc.endswith('_HUMAN'):
        return 'UniProtKB_AC-ID'
    #UniParc
    if acc.startswith('UPI'):
        return 'UniParc'
    #UniRef
    if acc.startswith('UniRef50'):
        return 'UniRef50'
    elif acc.startswith('UniRef90'):
        return 'UniRef90'
    elif acc.startswith('UniRef100'):
        return 'UniRef100'
    #CCDS 
    if acc.startswith('CCDS'):
        return 'CCDS'
    
    #ensembl based accessions
    if acc.startswith('ENSG'):
        return 'Ensembl'
    if acc.startswith('ENSP'):
        return 'Ensembl_Protein'
    if acc.startswith('ENST'):
        return 'Ensembl_Transcript'
    if acc.startswith('NM_') or acc.startswith('XM_'):
        return 'RefSeq_Nucleotide'
    if acc.startswith('NP_') or acc.startswith('XP_'):
        return 'RefSeq_Protein'
    if acc.startswith('gi|'):
        return 'GI_number'
    
    #uniprot id after checking others since it is a little harder to identify (remove isoform detail ('-1') before checking)
    if is_swissprot_accession(acc.split('-')[0]):
        return 'UniProtKB'
    
    return 'genename/unrecognized'

def identify_accession_type_in_row(row, accession_col, id_sep = None):
    """
    Given a row of a dataframe and the column containing the accessions, identify the type of accession for each accession in the cell and return a list of the unique accession types found in the cell. This is useful for identifying the types of accessions in a dataset that may have a mix of different types of accessions or when the type of accession is not specified.

    Parameters
    ----------
    row : pandas.Series
        Row of the dataframe containing the accession column
    accession_col : str
        Name of the column containing the protein accessions
    id_sep : str, optional
        Separator for multiple accessions in the same cell (default is None, meaning there is only one accession per cell)

    Returns
    -------
    str
        Most common type of accession found in the cell
    """
    if id_sep is not None:
        if id_sep in row[accession_col]:
            accessions = row[accession_col].split(id_sep)
            accession_type = identify_most_common_accession_type(accessions)
            return accession_type
        else:
            return identify_accession_type(row[accession_col])
    else:
        return identify_accession_type(row[accession_col])

def identify_most_common_accession_type(accessions):
    """
    Given a list of protein accessions, identify the most common type of accession
    
    Parameters
    ----------
    accessions : list
        List of protein accessions

    Returns
    -------
    str
        Most common type of accession. Options include:
        'UniProtKB_AC-ID' : UniProt Accession or Entry name
        'Ensembl Gene' : Ensembl Gene ID
        'Ensembl Protein' : Ensembl Protein ID
        'Ensembl Transcript' : Ensembl Transcript ID
        'RefSeq_Nucleotide' : RefSeq Nucleotide ID
        'RefSeq_Protein' : RefSeq Protein ID
        'CCDS' : CCDS ID
    """
    type_counts = {}
    for acc in accessions:
        acc_type = identify_accession_type(acc)
        if acc_type:
            if acc_type not in type_counts:
                type_counts[acc_type] = 0
            type_counts[acc_type] += 1
    
    
    #return the most common type
    most_common_type = max(type_counts, key=type_counts.get)

    if most_common_type == 'genename/unrecognized' and len(type_counts) > 1:
        #check if there is a second most common type that is not 'genename/unrecognized' and if so, return that type instead since it is likely that the accessions are actually that type but some gene names are mixed in or there is some error in the accessions
        second_most_common_type = sorted(type_counts, key=type_counts.get, reverse=True)[1]
        return second_most_common_type
    
    return most_common_type

def get_unique_uniprot_ids(id_string, id_sep = None):
    """
    Given a string of accessions separated by a separator, return a list of the unique UniProtKB accessions in the string, removing any isoform information (e.g. '-1') from the accessions. This is useful for ensuring that all accessions are in a consistent format for mapping to the KSTAR reference and there won't be redundant entries when mapping to the same UniProtKB accession with different isoform information.

    Parameters
    ----------
    id_string : str
        String of accessions separated by a separator
    id_sep : str
        Separator for multiple accessions in the same cell (default is ';')
    
    Returns
    -------
    string
        String of unique UniProtKB accessions separated by the same separator, with isoform information removed
    """
    if id_sep is not None:
        accessions = id_string.split(id_sep)
        unique_accessions = np.unique([acc.split('-')[0] for acc in accessions])
        return id_sep.join(unique_accessions)
    else:
        return id_string.split('-')[0]
