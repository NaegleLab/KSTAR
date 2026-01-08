import re
import time
import json
import zlib
from xml.etree import ElementTree
from urllib.parse import urlparse, parse_qs, urlencode
import requests
from requests.adapters import HTTPAdapter, Retry
import pandas as pd

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


def check_response(response):
    try:
        response.raise_for_status()
    except requests.HTTPError:
        print(response.json())
        raise


def submit_id_mapping(from_db, to_db, ids, taxonID=9606):
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

def convert_id(df, accession_col, from_id = 'UniProtKB_AC-ID', to_id = 'UniProtKB-Swiss-Prot', taxonID = 9606):
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
    """
    #check from_id 
    if from_id not in ['UniProtKB_AC-ID', 'Gene_Name']:
        raise ValueError(f'from_id {from_id} not supported. Supported options at this time include: UniProtKB_AC-ID, Gene_Name')
    
    if from_id == 'Gene_Name' and to_id != 'UniProtKB-Swiss-Prot':
        print('Warning: Mapping from Gene_Name may result in many-to-one mappings. It is recommended to map to "UniProtKB-Swiss-Prot" as `to_id` to ensure only reviewed entries are returned.')

    #check to make sure to_id is valid with from_id
    valid_to_ids = get_valid_to_ids(from_id)
    if to_id not in valid_to_ids:
        raise ValueError(f'to_id {to_id} not valid for from_id {from_id}. Valid to_id options for this from_id include: {valid_to_ids}')
    
    #remove rows without accession ids
    df = df.dropna(subset = [accession_col]).copy()
    #convert UniProt entry name to accession with uniprot mapping services
    ids = df[accession_col].unique().tolist()
    job_id = submit_id_mapping(
        from_db=from_id, to_db=to_id, ids=ids, taxonID=taxonID
    )

    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        results = get_id_mapping_results_search(link)

    #create mapping dictionary
    id_to_acc = {}

    for entry in results['results']:
        from_id = entry['from']
        to_id = entry['to']['primaryAccession']
        id_to_acc[from_id] = to_id


    df['Accession'] = df[accession_col].map(id_to_acc)
    num_rows_without_acc = df['Accession'].isna().sum()
    if num_rows_without_acc > 0:
        print(f'Warning: {num_rows_without_acc} rows could not be mapped to UniProtKB accessions and will be removed from dataset')
        missing_rows = df[df['Accession'].isna()]
        df = df.dropna(subset = ['Accession'])
    else:
        missing_rows = pd.DataFrame()

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
    # Swiss-Prot (UniProtKB/Swiss-Prot) accession regex
    SWISSPROT_REGEX = re.compile(r'^(?:[A-Z][0-9][A-Z0-9]{3}[0-9]|A0A[0-9][A-Z0-9]{6})$')
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
    
    #uniprot id after checking others since it is a little harder to identify
    if is_swissprot_accession(acc):
        return 'UniProtKB_AC-ID'
    
    return 'genename/unrecognized'

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
    return most_common_type