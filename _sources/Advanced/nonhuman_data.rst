KSTAR with Non-Human Datasets
=============================

By default, KSTAR can only work with data originating from human samples. This is largely due to limited kinase-substrate information for non-human proteins and the fact that prediction networks were generated using human data. 

However, it is possible to run KSTAR on non-human data by matching proteins/phosphosites from the original dataset (such as mouse) to homologous human proteins/phosphosites prior to running KSTAR. The easiest way to do this by using a recent tool that was published by another group, called `PTMoreR <https://doi.org/10.1016/j.crmeth.2024.100859>`_, for mapping phosphorylation sites between species. This tutorial will guide you how to run KSTAR on non-human data. In brief, you will follow these steps (* indicate new steps).


1. *Extract UniProt IDs and peptide sequence for upload to PTMoreR
2. *Homology map to homologous human proteins with PTMoreR
3. *Append human IDs and peptide sequence to original dataset
4. Process the dataset as normal (mapping to reference phosphoproteome, etc.)
5. Run through KSTAR
6. Analyze


Convert IDs/peptide sequence to human homologs with PTMoreR
--------------------------------------------------------------

The next step is to convert the original dataset (for example, mouse) to human homologs. One way this can be done is with the `PTMoreR tool <https://yslproteomics.shinyapps.io/PTMoreR/>`_.



1. In most cases, you will need to reformat the peptide sequence into a format expected by PTMoreR (needs to have a specific label for the modification, '#' by default). For example, if you have peptide sequences with phosphosites indicated by lowercase letters, you can use the following function to extract the necessary information and reformat the peptide sequence for PTMoreR, 

    .. code-block:: python

        def format_peptide_for_ptmorer(peptide):
            # Replace lowercase letters with uppercase and add '#' after the phosphorylated residue
            formatted_peptide = ''
            for char in peptide:
                if char.islower():
                    formatted_peptide += char.upper() + '#'
                else:
                    formatted_peptide += char
            return formatted_peptide

        def get_data_for_ptmorer(df, id_col='accession_id', peptide_col='peptide_sequence'):
            # Extract UniProt IDs and format peptide sequence for PTMoreR
            df['Uploaded.Peptides'] = df[peptide_col].apply(format_peptide_for_ptmorer)

            #rename id column to 'Uniprot.ID' for PTMoreR
            df.rename(columns={id_col: 'UniProt.ID'}, inplace=True)

            return df[['UniProt.ID', 'Uploaded.Peptides']]

        df = get_data_for_ptmorer(df)
        df.to_csv('data_for_ptmorer.csv', index=False)

2. Navigate to the `PTMoreR tool <https://yslproteomics.shinyapps.io/PTMoreR/>`_.
3. Follow the PTMoreR instructions until Step 3
    1. Upload the formatted dataset, indicating how modifications are labeled and which species the data is from
    2. On Step 2, click calculate. You don't need to check the box "Check if containing some regular sequence"
    3. Once Step 2 completes, move to Step 3. Specify the parameters for homology mapping/alignment. We recommend using the default parameters, except:
        - Set 'Central amino acid matching degree' to Exact Matching
        - Check "Whether setting BLOSUM50 score"
    4. Click calculate
    5. Once Step 3 completes, download the results table
5. Check how successful the mapping was:

    .. code-block:: python
        
        # Load the original dataset and the PTMoreR results
        data = pd.read_csv('original_dataset.csv')
        ptmorer_results = pd.read_csv('ptmorer_results.csv')

        print('Fraction of sites successfully mapped:', ptmorer_results.shape[0]/data.shape[0])


4. Append the human UniProt IDs and peptide sequences to the original dataset, matching on the original UniProt ID and peptide sequence. You can use the following code to do this:

    .. code-block:: python

        # Load the original dataset and the PTMoreR results
        data = pd.read_csv('original_dataset.csv')
        ptmorer_results = pd.read_csv('ptmorer_results.csv')

        #process the PTMoreR results
        #grab site information
        ptmorer_results['Human Site'] = ptmorer_results['Center.amino.acids.Other'] + ptmorer_results['PROindex.from.Other'].astype(str)

        #grab only the columns we need
        ptmorer_results = ptmorer_results[['PRO.from.Database',  'Pep.upload', 'PRO.from.Other', 'Human Site', 'Seqwindows.Other']]
        
        #rename to more informative column names
        ptmorer_results = ptmorer_results.rename(columns = {'PRO.from.Database':'UniProt.ID', 'Pep.upload':'Uploaded.Peptides', 'PRO.from.Other':'Human Accession', 'Seqwindows.Other':'Human Blasted Peptide'})

        
        #lowercase the modified residue for use with KSTAR
        ptmorer_results['Human Blasted Peptide'] = ptmorer_results['Human Blasted Peptide'].apply(lambda x: x[0:7] + x[7].lower() + x[8:])

        # Merge the datasets on the original UniProt ID and peptide sequence
        merged_df = pd.merge(data, ptmorer_results, left_on=['accession_id', 'peptide_sequence'], right_on=['UniProt.ID', 'Uploaded.Peptides'], how='left')

        # Now merged_df contains the original data along with the human UniProt IDs and peptide sequences from PTMoreR

Run KSTAR as normal
-------------------

Once human information has been appended to the original dataset, you can run KSTAR as normal, making sure to specify the appropriate columns containing the human UniProt IDs and peptide sequences (not the original species)










