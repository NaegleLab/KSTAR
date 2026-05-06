.. _data_prep:

Preparing your Data
===================

Requirements for Input Data
---------------------------

KSTAR can be run on any phosphoproteomics pipeline (DDA, DIA, targeted, etc.) or quantification type. As input, KSTAR uses processed phosphoproteomic data in the form of a pandas dataFrame with the following columns:

+-------------------------------------+------------------------------------------------------+----------------+
| Column(s) Description               | Required                                             |Example Value   |
+=====================================+======================================================+================+
| UniProt Accession                   | Yes                                                  | P00533         |
+-------------------------------------+------------------------------------------------------+----------------+
| Sample data/quantification          | Yes                                                  | 0.5            |
+-------------------------------------+------------------------------------------------------+----------------+
| Peptide Sequence                    | Recommended, required if site position not provided  | AENAEyLRVAPQS  |
+-------------------------------------+------------------------------------------------------+----------------+
| Site Position                       | If peptide sequence not provided                     | Y1197          | 
+-------------------------------------+------------------------------------------------------+----------------+

.. warning::
    KSTAR expects a specific format for the peptide sequence (lowercased residue indicating the phosphorylated residue) and site position (residue followed by position). If your data is not in this format, you can use some of KSTAR's built-in functions to reformat these into the expected format. See the next section for more details on how to prepare your data for KSTAR analysis.


Preparing your Data for KSTAR analysis
---------------------------------------

Before running KSTAR, go through the following checklist to see if your data is ready for KSTAR analysis

- [ ] Do all IDs convert to UniProt accessions? 
- [ ] Does your data have at least one column with peptide sequences or site positions?
- [ ] If using, do the peptide sequences look correct, with the phosphorylated residue shown as a lowercase letter?
- [ ] Does each row contain only one peptide sequence or site position?
- [ ] Do the site positions look correct, with the format of a single letter followed by a number (e.g., Y1197)?
- [ ] Do the data columns only contain numeric values (for example, log fold change, intensity, etc.)?

If the answer is yes to all these questions, your good to proceed to KSTAR activity calculation! If not, you can use some of KSTAR's built-in functions to help you prepare your data for KSTAR analysis. See the next sections for more details on how to prepare your data for KSTAR analysis.

Reformatting Peptide Sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KSTAR requires peptide sequences with the phosphorylated residue shown as a lowercase letter (e.g., ACDEsGHIK). If your peptides use a different format, use :func:`format_peptide_from_df()<kstar.dataset_processing.peptides.format_peptide_from_df>` to convert them. This function detects and reformats peptide sequences automatically.

Before running, ensure each row contains only one peptide. If not, split multiple peptides into separate rows first.

.. code-block:: python

    from kstar.dataset_processing.peptides import format_peptide_from_df

    #split multiple peptides in a single row into individual rows
    peptide_sep = ';' 
    df[peptide_column] = df[peptide_column].str.split(peptide_sep)
    df = df.explode(peptide_column).reset_index(drop=True)

    #format peptides
    df = format_peptide_from_df(df, peptide_column='peptide_sequence', separator=';')

This will add new columns to the dataframe with the reformatted peptide sequence ('Formatted Peptide'), as well as the detected format of the original peptide sequence ('Detected Peptide Format'). 


.. warning::
    Please compare the original and reformatted peptide sequences to ensure that the function is working as expected. If there are any issues with the reformatted peptide sequences, you may need to investigate the original peptide format further and reformat these manually.


Converting IDs to UniProt accessions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To convert from other ID types to UniProt accessions, use the :func:`automatic_id_conversion()<kstar.dataset_processing.accessions.automatic_id_conversion>` function. This function will automatically identify the most common ID type used in the dataset, and convert these to UniProt accessions. 

To run this function, run the below code, specifying the column containing the IDs and any separators used if multiple IDs are present in a single row

.. code-block:: python

    from kstar.dataset_processing.accessions import automatic_id_conversion

    df, missed = automatic_id_conversion(df, id_column='accession_id', separator=';')

This will return a new dataframe with the converted UniProt accessions, as well as a list of any IDs that could not be converted. 

The following types of IDs are currently supported for conversion:

- UniProt entry name (for example, "EGFR_HUMAN")
- Ensembl gene, transcript, or protein IDs (for example, "ENSG00000146648", "ENST00000275493", "ENSP00000275493")
- RefSeq mRNA or protein IDs (for example, "NM_005228", "NP_005219")
- GI numbers (for example, "gi|123456789")
- CCDS IDs (for example, "CCDS12345")
- Gene symbols (for example, "EGFR")
- UniParc IDs (for example, "UPI0000000001")
- UniRef IDs 

.. warning::
    If your dataset contains a large number of IDs that cannot be converted, you may want to investigate these further to see if they can be manually converted or if they indicate an issue with the dataset. KSTAR will not be able to make predictions for any sites associated with IDs that cannot be converted to UniProt accessions.



FAQ
-----------------

* What if I have a different type of data, such as protein-level data or non-phosphoproteomic data?

    KSTAR is designed for phosphoproteomic data, and the activity predictions are based on the presence of specific phosphosites. If you have protein-level data or non-phosphoproteomic data, KSTAR may not be the best tool for your analysis. 

* Can a peptide sequence have multiple phosphorylated residues?

    Yes, a peptide sequence can have multiple phosphorylated residues. However, for KSTAR analysis, each row should contain only one peptide sequence.

* What if I have nonhuman UniProt accessions?

    KSTAR is primarily designed for human data, and the activity predictions are based on human kinase-substrate relationships. If you have non-human UniProt data, you may need to convert these to their human homologs before running KSTAR analysis. See the :doc:`tutorial for converting non-human data for use with KSTAR <../Advanced/nonhuman_data>` for more details on how to do this.





