.. _data_prep:

Preparing your Data
===================

As input, KSTAR uses processed phosphoproteomic data in the form of a pandas DataFrame. KSTAR can be run on any phosphoproteomics pipeline (DDA, DIA, targeted, etc.) or quantification type. 

This dataset must contain the following information:

- **accession_id**: The UniProt accession ID of the protein containing the phosphosite.
- **data_columns**: One or more columns containing quantitative data for each phosphopeptide/site. These columns should be numeric and can represent any quantitative measure, such as intensity, fold change, or log2 fold change, binary, etc.

and must contain at least one of the following:

- **site_position**: The residue and position of the phosphorylated residue within the protein sequence (for example, Y1197)
- **peptide_sequence**: The peptide sequence containing the phosphorylated residue, with the phosphorylated residue indicated by a lowercase letter (for example, ACDEsGHIK for a phosphorylated serine)

Before proceeding, ensure that your data is in this format (correct ID type, correct site notation or peptide format, and numeric data columns). If your data is not in this format, please reformat it accordingly before running KSTAR.

