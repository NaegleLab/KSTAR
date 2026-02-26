Troubleshooting the KSTAR Tool
==============================

Critical Errors (no activities predicted)
-----------------------------------------

Memory Error/Killed Job
~~~~~~~~~~~~~~~~~~~~~~~
Large datasets may cause the tool to exceed memory limits for running on the Galaxy server, which will cause tool to be killed and fail. If you encounter this, we recommend splitting up your data columns into separate csv/tsv files and running in batch (data columns are run independently, so this will not effect results).

Parameter Error (Exit Code 2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This error is typically due to an issue with the parameters selected for the KSTAR run and will cause tool failure. The most common parameter issues are:
- Not selecting at least one data column
- Not selecting the appropriate columns for mapping (at least one of "site_column" or "peptide_column" must be selected

Usually, this error is an easy fix and just requires going back to the KSTAR tool form and adjusting the parameters. 

Dataset Error (Exit Code 3)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This error is caused by issues with the uploaded dataset and will usually cause tool failure. These can be broken into two categories:

1. Issues with the dataset format (e.g. peptides or sites not in the correct format)
2. Too few sites in the dataset to generate a KSTAR run. This is most common for tyrosine kinase runs (especially for global phosphoproteomics). If you encounter this error, we don't recommend running KSTAR on this dataset as the results will likely not be meaningful given the small number of sites. If you want to run KSTAR on this dataset anyway, this is possible to do so via python (see :doc:`Python Quickstart<python-quickstart>` for more details).


Mapping Error (Exit Code 4)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
This error is caused by issues during the mapping process, most commonly resulting from a large number of unmapped sites. This can be similar to the dataset error of too few sites, but in this case the issue is that the sites were not able to be mapped to the KSTAR reference phosphoproteome. This is most commonly due to incorrect accessions (not human SwissProt IDs) or incorrect peptide/site formats not caught during argument processing. If you encounter this error, we recommend inspecting the mapping stats file to see how many sites were mapped vs. lost during the mapping process, and the reasons for this. 

Threshold/Evidence Check Error
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you encounter an error related to thresholding, this is likely due to the fact that the there are either too few or too many sites used as evidence. Below shows the bounds for the number of sites used as evidence for the Galaxy tool:

    +------------------+--------------+-------------+
    | Modification     | Min. Sites   |  Max Sites  |
    +==================+==============+=============+
    | Tyrosine         | 50           | 1000        |
    +------------------+--------------+-------------+
    | Serine/Threonine | 500          | 20000       |
    +------------------+--------------+-------------+

If using a threshold, take a look at the threshold plot provided in the KSTAR outputs to see the impact of different thresholds. You want to choose a threshold that still maintains greater than 50 sites (for Y) and greater than 500 (for ST) to have enough information for activity calculation. 

If using either the Top N, Bottom N, or All sites option, make sure to select a top N value that is less than the maximum number of sites allowed (1000 for Y and 20000 for ST). 


Non-Critical Warnings (activities predicted, but potential issues with results))
--------------------------------------------------------------------------------

It's possible that activities are calculated, but some warnings are present in the error log. These warnings may be indicative of potential issues with the results, but not necessarily. We recommend inspecting these warnings to see if they are indicative of a larger issue with the results. 

Below are some things to be on the lookout for:

Poor Mapping to Reference Phosphoproteome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The tool will raise an warning if less than 70\% of sites/peptides are successfully mapped to the reference phosphoproteome (i.e. matching phosphosite was found in our reference), even if the tool doesn't fail. If you see warnings about poor mapping to the reference phosphoproteome, this may be indicative of an issue with your dataset (e.g. incorrect accessions or site formats). This is not necessarily a problem -- for example, if your data includes non-SwissProt UniProt IDs (i.e non-reviewed entries), these will be dropped from your data. This should not effect results and is fine. 

We recommend inspecting the mapping stats file to see how many sites were mapped vs. lost during the mapping process, and the reasons for this. If you have very poor mapping to the reference phosphoproteome, we don't recommend using these results as they will likely not be meaningful. 

Dropped Columns due to Improper Evidence Size
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In some cases, you may not have all of your data columns included in the final KSTAR run. If the number of sites used as evidence for a particular sample column is too small (less than 50 for Y and less than 500 for ST), that column will be dropped from the run and a warning will be raised. If this is a problem, you may consider adjusting your evidence selection (e.g. using a less stringent threshold or using the Top N sites option instead) to include more sites as evidence and avoid excluding your sample columns from the run.

*Note: if using automatic thresholding, KSTAR will prioritize using a threshold value that does not drop any columns, but in some cases this may not be possible due to few sites in that column. In this case, KSTAR will pick a threshold that achieves the desired evidence size, even this means dropping some columns.*


