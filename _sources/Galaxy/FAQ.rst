Frequently Asked Questions about Tool
=====================================

KSTAR Parameters
-----------------------------------------

What does it mean to use "All Sites" as evidence?
    This means that all sites with numeric, non-NaN, and non-zero values in the sample column will be used as evidence for activity. Importantly, this doesn't necessarily mean that it uses all rows in the dataset. This may be desired in single-sample contexts (like tissue-specific analysis) where you want to use all of the information available in the sample to get a general, perturbation independent view of kinase activity. You can also use this approach if you just have a set of sites for which you want to perform statistical enrichment analysis to see which kinases are most likely to be associated with these sites. 


KSTAR Outputs
-------------

The tool outputted a table of activities, but they don't include all of my sample columns. What happened?
    The KSTAR tool will only generate predictions for samples that contain a certain amount of phosphorylation sites (50-1000 sites for tyrosine kinases, 500-20000 sites for serine/threonine kinases). There are two main reasons for this restriction:
    1. We have observed these to the ranges for which KSTAR will generate the most reliable and meaningful predictions.
    2. The KSTAR tool uses a pregenerated random activity distribution for predictions, of which need to be comparable to the number of sites in the sample. If there are too few sites, the random activity distribution will not be comparable and predictions cannot be made. If there are too many sites, the random activity distribution will not be comparable and predictions cannot be made. For memory reasons, we have restricted this to a range that we felt was reasonable for most datasets.

    Should you wish to run KSTAR on samples with any number of sites, this is possible to do so via python (see :doc:`Python Quickstart<python-quickstart>` for more details). 

Other Questions
-----------------

My dataset contains non-human data, can I still run KSTAR?
    KSTAR is primarily designed for human data, and the activity predictions are based on human kinase-substrate relationships. If you have non-human data, you may need to convert these to their human homologs before running KSTAR analysis. See the :doc:`tutorial for converting non-human data for use with KSTAR <../Advanced/nonhuman_data>` for more details on how to do this. You will want to use the predictions with caution, as the predictions will be based on human kinase-substrate relationships and may not be as accurate for non-human data.

