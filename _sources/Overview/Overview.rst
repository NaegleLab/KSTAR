Overview
=============

KSTAR is a Python package for inferring kinase actiivities from phosphoproteomic data. It requires the following steps, given a phosphoproteomic dataset of interest.

1. :doc:`Network Generation <../Advanced/Pruning>`: This produces an ensemble of binary, heuristically pruned kinase-substrate graphs to be used in subsequent analyses. For most users, pre-generated networks based on NetworKIN can be installed (see :doc:`Getting Started<../Tutorial/quickstart>`).

2. :doc:`Mapping Dataset to Reference Phosphoproteome<../Tutorial/Map_Datasets>`: Given that the proteome is regularly updated, this step maps phosphopeptides identified in the dataset to the reference phosphoproteome to ensure site positions match the kinase-substrate networks.

3. :doc:`Kinase Activity Calculation<../Tutorial/Activity_Calculation>`: Given sites identified in an experiment, for each experiment (column) in a dataset, calculate the likelihood that those sites were pulled randomly from the kinase-substrate networks. 
4. :doc:`Analysis/Plotting<../Tutorial/plotting_tutorial>`: Various plotting and analysis functions are provided to visualize and interpret KSTAR results.

For more details about the algorithm, as well as best use cases, see the original publication: `KSTAR Paper`_



.. _KSTAR Paper: https://www.nature.com/articles/s41467-022-32017-5

