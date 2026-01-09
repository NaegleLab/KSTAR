KSTAR Reference
=============================
.. automodule:: kstar 

The "Config" Module
--------------------

.. automodule:: kstar.config
	:members: install_network_files, update_network_directory, check_configuration, update_configuration, find_available_networks, get_package_memory
	 
The "Prune" Module
------------------

The "Pruner" Class
~~~~~~~~~~~~~~~~~~
.. autoclass:: kstar.prune.Pruner
    :members:
	
Functions to Perform Pruning
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: kstar.prune
    :members: run_pruning, save_pruning, save_run_information
	
The "ExperimentMapper" class
----------------------------
.. autoclass:: kstar.mapping.ExperimentMapper
    :members:


Functions for Activity Calculation
----------------------------------

The "KinaseActivity" class
~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: kstar.calculate.KinaseActivity
    :members: check_data_columns, set_data_columns, test_threshold, test_threshold_range, get_allowable_threshold, recommend_threshold, create_binary_evidence, calculate_kinase_activities, get_random_activities, calculate_Mann_Whitney_sig, get_run_information_content, get_param_dict, make_summary_pdf, make_dotplot

Master Functions for Running KSTAR Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: kstar.calculate
    :members: enrichment_analysis, randomized_analysis, Mann_Whitney_analysis, run_kstar_analysis
	
Functions for Saving and Loading KSTAR results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
.. automodule:: kstar.calculate
    :members: save_kstar, from_kstar, from_kstar_nextflow

Plotting/Analysis Functions
---------------------------

The "DotPlot" class
~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: kstar.plot.DotPlot
    :members: 

The "KSTAR_PDF" class
~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: kstar.plot.KSTAR_PDF
    :members: 

Downstream Analysis Modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: kstar.analysis.interactions
    :members: 

.. automodule:: kstar.analysis.coverage
    :members:

.. automodule:: kstar.analysis.kinase_MI
    :members: kinase_mutual_information, plot_kinase_heatmap

Dataset Processing Functions
----------------------------
.. automodule:: kstar.dataset_processing.peptides
    :members: convert_ids, get_valid_from_ids, get_valid_to_ids, get_field_information, check_if_id_supported, identify_accession_type, identify_most_common_accession_type

.. automodule:: kstar.dataset_processing.accessions
    :members: detect_peptide_format, detect_most_common_format, format_peptide, format_peptide_list, format_peptide_from_df


Other Helper Functions
----------------------
.. automodule:: kstar.helpers
    :members:
	
