Updating Configuration
======================


If desired, you can alter a number of default options used by KSTAR. To update these, you can use the :func:`update_configuration<kstar.config.update_configuration>` function. Most of these are related to how to handle pregenerated random experiments (see :doc:`Pregenerated_Experiments<Advanced/Pregenerated_Experiments>`), which were introduced in KSTAR v1.0.4:

+------------------------------------------+-----------------------------------+---------------------------------------------------------------------------------+---------------+
| Config Parameter                         | Global Name                       | Description                                                                     | Default Value |
+==========================================+===================================+=================================================================================+===============+
| network_dir                              | NETWORK_DIR                       | where all network files/random experiments are stored                           | ./NETWORKS/   |
+------------------------------------------+-----------------------------------+---------------------------------------------------------------------------------+---------------+
| y_network_name                           | NETWORK_NAME['Y']                 | default name of tyrosine network within NETWORK_DIR to use                      | Default       |
+------------------------------------------+-----------------------------------+---------------------------------------------------------------------------------+---------------+
| st_network_name                          | NETWORK_NAME['ST']                | default name of serine/threonine network within NETWORK_DIR to use              | Default       |
+------------------------------------------+-----------------------------------+---------------------------------------------------------------------------------+---------------+
| use_pregenerated_random_activities       | USE_PREGENERATED_RANDOM_ACTIVITIES| boolean, whether to use pregenerated random activities                          | True          |
+------------------------------------------+-----------------------------------+---------------------------------------------------------------------------------+---------------+
| save_random_experiments                  | SAVE_RANDOM_EXPERIMENTS           | boolean, whether to save random experiments when generated                      | False         |
+------------------------------------------+-----------------------------------+---------------------------------------------------------------------------------+---------------+
| save_new_random_activities               | SAVE_NEW_RANDOM_ACTIVITIES        | boolean, whether to save random activities when generated from scratch          | False         |
+------------------------------------------+-----------------------------------+---------------------------------------------------------------------------------+---------------+
| custom_pregenerated_activities_dir       | CUSTOM_RANDOM_ACTIVITIES_DIR      | file path to where new random activities are stored and accessed in future runs | None          |
+------------------------------------------+-----------------------------------+---------------------------------------------------------------------------------+---------------+