# Getting Started

Here, we have provided a quick start guide that will allow you to get up and running quickly, especially if you wish to simply use the
same pruned networks as utilized in the KSTAR paper. 

## Installation

KSTAR can be installed via `pip`, `conda`, tarball, and directly from the Git repository. We recommend using conda forge to install the scientific packages `conda install conda-forge` before installing KSTAR.

#### Pip

To install via pip, execute `pip install kstar`.

#### Conda Forge
Please note: the conda-forge build of KSTAR is now required. Versions up to and including 1.0.3 are
deprecated. Install the latest package from conda-forge: [Conda-Forge KSTAR](https://anaconda.org/conda-forge/kstar).

To install via conda-forge, execute `conda install -c conda-forge kstar`.

#### Tarball

To install via a tarball, head over to the [Releases page](https://github.com/NaegleLab/KSTAR/releases) and download the latest stable tar release.

Afterwards, navigate to your downloads directory and execute the following commands, substituting <version> for the release's version number:
```
tar -xvf KSTAR-<version>.tar.gz
cd KSTAR-<version>
python setup.py install
```

#### Git

If you want to try out the latest commit, you can install directly from the Git repository by executing the following commands:
```
git clone https://github.com/NaegleLab/KSTAR
cd KSTAR
python setup.py install
```

## Configuring your KSTAR environment

After installing KSTAR, all necessary resource files (reference proteome and phosphoproteome) and networks (either downloaded from
FigShare or generated with the Pruner class) will need to be downloaded and configured so that KSTAR can find these files.

### Download Networks
#### for KSTAR version 1.0.4 and newer
```python
#For KSTAR versions 1.0.4 and newer, KSTAR will automatically download the network files using the code below.
#Please specifiy your working directory for the targer_dir parameter. Example directory shown below.  
config.install_network_files(target_dir="/Users/naeglelab/Documents/KSTAR")
```

#### for KSTAR version 0.5.3 and older
Note: If you are using KSTAR version 0.5.3 or older, please use the guide below to download the networks.

In addition to the resource files above, KSTAR also requires heuristically pruned kinase-substrate graphs used for activity calculation.
Pre-generated networks are available for download, which were generated based on NetworKIN. You may also generate your own networks if preferred.

If using the published networks:

1. Go to [Network Figshare](https://figshare.com/articles/dataset/NETWORKS/14944305)
2. Download the networks, decompress/unzip the files, and store in easily accessible folder.

#### If Using Self-generated Networks

1. Identify and download the weighted kinase-substrate prediction graph you would like to use for network generation. This should include weighted edges indicating the likelihood that a kinase phosphorylates a particular substrate, and should include predictions for all sites in the phosphoproteome.
2. Follow the tutorial for network generation to produce pruned networks from your weighted network (found in Tutorial section of this documentation).
3. Store the generated networks in an easily accessible location. Individual network files should be placed in a directory within the network directory in a folder titled 'INDIVIDUAL_NETWORKS'. 

### Configure KSTAR to point to the correct network directory
*skip if you used `install_network_files()` function

Once the networks have been downloaded/generated, the last step is to tell KSTAR where to find these networks. Use the *update_network_directory()* function in config.py to tell config where the network directory is located. On install, KSTAR is set to look in the './NETWORKS/NetworKIN'.

The python code used should look similar to below:
```python
from kstar import config

#update network directory: If KSTAR does not find this directory + necessary files, it will notify you
config.update_network_directory('/path/to/NETWORKS_folder/')
```

### Verify that KSTAR environment is ready

To check to make sure the previous steps all worked as desired, run *check_configuration()*:
```
from kstar import config

config.check_configuration()
```
This will indicate whether you are ready to generate kinase activity predictions or not. If you are not, it will tell you what steps still need to
be performed. 

## Follow the provided tutorial

That is all you need to do to set up your KSTAR environment (if working with large datasets where memory/time is a concern, see 'KSTAR in Parallel', as set up is slightly different). We recommend working through the tutorial in the following section to get an idea of the KSTAR workflow, either with the example dataset provided in our supplementary data figshare or with your own dataset of interest. You will only need to follow the 'Network Generation' section if you would like to use your own networks for analysis, otherwise go straight to 'Activity Prediction'.

## Updating Default Configuration

If desired, you can alter a number of default options used by KSTAR. To update these, you can use the `config.update_configuration` function. Most of these are related to how to handle pregenerated random experiments, which were introduced in KSTARv1.0.4:


| Config Parameter                        | Global Name                       | Description                                                                     | Default Value |
|-----------------------------------------|-----------------------------------|---------------------------------------------------------------------------------|---------------|
| network_dir                             | NETWORK_DIR                       | where all network files/random experiments are stored                           |./NETWORKS/    |
| y_network_name                          | NETWORK_NAME['Y']                 | default name of tyrosine network within NETWORK_DIR to use                      |Default        |
| st_network_name                         | NETWORK_NAME['ST']                | default name of serine/threonine network within NETWORK_DIR to use              |Default        |
| use_pregenerated_random_activities      | USE_PREGENERATED_RANDOM_ACTIVITIES| boolean, whether to use pregenerated random activities                          |True           |
| save_random_experiments                 | SAVE_RANDOM_EXPERIMENTS           | boolean, whether to save random experiments when generated                      |False          |
| save_new_random_activities              | SAVE_NEW_RANDOM_ACTIVITIES        | boolean, whether to save random activities when generated from scratch          |False          |
| custom_pregenerated_activities_dir      | CUSTOM_RANDOM_ACTIVITIES_DIR      | file path to where new random activities are stored and accessed in future runs |None           |