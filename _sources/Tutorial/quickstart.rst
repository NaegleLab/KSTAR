.. _python-quickstart:

Getting Started
===============

Here, we have provided a quick start guide that will allow you to get the KSTAR python package up and running quickly with the default settings.

Installation
------------

KSTAR can be installed via ``pip``, ``conda``, tarball, and directly from the Git repository. We recommend using conda forge to install the scientific packages::

    conda install conda-forge

before installing KSTAR.

Pip
~~~
To install via pip, execute::

    pip install kstar

Conda Forge
~~~~~~~~~~~
Please note: the conda-forge build of KSTAR is now required. Versions up to and including 1.0.3 are deprecated. Install the latest package from conda-forge: `Conda-Forge KSTAR <https://anaconda.org/conda-forge/kstar>`_.

To install via conda-forge, execute::

    conda install -c conda-forge kstar

Tarball
~~~~~~~
To install via a tarball, head over to the `Releases page <https://github.com/NaegleLab/KSTAR/releases>`_ and download the latest stable tar release.

Afterwards, navigate to your downloads directory and execute the following commands, substituting <version> for the release's version number::

    tar -xvf KSTAR-<version>.tar.gz
    cd KSTAR-<version>
    python setup.py install

Git
~~~
If you want to try out the latest commit, you can install directly from the Git repository by executing the following commands::

    git clone https://github.com/NaegleLab/KSTAR
    cd KSTAR
    python setup.py install

Configuring your KSTAR environment
----------------------------------

After installing KSTAR, all necessary resource files (reference proteome and phosphoproteome) and networks (either downloaded from FigShare or generated with the Pruner class) will need to be downloaded and configured so that KSTAR can find these files.

Download Networks
~~~~~~~~~~~~~~~~~

for KSTAR version 1.0.4 and newer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For KSTAR versions 1.0.4 and newer, KSTAR will automatically download the network files using the code below. Please specify your working directory for the ``target_dir`` parameter. Example directory shown below::

    config.install_network_files(target_dir="/Users/naeglelab/Documents/KSTAR")

for KSTAR version 0.5.3 and older
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Note: If you are using KSTAR version 0.5.3 or older, please use the guide below to download the networks.

In addition to the resource files above, KSTAR also requires heuristically pruned kinase-substrate graphs used for activity calculation. Pre-generated networks are available for download, which were generated based on NetworKIN. 

    1. Go to `Network Figshare <https://figshare.com/articles/dataset/NETWORKS/14944305>`_
    2. Download the networks, decompress/unzip the files, and store in an easily accessible folder.



Configure KSTAR to point to the correct network directory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*skip if you used ``install_network_files()`` function*

Once the networks have been downloaded/generated, the last step is to tell KSTAR where to find these networks. Use the :func:`update_network_directory()<kstar.config.update_network_directory>` function in config.py to tell config where the network directory is located. On install, KSTAR is set to look in the './NETWORKS/NetworKIN'.

The python code used should look similar to below::

    from kstar import config

    # update network directory: If KSTAR does not find this directory + necessary files, it will notify you
    config.update_network_directory('/path/to/NETWORKS_folder/')

Verify that KSTAR environment is ready
--------------------------------------

To check to make sure the previous steps all worked as desired, run :func:`check_configuration()<kstar.config.check_configuration>`::

    from kstar import config

    config.check_configuration()

This will indicate whether you are ready to generate kinase activity predictions or not. If you are not, it will tell you what steps still need to be performed.

Follow the provided tutorial
----------------------------

That is all you need to do to set up your KSTAR environment (if working with large datasets where memory/time is a concern, see 'KSTAR in Parallel', as set up is slightly different). We recommend working through the tutorial in the following section to get an idea of the KSTAR workflow, either with the example dataset provided in our supplementary data figshare or with your own dataset of interest. You will only need to follow the 'Network Generation' section if you would like to use your own networks for analysis, otherwise go straight to 'Activity Prediction'.

