# KSTAR instantiation for Galaxy

This directory contains all files needed to create a docker image of KSTAR for use with Galaxy. This largely uses the same code as base KSTAR, but omits certain modules that are not necessary (prune, pregenerate, etc.). In addition, it uses a unique config file that specifically points to the resource files within the docker image.

## Creating Docker image
To create this image, follow these steps:

1. Copy the kstar/ directory into the galaxy folder. Move the unique config.py file located in the galaxy directory into this directory.
2. Copy the most up-to-date version of the RESOURCE_FILES and NETWORKS directories into the galaxy folder. After these steps, your galaxy directory should look similar to this (**indicates unique to galaxy configuration):
```
KSTAR/galaxy
├───kstar
│   ├───<default_kstar_modules>
│   └───config.py**
├───NETWORKS
│   ├───ST
│   │   └───Default
│   │       ├───INDIVIDUAL_NETWORKS
│   │       └───RANDOM_ACTIVITIES
│   └───Y
│       └───Default
│           ├───INDIVIDUAL_NETWORKS
│           └───RANDOM_ACTIVITIES
└───RESOURCE_FILES
```
3. To make docker image, run  `docker build -t kstar-galaxy .` (if on Windows, it's easier to do this in a WSL2 environment)
4. Check that image was created and its size using ``docker images`
5. You can check to make sure the docker image successfully installed kstar and resources with a quick import check in the docker image:
```
docker run -it kstar-galaxy python -c "from kstar import config, calculate
```
You can also run a specific script with something like this, where your script is located in your current directory
```
docker run -it -v "$PWD":/app -w /app kstar-galaxy_trim:latest python calculate.py --input_file 'Asmussen_2013_Dasatinib.csv' --mapped 'No' --accession_column 1 --site_column 3 --peptide_column 2 --data_columns 4 --phospho_type 'Y' --threshold 0.2 --evidence_size 100
```
Change the `C:\path\to\your\script.py` to your path and update the script name. This mounts the script within the docker container and then runs it

## Setting up local Galaxy instantiation

To test Galaxy version of KSTAR (located at -------), you first need to download and setup Galaxy.


## Testing and updating Galaxy app

