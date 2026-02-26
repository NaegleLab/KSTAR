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

## Pushing docker image to Docker Hub

For galaxy to access the docker image, we need to make it accessible on docker hub. Check the name of the image you just created:
```
docker images
```
If it does not begin the name "naeglelab/kstar-galaxy" with the tag associated with the current KSTAR version (ex. "1.1.0"), you'll want to update the tag name to match the global repository. Do this below:
```
docker tag kstar-galaxy naeglelab/kstar-galaxy:1.1.0
```
Then push to docker hub:
```
docker push naeglelab/kstar-galaxy:1.1.0
```

## Setting up local Galaxy instantiation

To test Galaxy version of KSTAR (located at -------):
1. First need to download and setup Galaxy. See instructions for downloading and installing on Galaxy web site.
2. To make the KSTAR tool, make a new directory in the 'tools' folder of the local galaxy repository called KSTAR. Within this folder, add the .xml and .py files from [KSTAR_Galaxy repository](https://github.com/NaegleLab/KSTAR_Galaxy).
3. Navigate to the config folder. You should find a file called "tool_conf.xml.sample". Create a copy of this file and remove the ".sample" part of the file name. To this file, add the KSTAR tool after the toolbox tag:
```xml
  <section id='kstar' name='KSTAR Tools'>
    <tool file='KSTAR/calculate.xml' />
  </section>
```
4. Finally, since this tool uses docker, you will need to configure use of docker. Create a new file called 'job_conf.yml' and paste the following lines:
```yml
runners:
  local:
    load: galaxy.jobs.runners.local:LocalJobRunner
    workers: 4

execution:
  default: docker_dispatch
  environments:
    local:
      runner: local

    docker_local:
      runner: local
      docker_enabled: true
      docker_set_user:

      # InteractiveTools do need real hostnames or URLs to work - simply specifying IPs will not work.
      # If you develop interactive tools on your 'localhost' and don't have a proper domain name
      # you need to tell all Docker containers a hostname where Galaxy is running.
      # This can be done via the add-host parameter during the `docker run` command.
      # 'localhost' here is an arbitrary hostname that matches the IP address of your
      # Galaxy host. Make sure this hostname ('localhost') is also set in your galaxy.yml file, e.g.
      # `galaxy_infrastructure_url: http://localhost:8080`.
      # docker_run_extra_arguments: --add-host localhost:host-gateway --platform linux/amd64

    docker_dispatch:
      runner: dynamic
      type: docker_dispatch
      docker_destination_id: docker_local
      default_destination_id: local
```
This will let you use Docker images in the requirements tag



## Testing and updating Galaxy app

To test the galaxy app, set up tests with the test tag in the xml file, and add any input/output datasets for the test under a folder called 'test-data' in the same repo. You can then test the tool using `planemo test --docker calculate.xml`, replacing calculate.xml with whatever your xml file is called.

