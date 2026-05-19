# KSTAR instantiation for Galaxy

This directory contains all files needed to create a docker image of KSTAR for use with Galaxy. This largely uses the same code as base KSTAR, but omits certain modules that are not necessary (prune, pregenerate, etc.). In addition, it uses a unique config file that specifically points to the resource files within the docker image.

- [Current version of KSTAR tool](https://toolshed.g2.bx.psu.edu/view/naegle_lab/kstar_calculate/c1eae1b6178d)
- [Galaxy web server](https://usegalaxy.org/)
- [Galaxy XML syntax](https://docs.galaxyproject.org/en/latest/dev/schema.html#)

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

To test Galaxy version of KSTAR, follow the below instructions:
1. First need to download and setup Galaxy. See instructions for downloading and installing on Galaxy web site.
2. To make the KSTAR tool, make a new directory in the 'tools' folder of the local galaxy repository called KSTAR. Within this folder, add the .xml and .py files from [KSTAR_Galaxy repository](https://github.com/NaegleLab/KSTAR_Galaxy).
3. Navigate to the config folder. You should find a file called "tool_conf.xml.sample". Create a copy of this file and remove the ".sample" part of the file name. To this file, add the KSTAR tool after the toolbox tag:
```xml
  <section id='kstar' name='KSTAR Tools'>
    <tool file='KSTAR/calculate.xml' />
  </section>
```
4. Finally, since this tool uses docker, you will need to configure use of docker. Create a new file called 'job_conf.yml' in the config folder and paste the following lines:
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

## Running Galaxy and Testing KSTAR

To run your local galaxy instantiation, simply open the terminal (WSL distro if on windows) and navigate to galaxy repository. Run 
```
sh run.sh > galaxy.log >2&1
``` 
This will start galaxy on a localhost and send output information to a galaxy.log file in case you would like to inspect.

Once Galaxy is running, simply navigate to the tool and test out the tool (should be under KSTAR tools). If you don't see it, it's likely in issue with the underlying xml file or cheetah code in the command section that's causing the tool to fail to build. H

## Update Galaxy Tool and Pushing to Web Service

If you want to make any updates to the KSTAR tool (or make another tool), you only need to update the xml files and python script associated with the tool. Make sure to test with your local galaxy. Once you are happy with the tool, proceed to the next step to do some quick checks

### Running tests with Planemo

Galaxy uses a package called planemo for a lot of their services, and we can use it to test our tool to make sure it 1) fits formatting compatible with galaxy and 2) runs as expected. Planemo does not work in windows, so you will need to install it in your WSL environment if using windows.

#### Checking format of xml file

First, you can do a quick check to make sure the xml file. With planemo installed, run the following command to check your xml file
```
planemo lint my_tool.xml
```
This will provide any errors (not compatible with Galaxy) or warnings (suggested changes, often missing details or things recommended against). Make sure to fix any of the explicit errors, and consider following any suggestions from the warnings.

#### Testing the tool

To test the galaxy app, set up tests with the test tag in the xml file, and add any input/output datasets for the test under a folder called 'test-data' in the same repo. You can then test the tool using `planemo test --docker calculate.xml`, replacing calculate.xml with whatever your xml file is called.

#### Update tool meta data (.shed.yml file)

The tool meta data for each galaxy tool is stored in a .shed.yml file. This is what will be showed when people look at the tool in the toolshed. An example of one is below:

```yml
name: KSTAR Activity Prediction
owner: naegle_lab
long_description: |
  KSTAR (Kinase-Substrate Transfer to Activity Relationships) is a statistical analysis tool designed to infer kinase activities from phosphoproteomic data in an error- and bias-aware manner. KSTAR can handle most types of phosphoproteomic datasets, as long as they contain phosphorylation site information.
categories: 
  - Proteomics
  - Systems Biology
homepage_url: https://github.com/NaegleLab/KSTAR
remote_repository_url: https://github.com/NaegleLab/KSTAR_Galaxy/tree/main/tools/KSTAR/
auto_tool_repositories:
  name_template: "{{ tool_id }}"
  description_template: "Wrapper KSTAR tool: {{ tool_name }}"
```
Either create this file, or make any changes to an existing one


#### Pushing to the KSTAR_Galaxy repo

Add your updated xml, python, and .shed.yml files to the KSTAR Galaxy repository in the tools folder. Commit the changes and push to the repository. After each commit, checks will automatically be made and if all checks are passed, the updated tool will be deployed on the galaxy toolshed.

#### Creating a PR to add KSTAR tool to main server

The main toolset for usegalaxy servers can be found on https://github.com/galaxyproject/usegalaxy-tools

To add the KSTAR tool (or update), you need to create a pull request, following instructions in the README file of the usegalaxy-tools repo

