# KSTAR Nextflow pipeline

## parameters
- experiment_file : full filename of experiment
- name : name of expereiment
- phospho_event : phosphorylation event, choices: (Y, ST)
- data_columns : list of pre-mapped data columns to analyze. cannot contain : in string
- add_data_before : (boolean) whether data columns has data: in name
- resource_directory : directory of resource files
- network_directory : path relative to resource_directory of networks 
    - individual networks are used in pipeline
- num_random_experiments : number of random experiments to use in normalization and mann-whitney
- threshold : what threshold to use for generating binarized experiment
- activity_aggregate : aggregate for binarizing experiment. choices : (count, mean, median, max, min)
- greater : boolean of whether data above threshold are seen as present in sample
- fpr_alpha : fpr alpha to use in Mann-Whitney
- number_of_sig_trials : number of trials to run in Mann-Whitney
- chunk_size : chunk size to use in hypergeometric - useful for managing memory issues

 
## Running pipeline
- create .config file in /config to include specific parameters of pipeline run
- add config file to nextflow parameters
    - open nextflow.config
    - add name of profile and location to match
        ``` 
        profiles { 
            test { includeConfig 'config/test.config' } 
        }
        ```
### Local
- start new conda environment
- install nextflow into conda environment 
- run nextflow
    ```
    conda install nextflow
    nextflow run main.nf -profile test
    ```


### Rivanna
```
module load nextflow
nextflow run main.nf -profile test
```

## Editing Code in the pipeline
- all edited code must be in the /src directory
- Remake Docker image
    ```
    docker build -t naeglelab/kstar-nextflow:latest .
    docker push naeglelab/kstar-nextflow:latest
    ```
    - if using Rivanna Singularity image must be used - Docker does not work
    - if local then Docker will work, just enable through nextflow.config 
        ```
        docker.enabled = true
        process.container = 'naeglelab/kstar-nextflow:latest'
        ```
- Remake Singularity image
    - make singularity image (from Rivanna)
        ```
        module load singularity/3.5.2
        singularity pull docker://account/image
        ```
    - edit config file to use singularity
        ```
        process.container = '/path/to/singularity.img'
        singularity.enabled = true  
        ```
    

## Results
- Results are stored in output_directory/name/phospho_event
- Results are stored in folder corresponding to the proceess in the pipeline they were generated from 
- Additional pipeline info can be found in pipeline_info
    - svg of pipeline
    - execution timeline
    - execution report : contains info such as memory used, cpu time, etc
    - execution trace  

