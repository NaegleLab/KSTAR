# KSTAR Nextflow pipeline

## parameters
- experiment_file : full filename of experiment
- name : name of expereiment
- phospho_event : phosphorylation event, choices: (Y, ST)
- data_columns : list of pre-mapped data columns to analyze
    - cannot include parentheses, spaces, or commas in name
    - if blank then all columns that start with data: are analyzed
    - data columns are split by comma without spaces between
- add_data_before : (boolean) whether data columns has data: in name
    - if true then each column adds 'data:' before name
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
 - config Files to use 
    - Y : tyrosine config
    - ST : serine-threonine config
    - docker : use docker container
    - singularity : use singularity container
### Local
- start new conda environment
- install nextflow into conda environment 
    - `conda install nextflow`
- run nextflow
    ```
    nextflow run main.nf -profile Y,docker \
    --name example \
    --phospho_event Y \
    --experiment_file data/example_data_mapped.tsv \
    --data_columns data:EOE,data:PRE,data:HDP3,data:HDP6 \
    --num_random_experiments 150 \
    --outdir ./results \
    --resource_directory ../RESOURCE_FILES \
    --network_directoyr /NETWORKS/NetworKIN \
    --threshold 0.5 \
    --activity_aggregate mean \
    --fpr_alpha 0.05 \
    --number_of_sig_trials 100

    ```


### Rivanna
```
module load anaconda/2020.11-py3.8
conda create -n nextflow
source activate nextflow
conda install -c bioconda nextflow
```
- run nextflow
    ```
    nextflow run main.nf -profile Y,singularity \
    --name example \
    --phospho_event Y \
    --experiment_file data/example_data_mapped.tsv \
    --data_columns data:EOE,data:PRE,data:HDP3,data:HDP6 \
    --num_random_experiments 150 \
    --outdir ./results \
    --resource_directory ../RESOURCE_FILES \
    --network_directoyr /NETWORKS/NetworKIN \
    --threshold 0.5 \
    --activity_aggregate mean \
    --fpr_alpha 0.05 \
    --number_of_sig_trials 100

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

