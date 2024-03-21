nextflow.enable.dsl=2
include { binarize_experiment } from './modules/binarizeExperiment.nf'
include { hypergeometric_activity } from './modules/hypergeometricActivity.nf'
include { generate_random_experiments } from './modules/generateRandomExperiments.nf'
include { random_hypergeometric_activity } from './modules/randomHypergeometricActivity.nf'
include { summarize_random_aggregated } from './modules/summarizeRandomAggregated.nf'
include { generate_normalization_values } from './modules/generateNormalizationValues.nf'
include { normalize_activity } from './modules/normalizeActivity.nf'
include { summarize_normalized } from './modules/summarizeNormalized.nf'
include { mann_whitney } from './modules/mannWhitney.nf'
include { summarize_mann_whitney } from './modules/summarizeMannWhitney.nf'

/*--------------------------------------------------
Map Data Columns
---------------------------------------------------- 
-if data columns are provided
    -if add_data_before 
        add "data:" before each column name
-else
    -use all columns that start with "data:"
---------------------------------------------------*/

if(params.data_columns != false){
    split_data_columns=params.data_columns?.tokenize(",")
    mapped_data_columns = []
    if(params.add_data_before){
        for(col in split_data_columns){
            mapped_data_columns.add("data:"+col)
        }
    }
    else{
        for(col in split_data_columns){
            mapped_data_columns.add(col)
        }
    }
}
else{
    myFile = file(params.experiment_file)
    myReader = myFile.newReader()

    line = myReader.readLine()
    columns = line.split('\t')
    mapped_data_columns = []
    for(col in columns){
        if(col.startsWith("data:")){
            mapped_data_columns.add(col)
        }
    }
}

mapped_data_column_string = mapped_data_columns.join(" ")
print(mapped_data_column_string)

Channel
    .from(mapped_data_columns)
    .set { ch_mapped_data_columns}

Channel
    .from(mapped_data_column_string)
    .set { ch_mapped_data_column_string }

Channel
    .value(file(params.experiment_file))
    // .ifEmpty { error "Cannot find gtf file for parameter --experiment_file: ${params.experiment_file}" }
    .set { ch_experiment }  

Channel
    .value(file(params.resource_directory + params.compendia))
    // .ifEmpty { error "Cannot find gtf file for parameter --compendia: ${params.compendia}" }
    .set { ch_compendia }  

Channel
    .value(file(params.resource_directory + params.reference_fasta))
    // .ifEmpty { error "Cannot find gtf file for parameter --reference_fasta: ${params.reference_fasta}" }
    .set { ch_reference_fasta} 

Channel
    .fromPath(params.resource_directory + params.network_directory, type:'dir')
    .set{ch_network_directory}

    ch_network_directory.set { ch_network_directory_hypergeometric }
    ch_network_directory.set { ch_network_directory_random }

Channel
    .fromPath("${params.resource_directory}${params.network_directory}/${params.phospho_event}/INDIVIDUAL_NETWORKS/*.tsv")
    .count()
    .set {ch_number_of_networks}

workflow{

/*--------------------------------------------------
Binarize Experiment
---------------------------------------------------*/
BINARIZED = binarize_experiment(ch_experiment, ch_mapped_data_column_string)
ch_binary_experiment = BINARIZED
ch_binary_experiment_for_hypergeometric = ch_binary_experiment
ch_binary_experiment_for_random_experiments = ch_binary_experiment

/*--------------------------------------------------
Kinase Hypergeometric Activity 
---------------------------------------------------*/
HYPERGEOMETRIC_ACTIVITIES = hypergeometric_activity(ch_binary_experiment_for_hypergeometric, ch_network_directory_hypergeometric, ch_mapped_data_column_string)
ch_experiment_activites_files = HYPERGEOMETRIC_ACTIVITIES[0]
ch_experiment_activities = HYPERGEOMETRIC_ACTIVITIES[1]
ch_experiment_activities_list = HYPERGEOMETRIC_ACTIVITIES[2]
ch_experiment_activities_list_for_normalization = ch_experiment_activities_list
ch_experiment_activities_list_for_mann_whitney = ch_experiment_activities_list

/*--------------------------------------------------
Generate Random Experiments 
---------------------------------------------------*/
GENERATE_RANDOM_EXPERIMENTS = generate_random_experiments(ch_mapped_data_columns, ch_binary_experiment_for_random_experiments, ch_compendia)
ch_random_experiments = GENERATE_RANDOM_EXPERIMENTS

/*--------------------------------------------------
Hypergeometric Activity on Random Experiments
---------------------------------------------------*/
ch_random_experiments_and_network_directory = ch_random_experiments.combine(ch_network_directory_random)
RANDOM_HYPERGEOMETRIC_ACTIVITY = random_hypergeometric_activity(ch_random_experiments_and_network_directory)
ch_random_activity_files = RANDOM_HYPERGEOMETRIC_ACTIVITY[0]
ch_random_activity = RANDOM_HYPERGEOMETRIC_ACTIVITY[1]
ch_random_activity_list = RANDOM_HYPERGEOMETRIC_ACTIVITY[2]
ch_random_aggregated_activities_single = RANDOM_HYPERGEOMETRIC_ACTIVITY[3]

/*--------------------------------------------------
Concatenate Random Activity Files
---------------------------------------------------*/
ch_random_aggregated_activities_single
    .collectFile(
        storeDir:"${params.outdir}/${params.name}/${params.phospho_event}/random_hypergeometric_activity",
        name:"${params.name}_aggregated_activities.tsv",
        keepHeader:true,
        skip:1
    ) 
    .set { ch_collected_random_aggregate_activities} 
    
/*--------------------------------------------------
Summarize Hypergeometric Activity on Random Experiments
---------------------------------------------------*/
SUMMARIZE_RANDOM_AGGREGATED = summarize_random_aggregated(ch_collected_random_aggregate_activities)
summarize_random_aggregated_files = SUMMARIZE_RANDOM_AGGREGATED

/*--------------------------------------------------
Genrate Normalization Values
---------------------------------------------------*/
ch_random_activity_experiment = ch_random_activity.combine(ch_experiment_activities)
GENERATE_NORMALIZATION_VALUES = generate_normalization_values(ch_random_activity_experiment)
ch_normalization_values = GENERATE_NORMALIZATION_VALUES[0]
ch_experiment_and_normalizers = ch_normalization_values.combine(ch_experiment_activities_list_for_normalization)

/*--------------------------------------------------
Calculate Normalizated Activity
---------------------------------------------------*/
CALCULATE_NORMALIZED_ACTIVITY = normalize_activity(ch_experiment_and_normalizers)
ch_normalized_actviity_files = CALCULATE_NORMALIZED_ACTIVITY[0]
ch_normalized_agg_activity = CALCULATE_NORMALIZED_ACTIVITY[1]
ch_normalized_activity_path = CALCULATE_NORMALIZED_ACTIVITY[2]

/*--------------------------------------------------
Concatenate Normalizated Activity Files
---------------------------------------------------*/
ch_normalized_activity_path
    .collectFile(
        storeDir:"${params.outdir}/${params.name}/${params.phospho_event}/normalized_activity",
        name:"${params.name}_normalized_aggregate_activity.tsv",
        keepHeader:true,
        skip:1
    ) 
    .set { ch_normalized_combined } 


/*--------------------------------------------------
Summarize Normalizated Activity
---------------------------------------------------*/
SUMMARIZE_NORMALIZED_ACTIVITY = summarize_normalized(ch_normalized_combined)
summarize_normalized_activity_files = SUMMARIZE_NORMALIZED_ACTIVITY
ch_mann_whitney_input = ch_random_activity_list.combine(ch_normalized_agg_activity).combine(ch_experiment_activities_list_for_mann_whitney)

/*--------------------------------------------------
Mann Whitney on each experiment
---------------------------------------------------*/
MANN_WHITNEY = mann_whitney(ch_mann_whitney_input, ch_number_of_networks)
ch_mann_whitney = MANN_WHITNEY

/*--------------------------------------------------
Combine Results of Mann Whitney
---------------------------------------------------*/
ch_mann_whitney
    .collectFile(
        storeDir:"${params.outdir}/${params.name}/${params.phospho_event}/mann_whitney",
        name:"mann_whitney_combined.tsv",
        keepHeader:true,
        skip:1
    ) 
    .set { ch_mann_whitney_combined } 

/*--------------------------------------------------
Summarize Results of Mann Whitney
---------------------------------------------------*/
SUMMARIZE_MANN_WHITNEY = summarize_mann_whitney(ch_mann_whitney_combined)
summarize_mann_whitney_files = SUMMARIZE_MANN_WHITNEY
}