
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


ch_network_directory.into{
    ch_network_directory_hypergeometric
    ch_network_directory_random

}

Channel
    .fromPath("${params.resource_directory}${params.network_directory}/${params.phospho_event}/INDIVIDUAL_NETWORKS/*.tsv")
    .count()
    .set {ch_number_of_networks}


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
    mapped_data_columns = []
    if(params.add_data_before){
        for(col in params.data_columns){
            mapped_data_columns.add("data:"+col)
        }
    }
    else{
        for(col in params.data_columns){
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

// mapped_data_columns = params.data_columns
mapped_data_column_string = mapped_data_columns.join(" ")
print(mapped_data_column_string)

Channel
    .from(mapped_data_columns)
    .set { ch_mapped_data_columns}
    

// /*--------------------------------------------------
// Binarize Experiment
// ---------------------------------------------------*/
process binarize_experiment{
    publishDir "${params.outdir}/${params.name}/${params.phospho_event}/binary_experiment", mode: 'copy'

    input:
        file(experiment) from ch_experiment
    output:
        file("${params.name}_binarized_experiment.tsv") into ch_binary_experiment
    
    script:
        """
        binarize_experiment.py \
        --evidence $experiment \
        --data_columns $mapped_data_column_string \
        --name ${params.name} \
        --threshold ${params.threshold} \
        --agg ${params.activity_aggregate} \
        --greater ${params.greater}
        """
}

ch_binary_experiment.into {
    ch_binary_experiment_for_hypergeometric
    ch_binary_experiment_for_random_experiments
    }
/*--------------------------------------------------
Kinase Hypergeometric Activity 
---------------------------------------------------*/
process hypergeometric_activity{
    tag "${params.phospho_event}"
    publishDir "${params.outdir}/${params.name}/${params.phospho_event}/hypergeometric_activity", mode: 'copy' 
    label "all_experiments"
    label "big_memory"

    input:
        file(experiment) from ch_binary_experiment_for_hypergeometric
        path(network_directory) from ch_network_directory_hypergeometric

    output:
        file("*")
        file("${params.name}_activities.tsv") into ch_experiment_activities
        file("${params.name}_activities_list.tsv")  into ch_experiment_activities_list

    script:

        """
        hypergeometric_activity_binary.py \
        --evidence_file $experiment \
        --network_directory $network_directory/${params.phospho_event}/INDIVIDUAL_NETWORKS \
        --pevent ${params.phospho_event} \
        --name ${params.name} \
        --data_columns $mapped_data_column_string \
        --max_cpus ${task.cpus}
        """
}


ch_experiment_activities_list.into{
    ch_experiment_activities_list_for_normalization
    ch_experiment_activities_list_for_mann_whitney

}

/*--------------------------------------------------
Generate Random Experiments 
---------------------------------------------------*/
process generate_random_experiments {
  tag "${data_col}"
  cpus 1
  publishDir "${params.outdir}/${params.name}/${params.phospho_event}/individual_experiments/$data_col/random_experiments", mode: 'copy'
  label "big_memory"
  input:
    each data_col from ch_mapped_data_columns
    file(experiment) from ch_binary_experiment_for_random_experiments
    file(compendia) from ch_compendia

    output:
        tuple( val(data_col), file("*")) into ch_random_experiments
  
  
  script:
        """
        generate_random_experiments_single_data_column.py \
        --experiment_file $experiment \
        --phosphorylation_events ${params.phospho_event} \
        --data $data_col \
        --num_random_experiments ${params.num_random_experiments} \
        --compendia $compendia
        """
}

/*--------------------------------------------------
Hypergeometric Activity on Random Experiments
---------------------------------------------------*/
ch_random_experiments_and_network_directory = ch_random_experiments.combine(ch_network_directory_random)

process random_hypergeometric_activity{
    tag "${data_col}"
    publishDir "${params.outdir}/${params.name}/${params.phospho_event}/individual_experiments/$data_col/random_hypergeometric_activity", mode: 'copy'
    label "single_experiment"
    label "big_memory"
    // memory '16 GB'
    

    input:
        // file(network_y) from ch_network_y
        // file(network_st) from ch_network_st
        // path(network_y) from ch_network_y_path_random
        // path(network_st) from ch_network_st_path_random
        // path(network_directory) from ch_network_directory_random
        tuple val(data_col), file(experiment), path(network_directory) from ch_random_experiments_and_network_directory
    output:
        file("*")
        tuple( val(data_col), file("*_random_activities.tsv") ) into ch_random_activity
        tuple( val(data_col),file("*random_activities_list.tsv")) into ch_random_activity_list
        path("*_aggregated_activities.tsv") into ch_random_aggregated_activities_single


    script:
        """
        random_hypergeometric_activity.py \
        --evidence_file $experiment \
        --network_directory $network_directory/${params.phospho_event}/INDIVIDUAL_NETWORKS \
        --pevent ${params.phospho_event} \
        --name ${data_col}_random \
        --max_cpus ${task.cpus}
        """
}

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
process summarize_random_aggregated{
    publishDir "${params.outdir}/${params.name}/${params.phospho_event}/random_hypergeometric_activity", mode: 'copy'
    input: 
        file(combined_agg) from ch_collected_random_aggregate_activities
    output:
        file("*")
    
    script:
        """
        summarize_activities.py \
        --activities $combined_agg \
        --method median_activity \
        --name ${params.name}_random_activities
        """
}


/*--------------------------------------------------
Genrate Normalization Values
---------------------------------------------------*/

ch_random_activity_experiment = ch_random_activity.combine(ch_experiment_activities)

process generate_normalization_values{
    tag "$data_col"
    publishDir "${params.outdir}/${params.name}/${params.phospho_event}/individual_experiments/$data_col/normalizers", mode: 'copy'

    input:
        // file(experiment_activity) from ch_experiment_activities
        // tuple val(data_col), file(random_activity) from ch_random_activity
        // file(experiment_activity) from ch_experiment_activities
        tuple val(data_col), file(random_activity), file(experiment_activity) from ch_random_activity_experiment

    output:
        tuple( val(data_col), file("*_normalization.tsv") ) into ch_normalization_values
    script:
        """
        calculate_fpr.py \
        --data_column $data_col \
        --experiment_activity $experiment_activity \
        --random_activity $random_activity \
        --fpr_alpha ${params.fpr_alpha}

        """
}

ch_experiment_and_normalizers = ch_normalization_values.combine(ch_experiment_activities_list_for_normalization)
/*--------------------------------------------------
Calculate Normalizated Activity
---------------------------------------------------*/
process normalize_activity {
    tag "$data_col"
    publishDir "${params.outdir}/${params.name}/${params.phospho_event}/individual_experiments/$data_col/normalized_activity", mode: 'copy'

    input:
        // file(experiment_activity) from ch_experiment_activities_list_for_normalization
        // tuple val(data_col), file(normalizers) from ch_normalization_tuple

        tuple val(data_col), file(normalizers), file(experiment_activity) from ch_experiment_and_normalizers

        

    output:
        file("*")
        tuple( val(data_col), file("*_normalized_aggregate_activity.tsv") ) into ch_normalized_agg_activity
        path("*_normalized_aggregate_activity.tsv") into ch_normalized_activity_path
    script:
        """
        normalize.py \
        --data_column $data_col \
        --experiment_activity $experiment_activity \
        --normalizers $normalizers \
        """
}

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
process summarize_normalized{
    publishDir "${params.outdir}/${params.name}/${params.phospho_event}/normalized_activity", mode: 'copy'

    input: 
        file(normalized) from ch_normalized_combined
    output:
        file("*")
    
    script:
        """
        summarize_activities.py \
        --activities $normalized \
        --method median_normalized_activity \
        --name ${params.name}_normalized_activities
        """
}



ch_mann_whitney_input = ch_random_activity_list
                            .combine(ch_normalized_agg_activity)
                            .combine(ch_experiment_activities_list_for_mann_whitney)
/*--------------------------------------------------
Mann Whitney on each experiment
---------------------------------------------------*/
process mann_whitney {
    tag "$rand_exp"
    publishDir "${params.outdir}/${params.name}/${params.phospho_event}/$rand_exp/mann_whitney", mode: 'copy'
    label "single_experiment"
    label "big_memory"

    input:
        tuple val(rand_exp), file(random_activity), val(norm_exp), file(normalized_activity), file(activity) from ch_mann_whitney_input
        val(number_of_networks) from ch_number_of_networks
    when:
        rand_exp == norm_exp
    output:
        path("*_mann_whittney.tsv") into ch_mann_whitney

    script:
        """
        mann_whitney.py \
        --activity_list $activity \
        --random_activity_list $random_activity \
        --normalized_activities $normalized_activity \
        --num_networks $number_of_networks \
        --num_sig_trials ${params.number_of_sig_trials} \
        --num_random_experiments ${params.num_random_experiments} \
        --experiment_name $norm_exp \
        --max_cpus ${task.cpus}
        """
}
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
process summarize_mann_whitney{
    publishDir "${params.outdir}/${params.name}/${params.phospho_event}/mann_whitney", mode: 'copy'

    input: 
        file(mann_whitney) from ch_mann_whitney_combined
    output:
        file("*")
    
    script:
        """
        summarize_activities.py \
        --activities $mann_whitney \
        --method mann_whitney_activities \
        --name ${params.name}_mann_whitney_activities

        summarize_activities.py \
        --activities $mann_whitney \
        --method mann_whitney_fpr \
        --name ${params.name}_mann_whitney_fpr
        """
}