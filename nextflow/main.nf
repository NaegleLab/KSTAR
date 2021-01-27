
Channel
    .value(file(params.experiment_file))
    .ifEmpty { error "Cannot find gtf file for parameter --experiment_file: ${params.experiment_file}" }
    .set { ch_experiment }  

Channel
    .value(file(params.resource_directory + params.compendia))
    .ifEmpty { error "Cannot find gtf file for parameter --compendia: ${params.compendia}" }
    .set { ch_compendia }  

Channel
    .value(file(params.resource_directory + params.reference_fasta))
    .ifEmpty { error "Cannot find gtf file for parameter --reference_fasta: ${params.reference_fasta}" }
    .set { ch_reference_fasta} 

Channel
    .value(file(params.resource_directory + params.network_y))
    .ifEmpty { error "Cannot find gtf file for parameter --network_y: ${params.network_y}" }
    .set { ch_network_y} 

Channel
    .value(file(params.resource_directory + params.network_st))
    .ifEmpty { error "Cannot find gtf file for parameter --network_st: ${params.network_st}" }
    .set { ch_network_st} 


mapped_data_columns = []
    for(col in params.data_columns){
        mapped_data_columns.add("data:"+col)
    }
// mapped_data_columns = params.data_columns
mapped_data_column_string = mapped_data_columns.join(" ")

Channel
    .from(mapped_data_columns)
    .set { ch_mapped_data_columns}
    
/*--------------------------------------------------
Map Experiment
---------------------------------------------------*/

/*--------------------------------------------------
Binarize Experiment
---------------------------------------------------*/
process binarize_experiment{
    publishDir "${params.outdir}/${params.phospho_event}/binary_experiment", mode: 'copy'

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
/*--------------------------------------------------
Kinase Hypergeometric Activity 
---------------------------------------------------*/
process hypergeometric_activity{
    tag "${params.phospho_event}"
    publishDir "${params.outdir}/${params.phospho_event}/hypergeometric_activity", mode: 'copy' 
    label "all_experiments"

    input:
        file(experiment) from ch_binary_experiment
        file(network_y) from ch_network_y 
        file(network_st) from ch_network_st
        // tuple (data_columns) from params.data_columns
    output:
        file("*")
        file("${params.name}_activities.tsv") into ch_experiment_activities
        file("${params.name}_activities_list.tsv")  into ch_experiment_activities_list

    script:
        if (params.phospho_event =="ST")
            """
            hypergeometric_activity_binary_evidence.py \
            --experiment_file $experiment \
            --networks $network_st \
            --pevent ${params.phospho_event} \
            --name ${params.name} \
            --data_columns $mapped_data_column_string \
            """
        else if (params.phospho_event =="Y")
            """
            hypergeometric_activity_binary_evidence.py \
            --experiment_file $experiment \
            --networks $network_y \
            --pevent ${params.phospho_event} \
            --name ${params.name} \
            --data_columns $mapped_data_column_string \
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
  publishDir "${params.outdir}/${params.phospho_event}/$data_col/random_experiments", mode: 'copy'
  label "single_experiment"
  input:
    each data_col from ch_mapped_data_columns
    file(experiment) from ch_experiment
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
process random_hypergeometric_activity{
    tag "${data_col}"
    publishDir "${params.outdir}/${params.phospho_event}/$data_col/random_hypergeometric_activity", mode: 'copy'
    label "single_experiment"

    input:
        file(network_y) from ch_network_y
        file(network_st) from ch_network_st
        tuple val(data_col), file(experiment) from ch_random_experiments
    output:
        file("*")
        tuple( val(data_col), file("*_random_activities.tsv") ) into ch_random_activity
        tuple( val(data_col),file("*random_activities_list.tsv")) into ch_random_activity_list
        path("*_aggregated_activities.tsv") into ch_random_aggregated_activities_single


    script:
        if (params.phospho_event == "ST")
            """
            hypergeometric_activity_binary_evidence.py \
            --experiment_file $experiment \
            --networks $network_st \
            --pevent ${params.phospho_event} \
            --name ${data_col}_random_activity \

            """
        else if (params.phospho_event == "Y")
            """
            hypergeometric_activity_binary_evidence.py \
            --experiment_file $experiment \
            --networks $network_y \
            --pevent ${params.phospho_event} \
            --name ${data_col}_random \
            """
}

/*--------------------------------------------------
Concatenate Random Activity Files
---------------------------------------------------*/
ch_random_aggregated_activities_single
    .collectFile(
        storeDir:"${params.outdir}/${params.phospho_event}/random_hypergeometric_activity",
        name:"${params.name}_aggregated_activities.tsv",
        keepHeader:true,
        skip:1
    ) 
    .set { ch_collected_random_aggregate_activities} 
    
/*--------------------------------------------------
Summarize Hypergeometric Activity on Random Experiments
---------------------------------------------------*/
process summarize_random_aggregated{
    publishDir "${params.outdir}/${params.phospho_event}/random_hypergeometric_activity", mode: 'copy'

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

process generate_normalization_values{
    tag "${params.fpr_alpha}, $data_col"
    publishDir "${params.outdir}/${params.phospho_event}/$data_col/normalizers", mode: 'copy'
    label "single_experiment"

    input:
        file(experiment_activity) from ch_experiment_activities
        tuple val(data_col), file(random_activity) from ch_random_activity

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

/*--------------------------------------------------
Calculate Normalizated Activity
---------------------------------------------------*/
process normalize_activity {
    tag "${params.fpr_alpha}, $data_col"
    publishDir "${params.outdir}/${params.phospho_event}/$data_col/normalized_activity", mode: 'copy'
    label "single_experiment"

    input:
        file(experiment_activity) from ch_experiment_activities_list_for_normalization
        tuple val(data_col), file(normalizers) from ch_normalization_values

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
        storeDir:"${params.outdir}/${params.phospho_event}/normalized_activity",
        name:"${params.name}_normalized_aggregate_activity.tsv",
        keepHeader:true,
        skip:1
    ) 
    .set { ch_normalized_combined } 

/*--------------------------------------------------
Summarize Normalizated Activity
---------------------------------------------------*/
process summarize_normalized{
    publishDir "${params.outdir}/${params.phospho_event}/normalized_activity", mode: 'copy'

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

/*--------------------------------------------------
Mann Whitney on each experiment
---------------------------------------------------*/
process mann_whitney {
    tag "$rand_exp"
    publishDir "${params.outdir}/${params.phospho_event}/$rand_exp/mann_whitney", mode: 'copy'
    label "single_experiment"

    input:
        file(activity) from ch_experiment_activities_list_for_mann_whitney
        tuple val(rand_exp), file(random_activity) from ch_random_activity_list
        tuple val(norm_exp), file(normalized_activity) from ch_normalized_agg_activity
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
        --num_networks ${params.number_of_networks} \
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
        storeDir:"${params.outdir}/${params.phospho_event}/mann_whitney",
        name:"mann_whitney_combined.tsv",
        keepHeader:true,
        skip:1
    ) 
    .set { ch_mann_whitney_combined } 

/*--------------------------------------------------
Summarize Results of Mann Whitney
---------------------------------------------------*/
process summarize_mann_whitney{
    publishDir "${params.outdir}/${params.phospho_event}/mann_whitney", mode: 'copy'

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