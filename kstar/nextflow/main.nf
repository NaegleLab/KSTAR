
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
    // phosporylation_events=['Y','ST']
    // name='test'
    // data_columns='data:average:data:treated_to_untreated:EOE(drug washout)', 'data:average:data:treated_to_untreated:HDP3(3hrs post treatment)', 'data:average:data:treated_to_untreated:HDP6(6hrs post treatment)','data:average:data:treated_to_untreated:pre-treatment'
    // activity_agg='count'
    // threshold=1.0
    // num_random_experiments=10

mapped_data_columns = []
    for(col in params.data_columns){
        mapped_data_columns.add("data:"+col)
    }
mapped_data_column_string = mapped_data_columns.join(" ")

Channel
    .from(mapped_data_columns)
    .set { ch_mapped_data_columns}
    


        
/*--------------------------------------------------
Kinase Hypergeometric Activity 
---------------------------------------------------*/

process kinase_hypergeometric_activity{
    tag "${phospho_event}"
    publishDir "${params.outdir}/$phospho_event/kinase_hypergeometric_activity", mode: 'copy' 
    input:
        each phospho_event from params.phospho_events
        file(experiment) from ch_experiment
        file(network_y) from ch_network_y 
        file(network_st) from ch_network_st
        // tuple (data_columns) from params.data_columns
    output:
        file("*")
        file("${params.name}_binarized_experiment.tsv") into ch_binarized_experiment
    script:
        if (phospho_event =="ST")
            """
            hypergeometric_activity.py \
            --experiment_file $experiment \
            --networks $network_st \
            --pevent $phospho_event \
            --name ${params.name} \
            --activity_agg ${params.activity_aggregate} \
            --threshold ${params.threshold}
            """
        else if (phospho_event =="Y")
            """
            hypergeometric_activity.py \
            --experiment_file $experiment \
            --networks $network_y \
            --pevent $phospho_event \
            --name ${params.name} \
            --activity_agg ${params.activity_aggregate} \
            --threshold ${params.threshold} \
            --data_columns $mapped_data_column_string
            """
}


/*--------------------------------------------------
Generate Random Experiments 
---------------------------------------------------*/
process generate_random_experiments {
  tag "${phospho_event} ${data_col}"
  cpus 1
  publishDir "${params.outdir}/$phospho_event/random_experiments", mode: 'copy'

  input:
    each phospho_event from params.phospho_events
    each data_col from ch_mapped_data_columns
    file(experiment) from ch_binarized_experiment
    file(compendia) from ch_compendia

  
  
    output:
        tuple(val(phospho_event), val(data_col), file("*")) into ch_random_experiments
  
  
  script:
        """
        generate_random_experiments_single_data_column.py \
        --experiment_file $experiment \
        --phosphorylation_events $phospho_event \
        --data $data_col \
        --num_random_experiments ${params.num_random_experiments} \
        --compendia $compendia
        """
}

// /*--------------------------------------------------
// Run KSTAR on Random Experiments
// ---------------------------------------------------*/
process random_hypergeometric_activity{
    tag "${phospho_event} ${data_col}"
    publishDir "${params.outdir}/$phospho_event/random_hypergeometric_activity", mode: 'copy'

    input:
        file(network_y) from ch_network_y
        file(network_st) from ch_network_st
        tuple val(phospho_event), val(data_col), file(experiment) from ch_random_experiments
    output:
        file("*")
    script:
        if (phospho_event =="ST")
            """
            hypergeometric_activity.py \
            --experiment_file $experiment \
            --networks $network_st \
            --pevent $phospho_event \
            --name ${data_col}_random_activity 
            """
        else if (phospho_event =="Y")
            """
            hypergeometric_activity.py \
            --experiment_file $experiment \
            --networks $network_y \
            --pevent $phospho_event \
            --name ${data_col}_random_activity 
            """
}