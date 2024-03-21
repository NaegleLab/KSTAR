/*--------------------------------------------------
Hypergeometric Activity on Random Experiments
---------------------------------------------------*/
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
        tuple val(data_col), path(experiment), path(network_directory) //from ch_random_experiments_and_network_directory
    output:
        path("*")
        tuple( val(data_col), path("*_random_activities.tsv") ) //into ch_random_activity
        tuple( val(data_col),path("*random_activities_list.tsv")) //into ch_random_activity_list
        path("*_aggregated_activities.tsv") //ch_randointo ch_random_aggregated_activities_single

    script:
        """
        random_hypergeometric_activity.py \
        --experiment_file $experiment \
        --network_directory $network_directory/${params.phospho_event}/INDIVIDUAL_NETWORKS \
        --pevent ${params.phospho_event} \
        --name ${data_col}_random \
        --max_cpus ${task.cpus} \
        --chunk_size ${params.chunk_size}
        """
}