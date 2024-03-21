/*--------------------------------------------------
Kinase Hypergeometric Activity 
---------------------------------------------------*/
process hypergeometric_activity{
    tag "${params.name}"
    publishDir "${params.outdir}/${params.name}/${params.phospho_event}/hypergeometric_activity", mode: 'copy' 
    label "all_experiments"
    label "big_memory"

    input:
        path(experiment) //from ch_binary_experiment_for_hypergeometric
        path(network_directory) //from ch_network_directory_hypergeometric
        val(mapped_data_column_string) //from ch_mapped_data_columns

    output:
        path("*")
        path("${params.name}_activities.tsv") //from ch_experiment_activities
        path("${params.name}_activities_list.tsv") //from ch_experiment_activities_list

    script:

        """
        hypergeometric_activity_binary.py \
        --experiment_file $experiment \
        --network_directory $network_directory/${params.phospho_event}/INDIVIDUAL_NETWORKS \
        --pevent ${params.phospho_event} \
        --name ${params.name} \
        --data_columns $mapped_data_column_string \
        --max_cpus ${task.cpus} \
        --chunk_size ${params.chunk_size}
        """
}