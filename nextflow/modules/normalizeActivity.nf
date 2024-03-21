/*--------------------------------------------------
Calculate Normalizated Activity
---------------------------------------------------*/
process normalize_activity {
    tag "$data_col"
    publishDir "${params.outdir}/${params.name}/${params.phospho_event}/individual_experiments/$data_col/normalized_activity", mode: 'copy'

    input:
        // file(experiment_activity) from ch_experiment_activities_list_for_normalization
        // tuple val(data_col), file(normalizers) from ch_normalization_tuple
        tuple val(data_col), path(normalizers), path(experiment_activity) //from ch_experiment_and_normalizers

        

    output:
        path("*")
        tuple( val(data_col), path("*_normalized_aggregate_activity.tsv") ) //into ch_normalized_agg_activity
        path("*_normalized_aggregate_activity.tsv") //into ch_normalized_activity_path
    script:
        """
        normalize.py \
        --data_column $data_col \
        --experiment_activity $experiment_activity \
        --normalizers $normalizers \
        """
}
