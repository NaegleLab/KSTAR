/*--------------------------------------------------
Genrate Normalization Values
---------------------------------------------------*/
process generate_normalization_values{
    tag "$data_col"
    publishDir "${params.outdir}/${params.name}/${params.phospho_event}/individual_experiments/$data_col/normalizers", mode: 'copy'

    input:
        // file(experiment_activity) from ch_experiment_activities
        // tuple val(data_col), file(random_activity) from ch_random_activity
        // file(experiment_activity) from ch_experiment_activities
        tuple val(data_col), path(random_activity), path(experiment_activity) //from ch_random_activity_experiment

    output:
        tuple( val(data_col), path("*_normalization.tsv") ) //into ch_normalization_values
        path("*")
    script:
        """
        calculate_fpr.py \
        --data_column $data_col \
        --experiment_activity $experiment_activity \
        --random_activity $random_activity \
        --fpr_alpha ${params.fpr_alpha}

        """
}
