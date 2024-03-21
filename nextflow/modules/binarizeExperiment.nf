/*--------------------------------------------------
Binarize Experiment
---------------------------------------------------*/
process binarize_experiment {
    publishDir "${params.outdir}/${params.name}/${params.phospho_event}/binary_experiment", mode: 'copy'

    input:
        path(experiment) //from ch_experiment
        val(mapped_data_column_string)

    output:
        path("${params.name}_binarized_experiment.tsv") //into ch_binary_experiment

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
