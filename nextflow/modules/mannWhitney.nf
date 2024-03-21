/*--------------------------------------------------
Mann Whitney on each experiment
---------------------------------------------------*/
process mann_whitney {
    tag "$rand_exp"
    publishDir "${params.outdir}/${params.name}/${params.phospho_event}/individual_experiments/$rand_exp/mann_whitney", mode: 'copy'
    label "single_experiment"
    label "big_memory"

    input:
        tuple val(rand_exp), path(random_activity), val(norm_exp), path(normalized_activity), path(activity) //from ch_mann_whitney_input
        val(number_of_networks) //from ch_number_of_networks
    when:
        rand_exp == norm_exp
    output:
        path("*_mann_whittney.tsv") //into ch_mann_whitney

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
