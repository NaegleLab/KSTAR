/*--------------------------------------------------
Summarize Hypergeometric Activity on Random Experiments
---------------------------------------------------*/
process summarize_random_aggregated{
    tag "${params.name}"
    publishDir "${params.outdir}/${params.name}/${params.phospho_event}/random_hypergeometric_activity", mode: 'copy'
    input: 
        path(combined_agg) //from ch_collected_random_aggregate_activities
    output:
        path("*")
    
    script:
        """
        summarize_activities.py \
        --activities $combined_agg \
        --method median_activity \
        --name ${params.name}_random_activities
        """
}