/*--------------------------------------------------
Summarize Normalizated Activity
---------------------------------------------------*/
process summarize_normalized{
    tag "${params.name}"
    publishDir "${params.outdir}/${params.name}/${params.phospho_event}/normalized_activity", mode: 'copy'

    input: 
        path(normalized) //from ch_normalized_combined
    output:
        path("*")
    
    script:
        """
        summarize_activities.py \
        --activities $normalized \
        --method median_normalized_activity \
        --name ${params.name}_normalized_activities
        """
}