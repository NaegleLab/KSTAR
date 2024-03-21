/*--------------------------------------------------
Summarize Results of Mann Whitney
---------------------------------------------------*/
process summarize_mann_whitney{
    tag "${params.name}"
    publishDir "${params.outdir}/${params.name}/${params.phospho_event}/mann_whitney", mode: 'copy'

    input: 
        path(mann_whitney) //from ch_mann_whitney_combined
    output:
        path("*")
    
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