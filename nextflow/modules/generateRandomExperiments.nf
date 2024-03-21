/*--------------------------------------------------
Generate Random Experiments 
---------------------------------------------------*/
process generate_random_experiments {
  tag "${data_col}"
  cpus 1
  publishDir "${params.outdir}/${params.name}/${params.phospho_event}/individual_experiments/$data_col/random_experiments", mode: 'copy'
  label "big_memory"
  input:
    each data_col //from ch_mapped_data_columns
    path experiment //from ch_binary_experiment_for_random_experiments
    path compendia //from ch_compendia

    output:
        tuple( val(data_col), path("*")) //into ch_random_experiments
  
  
  script:
        """
        generate_random_experiments_single_data_column.py \
        --experiment_file $experiment \
        --phosphorylation_events ${params.phospho_event} \
        --data $data_col \
        --num_random_experiments ${params.num_random_experiments} \
        --compendia $compendia
        """
}
