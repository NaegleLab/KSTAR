python run_kstar.py \
--experiment_file /Users/bj8th/Documents/GitHub/KSTAR/analysis/BCR-ABL/PRE-MAPPED/experiment.1704.tsv \
--resource_directory /Users/bj8th/Documents/GitHub/KSTAR/RESOURCE_FILES \
--output_directory /Users/bj8th/Documents/GitHub/KSTAR/analysis/BCR-ABL \
--phosphorylation_events Y \
--name BCR-ABL \
--accession query_accession \
--site mod_sites \
--peptide aligned_peptides \
--window 7 \
--data_columns "average:data:treated_to_untreated:EOE(drug washout)" \
"average:data:treated_to_untreated:HDP3(3hrs post treatment)" \
"average:data:treated_to_untreated:HDP6(6hrs post treatment)" \
"average:data:treated_to_untreated:pre-treatment" \
--activity_agg mean \
--threshold 0.0 \
--normalize yes \
--num_random_experiments 10