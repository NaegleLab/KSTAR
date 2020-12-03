#!/bin/bash
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -o test_kstar.out
#SBATCH -p standard

python kstar_runner.py \
--experiment_file /home/bj8th/Documents/NaegleLab/KSTAR/analysis/BCR-ABL/PRE-MAPPED/experiment.1704.tsv \
--output_directory /Users/bj8th/Documents/GitHub/KSTAR/analysis/BCR-ABL \
--phospho_types Y \
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
--greater yes \
--normalize yes \
--num_random_experiments 10 \
--target_alpha 0.05
