
nextflow run main.nf -profile Y docker \
--name bcr-abl \
--phospho_event Y \
--experiment_file test_data/BCR-ABL_mapped.tsv \
--data_columns data:EOE,data:pre-treatment,data:HDP3,data:HDP6 \
--num_random_experiments 15 \
--outdir ./results \
--resource_directory ../RESOURCE_FILES \
--network_directoyr /NETWORKS/NetworKIN \
--threshold 0.5 \
--activity_aggregate mean \
--fpr_alpha 0.05 \
--number_of_sig_trials 10
