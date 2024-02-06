# KSTAR
Kinase-Substrate Translation to Activity Relationships 

Branch used to test various approaches to improve speed and memory usage in KSTAR. calculate.py contains the same code from the main branch (as of 2/6/2024), while calculate2.py contains updated functions for updating. Code used to test calculate.py and calculate2.py can be found in KinaseActivityApplications repository (private).

Current limitations of KSTAR speed (italicized limitations that are being addressed in calculate2.py)
1. *In calculate.randomized_analysis, KSTAR generates all 150 random experiments, stores them in unique variables, and then calculates activity. For large datasets, the memory burden of storing these random experiments can make it infeasible to run on most laptops. Because of this, almost all serine/threonine runs need to be done on Rivanna or another high performance computing environment.*
2. A unique set of random experiments is generated for each sample to best match the study bias of sites found in that sample. However, for datasets with many samples, this quickly balloons the time required for calculation. If study bias distribution is similar, it should obtain similar results to use the same random experiments for many samples.

Current differences in calculate2.py (update as changes are made)
1. After generating each random experiment, immediately calculate kinase activity for that experiment, save activities, and then delete experiment. Provide option to save random experiments to file. This reduces memory for semi-large datasets, but is slightly slower as is. Multiprocessing is not currently working and needs to be fixed.
