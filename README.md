# SNP_Prioritization_Pipeline

The SNP Pipeline leverages multilpe sources of developmental tissue data to identify key genetic variants associated with craniofacial disorders for further biological validation. 

Using the proprietary OroLab keratinocyte *in vitro* differentiation system, large amount of non-neural ectoderm tissue is produced. This developmental tissue offers new insights into the early origins of craniofacial disoders. 

<p align="center">
  <img src="https://github.com/OroLabStanford/SNP_pipeline/blob/master/images/non-neural_ectoderm.png">
</p>


The non-neural ectoderm differentiation is initiated by exposure to retinoid acid and BMP4 for 7 days. These cells give rise to craniofacial structures. 



<p align="center">
  <img src="https://github.com/OroLabStanford/SNP_pipeline/blob/master/images/SNP_Pipeline.png">
</p>


The SNP pipeline is summarized in the above diagram. In the black is a chromatin loop, which represents portions of the genome which are linearly far apart but come in close physical proximity due to DNA folding (HiChIP). 

Step 1 of the pipeline (pipeline_scripts/step1.py) is to find regulatory elements lying in either side of a chromatin loop. Step 2 of the pipeline (pipeline_scripts/step2.py) is to find SNPs lying within these identified regulatory elements. Step 3 of the pipeline (pipeline_scripts/step3.py) is to identify genes falling on the opposite end of chromatin loops as the matched SNPs. The final output represents a SNP, which falls in a region known to be important for gene expression regulation, that is coming in close physical proximity with a gene via chromatin looping. For convenience, all steps of the pipeline were combined in (pipeline_scripts/all_steps.py). Resulting trios can be filtered by a specific loci (pipeline_scripts/Trio_Loci_Filter.py).

Preliminary findings can be found in (pipeline_analysis/). These results are based off the use of in house data along with a dataset of SNPs from 16 genome-wide assocation studies in mature craniofacial tissue, where SNPs that appear significantly more often in the diseased vs. control population were pooled. This dataset included 262 statistically significant SNPs, along with ~8600 additional SNPs in linkage disequilibrium with these SNPs. 

In the future, an FDR (false discovery rate) will be established for each regulator element, SNP, gene trio using random simulation. These scripts can be found in (random_simulation_scripts/). 
