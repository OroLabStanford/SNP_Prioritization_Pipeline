# SNP_Prioritization_Pipeline

## Overview
The SNP Priortization Pipeline leverages multiple types of cell-specific data to prioritize genetic variants associated with a phenotype of interest. Variants associated with a phenotype of interest that reside within known important epigenetic regulatory regions in a cell-type interest are associated with their gene targets via three-dimensional chromatin looping within the cell-type of interest. The genetic variants that participate in the specified three-dimensional configuration are then prioritized for biological validation using betweeness centrality. This is calculated by building a bipartite graph of genetic variants and their associated gene targets and then calculating the number of times a node in the graph acts as a bridge along the shortest path between two other nodes. The higher the betweness centrality nodes is then used to rank sort the genetic variants. We recommend this pipeline be used to prioritize genetic variants for biological validation

<p align="center">
  <img src="https://github.com/OroLabStanford/SNP_pipeline/blob/master/images/SNP_Pipeline.png">
</p>

The SNP pipeline is summarized in the above diagram. In the black is a chromatin loop, which represents portions of the genome which are linearly far apart but come in close physical proximity due to DNA folding (HiC or HiChIP). 

## Input Files
All files need to be tab delimited in the same genome assembly of your choosing. Bed file refers to the 4 indexed format of chr start stop name.
- Bed file contiaining genetic variants significantly associated with your phenotype of interest. We recommend including all genetic variants that are in linkage disequilibrium (r=0.7) with variants that are siginicantly associated with your phenotype of interest via one or more Genome-Wide Association Study (GWAS).
- Bed file containing regulatory element positions in your cell-type of interest. These can include open chromatin regions (ATAC-seq or FAIRE-seq), histone modifications (ChIP-seq), or transcription factor binding sites (ChIP-seq) data.
- Bed file containing the genomic coordinates of genes.
- HiChIP file containing the three dimensional chromatin looping (e.g. HiCHIP or HiC) in the format: chr1 start1 stop1 chr2 start2 stop2 count fdr name. We recommend using high resolution chromatin looping data of 5 kb resolution or less.


## Script Files
- (identify_genetic_variant_gene_targets.py)[identify_genetic_variant_gene_targets.py] takes the input files listed above as input and outputs genetic variants that are in regulatory elements that are associated with a distal gene target via chromatin looping.
- (Centrality.ipynb)[Centrality.ipynb] is a jupyter notebook that conducts centrality on genetic variant and gene bipartite graphs and create a ranked list of genetic variants of interest. The input for this notebook is the output of identify_genetic_variant_gene_targets.py.
- (SNP_gene_all_RE_bipartite_graph.ipynb)[SNP_gene_all_RE_bipartite_graph.ipynb] generates visualizations of the connectivity between regulatory elements containing genetic variants of interest and their distal gene targets. The regulatory elements are colored squares whose label is the number of unique genetic variants of interest that are contained within the regulatory element and the color specifies the type of regulatory element. The input for this notebook is the output of identify_genetic_variant_gene_targets.py.
