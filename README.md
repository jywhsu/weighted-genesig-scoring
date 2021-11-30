# weighted-genesig-scoring


Developed using `R 3.6.3`

PURPOSE: To calculate a continuous score of likeness for a given transcriptome based on a weighted gene signature from training case-control data\

### 1. Data Requirements
* **Case-control "Training" Data**: Differential expression analysis output from edgeR.
* **Test data**: Z-score converted gene expression matrix. 

*Gene symbols must be the first column of both the training and test matrices, they must be labeled as "hgnc_symbol", and they must match (i.e. you cannot have gene symbols in the training dataset and ensembl ids in the test dataset).*

### 2. Content
#### R Script
`s-score_calculation.R`: S-scoring R script
* `pm.score.calc`: function to calculate weighted gene signature based on FDR and case-ctl fold change expression. 
* `SigScorePrep`: intermediate function to clean up training data and test data. 
* `SigScoreCalc`: function to calculate Signature or "S"-score

#### Example data
`casectl_tube_pax7_edgeR_analysis.txt`: Case-control data from [GSE129505](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129505).   
`six1kd_rnaseq_z-score_cpm_data.txt`: Test transcriptome data from [GSE173155](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173155)   

* For case-control data, cases were defined as iPSC-differentiated Myotubes; Controls were defined as PAX7+ putative muscle progenitors. 
* For more information about the details of the test data, please refer to our upcoming Cell Reports manuscript. This data will be made publicly available on April 29, 2022. 
