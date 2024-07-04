Here is the code and some related files used in processing and analysis of the spatial RNAseq data described in the paper 
**Spatial Transcriptomics and Network-Based Bioinformatics Differentiate Intestinal Phenotypes of Cardiac and Classical Necrotizing Enterocolitis** by
Kathryn Y. Burge et al.  
* **objSeurat_UMAPs.rda** includes the Seurat object with all the data and downstream analysis derived items (clusters, dimension reduction maps, sample metadata). The file can be loaded in R and the object inside accessed for further analysis using Seurat package functionality.
* **convolution.rda** contains the results of cell type deconvolution with SpatialDecon package from Bioconductor. It can also be loaded in R for further analysis.
* **chaabanProcessing_toSend** files are several version of the codding script used for data analysis. The RMD version can be loaded in R for replicating analysis.
* The md version, **[chaabanProcessing_toSend.md](https://github.com/ttgeo/codeShaabanPaper/blob/master/chaabanProcessing_toSend.md)**, is the friendlier markdown presentation of the code including, along with the code chuncks, text comments and some of the analysis results.
