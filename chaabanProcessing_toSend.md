Description Caaban Project
================

## Reading the DCC files

(following
**[GeomxSet_coercions.html](https://www.bioconductor.org/packages/release/bioc/vignettes/GeomxTools/inst/doc/GeomxSet_coercions.html)**
vignette of the Bioconductor package **GeomxTools**)

``` r
# load packages
library(GeomxTools)
library(Seurat)
library(SpatialExperiment)
library(SpatialDecon)
library(patchwork) # combine separate ggplots into the same graphic
library(RColorBrewer)
library(celldex)
library(SingleR)
library(dplyr)
dr="//derelict/hpc-prj/wren/chaaban_scRNAspatial/allDCC/"
```

``` r
## Read 

############ read All DCC files at once ##########

dr="//derelict/hpc-prj/wren/chaaban_scRNAspatial/allDCC/"
DCCFiles <- dir(dr, pattern=".dcc$", full.names=TRUE);
#PKCFiles <- unzip(zipfile = file.path(dr,"/Hs_R_NGS_WTA_v1.0.pkc_.zip"))
PKCFiles="./Hs_R_NGS_WTA_v1.0.pkc"
SampleAnnotationFile<-file.path(dr,"ROIKeyAnno_BothSet.xlsx") 

objData12=suppressWarnings(readNanoStringGeoMxSet(dccFiles=DCCFiles,pkcFiles=PKCFiles[1],phenoDataFile=SampleAnnotationFile,phenoDataSheet="29June2021_CTA_20210630T1934_La",phenoDataDccColName="Sample_ID",protocolDataColNames=c("aoi","NECType","roi","area","segment","Patient","segType")))

pData(objData12)$Set=c("Fr","Sc")[1*(sData(objData12)$Patient%in%c(5:7))+1]

# access objData12 components 
# objData12 # NanoStringGeoMxSet object ; assayData: 18815 features, 56 samples 
# mat=exprs(objData12);dim(mat) # 18815 56
# colData=sData(objData12); dim(colData); # 56 27
# colnames(colData); # slide name; scan name; panel; Set; FileVersion; SoftwareVersion; Date; SampleID; Plate_ID; Well; SeqSetId; Raw; Trimmed; Stitched; Aligned; umiQ30; rtsQ30; DeduplicatedReads; roi; segment; aoi; area; NECType; Patient; segType; NTC_ID; NTC
# table(colData$Patient,colData$segType)
# table(colData$Set,colData$Patient)

# QC filtering
objData12<- shiftCountsOne(objData12, useDALogic=TRUE) # Shift all 0 counts by 1  (package GeomxTools)
objData12<-setSegmentQCFlags(objData12, qcCutoffs = list(percentSaturation = 45)) # QCFlags data frame appended to protocolData ; 
objData12<-setBioProbeQCFlags(objData12) # Add probe QC flags to object feature data
lowSaturation <- which(protocolData(objData12)[["QCFlags"]]$LowSaturation) # low sequenced ROIs
lowQCprobes <- which(featureData(objData12)[["QCFlags"]]$LowProbeRatio |featureData(objData12)[["QCFlags"]]$GlobalGrubbsOutlier) # probes considered outliers 
#passedQC <- objData12[-lowQCprobes, -lowSaturation]
passedQC <- objData12 # no features/ROIs with low quality
assayDataElementNames(passedQC) # "exprs" "preLocalRemoval" "rawZero" # same as names(assayData(passedQC)) all 18815 56 matrices

# need to aggregated to Target level data before coercing.  (Changes row form probe ID to gene name; aggregate multiple probes per gene)
featureType(passedQC) # Probe
table(table(fData(passedQC)[,"TargetName"])) # there are multiple rows (probes) for each gene ("TargetName" column in fData data frame)
target_demoData <- aggregateCounts(passedQC) # use TargetName column of fData, aggregate on it
featureType(target_demoData) # "Target" (row name of exprs(target_demoData) become gene names)

# 3.    normalize, before coercing to Seurat, using a GeoMx specific model (recommended), rather than use Seurat normalization after)
norm_target_demoData <- normalize(target_demoData, norm_method="quant", desiredQuantile = .75, toElt = "q_norm") # normalized data in assayData slot “q_norm”; 
#norm_target_demoData <- normalize(target_demoData, norm_method="Housekeeping-Log2", toElt = "q_norm") # normalized data in assayData slot “q_norm”; 
assayDataElementNames(norm_target_demoData) # "exprs" "q_norm"

#   Seurat Coercion
#demoSeurat <- as.Seurat(norm_target_demoData, normData = "q_norm") # input need to be aggregated on Target; normalized (can force raw input with forceRaw=TRUE); correct name of norm Data specified;
objSeurat <- as.Seurat(norm_target_demoData, normData = "q_norm", ident = "segType")
```

## SCTransform normalization followed by UMAP and Clustering

``` r
#Initially I got Differentially expressed genes using voom/limma as if I had bulkRna data. First version of umap was done directly #using library(umap) on "vobjDream$E" matrix (see myCode.txt inside chaaban_scRNAspatial folder).

## SCTransform normalization followed by UMAP and Clustering

#Below there is a diffrent umap picture made after renormalizing the seurat object with SCT and going trough standard pipeline (it #merge better the two sets of DCC files, patients 1,2,3,4 respectively patients 5,6,7)

## SCTransform normalization followed by UMAP and Clustering
CDPreAdj=(objSeurat$segType=="CDPre")*(1-2*(objSeurat$Set=="Fr"))
CDPreAdj[CDPreAdj>0]=1/sum(CDPreAdj[CDPreAdj>0]);CDPreAdj[CDPreAdj<0]=1/sum(CDPreAdj[CDPreAdj<0])
PanPreAdj=(objSeurat$segType=="PanPre")*(1-2*(objSeurat$Set=="Fr"))
PanPreAdj[PanPreAdj>0]=1/sum(PanPreAdj[PanPreAdj>0]);PanPreAdj[PanPreAdj<0]=1/sum(PanPreAdj[PanPreAdj<0])
#CDPreAdj=(objSeurat$segType=="CDPre")*(objSeurat$Set=="Fr")
#PanPreAdj=(objSeurat$segType=="PanPre")*(objSeurat$Set=="Fr")
objSeurat$CDPreAdj=CDPreAdj
objSeurat$PanPreAdj=PanPreAdj
objSeurat<-SCTransform(objSeurat,method="glmGamPoi",vars.to.regress=c("CDPreAdj","PanPreAdj","nFeature_GeoMx"),verbose=FALSE,assay="GeoMx")

objSeurat <- RunPCA(objSeurat, verbose = FALSE,npcs=25)
objSeurat <- RunUMAP(objSeurat, dims = 1:24, verbose = FALSE,n.components=2) # use 30 PCs instead of 10 for downstream
objSeurat <- FindNeighbors(objSeurat, dims = 1:24, verbose = FALSE,k.param = 10) # parameter reduction = 'pca' (default) selects which dim reduction to use
objSeurat <- FindClusters(objSeurat, verbose = FALSE)
objSeurat[["umap"]]@cell.embeddings=objSeurat[["umap"]]@cell.embeddings[,2:1]; 
objSeurat[["umap"]]@cell.embeddings[,1]=-objSeurat[["umap"]]@cell.embeddings[,1]; 
colnames(objSeurat[["umap"]]@cell.embeddings)=colnames(objSeurat[["umap"]]@cell.embeddings)[2:1]

objSeurat$seurat_clusters=factor(c(2,3,1)[as.numeric(objSeurat$seurat_clusters)])
table(objSeurat$segType,objSeurat$seurat_clusters)
# separae PanPre and CDPre in 2 groups each
segTypeClst=objSeurat$segType;segTypeClst[segTypeClst=="PanPre"&objSeurat$seurat_clusters==2]="PanPre_1";segTypeClst[segTypeClst=="CDPre"&objSeurat$seurat_clusters==1]="CDPre_1";
objSeurat$segTypeClst=segTypeClst

save(objSeurat,file="objSeurat_UMAPs.rda")
```

#### Here is some brief data description

``` r
load("objSeurat_UMAPs.rda")
tt=table(objSeurat$Patient,objSeurat$segType)[c(1,3,2,4:7),c(1,3,2,4)]
```

``` r
objSeurat # Seurat object, 18677 features 56 samples 1 assay (Active assay: GeoMx)
```

    ## An object of class Seurat 
    ## 37354 features across 56 samples within 2 assays 
    ## Active assay: SCT (18677 features, 3000 variable features)
    ##  1 other assay present: GeoMx
    ##  2 dimensional reductions calculated: pca, umap

Data include expression counts for 18677 features (genes) on 56 ROIs.

Here is the distribution of ROIs per patient per disease/cell category.
The first four patients are from the first set (denoted “Fr” in the
pictures below), next three patients, 5,6,7 being from the second set
(“Sc”), as mentioned in your latest file describing the DCC samples from
the second set.

``` r
print(cbind(SET=c(rep("Fr",4),rep("Sc",3)),Patient=row.names(tt),tt),quote=F,row.names=F)
```

    ##   SET Patient CDCar PanCar CDPre PanPre
    ## 1 Fr  1       3     3      0     0     
    ## 3 Fr  3       4     4      0     0     
    ## 2 Fr  2       0     0      3     3     
    ## 4 Fr  4       0     0      2     2     
    ## 5 Sc  5       0     0      10    6     
    ## 6 Sc  6       0     0      4     6     
    ## 7 Sc  7       0     0      0     6

``` r
names(objSeurat) # "GeoMx" "SCT" "SCT_nn" "SCT_snn" "pca" "umap" 
slotNames(objSeurat) # assays; meta.data; active.assay; active.ident; graphs; neighbors; reductions; images; project.name; misc; version; commands; tools
#Assays(objSeurat) #  "GeoMx" # not working
names(objSeurat@reductions) # "pca"  "umap"
# access
names(objSeurat@assays) # "GeoMx" (before normalization) # "GeoMx" "SCT" after SCTnomalize
objSeurat@assays[["GeoMx"]][1:2,1:5] # 18677 56 (normalized values, with .)
slotNames(objSeurat@assays[["GeoMx"]]) #  counts (18677 56); data (18677 56); scale.data (0 0); key; assay.orig; var.features; meta.features; misc
slotNames(objSeurat[["GeoMx"]]) # same as slotNames(objSeurat@assays[["GeoMx"]])
slotNames(objSeurat@assays[["SCT"]]) # SCTModel.list; counts (18677 56); data (18677 56); scale.data (3000 56); key; assay.orig; var.features; meta.features; misc
head(objSeurat, 3) # most important ROI metadata 56 24 
colnames(head(objSeurat, 3)) # orig.ident; nCount_GeoMx; nFeature_GeoMx; slide.name; scan.name; panel; Set; NegGeoMean_Hs_R_NGS_WTA_v1.0; NegGeoSD_Hs_R_NGS_WTA_v1.0; q_norm_qFactors; SampleID; roi; segment; aoi; area; NECType; Patient; segType; CDPreAdj; PanPreAdj; nCount_SCT; nFeature_SCT; SCT_snn_res.0.8; seurat_clusters
dim(objSeurat@meta.data) #  56 24 # the file above (like pData) – ROI metadata
head(objSeurat@assays$GeoMx@meta.features) # gene metadata 18677 6; colnames: TargetName; Module; CodeClass; GeneID; SystematicName; Negative
```

### UMAPs

I had a second look at the way the data from the two sets were
integrated and come up with some updated version of the UMAP pictures,
hopefully mixing the ROIs from the two sets better this time.  
Here is the new UMAP picture showing the two sets in different colours,
I think the mixing seems a bit better than in the last UMAP I sent you.

``` r
#boxplot(objSeurat[["umap"]]@cell.embeddings[,"UMAP_1"]~objSeurat$Patient)
#boxplot(objSeurat[["umap"]]@cell.embeddings[,"UMAP_2"]~objSeurat$Patient)
load("objSeurat_UMAPs.rda")
#objSeurat$seurat_clusters=factor(c(2,3,1)[as.numeric(objSeurat$seurat_clusters)])
#table(objSeurat$segType,objSeurat$seurat_clusters)
## separae PanPre and CDPre in 2 groups each
#segTypeClst=objSeurat$segType;segTypeClst[segTypeClst=="PanPre"&objSeurat$seurat_clusters==2]="PanPre_1";segTypeClst[segTypeClst=="CDPre"&objSeurat$seurat_clusters==1]="CDPre_1";
#objSeurat$segTypeClst=segTypeClst
#table(objSeurat$segTypeClst)

DimPlot(objSeurat,label=F,group.by ="Set",dims = c(1, 2),pt.size=2)
```

![](chaabanProcessing_toSend_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

And here is the distribution and positioning on the UMAP of each Patient
Rois

``` r
DimPlot(objSeurat,label=F,group.by="Patient",dims=c(1,2),pt.size=2,shape.by="Set",cols=brewer.pal(length(unique(objSeurat$Patient)),"Accent"))
```

![](chaabanProcessing_toSend_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

One encouraging thing is the UMAP seems to sepparate pretty well the
four segment types annotated in the initial DSS file description.

``` r
DimPlot(objSeurat,label=F,group.by="segType",dims=c(1,2),pt.size=2,shape.by="Set")
```

![](chaabanProcessing_toSend_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

As it is customary for scRNA data I tried clustering the ROIs. The
software suggested three clusters as below:

``` r
DimPlot(objSeurat,label=F,group.by ="seurat_clusters",dims = c(1, 2),pt.size=2,shape.by="Set")
```

![](chaabanProcessing_toSend_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
#table(objSeurat$segType,objSeurat$seurat_clusters)
```

Here is how the clusters overlap the initial segment types:

``` r
DimPlot(objSeurat,label=F,group.by="segType",dims=c(1,2),pt.size=2,shape.by="seurat_clusters")
```

![](chaabanProcessing_toSend_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

And here is the number of ROIs per cluster per segment type

``` r
table(objSeurat$segType,objSeurat$seurat_clusters)
```

    ##         
    ##           1  2  3
    ##   CDCar   5  2  0
    ##   CDPre   4 15  0
    ##   PanCar  6  1  0
    ##   PanPre  0  6 17

### Deconvolution with cellDex

I also run again the cell type deconvoluton, this time with a diffrent
package (cellDex), and a very comprehensive reference set Human primary
cell atlas (HPCA). The cell types assigned to each ROI probably differ
substantially from what I sent you before (again I used a diffrent
package and a diffrence reference), but it seems they provide a good
separation of the epitelial versus immune cell.

``` r
#Deconvolution with cellDex as used in Susan project

# http://bioconductor.org/books/3.13/OSCA.basic/cell-type-annotation.html#assigning-cell-labels-from-reference-data
# https://bioconductor.org/packages/3.13/data/experiment/vignettes/celldex/inst/doc/userguide.html # lists options for ref data

cnts=GetAssayData(objSeurat,assay="SCT",slot="data")

library(celldex)
# ref <- BlueprintEncodeData() # Blueprint/ENCODE
ref <- HumanPrimaryCellAtlasData() # Human primary cell atlas (HPCA)
library(SingleR)
pred<-SingleR(test=cnts,ref=ref,labels=ref$label.main)
lbls=pred$labels
#baboon.combined.sct$CellType=lbls
predFine <- SingleR(test=cnts, ref=ref, labels=ref$label.fine)
lblsFine=predFine$labels

objSeurat$CellType=lbls
save(objSeurat,pred,file="objSeurat_UMAPs.rda")
```

Here is a heat map listing the cell type assigned per ROI. Each column
in the heat map is a ROI, each row a cell type. The highest color
intensity decides ROI assigned cell type as marked by the segments on
the top bar (note the yellow segment, epithelial cell, make a large
portion of the assignments) .

``` r
load("objSeurat_UMAPs.rda")

plotScoreHeatmap(pred)
```

![](chaabanProcessing_toSend_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

Here is the distribution of assigned cell types per segment type and
their position on the UMAP (only top 6 most abundant cell types are
included on the picture)

``` r
lbls=objSeurat$CellType
#table(objSeurat$CellType,objSeurat$seurat_clusters)
#table(objSeurat$CellType,objSeurat$segType)
#table(objSeurat$CellType,objSeurat$segTypeClst)
tt=sort(table(objSeurat$CellType),decreasing=T);clk=names(tt[tt>2]);jjj=which(lbls%in%clk);clr=rep("other",length(lbls));clr[jjj]=lbls[jjj];objSeurat$CellType6=clr
DimPlot(objSeurat,label=F,group.by ="CellType6",dims = c(1, 2),pt.size=2,shape.by="segType",cols = brewer.pal(length(clk)+1,"Accent"),order = rev(c(clk,"other")))
```

![](chaabanProcessing_toSend_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

Here are the number of ROIs assigned to each cell type per segment type

``` r
table(objSeurat$CellType,objSeurat$segType)
```

    ##                    
    ##                     CDCar CDPre PanCar PanPre
    ##   Astrocyte             0     1      0      0
    ##   B_cell                2     2      0      0
    ##   DC                    0     0      1      1
    ##   Epithelial_cells      0     0      5     21
    ##   Erythroblast          0     1      0      0
    ##   Gametocytes           0     1      0      0
    ##   iPS_cells             0     0      1      0
    ##   Macrophage            2     3      0      0
    ##   Monocyte              3     0      0      1
    ##   MSC                   0     1      0      0
    ##   Myelocyte             0     1      0      0
    ##   Neutrophils           0     3      0      0
    ##   Platelets             0     2      0      0
    ##   T_cells               0     3      0      0
    ##   Tissue_stem_cells     0     1      0      0

And here is the ROIs distribution per assigned cell type, per each of
the 3 clusters the software suggested

``` r
table(objSeurat$CellType,objSeurat$seurat_clusters)
```

    ##                    
    ##                      1  2  3
    ##   Astrocyte          0  1  0
    ##   B_cell             4  0  0
    ##   DC                 0  2  0
    ##   Epithelial_cells   5  4 17
    ##   Erythroblast       0  1  0
    ##   Gametocytes        0  1  0
    ##   iPS_cells          1  0  0
    ##   Macrophage         1  4  0
    ##   Monocyte           2  2  0
    ##   MSC                0  1  0
    ##   Myelocyte          0  1  0
    ##   Neutrophils        0  3  0
    ##   Platelets          0  2  0
    ##   T_cells            2  1  0
    ##   Tissue_stem_cells  0  1  0

``` r
## Convolution with SpatialDecon  
#Follow SpatialDecon_vignette.html from fost_Seurat/NanoString folder.  
#If we can add more than one matrix for convolution, the **Ileum_Wang** and **ImmuneTumor_safeTME** might be interesting—if not, the one you already selected  (**Gut_HCA**) looks good. (See HumanAdult_CellTypes.xlsx file, in chaaban_scRNAspatial folder, listing the options for matrices)


## NOT WORKING (but it did work last time, see how it was done)
library(Seurat)
library(SpatialDecon)
# download profile matrix from List: CellProfileLibrary GitHub Page (https://github.com/Nanostring-Biostats/CellProfileLibrary/tree/NewProfileMatrices) (use Gut_HCA.RData for your project)
# Elmentaite, R. et al. Cells of the human intestinal tract mapped across space and time. bioRxiv 2021.04.07.438755 (2021) doi:10.1101/2021.04.07.438755.
profileMat<-download_profile_matrix(species = "Human",age_group = "Adult",matrixname = "Gut_HCA") 
head(cellGroups) # list describing cell types; uploaded by download_profile_matrix function
metadata # general info about mousespleen (also uploaded during download)
heatmap(sweep(profileMat, 1, apply(profileMat, 1, max), "/"),labRow = NA, margins = c(10, 5), cexCol = 0.7)
## profileMatCollapsed=collapseCellTypes(fit=profileMat,matching=cellGroups) # safeTME.matches – list name group, cell types (see above)

profileMat_Gut=profileMat;
cellGroups_Gut=cellGroups;

profileMat<-download_profile_matrix(species="Human",age_group="Adult",matrixname="Ileum_Wang"); profileMat_IW=profileMat;
cellGroups_IW=cellGroups;

profileMat<-download_profile_matrix(species="Human",age_group="Adult",matrixname="ImmuneTumor_safeTME");profileMat_TME=profileMat;cellGroups_TME=cellGroups;

norm=objSeurat@assays[["SCT"]]@data
# get background matrix
#bg=derive_GeoMx_background(norm = norm, probepool = rep(1,nrow(norm)),negnames="NegProbe-WTX") # using package function
per.observation.mean.neg = norm["NegProbe-WTX", ]; 
bg = sweep(norm * 0, 2, per.observation.mean.neg, "+") 

#res1=spatialdecon(norm =as.matrix(norm), bg = bg, X = safeTME, cellmerges=safeTME.matches,raw=demoSeurat@assays[["GeoMx"]]@counts)
res1=spatialdecon(norm =as.matrix(norm), bg = bg, X =profileMat_Gut, cellmerges=cellGroups_Gut)
```

### Six groups clustering

Looking at all the pictures and tables above I decided to combine the
segment type classification with the clustering into a new, finer
clustering, including 6 groups, basically along the segment types, with
each CDPre and PanPre split in two.  
Here is the mapping of the 6 clusters onto UMAP representation:

``` r
## separae PanPre and CDPre in 2 groups each
#table(objSeurat$segTypeClst)
DimPlot(objSeurat,label=F,group.by="segTypeClst",dims=c(1,2),pt.size=2,shape.by="seurat_clusters")
```

![](chaabanProcessing_toSend_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

Here is the number of ROIs in each of the 6 clusters:

``` r
## separae PanPre and CDPre in 2 groups each
table(objSeurat$segTypeClst)
```

    ## 
    ##    CDCar    CDPre  CDPre_1   PanCar   PanPre PanPre_1 
    ##        7       15        4        7       17        6

``` r
#DimPlot(objSeurat,label=F,group.by="segTypeClst",dims=c(1,2),pt.size=2,shape.by="seurat_clusters")
```

Here is a depiction of the cell types assigned by deconvolution inside
the six clusters. The picture and the table below it should help in
giving a meaning to each cluster.

``` r
DimPlot(objSeurat,label=F,group.by ="CellType6",dims = c(1, 2),pt.size=2,shape.by="segTypeClst",cols = brewer.pal(length(clk)+1,"Accent"),order = rev(c(clk,"other")))
```

![](chaabanProcessing_toSend_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
table(objSeurat$CellType,objSeurat$segTypeClst)
```

    ##                    
    ##                     CDCar CDPre CDPre_1 PanCar PanPre PanPre_1
    ##   Astrocyte             0     1       0      0      0        0
    ##   B_cell                2     0       2      0      0        0
    ##   DC                    0     0       0      1      0        1
    ##   Epithelial_cells      0     0       0      5     17        4
    ##   Erythroblast          0     1       0      0      0        0
    ##   Gametocytes           0     1       0      0      0        0
    ##   iPS_cells             0     0       0      1      0        0
    ##   Macrophage            2     3       0      0      0        0
    ##   Monocyte              3     0       0      0      0        1
    ##   MSC                   0     1       0      0      0        0
    ##   Myelocyte             0     1       0      0      0        0
    ##   Neutrophils           0     3       0      0      0        0
    ##   Platelets             0     2       0      0      0        0
    ##   T_cells               0     1       2      0      0        0
    ##   Tissue_stem_cells     0     1       0      0      0        0

## Caracterize the Clusters

Also of some help in clusters characterization should be the lists
below, of top 6 genes per each cluster (I picked top 6 to fit in the
heatmap below, showing for each gene the expression profile across all
ROIs and the cluster expression profile separation.)  
Note that, for each cluster, the top genes are genes highly expressed in
that cluster compared to all other 5 cluster (Meaning if there is some
gene missing that you think it should be among these 6 on top, it might
be because that gene is also high in one of the other five clusters)

``` r
# load("objSeurat_UMAPs.rda")
Idents(objSeurat)="segTypeClst"
markersAll=FindAllMarkers(objSeurat,assay="SCT",slot="data",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
#markersAll$t=qnorm(markersAll$p_val)*(1-2*(markersAll$avg_log2FC>0))
#head(markersAll,n=3)
top10s=markersAll%>%group_by(cluster)%>%slice_min(n=6,order_by=p_val)
#print(top10s,n=3)
cbind(tapply(top10s$gene,top10s$cluster,function(x) paste0(x,collapse = ", ")))
```

    ##          [,1]                                           
    ## CDPre_1  "TRAC, JAK3, CORO1A, CXCL13, LTB, LDHB"        
    ## CDPre    "CSF3, SLCO1B7, PYGL, MCOLN2, F9, CATSPER3"    
    ## PanPre   "LLGL2, FXYD3, SLC44A4, CLDN3, UQCRC1, PIGR"   
    ## PanPre_1 "TAS2R45, TSNAX, MRGPRX1, GMPR, ELAVL4, KCNAB1"
    ## PanCar   "ACP1, ATP5MC2, FOXA1, SLC25A39, PERP, PKDCC"  
    ## CDCar    "MYCBP2, HLA-DQB1, PARVG, ATM, DR1, CSAD"

``` r
DoHeatmap(object=objSeurat,features=top10s$gene,assay="SCT",slot="data",size=2.5,draw.lines=T)+NoLegend() 
```

![](chaabanProcessing_toSend_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

Similarly here are the top 10 genes and the heatmap with their
expression profiles for the initial four segment types. Again these are
genes differential expressed between one group and all the other three,
and high in the reference group. These top gene information is supposed
to complement what I sent you before, not replace it. The genes
differential expressed between pairs of groups are those I already sent,
and also uploaded for interactive analysis in Ingenuity.

``` r
# load("objSeurat_UMAPs.rda")
Idents(objSeurat)="segType"
markersAll=FindAllMarkers(objSeurat,assay="SCT",slot="data",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
#markersAll$t=qnorm(markersAll$p_val)*(1-2*(markersAll$avg_log2FC>0))
#head(markersAll)
top10s=markersAll%>%group_by(cluster)%>%slice_min(n=10,order_by=p_val)
#print(top10s,n=3)
cbind(tapply(top10s$gene,top10s$cluster,function(x) paste0(x,collapse = ", ")))
```

    ##        [,1]                                                                       
    ## CDPre  "LRRC8C, ZNF224, TMC8, CD7, ZNF185, DGKA, IGHG2, TDO2, UCP2, MYH10"        
    ## PanPre "PIGR, FAM3D, DHCR24, ATP10B, AOC1, UGT2B15, ELF3, TSPAN13, SCAMP2, TJP3"  
    ## PanCar "ACP1, ATP5MC2, FOXA1, SLC25A39, PERP, PKDCC, AHCY, ARL6IP1, YKT6, MUC17"  
    ## CDCar  "MYCBP2, HLA-DQB1, PARVG, ATM, DR1, CSAD, DENND4B, HLA-DPB1, CHST15, VOPP1"

``` r
DoHeatmap(object=objSeurat,features=top10s$gene,assay="SCT",slot="data",size=2.5,draw.lines=T)+NoLegend() 
```

![](chaabanProcessing_toSend_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

Usually most scRNA papers try to identify some small group of cells that
behaves unexpectedly and try to characterize it (this was the main
reason for us to split the PanPre and CDPre in two). As an attempt to
characterize the two small subgroups below are heat maps of the top
differential expressed genes between PanPre_1 and the rest of PanPre,
respectivele CDPre_1 and the rest of CDPre. Here the selected genes
might be either over expressed or under expressed in the reference group
(unlike as above were it had to be over expressed)  
These lists of markers can be upload in Ingenuity for interactive
analysis and altered pathway identification, along with the others.

#### PanPre_1 markers

``` r
# load("objSeurat_UMAPs.rda")
Idents(objSeurat)="segTypeClst";
jjj=which(objSeurat$segTypeClst%in%c("PanPre_1","PanPre"))
#markers=FindMarkers(objSeurat[,jjj],assay="SCT",slot="data",only.pos = TRUE, min.pct=0.25, logfc.threshold=0.25,ident.1="PanPre_1",ident.2="PanPre" ) 
markers=FindMarkers(objSeurat[,jjj],assay="SCT",slot="data",only.pos = FALSE, min.pct=0.25, logfc.threshold=0.25,ident.1="PanPre_1",ident.2="PanPre" ) 
#markersAll$t=qnorm(markersAll$p_val)*(1-2*(markersAll$avg_log2FC>0))
#head(markersAll)
top10s=markers%>%slice_min(n=25,order_by=p_val)
#print(top10s,n=3)
paste0(row.names(top10s),collapse=", ")
```

    ## [1] "ZNF311, NFE2, UBQLN4, CCIN, PNOC, KDELR1, RAB2A, COX7C, TAGLN2, YWHAE, C19orf33, DDX39B, CST3, H2AC11, LLGL2, EIF3K, CTNND1, LGALS3, COX6C, ELOB, ARF1, VAMP8, HLA-E, HNRNPU, XBP1"

``` r
DoHeatmap(object=objSeurat[,jjj],features=row.names(top10s),assay="SCT",slot="data",size=2.5,draw.lines=T,label=F) 
```

![](chaabanProcessing_toSend_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

#### CDPre_1 markers

``` r
# load("objSeurat_UMAPs.rda")
Idents(objSeurat)="segTypeClst";
jjj=which(objSeurat$segTypeClst%in%c("CDPre_1","CDPre"))
#markers=FindMarkers(objSeurat[,jjj],assay="SCT",slot="data",only.pos = TRUE, min.pct=0.25, logfc.threshold=0.25,ident.1="CDPre_1",ident.2="CDPre" ) 
markers=FindMarkers(objSeurat[,jjj],assay="SCT",slot="data",only.pos =FALSE, min.pct=0.25, logfc.threshold=0.25,ident.1="CDPre_1",ident.2="CDPre" ) 
#markersAll$t=qnorm(markersAll$p_val)*(1-2*(markersAll$avg_log2FC>0))
#head(markersAll)
top10s=markers%>%slice_min(n=25,order_by=p_val)
#print(top10s,n=3)
paste0(row.names(top10s),collapse = ", ")
```

    ## [1] "SLC16A5, PRSS58, RAB6C, MOSPD2, ERV3-1, FER1L6, FAM222B, BTBD18, STPG3, PYGL, RELB, DLEC1, OTUD1, SLC38A5, GH1, MATN3, RBM24, CCDC89, OR10G6, CYP2B6, OR2F2, CTRC, CFAP161, PGLYRP4, EPPIN"

``` r
DoHeatmap(object=objSeurat[,jjj],features=row.names(top10s),assay="SCT",slot="data",size=2.5,draw.lines=T,label=F) 
```

![](chaabanProcessing_toSend_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

### Method Section

And here is a brief description of the analysis tools we used in your
paper.

**ScRNA-seq data preprocessing, quantification and analysis**  
GeoMX digital spatial profiling paraffin embedded tissues from 5 preterm
NEC and 2 cardiac NEC patients were processed and analyzed at OMRF
NanoString GeoMx service Imaging core using a combination of
fluorescently labeled antibodies and the GeoMx Whole Transcriptome Atlas
(WTA), targeting over 18,000 human protein coding genes. Roughly equal
number of immune-rich (CD45) and epithelial-rich (Pan-CK) regions of
interests (ROIs) were selected per patient, summing up to 7 ROIs per
condition in cardiac samples and about 21 ROIs per condition in preterm
samples. After selecting regions of interest (ROIs) on the GeoMx DSP,
the DSP barcodes were UV cleaved and collected. During library
preparation, the DSP barcodes were tagged with their ROI location then
sequenced on an Illumina sequencer. Sequenced oligonucleotides were
processed then imported back into the GeoMx DSP platform for integration
with the slide images and ROI selections for spatially-resolved RNA
expression. Using the Nanostring GeoMx Digital Spatial Profiler (DSP)
sequencing technology whole human transcriptome abundance for each
region of interest (ROI) was quantified. The DSS files with raw count
data were imported into R, QC filtered and normalized using GeomxTools
package from Bioconductor. RNA probe counts used in the analyses were
selected following a sequencing QC according to Nanostring protocols,
where counts from each area of interest are analyzed and under-sequenced
samples are dropped (field of view percentage of 75% and Binding density
from 0.1 to 2.25), and a probe QC, where mRNAs are targeted by multiple
probes and outlier probes are dropped from downstream data analysis
(positive spike-in normalization factor between 0.3 and 3). Then RNA
counts were normalized using a signal-based normalization, in which
individual counts are normalized against the 75th percentile of signal
from their own area of interest. Normalized data were then coerced into
a Seurat object and standard subsequent analyses was performed using the
integrated pipeline of Seurat R package. After selection of highly
variable features the integrated data over all 7 samples were subject to
PCA and nonlinear dimensionality reduction with UMAP/sSNE. 3 clusters of
ROIs were generated using the K-nearest neighbor (KNN) approach with SLM
modularity optimization implemented by FindClusters function (resolution
= 0.015) of the Seurat package. UMAP 2d representation was used to
visualize the result of clustering. Differential expression analysis
between ROI conditions and cancer types, and identification of conserved
cell type markers within each cluster, was performed with Seurat
FindMarkers function and expression profiles of the selected features
were visualized with DoHeatmap function. Adjusted p\<0.05 and
\|log2FC\|\>0.25 was used to defined significant DEGs.  
Functional analysis, identifying sets of genes sharing the same
functionality (GO, KEGG pathways), overrepresented among the lists of
differentially expressed features, was performed with specialized
packages from Bioconductor. Ingenuity Pathway Analysis (IPA, QIAGEN,
Redwood City CA) was used to identify and explore significantly altered
pathways, functional sets and gene networks interactively.  
Mixed cell deconvolution on each ROI was perform with celldex package
functionality. Using the gene expression data, abundance of each cell
type within each observation was estimated and the dominant category was
assigned as ROI cell type. Relative percentage distribution and
enrichment of cell types across clusters, disease type and Pan-CK/CD45
condition was studied and visualized using Correspondence Analysis
tools.
