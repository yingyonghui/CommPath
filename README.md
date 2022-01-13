# Commpath
Commpath is an R package for inference and analysis of ligand-receptor interactions from single cell RNA sequencing data
## Installation
Commpath R package can be easily installed from Github using devtools:
```
devtools::install_github("yingyonghui/Commpath")
library(Commpath)
```
### Dependencies
- [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html)
- [circlize](https://cran.r-project.org/web/packages/circlize/index.html)
- [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
- [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
- [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)
- [GSVA](https://www.bioconductor.org/packages/release/bioc/html/GSVA.html) (suggested)

## Tutorials
#### Load the built-in sample dataset in commpath
```
data("HCC.tumor.3k",package='Commpath')
```
***tumor.expr*** : expression matrix of gene * cell. Expression values are required to be first normalized by the library-size and log-transformed

***tumor.label*** : data frame containing the meta infomation of cells

***tumor.marker*** : data frame of marker genes for each identity class, usually calculated by FindAllMarkers from [Seurat](https://satijalab.org/seurat/)

***tumor.gsva*** : precomputed gsva scores of pathways for the example dataset
#### Identification of marker ligands and receptors
We start Commpath analysis by creating a Commpath object, which is a S4 object and consists of five slots including (i) data, a matrix containing the normalized expression values by gene * cell; (ii) meta.info, a list containing  the meta information about cells and some important parameters used during the  analysis; (iii) LR.marker, a data.frame containing the result of differential expression test of ligands and receptors; (iv) interact, a list containing the  information of LR interaction among clusters; (vi) pathway, a list containing the information of pathways related to the ligands and receptors.
```
### to create a Commpath object
object <- createCommpath(expr.mat=tumor.expr, 
		cell.info=tumor.label, 
		species='hsapiens')
		
# type ?createCommpath to get more information about each parameter
?createCommpath
```
Firstly we're supposed to identify marker ligands and receptors (ligands and receptors that are significantly highly expressed) for each identity class of cells in the expression matrix. commpath provide **findLRmarker** to identify these markers by *t.test* or *wilcox.test*.
```
# to identify marker ligands and receptors
# object <- findLRmarker(object, method='wilcox.test')
# to save time, we have pre-identified marker ligands and receptors
# and saved it in the varible sample.marker
object@LR.marker <- tumor.marker
```

#### Identification of ligand-receptor (L-R) associations
```
# find significant L-R pairs
logFC.thre = 0
p.thre = 0.05
object <- findLRpairs(object,
		logFC.thre=logFC.thre, 
		p.thre=p.thre)
```
The number of significant L-R pairs is then stored in object@interact[['InteractNumer']]，and the detailed information of each L-R pair is stored in object@interact[['InteractGeneUnfold']]

Then you can visualize the interaction through a circos plot:
```
# plot interaction for all cluster
circosPlot(object)
```
<img src="https://github.com/yingyonghui/Commpath/blob/main/pic/circosPlot.png" height=300, width=300>

```
# you may want to highlight the interaction of specific cluster
# here we take the endothelial cell as an example
ident='Endothelial'
circosPlot(object, ident=ident)
```
<img src="https://github.com/yingyonghui/Commpath/blob/main/pic/circosPlot-Endothelial.png" height=300, width=300>

For a specific cluster of interest, commpath also provides dot plots to investigate its upstream clusters which release specific ligands and its downstream clusters which expressed specific receptors: 
```
ident='Endothelial'
# to investigate the upstream clusters which release specific ligands to the interested cluster
dotPlot(object, receptor.ident=ident)
```
<img src="https://github.com/yingyonghui/Commpath/blob/main/pic/dotPlot-ligand.png" height=300, width=400>

```
# to investigate the downstream clusters which expressed specific receptors for the interested cluster
dotPlot(object, ligand.ident=ident)
```
<img src="https://github.com/yingyonghui/Commpath/blob/main/pic/dotPlot-receptor.png" height=300, width=400>

Also commpath provides function **findLigand** (**findReceptor**) to find the upstream (downstream) cluster and the corresponding ligand (receptor) for specific cluster and receptor (ligand) 
```
# For the selected cluster and selected receptor, find the upstream cluster
select.ident = 'Endothelial'
select.receptor = 'ITGB1'

ident.up.dat <- findLigand(object, 
    select.ident=select.ident, 
    select.receptor=select.receptor)
head(ident.up.dat)

# For the selected cluster and selected ligand, find the downstream cluster
select.ident = 'Endothelial'
select.ligand = 'CCL14'

ident.down.dat <- findReceptor(object, 
    select.ident=select.ident, 
    select.ligand=select.ligand)
head(ident.down.dat)
```
#### Pathway analysis
```
# find pathways in which genesets show overlap 
# with the ligands and receptors in the example dataset
object <- findLRpath(object, category='kegg')
```
Now genesets show overlap with the ligands and receptors in the exsample dataset are saved in object@interact[['pathwayLR']]

Scoring the pathways：
```
# to compute pathway activation score by the gsva algorithm or an average manner
# library(GSVA)
# object <- scorePath(object, method='gsva', min.size=10, parallel.sz=10)
# to save time, we have precomputed gsva score and saved it in the varible *gsva.mat*
object@pathway$acti.score <- gsva.mat
# normal.gsva <- object@pathway$acti.score
```
Pathway differential activation analysis：
```
# to find the different activated pathways for cells in the selected identity class 
# and the receptor and ligand in the pathway
ident.label = sample.label
select.ident.1 = 'Endothelial'
method = 't.test'
ident.path.dat <- diffPath(object, 
    select.ident.1=select.ident.1,
    method=method)

head(ident.path.dat)

# perform diffPath for all clusters
all.path.dat <- diffAllPath(object, method=method)
# get all significant pathways
all.path.dat <- subset(all.path.dat, p.val.adj < 0.05)

head(all.path.dat)
```
Columns ***mean.diff***, ***mean.1***, ***mean.2***, ***t***, ***df***, ***p.val***, ***p.val.adj*** show the statistic result; *description* shows the name of pathway; 

Columns ***cell.up*** and ***ligand.up*** show the upstream identity classes which would release specific ligands to interact with the receptors from the current identity class; 

Column ***receptor.in.path*** shows the marker receptors expressed by the current identity class and these receptors are included in the current pathway;

Column ***ligand.in.path*** shows the marker ligands released by the current identity class and these ligands are also included in the current pathway.

Then we use **pathHeatmap** to plot a heatmap of those differentially enriched pathways for each cluster to display the highly variable pathways:
```
pathHeatmap(object,
       all.path.dat,
       top.n.pathway = 10,
       sort = "p.val.adj")
```
<img src="https://github.com/yingyonghui/Commpath/blob/main/pic/pathHeatmap.png" height=350, width=600>

#### Cell-cell interaction flow via pathways

```
### first we identify differentially activated pathways associated with receptors in the selected ident
ident.label = sample.label
select.ident.1 = 'Endothelial'
method = 't.test'
ident.path.dat <- diffPath(object, 
    select.ident.1=select.ident.1,
    method=method)

# save.image('before.lineplot.hcc.RData')
# load('before.lineplot.hcc.RData')

# visualization of the identified pathways
# plot to identify receptors and the associated activated pathways
receptorPathPlot(object, 
    select.ident=select.ident.1, 
    ident.path.dat=ident.path.dat)
```
<img src="https://github.com/yingyonghui/Commpath/blob/main/pic/receptorPathPlot.png" height=300, width=380>

```
# plot to identify receptors, the associated activated pathways, and the downstream clusters
pathInterPlot(object, 
    select.ident=select.ident.1, 
    ident.path.dat=ident.path.dat)
```
<img src="https://github.com/yingyonghui/Commpath/blob/main/pic/pathInterPlot.png" height=300, width=600>

#### sessionInfo()
```
R version 4.0.3 (2020-10-10)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /home/luh/miniconda3/envs/seurat4/lib/libopenblasp-r0.3.17.so

locale:
 [1] LC_CTYPE=zh_CN.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=zh_CN.UTF-8        LC_COLLATE=zh_CN.UTF-8    
 [5] LC_MONETARY=zh_CN.UTF-8    LC_MESSAGES=zh_CN.UTF-8   
 [7] LC_PAPER=zh_CN.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=zh_CN.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] GSVA_1.38.2     ggplot2_3.3.5   dplyr_1.0.7     reshape2_1.4.4 
[5] circlize_0.4.13 commpath_0.1.0 

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.7                  lattice_0.20-44            
 [3] digest_0.6.27               assertthat_0.2.1           
 [5] utf8_1.2.2                  R6_2.5.1                   
 [7] GenomeInfoDb_1.26.7         plyr_1.8.6                 
 [9] stats4_4.0.3                RSQLite_2.2.8              
[11] httr_1.4.2                  pillar_1.6.2               
[13] zlibbioc_1.36.0             GlobalOptions_0.1.2        
[15] rlang_0.4.11                annotate_1.68.0            
[17] blob_1.2.2                  S4Vectors_0.28.1           
[19] Matrix_1.3-4                labeling_0.4.2             
[21] BiocParallel_1.24.1         stringr_1.4.0              
[23] RCurl_1.98-1.4              bit_4.0.4                  
[25] munsell_0.5.0               DelayedArray_0.16.3        
[27] compiler_4.0.3              pkgconfig_2.0.3            
[29] BiocGenerics_0.36.1         shape_1.4.6                
[31] tidyselect_1.1.1            SummarizedExperiment_1.20.0
[33] tibble_3.1.3                GenomeInfoDbData_1.2.4     
[35] IRanges_2.24.1              matrixStats_0.60.1         
[37] XML_3.99-0.7                fansi_0.5.0                
[39] crayon_1.4.1                withr_2.4.2                
[41] bitops_1.0-7                grid_4.0.3                 
[43] xtable_1.8-4                GSEABase_1.52.1            
[45] gtable_0.3.0                lifecycle_1.0.0            
[47] DBI_1.1.1                   magrittr_2.0.1             
[49] scales_1.1.1                graph_1.68.0               
[51] stringi_1.7.4               cachem_1.0.6               
[53] farver_2.1.0                XVector_0.30.0             
[55] ellipsis_0.3.2              generics_0.1.0             
[57] vctrs_0.3.8                 tools_4.0.3                
[59] bit64_4.0.5                 Biobase_2.50.0             
[61] glue_1.4.2                  purrr_0.3.4                
[63] MatrixGenerics_1.2.1        parallel_4.0.3             
[65] fastmap_1.1.0               AnnotationDbi_1.52.0       
[67] colorspace_2.0-2            GenomicRanges_1.42.0       
[69] memoise_2.0.0              
```
