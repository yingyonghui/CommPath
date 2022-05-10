# CommPath
CommPath is an R package for inference and analysis of ligand-receptor interactions from single cell RNA sequencing data.
## Installation
CommPath R package can be easily installed from Github using devtools:
```
devtools::install_github("yingyonghui/CommPath")
library(CommPath)
```
### Dependencies
- [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html)
- [circlize](https://cran.r-project.org/web/packages/circlize/index.html)
- [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
- [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
- [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)
- [GSVA](https://www.bioconductor.org/packages/release/bioc/html/GSVA.html) (suggested)

## Tutorials
In this vignette we show CommPath's steps and functionalities for inference and analysis of ligand-receptor interactions by applying it to a scRNA-seq data (GEO accession number: GSE156337) on cells from hepatocellular carcinoma (HCC) patients.
### Brief description of CommPath object
We start CommPath analysis by creating a CommPath object, which is a S4 object and consists of six slots including (i) ***data***, a matrix containing the normalized expression values by gene $\times$ cell; (ii) ***cell.info***, a dataframe contain the information of cells; (iii) ***meta.info***, a list containing some important parameters used during the analysis; (iv) ***LR.marker***, a dataframe containing the result of differential expression test of ligands and receptors; (v) ***interact***, a list containing the  information of LR interaction among clusters; (v) ***interact.filter***, a list containing the  information of filtered LR interaction among clusters; (vi) ***pathway***, a list containing the information of pathways related to the ligands and receptors.
### CommPath input
The expression matrix and cell indentity information are required for CommPath input. We downloaded the processed HCC scRNA-seq data from [Mendeley data](https://doi.org/10.17632/6wmzcskt6k.1). For a fast review and illustration of CommPath's functionalities, we randomly selected the expression data of 3000 cells across the top 5000 highly variable genes from the tumor and normal tissues, respectively. The example data are available in [figshare](https://figshare.com/articles/dataset/HCC_tumor_normal_3k_RData/19090553).
We here illustrate the CommPath steps for data from the tumor tissues. And analysis for data from the normal tissues would be roughly in the same manner.
```
# load(url("https://figshare.com/ndownloader/files/34914603"))
load("path_to_download/HCC.tumor.3k.RData")
```
This dataset consists of 3 varibles, 2 of which are required for CommPath input:
***tumor.expr*** : the expression matrix of gene $\times$ cell. Expression values are required to be first normalized by the library-size and log-transformed;
***tumor.label*** : a vector of lables indicating identity classes of cells in the expression matrix, and the order of lables should match the order of cells in the expression matrix; usrs may also provide a dataframe containing the meta infomation of cells with the row names matching the cells in the expression matrix and a column named ***Cluster*** must be included to indicate identity classes of cells.
The remaining one variable ***tumor.obj*** is the CommPath object created from ***tumor.expr*** and ***tumor.label*** and processed by CommPath standard procedures. This variable would be temporarily ignored since we will recreate it by the following steps. 
#### Identification of marker ligands and receptors
We start CommPath analysis by creating a CommPath object:
```
tumor.obj <- createCommPath(expr.mat = tumor.expr, 
		cell.info = tumor.label, 
		species = "hsapiens")
```
Firstly we're supposed to identify marker ligands and receptors (ligands and receptors that are significantly highly expressed) for each cluster of cells in the expression matrix. CommPath provide **findLRmarker** to identify these markers by *t.test* or *wilcox.test*.
```
tumor.obj <- findLRmarker(object = tumor.obj, method = "wilcox.test")
```

#### Statistical identification of potential ligand-receptor (LR) associations
```
# To find significant LR pairs
tumor.obj <- findLRpairs(object = tumor.obj,
		logFC.thre = 0, 
		p.thre = 0.05)
```
The counts of Statistically significant LR pairs and overall interaction intensity among cell clusters are then stored in ***tumor.obj@interact[['InteractNumer']]***，and the detailed information of each LR pair is stored in ***tumor.obj@interact[['InteractGene']]***. 

Then we can visualize all interactions through a circos plot:
```
# To show the counts of LR associations among all clusters
# Here we set the parameter "filter" as FALSE, which means that those LR interactions are identified only based on their expression profiles, not filtered by pathways in the receiver cells (as described in the later sections)
pdf('circosPlot-count.nonfiltered.pdf',height=6,width=6)
circosPlot(object = tumor.obj, filter=FALSE)
dev.off()
```
<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/circosPlot-count.nonfiltered.png" height=300, width=300>

In the above circos plot, the directions of lines indicate the associations from ligands to receptors, and the widths of lines represent the counts of LR pairs among clusters.

```
# To show the overall interaction intensity of LR interactions among all clusters
pdf('circosPlot-intensity.nonfiltered.pdf',height=6,width=6)
circosPlot(object = tumor.obj, plot="intensity", filter=FALSE)
dev.off()
```
<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/circosPlot-intensity.nonfiltered.png" height=300, width=300>

Now the widths of lines represent the overall interaction intensity among clusters.
#### Pathway enrichment analysis
CommPath conducts pathway analysis to identify  dysregulated signaling pathways containing the marker ligands and receptors for each cluster.
```
# To find pathways in which genesets show overlap with the marker ligands and receptors
# CommPath provides pathway annotations from KEGG pathways, WikiPathways, reactome pathways, and GO terms
# Here we take the KEGG pathways as an example
tumor.obj <- findLRpath(object = tumor.obj, category = "kegg")
```
Now genesets showing overlap with the marker ligands and receptors are stored  in ***tumor.obj@pathway[['pathwayLR']]***. Then we score the pathways to measure the activation levels for each pathway in each cell.
```
# To compute pathway activation score by the gsva algorithm or in an average manner
# For more information about the gsva algorithm, see the GSVA package (PMID23323831)
tumor.obj <- scorePath(object = tumor.obj, method = "gsva", min.size = 10)
```
After that CommPath provides **diffAllPath** to perform pathway differential activation analysis for each cluster and find the receptors and ligands in each pathway:
```
# To get significantly up-regulated pathways in each cluster
acti.path.dat <- diffAllPath(object = tumor.obj, only.posi = TRUE, only.sig = TRUE)
head(acti.path.dat)
```
There are several columns stored in the variable ***acti.path.dat***:
Columns ***mean.diff***, ***mean.1***, ***mean.2***, ***t***, ***df***, ***p.val*** and ***p.val.adj*** show the statistic results; ***description*** shows the name of pathway; 
Columns ***cell.up*** and ***ligand.up*** show the upstream clusters which would release specific ligands to interact with the receptors expressed by the current cluster; 
Column ***receptor.in.path*** and ***ligand.in.path*** show the marker receptors and ligands expressed by the current cluster and included in the current pathway;

Then we use **pathHeatmap** to present a heatmap of those differentially activated pathways for each cluster:
```
pdf('pathHeatmap.pdf',height=10,width=7)
pathHeatmap(object = tumor.obj,
       acti.path.dat = acti.path.dat,
       top.n.pathway = 10,
       cell.aver = TRUE)
dev.off()
```

<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/pathHeatmap.png" height=650, width=455>

#### Screening LR interactions associated with activated pathways
For each cell cluster, CommPath identifies LR associations involved in the activated pathways to screen functional LR interactions. Those pathway-filtered interacctions are considered to be more likely to trigger the corresponding molecular pathways in the receiver cells.
```
# To screen functional LR interactions
tumor.obj <- filterLR(object = tumor.obj, acti.path.dat = acti.path.dat)
```
Then counts of filtered significant LR pairs and the corresponding overall interaction intensity among cell clusters are then stored in ***tumor.obj@interact.filter[['InteractNumer']]***，and the detailed information of each LR pair is stored in ***tumor.obj@interact.filter[['InteractGene']]***.

Filtered interactions could be also visualized through the circos plot:
```
# To show the counts of filtered LR associations among all clusters
pdf('circosPlot-count.filtered.pdf',height=6,width=6)
circosPlot(object = tumor.obj, filter=TRUE)
dev.off()
```
<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/circosPlot-count.filtered.png" height=300, width=300>

```
# To show the overall interaction intensity of filtered LR interactions among all clusters
pdf('circosPlot-intensity.filtered.pdf',height=6,width=6)
circosPlot(object = tumor.obj, plot="intensity", filter=TRUE)
dev.off()
```
<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/circosPlot-intensity.filtered.png" height=300, width=300>

Users would highlight the interaction of specific clusters. Here we take the endothelial cell as an example:
```
select.ident = 'Endothelial'
pdf('circosPlot-Endothelial-count.filtered.pdf',height=6,width=6)
circosPlot(object = tumor.obj, select.ident = select.ident, filter=TRUE)
dev.off()
```
<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/circosPlot-Endothelial-count.filtered.png" height=300, width=300>

For a specific cluster of interest, CommPath provides function **findLigand** (**findReceptor**) to find the upstream (downstream) cluster and the corresponding ligand (receptor) interacting with the specific cluster: 
```
# To find upstream clusters and ligands for the selected cluster and receptor 
select.ident = "Endothelial"
select.receptor = "ACKR1"

ident.up.dat <- findLigand(object = tumor.obj, 
    select.ident = select.ident, 
    select.receptor = select.receptor)
head(ident.up.dat)

# To find downstream clusters and receptors for the selected cluster and ligand 
select.ident = "Endothelial"
select.ligand = "CXCL12"

ident.down.dat <- findReceptor(object = tumor.obj, 
    select.ident = select.ident, 
    select.ligand = select.ligand)
head(ident.down.dat)
```

There are 7 columns stored in the variable ***ident.up.dat***/***ident.down.dat***:
Columns ***cell.from***, ***cell.to***, ***ligand***, ***receptor*** show the upstream and dowstream clusters and the specific ligands and receptors for the LR associations;
Columns ***log2FC.LR***, ***P.val.LR***, ***P.val.adj.LR*** show the interaction intensity (measured by the product of *log2FCs* of ligands and receptors) and the corresponding original and adjusted ***P*** values for hypothesis tests of the LR pairs.

CommPath also provides dot plots to investigate the upstream clusters which release specific ligands and the downstream clusters which express specific receptors: 
```
# To investigate the upstream clusters which release ligands to the selected cluster
pdf('dotPlot-ligand.pdf',height=5,width=10)
dotPlot(object = tumor.obj, receptor.ident = select.ident)
dev.off()
```
<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/dotPlot-ligand.png" height=300, width=600>

```
# To investigate the downstream clusters which express receptors for the selected cluster
pdf('dotPlot-receptor.pdf',height=5,width=10.5)
dotPlot(object = tumor.obj, ligand.ident = select.ident)
dev.off()
```

<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/dotPlot-receptor.png" height=300, width=630>


Then CommPath provides network graph tools to visualize the pathways and associated functional LR interactions:
```
# First to integrate the statistics of activated pathways and their associated LR interactions 
tumor.obj <- pathNet(object = tumor.obj, acti.path.dat = acti.path.dat)
# To visualize the pathways and LR interactions in a network plot
pdf('LR.pathway.net.Endothelial.pdf',width=6,heigh=6)
set.seed(1234)
pathNetPlot(object = tumor.obj, select.ident =  select.ident,
    layout='layout.davidson.harel', vert.size.path.adj=8, vert.size.LR=3,
    LR.label=TRUE, pathway.label=TRUE, vertex.label.cex=0.25)
dev.off()
```
<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/LR.pathway.net.Endothelial.png" height=500, width=530>

In the above network plot, the pie charts represent the activated pathways in the selected  cells (here endothelial cells) and the scatter points represent the LR pairs of which the receptors are included in the genesets of the linked pathways. Colors of scatter points indicate the upstream clusters releasing the corresponding ligands. Sizes of pie charts indicate their total in-degree and the proportions indicate the in-degree from different upstream clusters.
#### Identification of pathway-mediated cell-cell communication chain
For a specific cell cluster, here named as B for demonstration, CommPath identifies the upstream cluster A sending signals to B, the downstream cluster C receiving signals from B, and the significantly activated pathways in B to mediate the A-B-C communication chain. More exactly, through LR and pathways analysis described above, CommPath is able to identify LR pairs between A and B, LR pairs between B and C, and pathways activated in B. Then CommPath screens for pathways in B which involve both the receptors to interact with A and ligands to interact with C.
```
# To identify receptors and the associated activated pathways for a specific cluster
select.ident = 'Endothelial'
pdf('pathPlot.pdf',height=6,width=10)
pathPlot(object = tumor.obj, 
    select.ident = select.ident, 
    acti.path.dat = acti.path.dat)
dev.off()
```

<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/pathPlot.png" height=300, width=390>

In the above line plot, the widths of lines between ***Upstream*** cluster and ***Receptor*** represent the overall interaction intensity between the upstream  cluster and endothelial cells via the specific receptors; the sizes and colors of dots in the ***Receptor*** column represent the average *log2FC* and *-log10(P)* of the expression of receptors in endothelial cells compared to all other cells; the lengths and colors of bars in the ***Pathway annotation*** column represent the mean difference and *-log10(P)* of pathway scores in endothelial cells compared to all other cells.

```
# To identify receptors, the associated activated pathways, and the downstream clusters
pdf('pathInterPlot.pdf',height=6,width=14)
pathInterPlot(object = tumor.obj, 
    select.ident = select.ident, 
    acti.path.dat = acti.path.dat)
dev.off()
```

<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/pathInterPlot.png" height=300, width=600>

The legend of the above line plot is generally the same to that of the previous plot from **pathPlot**.

#### Comparison of cell-cell communication between two conditions
CommPath also provides useful utilities to compare cell-cell interactions between two conditions such as disease and control. Here we, for example, use CommPath to compare the cell-cell interactions between cells from HCC tumor and normal tissues. We have pre-created the CommPath object for the normal samples following the above steps, and the data are also available in [figshare](https://figshare.com/articles/dataset/HCC_tumor_normal_3k_RData/19090553).
```
# load(url("https://figshare.com/ndownloader/files/34914606"))
load("path_to_download/HCC.normal.3k.RData")
```
This dataset consists of 3 varibles:
***normal.expr*** : expression matrix for cells from normal tissues;
***normal.label*** : indentity lables for cells from normal tissues;
***normal.obj*** : CommPath object created from ***normal.expr*** and ***normal.label***, and processed by CommPath steps described above.

To compare 2 CommPath objects, we shall first identify the differentially activated pathways between the same cluster of cells in the 2 objects.
```
select.ident <- 'Endothelial'
diff.path.dat <- comparePath(object.1 = tumor.obj, 
		object.2 = normal.obj, 
		select.ident = select.ident)
```
Then we compare the differentially activated pathways and the cell-cell communication chain mediated by those pathways.
```
# To compare differentially activated pathways and the involved receptors between the selected clusters in 2 CommPath objects
pdf('pathPlot-compare.pdf',height=6,width=10)
pathPlot.compare(object.1 = tumor.obj, 
		object.2 = normal.obj, 
		select.ident = select.ident, 
		diff.path.dat = diff.path.dat)
dev.off()
```

<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/pathPlot-compare.png" height=300, width=390>

In the above line plot, the widths of lines between ***Upstream*** cluster and ***Receptor*** represent the overall interaction intensity between the upstream clusters and endothelial cells via the specific receptors, and the colors indicate the interaction intensity is upregulated (red) or downregulated (blue) in tumor tissues (object.1) compared to that in normal tissues (object.2); the sizes and colors of dots in the ***Receptor*** column represent the average *log2FC* and *-log10(P)* of expression of receptors in endothelial cells compared to all other cells in tumor tissues; the lengths and colors of bars in the ***Pathway annotation*** column represent the mean difference and *-log10(P)* of pathway scores of endothelail cells in tumor tissues compared to that in normal tissues.

```
# To compare the pathway-mediated cell-cell communication chain for a specific cluster between 2 CommPath objects
pdf('pathInterPlot-compare.pdf',height=6,width=14)
pathInterPlot.compare(object.1 = tumor.obj, 
		object.2 = normal.obj, 
		select.ident = select.ident, 
		diff.path.dat = diff.path.dat)
dev.off()
```

<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/pathInterPlot-compare.png" height=300, width=600>

The legend of the above line plot is generally the same to that of the previous plot from **pathPlot.compare**.

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
[1] CommPath_1.0.0  igraph_1.2.6    Matrix_1.3-4    ggplot2_3.3.5
[5] dplyr_1.0.7     plyr_1.8.6      reshape2_1.4.4  circlize_0.4.13

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.7          pillar_1.6.2        compiler_4.0.3
 [4] tools_4.0.3         lifecycle_1.0.0     tibble_3.1.3
 [7] gtable_0.3.0        lattice_0.20-44     pkgconfig_2.0.3
[10] rlang_1.0.2         cli_3.0.1           DBI_1.1.1
[13] withr_2.4.2         stringr_1.4.0       generics_0.1.0
[16] GlobalOptions_0.1.2 vctrs_0.3.8         grid_4.0.3
[19] tidyselect_1.1.2    glue_1.4.2          R6_2.5.1
[22] fansi_0.5.0         purrr_0.3.4         magrittr_2.0.1
[25] scales_1.1.1        ellipsis_0.3.2      assertthat_0.2.1
[28] shape_1.4.6         colorspace_2.0-2    utf8_1.2.2
[31] stringi_1.7.4       munsell_0.5.0       crayon_1.4.1        
```
