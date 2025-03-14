# CommPath
CommPath is an R package for inference and analysis of ligand-receptor interactions from single cell RNA sequencing data.
## Installation
CommPath R package can be easily installed from Github using devtools:
```
devtools::install_github("yingyonghui/CommPath")
library(CommPath)
```
### Dependencies
- library([Matrix](https://cran.r-project.org/web/packages/Matrix/index.html))
- library([circlize](https://cran.r-project.org/web/packages/circlize/index.html))
- library([ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html))
- library([dplyr](https://cran.r-project.org/web/packages/dplyr/index.html))
- library([reshape2](https://cran.r-project.org/web/packages/reshape2/index.html))
- library([GSVA](https://www.bioconductor.org/packages/release/bioc/html/GSVA.html)) (suggested)

## Tutorials
In this vignette we show CommPath's steps and functionalities for inference and analysis of ligand-receptor interactions by applying it to a scRNA-seq data (GEO accession number: GSE156337) on cells from hepatocellular carcinoma (HCC) patients.
### Brief description of CommPath object
We start CommPath analysis by creating a CommPath object, which is a S4 object and consists of six slots including (i) ***data***, a matrix containing the **normalized expression values** by gene * cell; (ii) ***cell.info***, a dataframe contain the information of cells; (iii) ***meta.info***, a list containing some important parameters used during the analysis; (iv) ***LR.marker***, a dataframe containing the result of differential expression test of ligands and receptors; (v) ***interact***, a list containing the  information of LR interaction among clusters; (vi) ***interact.filter***, a list containing the  information of filtered LR interaction among clusters; (vii) ***pathway***, a list containing the information of pathways related to the ligands and receptors; (viii) ***pathway.net***, a list containing the integrated information of the statistics of LR interactions and associated pathways.
### CommPath input
The expression matrix and cell identity information are required for CommPath input. We downloaded the processed HCC scRNA-seq data from [Mendeley data](https://doi.org/10.17632/6wmzcskt6k.1). For a fast review and illustration of CommPath's functionalities, we randomly selected the expression data of 3000 cells across the top 5000 highly variable genes from the tumor and normal tissues, respectively. The example data are available in [figshare](https://figshare.com/articles/dataset/HCC_tumor_normal_3k_RData/19090553).
We here illustrate the CommPath steps for data from the tumor tissues. And analysis for data from the normal tissues would be roughly in the same manner.

```
# load(url("https://figshare.com/ndownloader/files/35185522"))
load("path_to_download/HCC.tumor.3k.RData")
```

This dataset consists of 3 varibles, 2 of which are required for CommPath input:
***tumor.expr*** : the expression matrix of gene * cell. Expression values are required to be first **normalized by the library-size and log-transformed**;
***tumor.label*** : a vector of labels indicating identity classes of cells in the expression matrix, and the order of labels should match the order of cells in the expression matrix.
The remaining one variable ***tumor.obj*** is the CommPath object created from ***tumor.expr*** and ***tumor.label*** and processed by CommPath standard procedures. This variable would be temporarily ignored since we will recreate it by the following steps. 
#### Identification of marker ligands and receptors
We start CommPath analysis by creating a CommPath object:

```
tumor.obj <- createCommPath(expr.mat = tumor.expr, 
		cell.info = tumor.label, 
		species = "hsapiens")
```

CommPath contains LR and pathway databases for human (hsapiens), mouse (mmusculus), rat (rnorvegicus), zebrafish (drerio), fruitfly (dmelanogaster), and worm (celegans) species.

Firstly we're supposed to identify marker ligands and receptors (ligands and receptors that are significantly highly expressed) for each cluster of cells in the expression matrix. CommPath provide **findLRmarker** to identify these markers by *t.test* or *wilcox.test*.

```
tumor.obj <- findLRmarker(object = tumor.obj, method = "wilcox.test")
```

#### Statistical identification of potential ligand-receptor (LR) associations
```
# To find significant LR pairs
tumor.obj <- findLRpairs(object = tumor.obj,
		logFC.thre = 0, 
		p.thre = 0.01)
```

The counts of Statistically significant LR pairs and overall interaction intensity among cell clusters are then stored in ***tumor.obj@interact[['InteractNumer']]***，and the detailed information of each LR pair is stored in ***tumor.obj@interact[['InteractGene']]***. 

Then we can visualize all interactions through a circos plot:
```
# To show the counts of LR associations among all clusters
# Here we set the parameter "filter" as FALSE, which means that those LR interactions are identified only based on their expression profiles, not filtered by pathways in the receiver cells (as described in the later sections)
pdf('circosPlot.count.nonfiltered.pdf',height=6,width=6)
circosPlot(object = tumor.obj, filter=FALSE)
dev.off()
```
<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/circosPlot.count.nonfiltered.png" height=300, width=300>

In the above circos plot, the directions of lines indicate the associations from ligands to receptors, and the widths of lines represent the counts of LR pairs among clusters.

```
# To show the overall interaction intensity of LR interactions among all clusters
pdf('circosPlot.intensity.nonfiltered.pdf',height=6,width=6)
circosPlot(object = tumor.obj, plot="intensity", filter=FALSE)
dev.off()
```
<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/circosPlot.intensity.nonfiltered.png" height=300, width=300>

Now the widths of lines represent the overall interaction intensity among clusters.
#### Pathway enrichment analysis
CommPath conducts pathway analysis to identify dysregulated signaling pathways containing the marker ligands and receptors for each cluster.

```
# To find pathways of which the genesets show overlap with the marker ligands and receptors
# CommPath provides pathway annotations from KEGG pathways, WikiPathways, reactome pathways, and GO terms
# Here we take the KEGG pathways as an example
tumor.obj <- findLRpath(object = tumor.obj, category = "kegg")
```

Now genesets showing overlap with the marker ligands and receptors are stored  in ***tumor.obj@pathway[['pathwayLR']]***. Then we score the pathways to measure the activation levels for each pathway in each cell.

```
# To compute pathway activation score by the gsva algorithm. For more information about the gsva algorithm, see the GSVA package (PMID: 23323831)
# Since the GSVA algorithm is time-consuming on large datasets, we also provide the UCell algorithm (PMID: 34285779) as an alternative. In addition, you can also assess the pathway activity in an average manner.
tumor.obj <- scorePath(object = tumor.obj, method = "gsva", min.size = 10)
```

After that CommPath provides **diffAllPath** to perform pathway differential activation analysis for each cluster and find the receptors and ligands in each pathway:

```
# To get significantly up-regulated pathways in each cluster
acti.path.dat <- diffAllPath(object = tumor.obj, only.posi = TRUE, only.sig = TRUE)
head(acti.path.dat)
```

There are several columns stored in the variable ***acti.path.dat***:

Columns ***mean.diff***, ***mean.1***, ***mean.2***, ***t***, ***df***, ***P.val*** and ***P.val.adj*** show the statistic results; ***description*** shows the name of pathway;

Column ***ligand.in.path*** and ***receptor.in.path*** show the marker ligands and receptors expressed by the current cluster and included in the current pathway;

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

Then counts of filtered significant LR pairs and the corresponding overall interaction intensity among cell clusters are stored in ***tumor.obj@interact.filter[['InteractNumer']]***，and the detailed information of each LR pair is stored in ***tumor.obj@interact.filter[['InteractGene']]***.

We next filter the pathways by removing those pathways containing only specific marker ligands which do not triger any pathway in the downstream clusters:

```
acti.path.filtered.dat <- filterPath(object = tumor.obj, acti.path.dat = acti.path.dat)
```

Filtered interactions could be also visualized through the circos plot:
```
# To show the counts of filtered LR associations among all clusters
pdf('circosPlot.count.filtered.pdf',height=6,width=6)
circosPlot(object = tumor.obj, filter=TRUE)
dev.off()
```
<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/circosPlot.count.filtered.png" height=300, width=300>

```
# To show the overall interaction intensity of filtered LR interactions among all clusters
pdf('circosPlot.intensity.filtered.pdf',height=6,width=6)
circosPlot(object = tumor.obj, plot="intensity", filter=TRUE)
dev.off()
```
<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/circosPlot.intensity.filtered.png" height=300, width=300>

Users would highlight the interaction of specific clusters. Here we take the Endothelial cells as an example:
```
select.ident = 'Endothelial'
pdf('circosPlot.Endothelial.count.filtered.pdf',height=6,width=6)
circosPlot(object = tumor.obj, select.ident = select.ident, filter=TRUE)
dev.off()
```
<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/circosPlot.Endothelial.count.filtered.png" height=300, width=300>

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

For the cluster of interest, CommPath also provides dot plots to investigate its upstream clusters which release specific ligands and its downstream clusters which express specific receptors: 
```
# To investigate the upstream clusters which release ligands to the selected cluster
pdf('dotPlot.ligand.Endothelial.pdf',height=5,width=10)
dotPlot.LR(object = tumor.obj, receptor.ident = select.ident)
dev.off()
```
<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/dotPlot.ligand.Endothelial.png" height=300, width=600>

```
# To investigate the downstream clusters which express receptors for the selected cluster
pdf('dotPlot.receptor.Endothelial.pdf',height=5,width=10.5)
dotPlot.LR(object = tumor.obj, ligand.ident = select.ident)
dev.off()
```
<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/dotPlot.receptor.Endothelial.png" height=300, width=630>

Then CommPath provides network graph tools to visualize the pathways and associated functional LR interactions:

```
# First to integrate the statistics of filtered activated pathways and their associated LR interactions 
tumor.obj <- pathNet(object = tumor.obj, acti.path.filtered.dat = acti.path.filtered.dat)
# To visualize the pathways and the associated upstream LR interactions in a network plot
pdf('pathNet.upstream.Endothelial.1.pdf',width=6,heigh=6)
set.seed(1234)
pathNetPlot(object = tumor.obj, select.ident =  select.ident, plot = "upstream",
    layout = 'layout.davidson.harel',
    vert.size.LR = 3, vert.size.path.adj = 10, 
    LR.label = 'R', vertex.label.cex.LR=0.25, vertex.label.cex.path=0.3)
dev.off()
```
<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/pathNet.upstream.Endothelial.png" height=600, width=780>

In the above network graph, the pie charts represent the activated pathways in the selected  cells (here Endothelial cells) and the scatter points represent the LR pairs of which the receptors are included in the genesets of the linked pathways. Colors of scatter points indicate the upstream clusters releasing the corresponding ligands. Sizes of pie charts indicate their total in-degree and the proportions indicate the in-degree from different upstream clusters. This plot would provide lots of important and detailed information of signaling pathways and the relevant LR interactions along the communication chains. To ensure the flexibility and interpretability, CommPath functions provide multiple optional arguments for selecting and limiting pathways and LR pairs to depict. 
```
To select the top 5 significantly upregulated pathways in Endothelial cells
pdf('pathNet.upstream.Endothelial.top5.pdf',width=6,heigh=6)
set.seed(1234)
pathNetPlot(object = tumor.obj, select.ident =  select.ident, plot = "upstream",
    top.n.path=5,
    layout = 'layout.davidson.harel',
    vert.size.LR = 3, vert.size.path.adj = 10, 
    LR.label = 'R', vertex.label.cex.LR=0.25, vertex.label.cex.path=0.3)
dev.off()
```

<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/pathNet.upstream.Endothelial.top5.png" height=416, width=550>

```
To select and show pathways of interest in the diagram by names of pathways directly
pdf('pathNet.upstream.Endothelial.selectpathway.pdf',width=6,heigh=6)
set.seed(1234)
pathNetPlot(object = tumor.obj, select.ident =  select.ident, plot = "upstream",
    select.path = c("JAK-STAT signaling pathway", "MAPK signaling pathway", "Ras signaling pathway", "PI3K-Akt signaling pathway", "AMPK signaling pathway"),
    layout = 'layout.davidson.harel',
    vert.size.LR = 3, vert.size.path.adj = 10, 
    LR.label = 'R', vertex.label.cex.LR=0.25, vertex.label.cex.path=0.3)
dev.off()
```

<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/pathNet.upstream.Endothelial.selectpathway.png" height=475, width=550>

```
# Also to visualize the pathways and the associated downstream LR interactions in a network plot
pdf('pathNet.downstream.Endothelial.pdf',width=6,heigh=6)
set.seed(1234)
pathNetPlot(object = tumor.obj, select.ident =  select.ident, plot = "downstream",
    layout = 'layout.davidson.harel',
    vert.size.LR = 3, vert.size.path.adj = 10, 
    LR.label = 'L', vertex.label.cex.LR=0.25, vertex.label.cex.path=0.3)
dev.off()
```
<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/pathNet.downstream.Endothelial.png" height=600, width=750>

The legend of the above network graph is generally the same to that of the previous network plot, except that: (i) the scatter points represent the LR pairs of which the ligands are included in the genesets of the linked pathways; (ii) colors of scatter points indicate the downstream clusters expressing the corresponding receptors; (iii) sizes of pie charts indicate their total out-degree and the proportions indicate the out-degree to different downstream clusters.

CommPath also provides dot plot to investigate the upstream and downstream LR pairs involved in the specific pathways in the selected clusters:
```
pathway = "Ras signaling pathway"
# To visualize the upstream LR pairs of which the receptors are expressed by Endothelial cells and are included in the pathway "Ras signaling pathway"
pdf('dotPlot.Ras.pathway.upstream.Endothelial.pdf',height=5,width=7.5)
dotPlot.pathway(object = tumor.obj, pathway = pathway, acti.path.filtered.dat = acti.path.filtered.dat, receptor.ident = select.ident, top.n.inter = 10)
dev.off()
```
<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/dotPlot.Ras.pathway.upstream.Endothelial.png" height=300, width=430>

```
# To visualize the downstream LR pairs of which the ligands are released by Endothelial cells and are included in the pathway "Ras signaling pathway"
pdf('dotPlot.Ras.pathway.downstream.Endothelial.pdf',height=5,width=7.5)
dotPlot.pathway(object = tumor.obj, pathway = pathway, acti.path.filtered.dat = acti.path.filtered.dat, ligand.ident = select.ident, top.n.inter = 10)
dev.off()
```
<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/dotPlot.Ras.pathway.downstream.Endothelial.png" height=300, width=430>

#### Identification of pathway-mediated cell-cell communication chain
For a specific cell cluster, here named as B for demonstration, CommPath identifies the upstream cluster A sending signals to B, the downstream cluster C receiving signals from B, and the significantly activated pathways in B to mediate the A-B-C communication chain. More exactly, through LR and pathways analysis described above, CommPath is able to identify LR pairs between A and B, LR pairs between B and C, and pathways activated in B. Then CommPath screens for pathways in B which involve both the receptors to interact with A and ligands to interact with C.
```
# To investigate the activated patways and the associated receptors for a specific cluster
select.ident = 'Endothelial'
pdf('pathPlot.Endothelial.pdf',height=6,width=10)
pathPlot(object = tumor.obj, 
    select.ident = select.ident, 
    acti.path.dat = acti.path.filtered.dat)
dev.off()
```

<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/pathPlot.Endothelial.png" height=300, width=400>

In the above line plot, the widths of lines between ***Upstream*** cluster and ***Receptor*** represent the overall interaction intensity between the upstream  cluster and Endothelial cells via the specific receptors; the sizes and colors of dots in the ***Receptor*** column represent the average *log2FC* and *-log10(P)* from differential expression tests comparing the receptor  expression in Endothelial cells to that in all other cells; the lengths and colors of bars in the ***Pathway annotation*** column represent the mean difference and *-log10(P)* form differential activation tests comparing the pathway scores in Endothelial cells to those in all other cells.
```
# To select self-defind cluster for the  central cell population
# And select pathways of interest by names
# Also select the top N receptors by interaction intensity 
select.ident = 'Endothelial'
pdf('pathPlot.Endothelial.selectpathway.pdf',height=4,width=10)
pathPlot(object = tumor.obj, 
    select.ident = select.ident, 
    up.ident = c("Treg","Myeloid","Endothelial","Fibroblast","Hepatocyte"),
    select.path = c("JAK-STAT signaling pathway", "MAPK signaling pathway", "Ras signaling pathway", "PI3K-Akt signaling pathway", "AMPK signaling pathway"),
    top.n.receptor = 5,
    acti.path.dat = acti.path.filtered.dat)
dev.off()
```

<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/pathPlot.Endothelial.selectpathway.png" height=200, width=400>


```
# To investigate the activated patways, the associated receptors and ligands for a specific cluster
pdf('pathChainPlot.Endothelial.pdf',height=6,width=14)
pathChainPlot(object = tumor.obj, 
    select.ident = select.ident, 
    up.ident = c("Treg","Myeloid","Endothelial","Fibroblast","Hepatocyte"),
    down.ident = c("B","Myeloid","Endothelial","Fibroblast","Hepatocyte"),
    select.path = c("JAK-STAT signaling pathway", "MAPK signaling pathway", "Ras signaling pathway", "PI3K-Akt signaling pathway", "AMPK signaling pathway"),
    acti.path.dat = acti.path.filtered.dat)
dev.off()
```

<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/pathChainPlot.Endothelial.png" height=300, width=600>

The legend of the above line plot is generally the same to that of the previous plot from **pathPlot**.

#### Comparison of cell-cell communication between two conditions
CommPath also provides useful utilities to compare cell-cell communication between two conditions such as disease and control. Here we, for example, use CommPath to compare the cell-cell communication between cells from HCC tumor and normal tissues. We have pre-created the CommPath object for the normal samples following the above steps, and the data are also available in [figshare](https://figshare.com/articles/dataset/HCC_tumor_normal_3k_RData/19090553).

```
# load(url("https://figshare.com/ndownloader/files/35185525"))
load("path_to_download/HCC.normal.3k.RData")
```

This dataset consists of 3 varibles:
***normal.expr*** : expression matrix for cells from normal tissues;
***normal.label*** : identity labels for cells from normal tissues;
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
pdf('pathPlot.compare.Endothelial.pdf',height=6,width=10)
pathPlot.compare(object.1 = tumor.obj, 
		object.2 = normal.obj, 
		select.ident = select.ident,
		up.ident = c("Treg","Myeloid","Endothelial","Fibroblast","Hepatocyte"),
		diff.path.dat = diff.path.dat)
dev.off()
```

<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/pathPlot.compare.Endothelial.png" height=300, width=400>

In the above line plot, the widths of lines between ***Upstream*** cluster and ***Receptor*** represent the overall interaction intensity between the upstream clusters and Endothelial cells via the specific receptors, and the colors indicate the interaction intensity is upregulated (red) or downregulated (blue) in tumor tissues (object.1) compared to that in normal tissues (object.2); the sizes and colors of dots in the ***Receptor*** column represent the average *log2FC* and *-log10(P)* of expression of receptors in Endothelial cells compared to all other cells in tumor tissues; the lengths and colors of bars in the ***Pathway annotation*** column represent the mean difference and *-log10(P)* of pathway scores of Endothelial cells in tumor tissues compared to that in normal tissues.

```
# To compare the pathway-mediated cell-cell communication chain for a specific cluster between 2 CommPath objects
pdf('pathChainPlot.compare.Endothelial.pdf',height=6,width=14)
pathChainPlot.compare(object.1 = tumor.obj, 
		object.2 = normal.obj, 
		select.ident = select.ident, 
		up.ident = c("Myeloid","Endothelial","Fibroblast","Hepatocyte"),
    down.ident = c("CD8 T","CD4 T","Treg","NK","Endothelial","Fibroblast","Hepatocyte"),
    select.path = c("Gastric cancer", "Relaxin signaling pathway", "Axon guidance"),
		diff.path.dat = diff.path.dat)
dev.off()
```

<img src="https://github.com/yingyonghui/SupplementaryData/blob/main/CommPath/tutorial_pic/pathChainPlot.compare.Endothelial.png" height=300, width=600>

The legend of the above line plot is generally the same to that of the previous plot from **pathPlot.compare**.

#### sessionInfo()
```
R version 4.4.3 (2025-02-28)
Platform: x86_64-conda-linux-gnu
Running under: Ubuntu 22.04.4 LTS

Matrix products: default

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

time zone: Asia/Shanghai
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
 [1] GSVA_2.0.5      UCell_2.10.1    CommPath_1.1.0  igraph_2.1.4
 [5] Matrix_1.7-2    ggplot2_3.5.1   dplyr_1.1.4     plyr_1.8.9
 [9] reshape2_1.4.4  circlize_0.4.16

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.1            blob_1.2.4
 [3] Biostrings_2.74.1           fastmap_1.2.0
 [5] SingleCellExperiment_1.28.1 XML_3.99-0.17
 [7] rsvd_1.0.5                  lifecycle_1.0.4
 [9] KEGGREST_1.46.0             RSQLite_2.3.9
[11] magrittr_2.0.3              compiler_4.4.3
[13] rlang_1.1.5                 tools_4.4.3
[15] data.table_1.17.0           S4Arrays_1.6.0
[17] bit_4.5.0.1                 DelayedArray_0.32.0
[19] abind_1.4-8                 BiocParallel_1.40.0
[21] HDF5Array_1.34.0            withr_3.0.2
[23] BiocGenerics_0.52.0         grid_4.4.3
[25] stats4_4.4.3                beachmat_2.22.0
[27] xtable_1.8-4                colorspace_2.1-1
[29] Rhdf5lib_1.28.0             scales_1.3.0
[31] SummarizedExperiment_1.36.0 cli_3.6.4
[33] crayon_1.5.3                generics_0.1.3
[35] rjson_0.2.23                httr_1.4.7
[37] DBI_1.2.3                   cachem_1.1.0
[39] rhdf5_2.50.2                stringr_1.5.1
[41] zlibbioc_1.52.0             parallel_4.4.3
[43] AnnotationDbi_1.68.0        XVector_0.46.0
[45] matrixStats_1.5.0           vctrs_0.6.5
[47] jsonlite_1.9.1              BiocSingular_1.22.0
[49] BiocNeighbors_2.0.1         IRanges_2.40.1
[51] S4Vectors_0.44.0            bit64_4.6.0-1
[53] irlba_2.3.5.1               GSEABase_1.68.0
[55] magick_2.8.5                annotate_1.84.0
[57] glue_1.8.0                  codetools_0.2-20
[59] stringi_1.8.4               shape_1.4.6.1
[61] gtable_0.3.6                GenomeInfoDb_1.42.3
[63] GenomicRanges_1.58.0        UCSC.utils_1.2.0
[65] ScaledMatrix_1.14.0         munsell_0.5.1
[67] tibble_3.2.1                pillar_1.10.1
[69] rhdf5filters_1.18.0         graph_1.84.1
[71] GenomeInfoDbData_1.2.13     R6_2.6.1
[73] sparseMatrixStats_1.18.0    lattice_0.22-6
[75] Biobase_2.66.0              png_0.1-8
[77] SpatialExperiment_1.16.0    memoise_2.0.1
[79] Rcpp_1.0.14                 SparseArray_1.6.2
[81] MatrixGenerics_1.18.1       pkgconfig_2.0.3
[83] GlobalOptions_0.1.2     
```
## Citation
[Lu H, Ping J, Zhou G, et al. CommPath: An R package for inference and analysis of pathway-mediated cell-cell communication chain from single-cell transcriptomics. Comput Struct Biotechnol J. 2022;20:5978-5983.](https://www.sciencedirect.com/science/article/pii/S2001037022004755)