#data('CommPathData', package=('CommPath'))

#' To identify marker ligands and marker receptors in the expression matrix
#' @param object CommPath object
#' @param method Method used for differential expression test, either 'wilcox.test' or 't.test'
#' @param p.adjust Method used for p value correction for multiple differential expression test; see p.adjust function for more information
#' @importFrom stats wilcox.test t.test p.adjust
#' @return Data frame containing the differential expression test
#' @export
findLRmarker <- function(object, method='wilcox.test', p.adjust='BH'){

	species <- object@meta.info$species
	expr.mat <- as.matrix(object@data)
	label <- object@cell.info[,'Cluster']
	
	lr.pair.dat <- CommPathData$DataLR[[species]]
	all.lig.reps <- unique(c(lr.pair.dat$L, lr.pair.dat$R))

	expr.mat <- expr.mat[which(rownames(expr.mat) %in% all.lig.reps), ]
	if (nrow(expr.mat)==0){
		stop("There is no ligand or receptor detected in the expression matrix!")
	}

	if (!is.factor(label)){ label <- factor(label) }
	ident.level <- levels(label)
	
	label <- as.character(label)
	### p.value calculated
	if (method!='wilcox.test' & method!='t.test'){
		stop("Select t.test or wilcox.test to conduct differential analysis")
	}

	if(method=='wilcox.test'){
		test.all.res <- lapply(ident.level, function(each.level) {
			message(paste0('Identifying marker genes for cluster ',each.level,' ...'))
			cell.ident <- which(label==each.level)
			cell.other <- which(label!=each.level)

			test.res <- apply(expr.mat, 1, function(row.expr){
				logFC <- log(mean(expm1(row.expr[cell.ident])) +1, base=2) - log(mean(expm1(row.expr[cell.other])) +1, base=2)
				p.value <- wilcox.test(x=row.expr[cell.ident], y=row.expr[cell.other])$p.value

				pct.1 <- length(which(row.expr[cell.ident] > 0))/length(cell.ident)
				pct.2 <- length(which(row.expr[cell.other] > 0))/length(cell.other)
				c(logFC, p.value, pct.1, pct.2)
			})

			test.res <- as.data.frame(t(test.res))
			colnames(test.res) <- c('avg_log2FC','p_val', 'pct.1', 'pct.2')
			test.res <- subset(test.res, pct.1>0 | pct.2>0)
			test.res$p_val_adj <- p.adjust(test.res$p_val, method=p.adjust)
			test.res$cluster <- each.level
			test.res$gene <- rownames(test.res)
			return(test.res)

		})
		test.all.res <- do.call(rbind, test.all.res)
	}else{
		test.all.res <- lapply(ident.level, function(each.level) {
			message(paste0('Identify marker genes for ',each.level,' ...'))
			cell.ident <- which(label==each.level)
			cell.other <- which(label!=each.level)

			test.res <- apply(expr.mat, 1, function(row.expr){
				logFC <- log(mean(expm1(row.expr[cell.ident])) +1, base=2) - log(mean(expm1(row.expr[cell.other])) +1, base=2)
				p.value <- t.test(x=row.expr[cell.ident], y=row.expr[cell.other])$p.value

				pct.1 <- length(which(row.expr[cell.ident] > 0))/length(cell.ident)
				pct.2 <- length(which(row.expr[cell.other] > 0))/length(cell.other)
				c(logFC, p.value, pct.1, pct.2)
			})

			test.res <- as.data.frame(t(test.res))
			colnames(test.res) <- c('avg_log2FC','p_val', 'pct.1', 'pct.2')
			test.res <- subset(test.res, pct.1>0 | pct.2>0)
			test.res$p_val_adj <- p.adjust(test.res$p_val, method=p.adjust)
			test.res$cluster <- each.level
			test.res$gene <- rownames(test.res)
			return(test.res)

		})
		test.all.res <- do.call(rbind, test.all.res)
	}

	object@LR.marker <- test.all.res
	return(object)
}

#' To find marker ligands and marker receptors
#' @param object CommPath object
#' @param logFC.thre logFC threshold, marker genes with a logFC > logFC.thre will be considered
#' @param p.thre p threshold, marker genes with a adjust p value < p.thre will be considered
#' @importFrom reshape2 melt
#' @return List containing the ligand-receptor interaction information
#' @export
findLRpairs <- function(object, logFC.thre=0, p.thre=0.05){
	options(stringsAsFactors=F)
	marker.dat <- object@LR.marker
	species <- object@meta.info$species

	lr.pair.dat <- CommPathData$DataLR[[species]]
	ligs <- lr.pair.dat$L
	reps <- lr.pair.dat$R

	cluster.level <- unique(marker.dat$cluster)
	num.cluster <- length(cluster.level)
	Interact.num.mat <- matrix(0,num.cluster,num.cluster)
	Interact.gene.mat <- matrix(NA,num.cluster,num.cluster)
	Interact.lig.mat <- matrix(NA,num.cluster,num.cluster)
	Interact.rep.mat <- matrix(NA,num.cluster,num.cluster)
	marker.lig.dat <- as.data.frame(matrix(NA,0,ncol(marker.dat)))
	marker.rep.dat <- as.data.frame(matrix(NA,0,ncol(marker.dat)))
	for (lig.idx in 1:num.cluster){
		cluster.l <- cluster.level[lig.idx]
		markers.l <- marker.dat[marker.dat$cluster == cluster.l, ]
		ligands <- markers.l[(markers.l$avg_log2FC > logFC.thre) & (markers.l$p_val_adj < p.thre), 'gene']
		for (rep.idx in 1:num.cluster){
			cluster.r <- cluster.level[rep.idx]
			markers.r <- marker.dat[marker.dat$cluster == cluster.r, ]
			receptors <- markers.r[(markers.r$avg_log2FC>logFC.thre) & (markers.r$p_val_adj<p.thre), 'gene']
			pair.valid <- which((ligs %in% ligands) & (reps %in% receptors))
			if (length(pair.valid) > 0){
				Interact.gene.mat[lig.idx,rep.idx] <- paste(paste(ligs[pair.valid],reps[pair.valid],sep='--'),collapse=';')
				Interact.lig.mat[lig.idx,rep.idx] <- paste(ligs[pair.valid],collapse=';')
				Interact.rep.mat[lig.idx,rep.idx] <- paste(reps[pair.valid],collapse=';')
				marker.lig.dat <- rbind(marker.lig.dat,markers.l[markers.l$gene %in% ligs[pair.valid],])
				marker.rep.dat <- rbind(marker.rep.dat,markers.r[markers.r$gene %in% reps[pair.valid],])
			}
			Interact.num.mat[lig.idx,rep.idx] <- length(pair.valid)
		}
	}

	if (sum(Interact.num.mat)==0){
		stop('No LR pairs detected! Select a lower logFC.thre or higher p.thre and try again')
	}
	rownames(Interact.num.mat) <- cluster.level
	colnames(Interact.num.mat) <- cluster.level
	Interact.num.dat <- melt(Interact.num.mat, varnames=c('Cell.From','Cell.To'),value.name="LR.Count",  na.rm=TRUE)
	Interact.num.dat <- factor.to.character(Interact.num.dat)
	rownames(Interact.num.dat) <- paste(Interact.num.dat$Cell.From,Interact.num.dat$Cell.To, sep='--')

	rownames(Interact.gene.mat) <- cluster.level
	colnames(Interact.gene.mat) <- cluster.level
	Interact.gene.dat <- melt(Interact.gene.mat,varnames=c('Cell.From','Cell.To'),value.name="LR.Info", na.rm = TRUE)
	Interact.gene.dat <- factor.to.character(Interact.gene.dat)
	rownames(Interact.gene.dat) <- paste(Interact.gene.dat$Cell.From,Interact.gene.dat$Cell.To, sep='--')
	
	rownames(Interact.lig.mat) <- cluster.level
	colnames(Interact.lig.mat) <- cluster.level
	Interact.lig.dat <- melt(Interact.lig.mat,varnames=c('Cell.From','Cell.To'),value.name="L.Info", na.rm = TRUE)
	Interact.lig.dat <- factor.to.character(Interact.lig.dat)

	rownames(Interact.rep.mat) <- cluster.level
	colnames(Interact.rep.mat) <- cluster.level
	Interact.rep.dat <- melt(Interact.rep.mat,varnames=c('Cell.From','Cell.To'),value.name="R.Info", na.rm = TRUE)
	Interact.rep.dat <- factor.to.character(Interact.rep.dat)
	
	Interact.gene.dat$L.Info <- Interact.lig.dat$L.Info
	Interact.gene.dat$R.Info <- Interact.rep.dat$R.Info
	rownames(Interact.gene.dat) <- 1:nrow(Interact.gene.dat)

	LRpair <- Interact.gene.dat$LR.Info
	names(LRpair) <- paste(Interact.gene.dat$Cell.From,Interact.gene.dat$Cell.To,sep='--')
	lr.split.list <- sapply(LRpair, function(LR){
		strsplit(LR,split=';')
	})
	ident.pair <- names(lr.split.list)
	for (each.pair in ident.pair){
		lr.split.list[[each.pair]] <- paste(each.pair, lr.split.list[[each.pair]],sep='--')
	}
	lr.split.list <- sapply((unlist(lr.split.list)), function(ident.LR.info){
		strsplit(ident.LR.info, split='--')
		})
	lr.unfold.dat <- as.data.frame(t(as.data.frame(lr.split.list)))
	colnames(lr.unfold.dat) <- c('Cell.From','Cell.To','Ligand','Receptor')
	rownames(lr.unfold.dat) <- paste(lr.unfold.dat$Cell.From,lr.unfold.dat$Cell.To,lr.unfold.dat$Ligand,lr.unfold.dat$Receptor, sep='--')

	marker.lig.dat <- unique(marker.lig.dat)
	marker.rep.dat <- unique(marker.rep.dat)

	fc.lig <- as.vector(apply(lr.unfold.dat, 1, function(x){ subset(marker.lig.dat, gene == x['Ligand'] & cluster == x['Cell.From'])[1,'avg_log2FC'] }))
	pval.lig <- as.vector(apply(lr.unfold.dat, 1, function(x){ subset(marker.lig.dat, gene == x['Ligand'] & cluster == x['Cell.From'])[1,'p_val_adj'] }))
	fc.rep <- as.vector(apply(lr.unfold.dat, 1, function(x){ subset(marker.rep.dat, gene == x['Receptor'] & cluster == x['Cell.To'])[1,'avg_log2FC'] }))
	pval.rep <- as.vector(apply(lr.unfold.dat, 1, function(x){ subset(marker.rep.dat, gene == x['Receptor'] & cluster == x['Cell.To'])[1,'p_val_adj'] }))

	lr.unfold.dat$Log2FC.LR <- fc.lig * fc.rep
	lr.unfold.dat$P.val.LR <- 1-(1-pval.lig)*(1-pval.rep)
	lr.unfold.dat$P.val.adj.LR <- p.adjust(lr.unfold.dat$P.val.LR, method='BH')
	
	lr.inten.all.name <- paste(lr.unfold.dat$Cell.From, lr.unfold.dat$Cell.To, sep='--')
	lr.inten.vec <- by(data=lr.unfold.dat$Log2FC.LR, INDICES=lr.inten.all.name, FUN=sum)
	lr.inten.vec.name <- names(lr.inten.vec)
	lr.inten.vec <- as.vector(lr.inten.vec)

	Interact.num.dat$Intensity <- lr.inten.vec[match(rownames(Interact.num.dat), lr.inten.vec.name)]
	Interact.num.dat$Intensity[is.na(Interact.num.dat$Intensity)] <- 0

	Interact <- list(InteractNumer=Interact.num.dat,InteractGene=Interact.gene.dat, InteractGeneUnfold=lr.unfold.dat, markerL=marker.lig.dat, markerR=marker.rep.dat)

	object@interact <- Interact

	object@meta.info$logFC.thre <- logFC.thre
	object@meta.info$p.thre <- p.thre
	return(object)
}


#' To find those pathways in which the genesets show overlap with the marker ligand and receptor genes in our dataset
#' @param object CommPath object
#' @param category Character to indicate which pathway to investigate; one of "go" (GO terms), "kegg" (for KEGG pathways), 'wiki' (for WikiPathways), and "reactome" (for reactome pathways), or "all" for all pathways
#' @return CommPath object containing the ligand-receptor interaction information and the pathways showing overlap with the marker ligand and receptor genes in the dataset
#' @export
findLRpath <- function(object, category='all'){
	options(stringsAsFactors=F)
	if (category!='all' & category!='go' & category!='kegg' & category!='wiki' & category!='reactome'){
		stop("Wrong category selected. Select one of 'go', 'kegg', 'wiki', and 'reactome', or 'all' for all pathways")
	}
	species <- object@meta.info$species

	### select one category
	if (category=='all'){
		path.go.list <- CommPathData$DataGOterm[[species]]
		path.kegg.list <- CommPathData$DataKEGGpathway[[species]]
		path.wiki.list <- CommPathData$DataWikiPathway[[species]]
		path.reac.list <- CommPathData$DataReactome[[species]]
		path.list <- c(path.go.list, path.kegg.list, path.wiki.list, path.reac.list)
		path.list <- path.list[which(!duplicated(names(path.list)))]

	}else if (category=='go'){
		path.list <- CommPathData$DataGOterm[[species]]
	}else if(category=='kegg'){
		path.list <- CommPathData$DataKEGGpathway[[species]]
	}else if(category=='wiki'){
		path.list <- CommPathData$DataWikiPathway[[species]]
	}else if(category=='reactome'){
		path.list <- CommPathData$DataReactome[[species]]
	}

	marker.lig.dat <- object@interact$markerL
	marker.rep.dat <- object@interact$markerR
	lr.gene <- unique(c(marker.lig.dat$gene,marker.rep.dat$gene))
	
	which.overlap.list <- unlist(
		lapply(path.list, function(x){ 
			if(any(x %in% lr.gene)){ return(TRUE) }else{return(FALSE)}
		}))
	path.lr.list <- path.list[which(which.overlap.list)]
	object@pathway$pathwayLR <- path.lr.list
	return(object)
}

#' To find different enriched pathway between two group cells
#' @param object CommPath object
#' @param method Method used for scoring the pathways, either 'gsva' of 'average'
#' @param min.size Minimum size of overlaping genes between candidate pathways and the expression matrix
#' @param ... Extra parameters passed to gsva
#' @return CommPath object with pathways activation scores stored in the slot pathway
#' @export
scorePath <- function(object, method='gsva', min.size=10, ...){
	expr.mat <- as.matrix(object@data)
	path.list <- object@pathway$pathwayLR
	if (is.null(path.list)){
		stop("No pathway detected, run findLRpath befor scorePath")
	}

	if (method=='gsva'){
		acti.score <- GSVA::gsva(expr.mat, path.list, min.sz=min.size, ...)
	}else if(method=='average'){
		acti.score <- t(as.data.frame(lapply(path.list, function(eachPath){
			overlap.gene <- intersect(eachPath, rownames(expr.mat))
			if (length(overlap.gene) < 10){
				return(rep(NA, ncol(expr.mat)))
			}else{
				return(colMeans(expr.mat[overlap.gene, ]))
			}
		})))
		rownames(acti.score) <- names(path.list)
		acti.score <- acti.score[which(!is.na(acti.score[,1])), ]
	}else{
		stop('Select "gsva" or "average" for pathway scoring')
	}
	object@pathway$acti.score <- acti.score
	object@pathway$method <- method
	return(object)
}


#' To find different enriched pathway between two group cells
#' @param object CommPath object
#' @param select.ident.1 Identity class to define cells for group 1
#' @param select.ident.2 Identity class to define cells for group 2 for comparison; if 'NULL', use all other cells for comparison
#' @param method Method used for differential enrichment analysis, either 't.test' of 'wilcox.test'
#' @param only.posi only save the information of pathways showing up regulation
#' @param only.sig only save the information of pathways showing significant differences
#' @return Data frame including the statistic result comparing the pathway enrichment sorces between group 1 and group 2, the significant recetor and ligand of group 1 in the pathways, and the corresponding up stream identity class which interact with group 1 by releasing specific ligand
#' @export
diffPath <- function(object, select.ident.1, select.ident.2=NULL, method='t.test', only.posi=FALSE, only.sig=FALSE){
	options(stringsAsFactors=F)
	acti.score <- object@pathway$acti.score
	ident.label <- object@cell.info$Cluster
	
	### input parameter check
	if (method!='t.test' & method!='wilcox.test'){
		stop("Select t.test or wilcox.test to conduct differential analysis")
	}else if(method=='t.test'){
		test.res.dat <- data.frame(matrix(NA,0,12))
		colnames(test.res.dat) <- c('mean.diff','mean.1','mean.2','t','df','p.val','p.val.adj','description','cell.up','ligand.up','receptor.in.path','ligand.in.path')
	}else{
		test.res.dat <- data.frame(matrix(NA,0,11))
		colnames(test.res.dat) <- c('median.diff','median.1','median.2','W','p.val','p.val.adj','description','cell.up','ligand.up','receptor.in.path','ligand.in.path')
	}

	### find those pathways showing overlap with the ligands and receptors in select ident
	path.lr.list <- object@pathway$pathwayLR
	if (is.null(path.lr.list)){
		stop("No pathway detected, run findLRpath befor diffPath")
	}

	marker.lig.rep <- subset(object@interact$markerL, cluster %in% select.ident.1)[, 'gene']
	marker.lig.rep <- unique(c(marker.lig.rep, subset(object@interact$markerR, cluster %in% select.ident.1)[, 'gene']))

	which.overlap.list <- unlist(
		lapply(path.lr.list, function(x){ 
			if(any(x %in% marker.lig.rep)){ return(TRUE) }else{ return(FALSE) }
		}))
	path.overlap.list <- path.lr.list[which(which.overlap.list)]
	if (length((path.overlap.list)) == 0){
		warning(paste0('There is no pathway showing overlap with the marker ligands and receptors of cluster ',select.ident.1))
		return(test.res.dat)
	}
	### test for overlaping pathways
	### since we set min.sz in the gsva function, there are some pathways not calculated in the gsva process, we shall remove those pathways
	path.overlap.list <- path.overlap.list[names(path.overlap.list) %in% rownames(acti.score)]
	acti.ident.score <- acti.score[names(path.overlap.list),]

	### if only 1 name in path.overlap.list, then acti.ident.score would be a vector
	### convert it to a matrix
	if (is.null(dim(acti.ident.score))){
		acti.ident.score <- t(as.matrix(acti.ident.score))
		rownames(acti.ident.score) <- names(path.overlap.list)
	}
	
	group <- as.character(ident.label)
	test.res.dat <- pathTest(acti.ident.score, group, select.ident.1, select.ident.2, method, only.posi=only.posi, only.sig=only.sig)

	### to find the ligand in the same pathway
	ident.lig.vec <- subset(object@interact$markerL, cluster %in% select.ident.1)[, 'gene']
	### for each pathway in the DEG result, find which pathway show overlap with marker ligands of the selected ident
	ident.lig.in.path <- sapply(test.res.dat$description, function(eachPath){
		each.set <- path.overlap.list[[eachPath]]
		if (any(ident.lig.vec %in% each.set)){
			return(paste(ident.lig.vec[ident.lig.vec %in% each.set],collapse=';'))
		}else{
			return(NA)
		}
	})
	test.res.dat$ligand.in.path <- unlist(ident.lig.in.path)

	### to find the receptor in the same pathway
	ident.rep.vec <- subset(object@interact$markerR, cluster %in% select.ident.1)[, 'gene']
	### for each pathway in the DEG result, find which pathway show overlap with marker receptors of the selected ident
	ident.rep.in.path <- sapply(test.res.dat$description, function(eachPath){
		each.set <- path.overlap.list[[eachPath]]
		if (any(ident.rep.vec %in% each.set)){
			return(paste(ident.rep.vec[ident.rep.vec %in% each.set],collapse=';'))
		}else{
			return(NA)
		}
	})
	test.res.dat$receptor.in.path <- unlist(ident.rep.in.path)

	return(test.res.dat)
}


#' To find differentially activated pathways in each identity class 
#' @param object CommPath object
#' @param method Method used for differential enrichment analysis, either 't.test' of 'wilcox.test'
#' @param only.posi only save the information of pathways showing up regulation
#' @param only.sig only save the information of pathways showing significant differences
#' @return Data frame including the statistic result comparing the pathway enrichment sorces between cells in each cluster and all other clusters, the significant recetor and ligand in the pathways, and the corresponding up stream identity class and ligand
#' @export
diffAllPath <- function(object, method='t.test', only.posi=FALSE, only.sig=FALSE){
	if (method=='t.test'){
		all.test.dat <- data.frame(matrix(NA,0,12))
	}else{
		all.test.dat <- data.frame(matrix(NA,0,11))
	}

	all.ident <- object@cell.info$Cluster
	if (!is.factor(all.ident)){ all.ident <- factor(all.ident) }
	unique.label <- levels(all.ident)
	for (each.ident in unique.label){
		message(paste0('Identifying pathways for cluster ',each.ident,'...'))
		test.res.dat <- diffPath(object, select.ident.1=each.ident, select.ident.2=NULL, method=method, only.posi=only.posi, only.sig=only.sig)
		if (nrow(test.res.dat)==0){ 
			next 
		}
		test.res.dat$cluster <- each.ident
		all.test.dat <- rbind(all.test.dat, test.res.dat)
	}

	return(all.test.dat)
}


#' To find the downstream identity class of specific ligand released by specific upstream identity class
#' @param object CommPath object
#' @param select.ident Upstream identity class; if 'NULL', use all identity classes
#' @param select.ligand Ligand released by upstream identity class; if 'NULL', use all ligands that are markers for the selected upstream identity class
#' @return Data frame including the interaction information
#' @export
findReceptor <- function(object, select.ident=NULL, select.ligand=NULL){
	options(stringsAsFactors=F)

	if (is.null(select.ident) & is.null(select.ligand)){
		stop("Either a select.ident or a select.ligand need to be asigned")
	}
	ident.down.dat <- object@interact$InteractGeneUnfold

	if (!is.null(select.ligand)){
		ident.down.dat <- subset(ident.down.dat, Ligand %in% select.ligand)
	}
	if (!is.null(select.ident)){
		ident.down.dat <- subset(ident.down.dat, Cell.From %in% select.ident)
	}

	if (nrow(ident.down.dat)==0){
		warning("No downstream ident found for the selected ident and ligand")
	}
	return(ident.down.dat)
}


#' To find the upstream identity classes and ligands of specific receptor expressed by specific downstream identity class
#' @param object CommPath object
#' @param select.ident Downstream identity class; if 'NULL', use all identity classes
#' @param select.receptor Receptor expressed by downstream identity class; if 'NULL', use all receptors that are markers for the selected downstream identity class
#' @return Data frame including the interaction information
#' @export
findLigand <- function(object, select.ident=NULL, select.receptor=NULL){
	options(stringsAsFactors=F)

	if (is.null(select.ident) & is.null(select.receptor)){
		stop("Either a select.ident or a select.receptor need to be asigned")
	}

	ident.up.dat <- object@interact$InteractGeneUnfold
	if (!is.null(select.receptor)){
		ident.up.dat <- subset(ident.up.dat, Receptor %in% select.receptor)
	}

	if (!is.null(select.ident)){
		ident.up.dat <- subset(ident.up.dat, Cell.To %in% select.ident)
	}

	if (nrow(ident.up.dat)==0){
		warning("No upstream ident found for the selected ident and receptor")
	}
	return(ident.up.dat)
}


#' Differential enrichment analysis by t.test or wilcox.txt
#' @param acti.ident.score Matrix of pathway scores, pathway * cell
#' @param group Vector of group labels of cells
#' @param select.ident.1 Identity class 1
#' @param select.ident.2 Identity class 2 for comparison
#' @param method Method for hypothesis test, either 't.test' or 'wilcox.test'
#' @param only.posi only save the information of pathways showing up regulation
#' @param only.sig only save the information of pathways showing significant differences
#' @importFrom stats median
#' @return Data frame including the statistic result
#' @export
pathTest <- function(acti.ident.score, group, select.ident.1, select.ident.2=NULL, method='t.test', only.posi=FALSE, only.sig=FALSE){
	if(method=='t.test'){
		if (is.null(select.ident.2)){
			t.result <- apply(acti.ident.score,1,function(geneExpr){
				t.test(x=geneExpr[group %in% select.ident.1],y=geneExpr[!(group %in% select.ident.1)])})
		}else{
			t.result <- apply(acti.ident.score,1,function(geneExpr){
				t.test(x=geneExpr[group %in% select.ident.1],y=geneExpr[group %in% select.ident.2])})
		}
		# t = testRes$statistic, df = testRes$parameter, mean.1 = testRes$estimate[1], mean.2 = testRes$estimate[2], pVal = testRes$p.value
		test.res.dat <- as.data.frame(lapply(t.result,function(testRes){
			c(testRes$estimate[1]-testRes$estimate[2],testRes$estimate[1],testRes$estimate[2],testRes$statistic,testRes$parameter,testRes$p.value)
		}))
		test.res.dat <- as.data.frame(t(test.res.dat))
		colnames(test.res.dat) <- c('mean.diff','mean.1','mean.2','t','df','p.val')
	}else{
		if (is.null(select.ident.2)){
			wil.result <- apply(acti.ident.score,1,function(geneExpr){
				wilcox.test(x=geneExpr[group %in% select.ident.1],y=geneExpr[!(group %in% select.ident.1)])})

			wil.median <- apply(acti.ident.score, 1, function(geneExpr){
				median.1 <- median(geneExpr[group %in% select.ident.1])
				median.2 <- median(geneExpr[!(group %in% select.ident.1)])
				median.diff <- median.1 - median.2
				return(c(median.diff, median.1, median.2))
			})
			wil.median <- as.data.frame(t(wil.median))
			colnames(wil.median) <- c('median.diff','median.1','median.2')

		}else{
			wil.result <- apply(acti.ident.score,1,function(geneExpr){
				wilcox.test(x=geneExpr[group %in% select.ident.1],y=geneExpr[group %in% select.ident.2])})
			wil.median <- apply(acti.ident.score, 1, function(geneExpr){
				median.1 <- median(geneExpr[group %in% select.ident.1])
				median.2 <- median(geneExpr[group %in% select.ident.2])
				median.diff <- median.1 - median.2
				return(c(median.diff, median.1, median.2))
			})
			wil.median <- as.data.frame(t(wil.median))
			colnames(wil.median) <- c('median.diff','median.1','median.2')
		}
		# W = testRes$statistic, pVal = testRes$p.value
		test.res.dat <- as.data.frame(lapply(wil.result,function(testRes){
			c(testRes$statistic,testRes$p.value)
		}))
		test.res.dat <- as.data.frame(t(test.res.dat))
		colnames(test.res.dat) <- c('W','p.val')

		test.res.dat <- cbind(wil.median,test.res.dat)
	}
	
	test.res.dat$p.val.adj <- p.adjust(test.res.dat$p.val, method='BH')
	test.res.dat$description <- rownames(acti.ident.score)

	if (only.posi){
		if (method=='t.test'){
			test.res.dat <- subset(test.res.dat, mean.diff > 0)
		}else{
			test.res.dat <- subset(test.res.dat, median.diff > 0)
		}
	}

	if (only.sig){
		test.res.dat <- subset(test.res.dat, p.val.adj < 0.05)
	}
	return(test.res.dat)
}


# #' To conduct defferential expression test in select.ident between  different conditions in one CommPath objects
# #' @param object The first CommPath object
# #' @param select.ident Identity class of interest
# #' @param compare compare
# #' @param method Method used for differential expression test, either 'wilcox.test' or 't.test'
# #' @param p.adjust Method used for p value correction for multiple differential expression test; see p.adjust function for more information
# #' @param only.posi only logFC > 0
# #' @param only.sig only p_val_adj < 0.05
# #' @return data.frame of differentially expressed geens between the same clusters in two CommPath object
# compCommPathMarker <- function(object, select.ident, compare, method='wilcox.test', p.adjust='BH', only.posi=FALSE, only.sig=TRUE){
# 	select.var <- compare[1]
# 	select.c1 <- compare[2]
# 	select.c2 <- compare[3]
# 	cell.info <- object@cell.info

# 	select.cell <- rownames(cell.info)[(cell.info[, 'Cluster'] == select.ident) & (cell.info[,select.var] %in% c(select.c1, select.c2))]


# 	expr.mat <- as.matrix(object@data[, select.cell])

# 	lr.pair.dat <- CommPathData$DataLR[[object@meta.info$species]]
# 	all.lig.reps <- unique(c(lr.pair.dat$L, lr.pair.dat$R))
# 	expr.mat <- expr.mat[which(rownames(expr.mat) %in% all.lig.reps), ]
# 	if (nrow(expr.mat)==0){
# 		stop(paste0("There is no ligand or receptor detected in the expression matrix of cluster ",select.ident))
# 	}

# 	group <- cell.info[select.cell, select.var]
# 	cell.ident <- which(group==select.c1)
# 	cell.other <- which(group==select.c2)
# 	if (length(cell.ident) < 3){
# 		stop(paste0('There is(are) ',length(cell.ident),' cell(s) in group ',select.c1,'\nselect other one group and try again'))
# 	}
# 	if (length(cell.other) < 3){
# 		stop(paste0('There is(are) ',length(cell.other),' cell(s) in group ',select.c2,'\nselect other one group and try again'))
# 	}

# 	if (method!='wilcox.test' & method!='t.test'){
# 		stop("Select t.test or wilcox.test to conduct differential analysis")
# 	}
# 	if(method=='wilcox.test'){
# 		test.res <- apply(expr.mat, 1, function(row.expr){
# 			logFC <- log(mean(expm1(row.expr[cell.ident])) +1, base=2) - log(mean(expm1(row.expr[cell.other])) +1, base=2)
# 			p.value <- wilcox.test(x=row.expr[cell.ident], y=row.expr[cell.other])$p.value
# 			pct.1 <- length(which(row.expr[cell.ident] > 0))/length(cell.ident)
# 			pct.2 <- length(which(row.expr[cell.other] > 0))/length(cell.other)
# 			c(logFC, p.value, pct.1, pct.2)
# 		})
# 	}else{
# 		test.res <- apply(expr.mat, 1, function(row.expr){
# 			logFC <- log(mean(expm1(row.expr[cell.ident])) +1, base=2) - log(mean(expm1(row.expr[cell.other])) +1, base=2)
# 			p.value <- t.test(x=row.expr[cell.ident], y=row.expr[cell.other])$p.value

# 			pct.1 <- length(which(row.expr[cell.ident] > 0))/length(cell.ident)
# 			pct.2 <- length(which(row.expr[cell.other] > 0))/length(cell.other)
# 			c(logFC, p.value, pct.1, pct.2)
# 	})
# 	}

# 	test.res <- as.data.frame(t(test.res))
# 	colnames(test.res) <- c('avg_log2FC','p_val', 'pct.1', 'pct.2')

# 	test.res <- subset(test.res, pct.1>0 | pct.2>0)
# 	test.res$p_val_adj <- p.adjust(test.res$p_val, method=p.adjust)
# 	test.res$gene <- rownames(test.res)

# 	if (only.posi){
# 		test.res <- subset(test.res, avg_log2FC > 0)
# 	}
# 	if (only.sig){
# 		test.res <- subset(test.res, p_val_adj < 0.05)
# 	}
	
# 	return(test.res)
# }

# #' To conduct defferential activation test in select.ident between different conditions in one CommPath objects
# #' @param object The first CommPath object
# #' @param select.ident Identity class of interest
# #' @param compare compare
# #' @param method Method used for differential expression test, either 'wilcox.test' or 't.test'
# #' @param p.adjust Method used for p value correction for multiple differential expression test; see p.adjust function for more information
# #' @param only.posi only logFC > 0
# #' @param only.sig only p_val_adj < 0.05
# #' @return Data frame of differentially expressed geens between the same clusters in two CommPath object
# compCommPathPath <- function(object, select.ident, compare, method='wilcox.test', p.adjust='BH', only.posi=FALSE, only.sig=TRUE){
# 	select.var <- compare[1]
# 	select.c1 <- compare[2]
# 	select.c2 <- compare[3]
# 	cell.info <- object@cell.info

# 	select.cell <- rownames(cell.info)[(cell.info[, 'Cluster'] == select.ident) & (cell.info[,select.var] %in% c(select.c1, select.c2))]

# 	expr.mat <- object@pathway$acti.score[, select.cell]
# 	group <- cell.info[select.cell, select.var]
# 	cell.ident <- which(group==select.c1)
# 	cell.other <- which(group==select.c2)
# 	if (length(cell.ident) < 3){
# 		stop(paste0('There is(are) ',length(cell.ident),' cell(s) in group ',select.c1,'\nselect other one group and try again'))
# 	}
# 	if (length(cell.other) < 3){
# 		stop(paste0('There is(are) ',length(cell.other),' cell(s) in group ',select.c2,'\nselect other one group and try again'))
# 	}

# 	if (method!='wilcox.test' & method!='t.test'){
# 		stop("Select t.test or wilcox.test to conduct differential analysis")
# 	}
# 	if(method=='wilcox.test'){
# 		wil.result <- apply(expr.mat,1,function(geneExpr){
# 			wilcox.test(x=geneExpr[cell.ident],y=geneExpr[cell.other])
# 		})
# 		test.res.dat <- as.data.frame(lapply(wil.result,function(testRes){
# 			c(testRes$statistic,testRes$p.value)
# 		}))
# 		test.res.dat <- as.data.frame(t(test.res.dat))
# 		colnames(test.res.dat) <- c('W','p.val')

# 		wil.median <- apply(expr.mat, 1, function(geneExpr){
# 			median.1 <- median(geneExpr[cell.ident])
# 			median.2 <- median(geneExpr[cell.other])
# 			median.diff <- median.1 - median.2
# 			return(c(median.diff, median.1, median.2))
# 		})
# 		wil.median <- as.data.frame(t(wil.median))
# 		colnames(wil.median) <- c('median.diff','median.1','median.2')
		
# 		test.res.dat <- cbind(wil.median,test.res.dat)
# 	}else{
		
# 		t.result <- apply(expr.mat,1,function(geneExpr){
# 			t.test(x=geneExpr[cell.ident],y=geneExpr[cell.other])
# 		})
# 		test.res.dat <- as.data.frame(lapply(t.result,function(testRes){
# 			return(c(testRes$estimate[1]-testRes$estimate[2],testRes$estimate[1],testRes$estimate[2],testRes$statistic,testRes$parameter,testRes$p.value))
# 		}))
# 		test.res.dat <- as.data.frame(t(test.res.dat))
# 		colnames(test.res.dat) <- c('mean.diff','mean.1','mean.2','t','df','p.val')
# 	}

# 	test.res.dat$p.val.adj <- p.adjust(test.res.dat$p.val, method=p.adjust)
# 	test.res.dat$description <- rownames(expr.mat)

# 	if (only.posi){
# 		if (method=='t.test'){
# 			test.res.dat <- subset(test.res.dat, mean.diff > 0)
# 		}else{
# 			test.res.dat <- subset(test.res.dat, median.diff > 0)
# 		}
# 	}

# 	if (only.sig){
# 		test.res.dat <- subset(test.res.dat, p.val.adj < 0.05)
# 	}
# 	return(test.res.dat)
# }


#' To conduct defferential expression test in select.ident between two CommPath objects
#' @param object.1 CommPath object 1
#' @param object.2 CommPath object 2 for comparison
#' @param select.ident Identity class of interest
#' @param method Method used for differential expression test, either 'wilcox.test' or 't.test'
#' @param p.adjust Method used for p value correction for multiple differential expression test; see p.adjust function for more information
#' @param only.posi only logFC > 0
#' @param only.sig only p_val_adj < 0.05
#' @return Data frame of differentially expressed geens between the same clusters in two CommPath object
#' @export
diffCommPathMarker <- function(object.1, object.2, select.ident, method='wilcox.test', p.adjust='BH', only.posi=FALSE, only.sig=TRUE){
	lr.pair.dat <- CommPathData$DataLR[[object.1@meta.info$species]]
	all.lig.reps <- unique(c(lr.pair.dat$L, lr.pair.dat$R))

	obj.1 <- subsetCommPath(object.1, ident.keep=select.ident)
	expr.mat.1 <- as.matrix(obj.1@data)
	expr.mat.1 <- expr.mat.1[which(rownames(expr.mat.1) %in% all.lig.reps), ]
	
	obj.2 <- subsetCommPath(object.2, ident.keep=select.ident)
	expr.mat.2 <- as.matrix(obj.2@data)
	expr.mat.2 <- expr.mat.2[which(rownames(expr.mat.2) %in% all.lig.reps), ]
	
	if (ncol(expr.mat.1) < 3){
		stop(paste0('There is(are) ',ncol(expr.mat.1),' cell(s) in object.1\nselect other one ident and try again'))
	}
	if (ncol(expr.mat.2) < 3){
		stop(paste0('There is(are) ',ncol(expr.mat.2),' cell(s) in object.2\nselect other one ident and try again'))
	}
	gene.keep <- intersect(rownames(expr.mat.1), rownames(expr.mat.2))
	expr.mat.1 <- expr.mat.1[gene.keep, ]
	cell.ident <- colnames(expr.mat.1)
	expr.mat.2 <- expr.mat.2[gene.keep, ]
	cell.other <- colnames(expr.mat.2)

	if (method!='wilcox.test' & method!='t.test'){
		stop("Select t.test or wilcox.test to conduct differential analysis")
	}
	if(method=='wilcox.test'){
		test.res <- apply(cbind(expr.mat.1, expr.mat.2), 1, function(row.expr){
			logFC <- log(mean(expm1(row.expr[cell.ident])) +1, base=2) - log(mean(expm1(row.expr[cell.other])) +1, base=2)
			p.value <- wilcox.test(x=row.expr[cell.ident], y=row.expr[cell.other])$p.value
			pct.1 <- length(which(row.expr[cell.ident] > 0))/length(cell.ident)
			pct.2 <- length(which(row.expr[cell.other] > 0))/length(cell.other)
			c(logFC, p.value, pct.1, pct.2)
		})
	}else{
		test.res <- apply(cbind(expr.mat.1, expr.mat.2), 1, function(row.expr){
			logFC <- log(mean(expm1(row.expr[cell.ident])) +1, base=2) - log(mean(expm1(row.expr[cell.other])) +1, base=2)
			p.value <- t.test(x=row.expr[cell.ident], y=row.expr[cell.other])$p.value

			pct.1 <- length(which(row.expr[cell.ident] > 0))/length(cell.ident)
			pct.2 <- length(which(row.expr[cell.other] > 0))/length(cell.other)
			c(logFC, p.value, pct.1, pct.2)
	})
	}

	test.res <- as.data.frame(t(test.res))
	colnames(test.res) <- c('avg_log2FC','p_val', 'pct.1', 'pct.2')
	test.res <- subset(test.res, pct.1>0 | pct.2>0)
	test.res$p_val_adj <- p.adjust(test.res$p_val, method=p.adjust)
	test.res$gene <- rownames(test.res)
	if (only.posi){
		test.res <- subset(test.res, avg_log2FC > 0)
	}
	if (only.sig){
		test.res <- subset(test.res, p_val_adj < 0.05)
	}
	return(test.res)
}

#' To conduct defferential activation test in select.ident between two CommPath objects
#' @param object.1 CommPath object 1
#' @param object.2 CommPath object 2 for comparison
#' @param select.ident Identity class of interest
#' @param method Method used for differential expression test, either 'wilcox.test' or 't.test'
#' @param p.adjust Method used for p value correction for multiple differential expression test; see p.adjust function for more information
#' @param min.size Minimum size of overlaping genes between candidate pathways and the expression matrix
#' @param only.posi only logFC > 0
#' @param only.sig only p_val_adj < 0.05
#' @param ... Extra parameters passed to gsva
#' @return Data frame of differentially activated pathways between the same clusters in two CommPath object
#' @export
diffCommPathPath <- function(object.1, object.2, select.ident, method='wilcox.test', p.adjust='BH', min.size=10, only.posi=FALSE, only.sig=TRUE, ...){
	obj.1 <- subsetCommPath(object.1, ident.keep=select.ident)
	obj.2 <- subsetCommPath(object.2, ident.keep=select.ident)
	score.mat.1 <- obj.1@pathway$acti.score
	# find pathways presenting in obj.1
	uniq.path.name <- rownames(score.mat.1)

	# add gsva analysis in obj.2 for the uniq.path.name pathways
	uniq.path.set <- obj.1@pathway$pathwayLR[uniq.path.name]
	obj.2.expr <- as.matrix(obj.2@data)
	if (obj.1@pathway$method=='gsva'){
		score.mat.2 <- GSVA::gsva(obj.2.expr, uniq.path.set, min.sz=min.size, ...)
	}else{
		score.mat.2 <- t(as.data.frame(lapply(uniq.path.set, function(eachPath){
			overlap.gene <- intersect(eachPath, rownames(obj.2.expr))
			if (length(overlap.gene) < min.size){
				return(rep(NA, ncol(obj.2.expr)))
			}else{
				return(colMeans(obj.2.expr[overlap.gene, ]))
			}
		})))
		rownames(score.mat.2) <- names(uniq.path.set)
		score.mat.2 <- score.mat.2[which(!is.na(score.mat.2[,1])), ]
	}

	path.keep <- intersect(rownames(score.mat.1), rownames(score.mat.2))
	score.mat.1 <- score.mat.1[path.keep, ]
	score.mat.2 <- score.mat.2[path.keep, ]

	### differential analysis
	if (ncol(score.mat.1) < 3){
		stop(paste0('There is(are) ',ncol(score.mat.1),' cell(s) in object.1\nselect other one ident and try again'))
	}
	if (ncol(score.mat.2) < 3){
		stop(paste0('There is(are) ',ncol(score.mat.1),' cell(s) in object.2\nselect other one ident and try again'))
	}

	expr.mat <- cbind(score.mat.1,score.mat.2)
	cell.ident <- c(1:ncol(score.mat.1))
	cell.other <- c(1:ncol(score.mat.2)) + ncol(score.mat.1)


	if (method!='wilcox.test' & method!='t.test'){
		stop("Select t.test or wilcox.test to conduct differential analysis")
	}
	if(method=='wilcox.test'){
		wil.result <- apply(expr.mat,1,function(geneExpr){
			wilcox.test(x=geneExpr[cell.ident],y=geneExpr[cell.other])
		})
		test.res.dat <- as.data.frame(lapply(wil.result,function(testRes){
			c(testRes$statistic,testRes$p.value)
		}))
		test.res.dat <- as.data.frame(t(test.res.dat))
		colnames(test.res.dat) <- c('W','p.val')

		wil.median <- apply(expr.mat, 1, function(geneExpr){
			median.1 <- median(geneExpr[cell.ident])
			median.2 <- median(geneExpr[cell.other])
			median.diff <- median.1 - median.2
			return(c(median.diff, median.1, median.2))
		})
		wil.median <- as.data.frame(t(wil.median))
		colnames(wil.median) <- c('median.diff','median.1','median.2')
		
		test.res.dat <- cbind(wil.median,test.res.dat)
	}else{
		t.result <- apply(expr.mat,1,function(geneExpr){
			t.test(x=geneExpr[cell.ident],y=geneExpr[cell.other])
		})
		test.res.dat <- as.data.frame(lapply(t.result,function(testRes){
			return(c(testRes$estimate[1]-testRes$estimate[2],testRes$estimate[1],testRes$estimate[2],testRes$statistic,testRes$parameter,testRes$p.value))
		}))
		test.res.dat <- as.data.frame(t(test.res.dat))
		colnames(test.res.dat) <- c('mean.diff','mean.1','mean.2','t','df','p.val')
	}

	test.res.dat$p.val.adj <- p.adjust(test.res.dat$p.val, method=p.adjust)
	test.res.dat$description <- rownames(expr.mat)

	ident.path.dat <- diffPath(object=object.1, select.ident.1=select.ident, only.posi=FALSE, only.sig=FALSE)
	test.res.dat <- subset(test.res.dat, description %in% ident.path.dat$description)
	ident.path.dat <- ident.path.dat[match(test.res.dat$description, ident.path.dat$description), c('ligand.in.path', 'receptor.in.path')]
	test.res.dat <- cbind(test.res.dat, ident.path.dat)

	if (only.posi){
		if (method=='t.test'){
			test.res.dat <- subset(test.res.dat, mean.diff > 0)
		}else{
			test.res.dat <- subset(test.res.dat, median.diff > 0)
		}
	}

	if (only.sig){
		test.res.dat <- subset(test.res.dat, p.val.adj < 0.05)
	}
	return(test.res.dat)

}
