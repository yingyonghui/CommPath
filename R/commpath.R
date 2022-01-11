data('CommpathData', package=('Commpath'))

#' To identify marker ligands and marker receptors in the expression matrix
#' @param object Commpath object
#' @param method Method used for differential expression test, either 'wilcox.test' or 't.test'
#' @param p.adjust Method used for p value correction for multiple differential expression test; see p.adjust function for more information
#' @return Data frame containing the differential expression test
#' @export
findLRmarker <- function(object, method='wilcox.test', p.adjust='BH'){

	species <- object@meta.info$species
	expr.mat <- as.matrix(object@data)
	label <- object@meta.info$cell.info[,'Cluster']
	
	lr.pair.dat <- CommpathData$DataLR[[species]]
	all.lig.reps <- unique(c(lr.pair.dat$L, lr.pair.dat$R))

	expr.mat <- expr.mat[which(rownames(expr.mat) %in% all.lig.reps), ]
	if (nrow(expr.mat)==0){
		stop("there is no ligand or receptor detected in the expression matrix!")
	}

	if (!is.factor(label)){ label <- factor(label) }
	ident.level <- levels(label)
	
	label <- as.character(label)
	### p.value calculated
	if (method!='wilcox.test' & method!='t.test'){
		stop("select t.test or wilcox.test to conduct differential analysis")
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
#' @param object Commpath object
#' @param logFC.thre logFC threshold, marker genes with a logFC > logFC.thre will be considered
#' @param p.thre p threshold, marker genes with a adjust p value < p.thre will be considered
#' @return List containing the ligand-receptor interaction information
#' @export
findLRpairs <- function(object, logFC.thre=0, p.thre=0.05){
	options(stringsAsFactors=F)
	marker.dat <- object@LR.marker
	species <- object@meta.info$species

	lr.pair.dat <- CommpathData$DataLR[[species]]
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
	rownames(Interact.num.mat) <- cluster.level
	colnames(Interact.num.mat) <- cluster.level
	Interact.num.dat <- melt(Interact.num.mat, varnames=c('Cell.From','Cell.To'),value.name="LR.Number",  na.rm=TRUE)
	Interact.num.dat <- factor.to.character(Interact.num.dat)

	rownames(Interact.gene.mat) <- cluster.level
	colnames(Interact.gene.mat) <- cluster.level
	Interact.gene.dat <- melt(Interact.gene.mat,varnames=c('Cell.From','Cell.To'),value.name="LR.Info", na.rm = TRUE)
	Interact.gene.dat <- factor.to.character(Interact.gene.dat)
	
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

	Interact <- list(InteractNumer=Interact.num.dat,InteractGene=Interact.gene.dat, InteractGeneUnfold=lr.unfold.dat, markerL=marker.lig.dat, markerR=marker.rep.dat)

	object@interact <- Interact

	object@meta.info$logFC.thre <- logFC.thre
	object@meta.info$p.thre <- p.thre
	return(object)
}


#' To find those pathways in which the genesets show overlap with the marker ligand and receptor genes in our dataset
#' @param object Commpath object
#' @param category Character to indicate which pathway to investigate; one of "go" (GO terms), "kegg" (for KEGG pathways), 'wiki' (for WikiPathways), and "reactome" (for reactome pathways), or "all" for all pathways
#' @return Interact list containing the ligand-receptor interaction information and the pathways showing overlap with the marker ligand and receptor genes in the dataset
#' @export
findLRpath <- function(object, category='all'){
	options(stringsAsFactors=F)
	if (category!='all' & category!='go' & category!='kegg' & category!='wiki' & category!='reactome'){
		stop("wrong category selected. Select one of 'go', 'kegg', 'wiki', and 'reactome', or 'all' for all pathways")
	}
	species <- object@meta.info$species

	### select one category
	if (category=='all'){
		path.go.list <- CommpathData$DataGOterm[[species]]
		path.kegg.list <- CommpathData$DataKEGGpathway[[species]]
		path.wiki.list <- CommpathData$DataWikiPathway[[species]]
		path.reac.list <- CommpathData$DataReactome[[species]]
		path.list <- c(path.go.list, path.kegg.list, path.wiki.list, path.reac.list)
		path.list <- path.list[which(!duplicated(names(path.list)))]

	}else if (category=='go'){
		path.list <- CommpathData$DataGOterm[[species]]
	}else if(category=='kegg'){
		path.list <- CommpathData$DataKEGGpathway[[species]]
	}else if(category=='wiki'){
		path.list <- CommpathData$DataWikiPathway[[species]]
	}else if(category=='reactome'){
		path.list <- CommpathData$DataReactome[[species]]
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
#' @param object Commpath object
#' @param method Method used for scoring the pathways, either 'gsva' of 'average'
#' @param min.size Minimum size of overlaping genes between candidate pathways and the expression matrix
#' @param ... Other parameters passed to gsva
#' @return Commpath object with pathways activation scores stored in the slot pathway
#' @export
scorePath <- function(object, method='gsva', min.size=10, ...){
	sample.expr <- as.matrix(object@data)
	path.list <- object@pathway$pathwayLR
	
	if (method=='gsva'){
		gsva.mat <- GSVA::gsva(sample.expr, path.list, min.sz=min.size, ...)
	}
	object@pathway$acti.score <- gsva.mat

	return(object)
}


#' To find different enriched pathway between two group cells
#' @param object Commpath object
#' @param select.ident.1 Identity class to define cells for group 1
#' @param select.ident.2 Identity class to define cells for group 2 for comparison; if 'NULL', use all other cells for comparison
#' @param method Method used for differential enrichment analysis, either 't.test' of 'wilcox.test'
#' @return Dataframe including the statistic result comparing the pathway enrichment sorces between group 1 and group 2, the significant recetor and ligand of group 1 in the pathways, and the corresponding up stream identity class which interact with group 1 by releasing specific ligand
#' @export
diffPath <- function(object, select.ident.1, select.ident.2=NULL, method='t.test'){
	options(stringsAsFactors=F)
	gsva.mat <- object@pathway$acti.score
	ident.label <- object@meta.info$cell.info$Cluster
	if (method!='t.test' & method!='wilcox.test'){
		stop("select t.test or wilcox.test to conduct differential analysis")
	}else if(method=='t.test'){
		test.res.dat <- data.frame(matrix(NA,0,12))
		colnames(test.res.dat) <- c('mean.diff','mean.1','mean.2','t','df','p.val','p.val.adj','description','cell.up','ligand.up','receptor.in.path','ligand.in.path')
	}else{
		test.res.dat <- data.frame(matrix(NA,0,11))
		colnames(test.res.dat) <- c('median.diff','median.1','median.2','W','p.val','p.val.adj','description','cell.up','ligand.up','receptor.in.path','ligand.in.path')
	}

	path.lr.list <- object@pathway$pathwayLR
	if (is.null(path.lr.list)){
		stop("no pathway detected, run findLRpath befor diffPath")
	}
	# marker.lig.dat <- Interact$markerL
	# marker.lig.dat <- marker.lig.dat[marker.lig.dat$cluster==select.ident,]
	marker.ident1.rep.dat <- subset(object@interact$InteractGeneUnfold, Cell.To %in% select.ident.1)

	if (nrow(marker.ident1.rep.dat)==0){
		warning(paste0("there is no significant receptor for cluster ", select.ident.1))
		return(test.res.dat)
	}

	### find those pathways in which genesets hava overlap with the marker receptors of the selected cluster
	ident.rep.gene <- marker.ident1.rep.dat$Receptor
	which.overlap.list <-lapply(path.lr.list, function(x){
		overlap.idx <-  x %in% ident.rep.gene
			if(any(overlap.idx)){ 
				return(paste(x[overlap.idx],collapse=',')) 
			}else{
				return(FALSE)
		}
	})
	if (all(which.overlap.list=='FALSE')){
		warning(paste0("there is no pathway showing overlap with the receptors in cluster ",select.ident.1))
		return(test.res.dat)
	}
	overlap.rep.list <- which.overlap.list[which(which.overlap.list!='FALSE')]
	### since we set min.sz in the gsva function, there are some pathways not calculated in the gsva process, we shall remove those pathways
	overlap.rep.list <- overlap.rep.list[names(overlap.rep.list) %in% rownames(gsva.mat)]

	### t.test or wilcox.text for the pathways for the selected cluster
	gsva.ident.mat <- gsva.mat[names(overlap.rep.list),]
	group <- as.character(ident.label)
	test.res.dat <- pathTest(gsva.ident.mat, group, select.ident.1, select.ident.2, method)

	
	### to find the upstream ident and ligand the ident.1 recieved
	test.res.dat$cell.up <- NA
	test.res.dat$ligand.up <- NA
	test.res.dat$receptor.in.path <- unlist(overlap.rep.list)
	for (each.row in 1:nrow(test.res.dat)){
		each.rep <- test.res.dat[each.row,'receptor.in.path']
		rep.vec <- strsplit(each.rep, split=',')[[1]]
		
		up.cell.inte <- c()
		up.lig.inte <- c()
		cur.rep.inte <- c()
		for (each.rep in rep.vec){
			each.unfold.rep <- subset(marker.ident1.rep.dat, Receptor==each.rep)
			up.cell.inte <- c(up.cell.inte,paste(each.unfold.rep$Cell.From, collapse=';'))
			up.lig.inte <- c(up.lig.inte,paste(each.unfold.rep$Ligand, collapse=';'))		
			cur.rep.inte <- c(cur.rep.inte,paste(each.unfold.rep$Receptor, collapse=';'))		
		}
		test.res.dat[each.row,'cell.up'] <- paste(up.cell.inte, collapse=';')
		test.res.dat[each.row,'ligand.up'] <- paste(up.lig.inte, collapse=';')
		test.res.dat[each.row,'receptor.in.path'] <- paste(cur.rep.inte, collapse=';')
	}

	### to find the ligand in the same pathway
	marker.lig.dat <- object@interact$markerL
	ident.lig.vec <- marker.lig.dat[marker.lig.dat$cluster %in% select.ident.1,'gene']
	### for each pathway in the DEG result, find which pathway show overlap with marker ligands of the selected ident
	ident.lig.in.path <- sapply(test.res.dat$description, function(eachPath){
		each.set <- path.lr.list[[eachPath]]
		if (any(ident.lig.vec %in% each.set)){
			return(paste(ident.lig.vec[ident.lig.vec %in% each.set],collapse=';'))
		}else{
			return(NA)
		}
	})
	test.res.dat$ligand.in.path <- unlist(ident.lig.in.path)

	return(test.res.dat)
}


#' To find different enriched pathways in each identity class 
#' @param object Commpath object
#' @param method Method used for differential enrichment analysis, either 't.test' of 'wilcox.test'
#' @return Dataframe including the statistic result comparing the pathway enrichment sorces between cells in each cluster and all other clusters, the significant recetor and ligand in the pathways, and the corresponding up stream identity class and ligand
#' @export
diffAllPath <- function(object, method='t.test'){
	if (method=='t.test'){
		all.test.dat <- data.frame(matrix(NA,0,12))
	}else{
		all.test.dat <- data.frame(matrix(NA,0,11))
	}
	all.rep.ident <- unique(object@interact$InteractGene$Cell.To)
	for (each.ident in all.rep.ident){
		message(paste0('Identifying pathways for cluster ',each.ident,'...'))
		test.res.dat <- diffPath(object, select.ident.1=each.ident, select.ident.2=NULL, method=method)
		if (nrow(test.res.dat)==0){ next }
		test.res.dat$cluster <- each.ident
		all.test.dat <- rbind(all.test.dat, test.res.dat)
	}

	return(all.test.dat)
}


#' To find the downstream identity class of specific ligand released by specific upstream identity class
#' @param object Commpath object
#' @param select.ident Upstream identity class; if 'NULL', use all identity classes
#' @param select.ligand Ligand released by upstream identity class; if 'NULL', use all ligands that are markers for the selected upstream identity class
#' @return Dataframe including the interaction information
#' @export
findReceptor <- function(object, select.ident=NULL, select.ligand=NULL){
	options(stringsAsFactors=F)

	if (is.null(select.ident) & is.null(select.ligand)){
		stop("either a select.ident or a select.ligand need to be asigned")
	}
	InteractGeneUnfold <- object@interact$InteractGeneUnfold
	if (is.null(select.ligand)){
		ident.down.dat <- subset(InteractGeneUnfold, Cell.From==select.ident)
	}else if(is.null(select.ident)){
		ident.down.dat <- subset(InteractGeneUnfold, Ligand==select.ligand)
	}else{
		ident.down.dat <- subset(InteractGeneUnfold, Cell.From==select.ident & Ligand==select.ligand)
	}

	if (nrow(ident.down.dat)==0){
		warning("no downstream ident found for the selected ident and ligand")
	}
	return(ident.down.dat)
}


#' To find the upstream identity classes and ligands of specific receptor expressed by specific downstream identity class
#' @param object Commpath object
#' @param select.ident Downstream identity class; if 'NULL', use all identity classes
#' @param select.receptor Receptor expressed by downstream identity class; if 'NULL', use all receptors that are markers for the selected downstream identity class
#' @return Dataframe including the interaction information
#' @export
findLigand <- function(object, select.ident=NULL, select.receptor=NULL){
	options(stringsAsFactors=F)

	if (is.null(select.ident) & is.null(select.receptor)){
		stop("either a select.ident or a select.receptor need to be asigned")
	}

	InteractGeneUnfold <- object@interact$InteractGeneUnfold
	if (is.null(select.receptor)){
		ident.up.dat <- subset(InteractGeneUnfold, Cell.To==select.ident)
	}else if(is.null(select.ident)){
		ident.up.dat <- subset(InteractGeneUnfold, Receptor==select.receptor)
	}else{
		ident.up.dat <- subset(InteractGeneUnfold, Cell.To==select.ident & Receptor==select.receptor)
	}

	if (nrow(ident.up.dat)==0){
		warning("no upstream ident found for the selected ident and receptor")
	}
	return(ident.up.dat)
}


#' Differential enrichment analysis by t.test or wilcox.txt
#' @param gsva.ident.mat Matrix of pathway scores, pathway * cell
#' @param group Vector of group labels of cells
#' @param select.ident.1 Identity class 1
#' @param select.ident.2 Identity class 2 for comparison
#' @param method Method for hypothesis test, either 't.test' or 'wilcox.test'
#' @return Dataframe including the statistic result
#' @export
pathTest <- function(gsva.ident.mat, group, select.ident.1, select.ident.2=NULL, method='t.test'){
	if(method=='t.test'){
		if (is.null(select.ident.2)){
			t.result <- apply(gsva.ident.mat,1,function(geneExpr){
				t.test(x=geneExpr[group %in% select.ident.1],y=geneExpr[!(group %in% select.ident.1)])})
		}else{
			t.result <- apply(gsva.ident.mat,1,function(geneExpr){
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
			wil.result <- apply(gsva.ident.mat,1,function(geneExpr){
				wilcox.test(x=geneExpr[group %in% select.ident.1],y=geneExpr[!(group %in% select.ident.1)])})

			wil.median <- apply(gsva.ident.mat, 1, function(geneExpr){
				median.1 <- median(geneExpr[group %in% select.ident.1])
				median.2 <- median(geneExpr[!(group %in% select.ident.1)])
				median.diff <- median.1 - median.2
				return(c(median.diff, median.1, median.2))
			})
			wil.median <- as.data.frame(t(wil.median))
			colnames(wil.median) <- c('median.diff','median.1','median.2')

		}else{
			wil.result <- apply(gsva.ident.mat,1,function(geneExpr){
				wilcox.test(x=geneExpr[group %in% select.ident.1],y=geneExpr[group %in% select.ident.2])})
			wil.median <- apply(gsva.ident.mat, 1, function(geneExpr){
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
	test.res.dat$description <- rownames(gsva.ident.mat)
	return(test.res.dat)
}


#' To conduct defferential expression test in select.ident between two Commpath objects
#' @param object.1 The first Commpath object
#' @param object.2 The second Commpath object for comparison
#' @param select.ident Identity class of interest
#' @param method Method used for differential expression test, either 'wilcox.test' or 't.test'
#' @param p.adjust Method used for p value correction for multiple differential expression test; see p.adjust function for more information
#' @param only.posi only logFC > 0
#' @param only.sig only p_val_adj < 0.05
#' @return data.frame of differentially expressed geens between the same clusters in two Commpath object
#' @export
diffCommpathMarker <- function(object.1, object.2, select.ident, method='wilcox.test', p.adjust='BH', only.posi=FALSE, only.sig=TRUE){
	obj.1 <- subsetCommpath(object.1, ident.keep=select.ident)
	obj.2 <- subsetCommpath(object.2, ident.keep=select.ident)
	expr.mat.1 <- obj.1@data
	expr.mat.2 <- obj.2@data
	if (ncol(expr.mat.1) < 3){
		stop(paste0('there is(are) ',ncol(expr.mat.1),' cell(s) in object.1\nselect other one ident and try again'))
	}
	if (ncol(expr.mat.2) < 3){
		stop(paste0('there is(are) ',ncol(expr.mat.2),' cell(s) in object.2\nselect other one ident and try again'))
	}
	gene.keep <- intersect(rownames(expr.mat.1), rownames(expr.mat.2))
	expr.mat.1 <- expr.mat.1[gene.keep, ]
	cell.ident <- colnames(expr.mat.1)
	expr.mat.2 <- expr.mat.2[gene.keep, ]
	cell.other <- colnames(expr.mat.2)

	if (method!='wilcox.test' & method!='t.test'){
		stop("select t.test or wilcox.test to conduct differential analysis")
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
	if (only.posi){
		test.res <- subset(test.res, avg_log2FC > 0)
	}
	if (only.sig){
		test.res <- subset(test.res, p_val_adj < 0.05)
	}
	return(test.res)
}


#' To compare two Commpath objects
#' @param object.1 The first Commpath object
#' @param object.2 The second Commpath object for comparison
#' @param select.ident Identity class of interest
#' @param diff.obj.marker diff.obj.marker computed by diffCommpathMarker
#' @return data.frame of differentially expressed geens between the same clusters in two Commpath object
#' @export
diffCommpathPlot <- function(object.1, object.2, select.ident, diff.obj.marker=NULL, top.n.receptor=5, order=NULL, top.n.path=10, p.thre=0.05, dot.ident.col=NULL, dot.receptor.col=NULL, bar.pathway.col=NULL, bar.pathway.width=10, dot.ident.size=1, dot.receptor.size=1, label.text.size=1, label.title.size=1, line.ident.width=1, line.path.width=1){
	diff.obj.marker <- diffCommpathMarker(object.1, object.2, select.ident, method='wilcox.test', p.adjust='BH', only.posi=FALSE, only.sig=TRUE)

	### check the input and select significant pathways
	if (is.null(ident.path.dat)){
		ident.path.dat <- diffPath(object, select.ident.1=select.ident.1)
	}
	if ('t' %in% colnames(ident.path.dat)){
		ident.path.dat <- subset(ident.path.dat, p.val.adj < p.thre & t > 0)
	}else if ('W' %in% colnames(ident.path.dat)){
		ident.path.dat <- subset(ident.path.dat, p.val.adj < p.thre & median.diff > 0)
	}else{
		stop('please input the integrate ident.path.dat computed from diffPath')
	}
	ident.path.dat <- ident.path.dat[order(ident.path.dat$p.val.adj, decreasing=TRUE),]
	all.sig.path <- ident.path.dat$description

	### preprocess and extract useful information
	plot.dat <- extract.info(ident.path.dat)
	# subset the plot.dat to exclude those markers not differentially expressed between two objects
	markerR.dat <- subset(object.1@interact$markerR, subset=(cluster==select.ident & gene %in% plot.dat$cur.rep))
	if (nrow(markerR.dat)==0){
		stop('there is no marker receptor for the selected ident compared to other clusters')
	}
	plot.dat <- subset(plot.dat, cur.rep %in% rownames(diff.obj.marker))
	markerR.dat <- subset(object.1@interact$markerR, subset=(cluster==select.ident & gene %in% plot.dat$cur.rep))
	if (nrow(markerR.dat)==0){
		stop('there is no marker receptor for the selected ident between two objects')
	}

	markerR.dat <- markerR.dat[order(abs(markerR.dat$avg_log2FC), decreasing=TRUE), ]
	if (top.n.receptor > nrow(markerR.dat)){
		warning(paste0('there is(are) ', nrow(markerR.dat),' marker receptor(s) in the selected ident between two objects, and the input top.n.receptor is ', top.n.receptor))
		top.n.receptor <- nrow(markerR.dat)
	}
	top.receptor.dat <- markerR.dat[1:top.n.receptor,]
	plot.dat <- subset(plot.dat, cur.rep %in% top.receptor.dat$gene)
	

	up.ident <- plot.dat$up.ident
	cur.rep <- plot.dat$cur.rep
	return(diff.obj.marker)
}

