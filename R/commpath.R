# mycolors <- hue_pal(c=100)(25)
data('commpathData', package=('commpath'))

#' To present a circos plot
#' @param Interact Interact list returned by findLRpairs
#' @param order Vector of orders of the identities arround the circos plot
#' @param col Vector of colors of each identity; names of the col vector are supposed to be assigned to indicate each color for each identity.
#' @param ident To highlight the interaction between a specific identity class and others; if 'NULL', plot interaction for all identity classes
#' @return Circos plot showing the ligand-receptor interaction
#' @export
circosPlot <- function(Interact, order=NULL, col=NULL, ident=NULL){
	options(stringsAsFactors=F)
	Interact.num.dat <- Interact$InteractNumer
	Interact.num.dat = Interact.num.dat[Interact.num.dat$LR.Number!=0,]
	all.ident <- unique(c(Interact.num.dat$Cell.From,Interact.num.dat$Cell.To))

	### to check the order parameter
	if(!is.null(order)){
		all.ident <- orderCheck(all.ident, order)
	}
	all.ident <- all.ident[order(all.ident)]
	
	### to check the col parameter
	if (is.null(col)){ 
		col <- scales::hue_pal(c=100)(length(all.ident))
		names(col) <- all.ident
	}else{
		if (is.null(names(col))){
			stop('the cols should be named according to the identity class')
		}else{
			colo.name <- names(col)
			ident.missed <- all.ident[!(all.ident %in% colo.name)]
			ident.missed <- pasteIdent(ident.missed)
			stop(paste0('the ident class ',ident.missed,' may be missed in the input color'))
		}
	}

	circos.clear()
	circos.par(start.degree=90, clock.wise=F)

	### plot interaction for all identity classes
	if (is.null(ident)){
		chordDiagram(Interact.num.dat, order=all.ident, grid.col=col, annotationTrack=c("name","grid"), transparency=0.1, directional=1, direction.type='arrows', link.arr.type = "big.arrow", preAllocateTracks = list(track.height=max(strwidth(all.ident))), annotationTrackHeight=convert_height(c(1, 2), "mm"))

	}else{
		### highlight the interaction for a specific identity class
		### check the ident parameter
		if (!(all(ident %in% all.ident)) | length(ident)>1){
			stop(paste0('select one identity class from ', pasteIdent(all.ident)))
		}

		line.col <- rep('gray',nrow(Interact.num.dat))
		line.col[(Interact.num.dat$Cell.From==ident)] <- col[as.character(ident)]
		ident.lig = Interact.num.dat[Interact.num.dat$Cell.To==ident,'Cell.From']
		line.col[(Interact.num.dat$Cell.To==ident)] <- col[as.character(ident.lig)]
		
		link.rank <- 1:nrow(Interact.num.dat)
		ident.line <- (Interact.num.dat$Cell.From==ident) | (Interact.num.dat$Cell.To==ident)
		link.rank[!ident.line] <- rank(-(Interact.num.dat[!ident.line,'LR.Number']))
		link.rank[ident.line] <- rank(-(Interact.num.dat[ident.line,'LR.Number']))+length(which(!ident.line))
		chordDiagram(Interact.num.dat, order=all.ident, grid.col=col, annotationTrack=c("name","grid"), transparency=0.1, directional=1, direction.type='arrows', link.arr.type = "big.arrow", link.rank=link.rank, col=line.col, preAllocateTracks = list(track.height=max(strwidth(all.ident))), annotationTrackHeight=convert_height(c(1, 2), "mm"))
	}
}

#' To identify marker ligands and marker receptors in the expression matrix
#' @param expr.mat Matrix or data frame of expression matrix, with genes in rows and cells in columns
#' @param label  Vector of identity classes of cells in the expression matrix
#' @param species  Species, either 'hsapiens', 'mmusculus', or 'rnorvegicus' 
#' @param method Method used for differential expression test, either 'wilcox.test' or 't.test'
#' @param p.adjust Method used for p value correction for multiple differential expression test; see p.adjust function for more information
#' @return Data frame containing the differential expression test
#' @export
findLRmarker <- function(expr.mat, label, species, method='wilcox.test', p.adjust='BH'){
	if (length(species)>1){
		stop("select one species once")
	}
	if (!(species %in% c('hsapiens', 'mmusculus', 'rnorvegicus'))){
		stop("select one species from 'hsapiens', 'mmusculus', and 'rnorvegicus'")
	}

	lr.pair.dat <- commpathData$DataLR[[species]]
	all.lig.reps <- unique(c(lr.pair.dat$L, lr.pair.dat$R))
	expr.mat <- expr.mat[which(rownames(expr.mat) %in% all.lig.reps), ]
	if (nrow(expr.mat)==0){
		stop("there is no ligand or receptor detected in the expression matrix!")
	}

	if (!is.factor(label)){ label <- as.factor(label) }
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
				c(logFC, p.value)
			})

			test.res <- as.data.frame(t(test.res))
			colnames(test.res) <- c('avg_log2FC','p_val')
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
				c(logFC, p.value)
			})

			test.res <- as.data.frame(t(test.res))
			colnames(test.res) <- c('avg_log2FC','p_val')
			test.res$p_val_adj <- p.adjust(test.res$p_val, method=p.adjust)
			test.res$cluster <- each.level
			test.res$gene <- rownames(test.res)
			return(test.res)

		})
		test.all.res <- do.call(rbind, test.all.res)
	}

	return(test.all.res)
}

#' To find marker ligands and marker receptors
#' @param marker.dat Data frame containing information of marker genes
#' @param species Species, either 'hsapiens', 'mmusculus', or 'rnorvegicus' 
#' @param logFC.thre logFC threshold, marker genes with a logFC > logFC.thre will be considered
#' @param p.thre p threshold, marker genes with a adjust p value < p.thre will be considered
#' @return List containing the ligand-receptor interaction information
#' @export
findLRpairs <- function(marker.dat, species, logFC.thre=0, p.thre=0.05){
	options(stringsAsFactors=F)

	lr.pair.dat <- commpathData$DataLR[[species]]
	ligs <- lr.pair.dat$L
	reps <- lr.pair.dat$R

	cluster.level <- as.factor(unique(marker.dat$cluster))
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

	rownames(Interact.gene.mat) <- cluster.level
	colnames(Interact.gene.mat) <- cluster.level
	Interact.gene.dat <- melt(Interact.gene.mat,varnames=c('Cell.From','Cell.To'),value.name="LR.Info", na.rm = TRUE)
	
	rownames(Interact.lig.mat) <- cluster.level
	colnames(Interact.lig.mat) <- cluster.level
	Interact.lig.dat <- melt(Interact.lig.mat,varnames=c('Cell.From','Cell.To'),value.name="L.Info", na.rm = TRUE)

	rownames(Interact.rep.mat) <- cluster.level
	colnames(Interact.rep.mat) <- cluster.level
	Interact.rep.dat <- melt(Interact.rep.mat,varnames=c('Cell.From','Cell.To'),value.name="R.Info", na.rm = TRUE)
	
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
	Interact <- list(InteractNumer=Interact.num.dat,InteractGene=Interact.gene.dat, InteractGeneUnfold=lr.unfold.dat, markerL=marker.lig.dat, markerR=marker.rep.dat, logFC.thre=logFC.thre, p.thre=p.thre, species=species)
	return(Interact)
}

#' To present a dot plot for specific ligand-receptor pairs in specific clusters
#' @param Interact Interact list returned by findLRpairs
#' @param ligand.ident Vector containing the ligand ident
#' @param receptor.ident Vector containing the receptor ident
#' @param ident.levels Vector of levels of the identities 
#' @param return.data Logical value indicating whether to return the data for the plot or not
#' @return Dotplot showing the ligand-receptor interaction between the selected ligand.ident and receptor.ident
#' @export
dotPlot <- function(Interact, ligand.ident=NULL, receptor.ident=NULL, ident.levels=NULL, return.data=FALSE){
	options(stringsAsFactors=F)
	#colorRampPalette(c("#440154" ,"#21908C", "#FDE725"))(100)
	if (is.null(ligand.ident) & is.null(receptor.ident)){
		stop("either ligand.ident or ligand.ident need to be asigned")
	}
	if (length(ligand.ident)>1 & length(receptor.ident)>1){
		stop("specify one cluster for ligand or receptor analysis")
	}

	# get the InteractGene dataframe
	inter.gene.dat <- Interact$InteractGene

	if (length(ligand.ident)==1){
		inter.gene.dat$Xaxis <- inter.gene.dat$Cell.To
		shape <- 16
	}else{
		inter.gene.dat$Xaxis <- inter.gene.dat$Cell.From
		shape <- 17
	}

	logFC.thre <- Interact$logFC.thre
	p.thre <- Interact$p.thre

	inter.ident.dat <- inter.gene.dat
	if (!is.null(receptor.ident)){
		inter.ident.dat <- inter.ident.dat[inter.ident.dat$Cell.To %in% receptor.ident,]
	}
	if (!is.null(ligand.ident)){
		inter.ident.dat <- inter.ident.dat[inter.ident.dat$Cell.From %in% ligand.ident,]
	}

	lr.ident.pair <- inter.ident.dat$LR.Info
	lr.ident.pair <- unlist(sapply(lr.ident.pair,function(x){strsplit(x,split=';')}))
	names(lr.ident.pair) <- NULL
	lr.ident.pair <- unique(lr.ident.pair)

	lr.ident.split.pair <- sapply(lr.ident.pair, function(x){strsplit(x,split='--')})
	ident.ligs <- unlist(lapply(lr.ident.split.pair,function(x){x[1]}))
	ident.reps <- unlist(lapply(lr.ident.split.pair,function(x){x[2]}))

	inter.ident.unfold.dat <- bind_rows(replicate(length(lr.ident.split.pair), inter.ident.dat, simplify = FALSE))
	inter.ident.unfold.dat$Lig <- rep(ident.ligs, each=nrow(inter.ident.dat))
	inter.ident.unfold.dat$Rep <- rep(ident.reps, each=nrow(inter.ident.dat))
	inter.ident.unfold.dat$LR.Info <- paste0(inter.ident.unfold.dat$Lig,' --> ',inter.ident.unfold.dat$Rep)
	### fc.lr and p.lr are to save the measured FC and pval of each LR pair
	markerL.dat <- Interact$markerL 
	markerR.dat <- Interact$markerR 
	fc.lr <- c()
	p.lr <- c()
	for (each.row in 1:nrow(inter.ident.unfold.dat)){
		current.from <- inter.ident.unfold.dat[each.row,'Cell.From']
		current.to <- inter.ident.unfold.dat[each.row,'Cell.To']
		current.lig <- inter.ident.unfold.dat[each.row,'Lig']
		current.rep <- inter.ident.unfold.dat[each.row,'Rep']

		fc.lig <- subset(markerL.dat, cluster==current.from & gene==current.lig)[,'avg_log2FC']
		p.lig <- subset(markerL.dat, cluster==current.from & gene==current.lig)[,'p_val_adj']

		fc.rep <- subset(markerR.dat, cluster==current.to & gene==current.rep)[,'avg_log2FC']
		p.rep <- subset(markerR.dat, cluster==current.to & gene==current.rep)[,'p_val_adj']

		### if the ligand have a FC > logFC.thre and a p.adj < p.thre
		### and if the receptor have a FC > logFC.thre and a p.adj < p.thre
		if ((length(fc.lig)>0) & (length(fc.rep)>0)){
			if ((fc.lig > logFC.thre) &  (fc.rep > logFC.thre) & (p.lig < p.thre) & (p.rep < p.thre)){
				fc.lr <- c(fc.lr, fc.lig*fc.rep)
				p.lr <- c(p.lr, 1-(1-p.lig)*(1-p.rep))
			}
		}else{
			fc.lr <- c(fc.lr, NA)
			p.lr <- c(p.lr, NA)
		}
		
	}

	inter.ident.unfold.dat$Log2FC_LR <- fc.lr
	inter.ident.unfold.dat$P_LR <- p.lr
	inter.ident.unfold.dat$Log10_P_adj <- -log10(p.adjust(p.lr, method='BH'))

	if (return.data){
		return(inter.ident.unfold.dat[,c('Cell.From', 'Cell.To', 'LR.Info', 'Log2FC_LR', 'P_LR', 'Log10_P_adj')])
	}
	inter.ident.unfold.dat[which(inter.ident.unfold.dat$Log10_P_adj > 30), 'Log10_P_adj'] <- 30 


	if (length(ligand.ident)==1){
		x.title <- paste0('Receptor clusters for cluster ',ligand.ident)
	}else{
		x.title <- paste0('Ligand clusters to cluster ',receptor.ident)
	}

	### you may want to adjust the order of the clusters in the x axis
	inter.ident.unfold.dat$Xaxis <- as.character(inter.ident.unfold.dat$Xaxis)
	if (!is.null(ident.levels)){
		x.levels <- ident.levels[ident.levels %in% inter.ident.unfold.dat$Xaxis]
		if (all(inter.ident.unfold.dat$Xaxis %in% x.levels)){
			inter.ident.unfold.dat$Xaxis <- factor(inter.ident.unfold.dat$Xaxis, levels=x.levels)
		}else{
			ident.missed <- unique(inter.ident.unfold.dat$Xaxis[!(inter.ident.unfold.dat$Xaxis %in% x.levels)])
			ident.missed <- pasteIdent(ident.missed)
			stop(paste0('the ident class ',ident.missed,' may be missed in the input ident.levels'))
		}
	}

	plot <- ggplot(inter.ident.unfold.dat,aes(Xaxis,LR.Info)) + 
		geom_point(aes(size=Log2FC_LR,col=Log10_P_adj),shape=shape) +
		#scale_colour_gradient(low="green",high="red") + 
		scale_color_gradientn(values = seq(0,1,0.2),colours=c("#0570b0", "grey", "#d7301f")) + 
		labs(color='-log10(p.adjust)',size='Log2FC',x=x.title,y="") + 
		theme(axis.text=element_text(colour='black',size=12),
		axis.title=element_text(colour='black',size=12),
		panel.background=element_rect(fill="white",color="black"),
		panel.grid=element_line(size=0.5,colour='gray')
		)
	return(plot)
}


#' To find those pathways in which the genesets show overlap with the marker ligand and receptor genes in our dataset
#' @param Interact Interact list returned by findLRpairs
#' @param category Character to indicate which pathway to investigate; one of "go", "kegg", 'wikipathway', and "reactome", or "all" for all pathways
#' @return Interact list containing the ligand-receptor interaction information and the pathways showing overlap with the marker ligand and receptor genes in the dataset
#' @export
findLRpath <- function(Interact, category='all'){
	options(stringsAsFactors=F)
	if (category!='all' & category!='go' & category!='kegg' & category!='wikipathway' & category!='reactome'){
		stop("wrong category selected. Select one of 'go', 'kegg', 'wikipathway', and 'reactome', or 'all' for all pathways")
	}
	if (category=='all'){
		path.list <- commpathData$DataPathway[[Interact$species]]
	}else if(category=='wikipathway'){
		path.list <- commpathData$DataWikiPathway[[Interact$species]]
	}
	marker.lig.dat <- Interact$markerL
	marker.rep.dat <- Interact$markerR
	lr.gene <- unique(c(marker.lig.dat$gene,marker.rep.dat$gene))
	
	which.overlap.list <- unlist(
		lapply(path.list, function(x){ 
			if(any(x %in% lr.gene)){ return(TRUE) }else{return(FALSE)}
		}))
	path.lr.list <- path.list[which(which.overlap.list)]
	Interact$pathwayLR <- path.lr.list
	return(Interact)
}

#' To find different enriched pathway between two group cells
#' @param Interact Interact list returned by findLRpath
#' @param gsva.mat Matrix containing the pathway enrichment sorces, with rows representing pathways and columns representing cells. Pathway scores are usually computed from gsva, or other methods aiming to measure the pathway enrichment in cells
#' @param ident.label Vector indicating the identity labels of cells, and the order of labels are required to match order of cells (columns) in the gsva.mat
#' @param select.ident.1 Identity class to define cells for group 1
#' @param select.ident.2 Identity class to define cells for group 2 for comparison; if 'NULL', use all other cells for comparison
#' @param method Method used for differential enrichment analysis, either 't.test' of 'wilcox.test'
#' @return Dataframe including the statistic result comparing the pathway enrichment sorces between group 1 and group 2, the significant recetor and ligand of group 1 in the pathways, and the corresponding up stream identity class which interact with group 1 by releasing specific ligand
#' @export
diffPath <- function(Interact, gsva.mat, ident.label, select.ident.1, select.ident.2=NULL, method='t.test'){
	options(stringsAsFactors=F)

	if (method!='t.test' & method!='wilcox.test'){
		stop("select t.test or wilcox.test to conduct differential analysis")
	}else if(method=='t.test'){
		test.res.dat <- data.frame(matrix(NA,0,12))
		colnames(test.res.dat) <- c('mean.diff','mean.1','mean.2','t','df','p.val','p.val.adj','description','cell.up','ligand.up','receptor.in.path','ligand.in.path')
	}else{
		test.res.dat <- data.frame(matrix(NA,0,11))
		colnames(test.res.dat) <- c('median.diff','median.1','median.2','W','p.val','p.val.adj','description','cell.up','ligand.up','receptor.in.path','ligand.in.path')
	}

	path.lr.list <- Interact$pathwayLR
	if (is.null(path.lr.list)){
		stop("no pathway detected, run findLRpath befor diffPath")
	}
	# marker.lig.dat <- Interact$markerL
	# marker.lig.dat <- marker.lig.dat[marker.lig.dat$cluster==select.ident,]
	marker.ident1.rep.dat <- subset(Interact$InteractGeneUnfold, Cell.To %in% select.ident.1)

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
	marker.lig.dat <- Interact$markerL
	ident.lig.vec <- marker.lig.dat[marker.lig.dat$cluster %in% select.ident.1,'gene']
	### for each pathway in the DEG result, find which pathway show overlap with marker ligands of the selected ident
	ident.lig.in.path <- sapply(test.res.dat$description, function(eachPath){
		each.set <- path.lr.list[[eachPath]]
		if (any(ident.lig.vec %in% each.set)){
			return(paste(ident.lig.vec[ident.lig.vec %in% each.set],collapse=','))
		}else{
			return(NA)
		}
	})
	test.res.dat$ligand.in.path <- unlist(ident.lig.in.path)

	return(test.res.dat)
}


#' To find different enriched pathways in each identity class 
#' @param Interact Interact list returned by findLRpath
#' @param gsva.mat Matrix containing the pathway enrichment sorces, with rows representing pathways and columns representing cells. Pathway scores are usually computed from gsva, or other methods aiming to measure the pathway enrichment in cells
#' @param ident.label Vector indicating the identity labels of cells, and the order of labels are required to match order of cells (columns) in the gsva.mat
#' @param method Method used for differential enrichment analysis, either 't.test' of 'wilcox.test'
#' @return Dataframe including the statistic result comparing the pathway enrichment sorces between cells in each cluster and all other clusters, the significant recetor and ligand in the pathways, and the corresponding up stream identity class and ligand
#' @export
diffAllPath <- function(Interact, gsva.mat, ident.label, method='t.test'){
	if (method=='t.test'){
		all.test.dat <- data.frame(matrix(NA,0,12))
	}else{
		all.test.dat <- data.frame(matrix(NA,0,11))
	}

	all.rep.ident <- unique(Interact$InteractGene$Cell.To)
	for (each.ident in all.rep.ident){
		message(paste0('Identifying pathways for cluster ',each.ident,'...'))
		test.res.dat <- diffPath(Interact, gsva.mat, ident.label, select.ident.1=each.ident, select.ident.2=NULL, method=method)
		if (nrow(test.res.dat)==0){ next }
		test.res.dat$cluster <- each.ident
		all.test.dat <- rbind(all.test.dat, test.res.dat)
	}

	return(all.test.dat)
}

#' To find the downstream identity class of specific ligand released by specific upstream identity class
#' @param Interact Interact list returned by findLRpairs
#' @param select.ident Upstream identity class; if 'NULL', use all identity classes
#' @param select.ligand Ligand released by upstream identity class; if 'NULL', use all ligands that are markers for the selected upstream identity class
#' @return Dataframe including the interaction information
#' @export
findReceptor <- function(Interact, select.ident=NULL, select.ligand=NULL){
	options(stringsAsFactors=F)

	if (is.null(select.ident) & is.null(select.ligand)){
		stop("either a select.ident or a select.ligand need to be asigned")
	}
	if (is.null(select.ligand)){
		ident.down.dat <- subset(Interact$InteractGeneUnfold, Cell.From==select.ident)
	}else if(is.null(select.ident)){
		ident.down.dat <- subset(Interact$InteractGeneUnfold, Ligand==select.ligand)
	}else{
		ident.down.dat <- subset(Interact$InteractGeneUnfold, Cell.From==select.ident & Ligand==select.ligand)
	}

	if (nrow(ident.down.dat)==0){
		warning("no downstream ident found for the selected ident and ligand")
	}
	return(ident.down.dat)
}

#' To find the upstream identity classes and ligands of specific receptor expressed by specific downstream identity class
#' @param Interact Interact list returned by findLRpairs
#' @param select.ident Downstream identity class; if 'NULL', use all identity classes
#' @param select.receptor Receptor expressed by downstream identity class; if 'NULL', use all receptors that are markers for the selected downstream identity class
#' @return Dataframe including the interaction information
#' @export
findReceptor <- function(Interact, select.ident=NULL, select.receptor=NULL){
	options(stringsAsFactors=F)

	if (is.null(select.ident) & is.null(select.receptor)){
		stop("either a select.ident or a select.receptor need to be asigned")
	}
	if (is.null(select.receptor)){
		ident.up.dat <- subset(Interact$InteractGeneUnfold, Cell.To==select.ident)
	}else if(is.null(select.ident)){
		ident.up.dat <- subset(Interact$InteractGeneUnfold, Receptor==select.receptor)
	}else{
		ident.up.dat <- subset(Interact$InteractGeneUnfold, Cell.To==select.ident & Receptor==select.receptor)
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

#' This is a plug-in function, aimming to paste a vector of idents into a ',' separated string
#' @param ident.missed Vector of idents
#' @return String with idents pasted
pasteIdent <- function(ident.missed){
	if (length(ident.missed) > 1){ 
		ident.missed <- paste0(paste(ident.missed[1:(length(ident.missed)-1)], collapse=', '), ' and ', ident.missed[length(ident.missed)])
	}
	return(ident.missed)
}

#' This is a plug-in function, aimming to check whether the specified order of idents by users contain all idents that present in the dataset
#' @param all.ident Vector of all idents that present in the dataset
#' @param order Vector of the specified order of idents by users
#' @return if all idents are not contained, stop the procedure 
orderCheck <- function(all.ident, order){
	if (all(all.ident %in% order)){
		order.overlap <- order[order %in% all.ident]
		all.ident <- factor(all.ident, levels=order.overlap)
		return(all.ident)
	}else{
		ident.missed <- all.ident[!(all.ident %in% order)]
		ident.missed <- pasteIdent(ident.missed)
		stop(paste0('the ident class ',ident.missed,' may be missed in the input order'))
	}
}

#' This is a plug-in function, aimming to find the highly variable pathways in the gsva score matrix
#' @param gsva.mat Matrix containing the pathway enrichment sorces, with rows representing pathways and columns representing cells. Pathway scores are usually computed from gsva, or other methods aiming to measure the pathway enrichment in cells
#' @param select.path selected pathways
#' @param n n
#' @return Vector of top n variable pathways
#' @export
variPath <- function(gsva.mat, select.path=NULL, n=10){
	if (!is.null(select.path)){
		gsva.mat <- gsva.mat[rownames(gsva.mat) %in% select.path,]
	}
	gsva.cv.vec <- apply(gsva.mat,1,sd)
	vari.path <- names(gsva.cv.vec[order(-gsva.cv.vec)][1:n])

	return(vari.path)
}

#' To present the interact and variable pathways in a line plot
#' @param Interact Interact list returned by findLRpairs
#' @param gsva.mat pathway score
#' @param all.test.dat all.test.dat
#' @param order Order of identity classes in the plot
#' @param n Top n variable pathways
#' @param hide.isol Whether to hide the isolated identity (those identity class in which cells show no interact with cells in other ientity class) or not; defult is TRUE
#' @param line.width line.width
#' @return String with idents pasted
#' @export
pathPlot <- function(Interact, gsva.mat, all.test.dat, order=NULL, n=10, hide.isol=TRUE, ident.size=1, label.size=1, path.size=1, ident.line.width=1, ident.line.alpha=0.8, path.line.alpha=0.8, path.line.width=1, label.dist=0.2){
	
	#order=paste0('C',1:20); n=10; hide.isol=TRUE; line.width=1; label.dist=0.2

	vari.path.name <- variPath(gsva.mat, all.test.dat$description, n=n)
	vari.path.name <- rev(vari.path.name)

	vari.path.dat <- subset(all.test.dat, description %in% vari.path.name)
	#vari.path.dat$cluster <- paste0('C',vari.path.dat$cluster)
	Interact.num.dat <- Interact$InteractNumer
	Interact.num.dezero.dat <- subset(Interact.num.dat, LR.Number!=0)
	#Interact.num.dat$Cell.From <- paste0('C',Interact.num.dat$Cell.From)
	#Interact.num.dat$Cell.To <- paste0('C',Interact.num.dat$Cell.To)

	if (hide.isol){ Interact.num.dat <- Interact.num.dezero.dat }
	# Interact.num.dat <- subset(Interact.num.dat, Cell.From %in% c(1:6) & Cell.To %in% c(1:6) )

	all.ident <- unique(c(Interact.num.dat$Cell.From,Interact.num.dat$Cell.To))
	if (!is.null(order)){
		order.overlap <- order[order %in% all.ident]
		all.ident <- factor(all.ident, levels=order.overlap)
	}
	all.ident <- all.ident[order(all.ident, decreasing=TRUE)]

	n.ident <- length(all.ident)
	col.ident <- rev(scales::hue_pal(c=100)(n.ident))

	lig.ident <- unique(Interact.num.dat$Cell.From)
	lig.coor <- match(lig.ident, all.ident)
	lig.col <- col.ident[lig.coor]

	rep.ident <- unique(Interact.num.dat$Cell.To)
	rep.coor <- match(rep.ident, all.ident)
	rep.col <- col.ident[rep.coor]

	line.y.coor <- match(Interact.num.dezero.dat$Cell.From, all.ident)
	line.yend.coor <- match(Interact.num.dezero.dat$Cell.To, all.ident)
	line.col <- col.ident[line.y.coor]
	line.size <- 0.1 * ident.line.width * Interact.num.dezero.dat$LR.Number

	pathway.coor <- seq(from=1, to=n.ident, length.out=length(vari.path.name))
	path.line.y.coor <- match(vari.path.dat$cluster, all.ident)
	path.line.yend.coor <- pathway.coor[match(vari.path.dat$description, vari.path.name)]
	path.line.col <- col.ident[path.line.y.coor]

	pathplot.theme <- theme(axis.title=element_blank(), axis.text=element_blank(),axis.ticks=element_blank(), axis.line=element_blank(), panel.background=element_rect(fill="white"))

	ggplot() + 
	annotate('text', hjust=0, x=1-label.dist, y=1:n.ident,label=all.ident, size=5*label.size) +
	geom_segment(aes(x=1, xend=2, y=line.y.coor,yend=line.yend.coor, size=line.size), color=line.col, alpha=ident.line.alpha, show.legend=F) +
	geom_point(aes(x=1,y=lig.coor),size=8*ident.size,color=lig.col) + 
	geom_point(aes(x=2,y=rep.coor),size=8*ident.size,color=rep.col) + 
	annotate('label', hjust=0, x=3, y=pathway.coor,label=vari.path.name, size=5*path.size) +
	geom_segment(aes(x=2, xend=3, y=path.line.y.coor,yend=path.line.yend.coor, size=0.1*path.line.width), color=path.line.col, alpha=path.line.alpha, show.legend=F) + 
	#scale_x_continuous(limits=c(0.8,5)) + 
	pathplot.theme

}

#' To present the interact and variable pathways in a line plot, for a specific ident
#' @param Interact Interact list returned by findLRpairs
#' @param select.ident select.ident
#' @param ident.path.dat ident.path.dat for the selected ident
#' @param top.n.receptor top.n.receptor
#' @param order order
#' @param p.thre p.thre
#' @param ident.size ident.size
#' @return String with idents pasted
#' @export
receptorPathPlot <- function(Interact, select.ident, ident.path.dat, top.n.receptor=5, top.n.path=10, order=NULL, p.thre = 0.05, ident.size=1, receptor.size=1, label.dist=0.4, label.size=1, line.size=1, ident.line.alpha=1, ident.line.width=1, path.line.width=1, path.line.alpha=0.8, path.size=1){
	options(stringsAsFactors=F)

	if ('t' %in% colnames(ident.path.dat)){
		ident.path.dat <- subset(ident.path.dat, p.val.adj < p.thre & t > 0)
	}else if ('W' %in% colnames(ident.path.dat)){
		ident.path.dat <- subset(ident.path.dat, p.val.adj < p.thre & median.diff > 0)
	}else{
		stop('please input the integrate ident.path.dat computed from diffPath')
	}
	ident.path.dat <- ident.path.dat[order(ident.path.dat$p.val.adj, decreasing=TRUE),]

	all.sig.path <- ident.path.dat$description
	#ident.path.dat <- ident.path.dat[1:n, ]

	### preprocess and extract useful information
	up.ident <- as.vector(unlist(sapply(ident.path.dat$cell.up, function(x){ strsplit(x, split=';') })))
	cur.rep <- as.vector(unlist(sapply(ident.path.dat$receptor.in.path, function(x){ strsplit(x, split=';') })))

	n.elem.each.row <- as.vector(unlist(lapply(sapply(ident.path.dat$receptor.in.path, function(x){ strsplit(x, split=';') }),length)))
	path.all <- rep(ident.path.dat$description, times=n.elem.each.row)
	plot.dat <- data.frame(up.ident=up.ident, cur.rep=cur.rep, path.name=path.all)
	
	### select the highly expressed receptors
	markerR.dat <- Interact$markerR
	#markerR.dat <- markerR.dat[markerR.dat$cluster==select.ident, ]
	markerR.dat <- subset(markerR.dat, subset=cluster==select.ident)
	if (nrow(markerR.dat)==0){
		stop('there is no marker receptor for the selected ident')
	}

	markerR.dat <- markerR.dat[order(markerR.dat$avg_log2FC, decreasing=TRUE), ]
	if (top.n.receptor > nrow(markerR.dat)){
		warning(paste0('there is(are) ', nrow(markerR.dat),' marker receptor(s) in the selected ident, and the input top.n.receptor is ', top.n.receptor))
		top.n.receptor <- nrow(markerR.dat)
	}
	top.receptor.dat <- markerR.dat[1:top.n.receptor,]
	plot.dat <- subset(plot.dat, cur.rep %in% top.receptor.dat$gene)
	
	### to select the top n pathways associated with the selected receptors
	path.uniq.name <- unique(plot.dat$path.name)
	### those pathways with a larger position index in the all.sig.path are more significant
	### the following codes is to get the pathways with the top n largest position index
	path.position <- match(path.uniq.name, all.sig.path)
	path.position <- path.position[order(path.position, decreasing=TRUE)]
	if (top.n.path > length(path.position)){
		warning(paste0('there is(are) ', length(path.position),' significant pathway(s) for the selected ident, and the input top.n.path is ', top.n.path))
		top.n.path <- length(path.position)
	}
	path.position <- path.position[1:top.n.path]
	top.path <- all.sig.path[path.position]
	plot.dat$path.name[which(!(plot.dat$path.name %in% top.path))] <- NA

	### data for plot of idents and receptors
	plot.ident.dat <- unique(plot.dat[,c('up.ident','cur.rep')])
	up.ident <- plot.ident.dat$up.ident
	cur.rep <- plot.ident.dat$cur.rep
	### the idents which release ligands
	up.uniq.ident <- unique(up.ident)
	### to check the order parameter
	if (!is.null(order)){
		up.uniq.ident <- orderCheck(up.uniq.ident, order)
	}
	up.uniq.ident <- up.uniq.ident[order(up.uniq.ident, decreasing=TRUE)]
	n.ident <- length(up.uniq.ident)
	col.ident <- rev(scales::hue_pal(c=100)(n.ident))
	lig.coor <- 1:n.ident

	### the receptor and their coordinate
	cur.uniq.rep <- factor(unique(cur.rep))
	cur.uniq.rep <- cur.uniq.rep[order(cur.uniq.rep, decreasing=TRUE)]
	n.rep <- length(cur.uniq.rep)
	rep.coor <- seq(from=1, to=n.ident, length.out=n.rep)
	rep.col <- rev(scales::hue_pal(c=100)(n.rep))

	logfc <- top.receptor.dat[match(cur.uniq.rep, top.receptor.dat$gene),'avg_log2FC']
	rep.size <- 5 * receptor.size * logfc
	
	### lines between ligand idents and receptors
	line.y.coor <- match(up.ident, up.uniq.ident)
	line.yend.coor <- rep.coor[match(cur.rep, cur.uniq.rep)]
	line.col <- col.ident[match(up.ident, up.uniq.ident)]

	### pathways coordinate and color
	plot.path.dat <- subset(plot.dat, subset=!(is.na(path.name)), select=c('cur.rep','path.name') )
	plot.path.dat <- unique(plot.path.dat)
	path.name <- plot.path.dat$path.name
	path.uniq.name <- unique(path.name)
	path.uniq.name <- factor(path.uniq.name, levels=top.path)
	path.uniq.name <- path.uniq.name[order(path.uniq.name, decreasing=TRUE)]
	width <- max(sapply(as.vector(path.uniq.name),nchar))/10
	pathway.coor <- seq(from=1, to=n.ident, length.out=length(path.uniq.name))
	
	### lines between receptors and pathways
	cur.rep <- plot.path.dat$cur.rep
	path.line.y.coor <- rep.coor[match(cur.rep, cur.uniq.rep)]
	path.line.yend.coor <- pathway.coor[match(path.name, path.uniq.name)]
	path.line.col <- rep.col[match(cur.rep, cur.uniq.rep)]

	pathplot.theme <- theme(axis.title=element_blank(), axis.text=element_blank(),axis.ticks=element_blank(), axis.line=element_blank(), panel.background=element_rect(fill="white"))

	ggplot() + 
	scale_x_continuous(limits=c(0.5,width)) + 
	geom_point(aes(x=1,y=lig.coor),size=8*ident.size,color=col.ident) + 
	annotate('text', hjust=0, x=1-label.dist, y=1:n.ident,label=up.uniq.ident, size=5*label.size) +
	geom_segment(aes(x=1, xend=2, y=line.y.coor,yend=line.yend.coor), size=ident.line.width, color=line.col, alpha=ident.line.alpha, show.legend=F) +
	annotate('label', hjust=0, x=3, y=pathway.coor,label=path.uniq.name, size=5*path.size) +
	geom_segment(aes(x=2, xend=3, y=path.line.y.coor,yend=path.line.yend.coor),size=path.line.width, color=path.line.col, alpha=path.line.alpha, show.legend=F) + 
	geom_point(aes(x=2,y=rep.coor),size=rep.size,color=rep.col) + 
	annotate('text', hjust=0.5, x=2, y=rep.coor-0.4,label=cur.uniq.rep, size=5*label.size) +
	annotate('text', hjust=0.5, x=1:2, y=max(lig.coor)+0.6,label=c('Source', 'Receptor'), size=5*label.size) +
	annotate('text', hjust=0, x=3, y=max(lig.coor)+0.6,label=c('Pathway'), size=5*label.size) +
	pathplot.theme

}