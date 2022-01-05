data('commpathData', package=('commpath'))

#' To present a circos plot
#' @param Interact Interact list returned by findLRpairs
#' @param order Vector of orders of the identities arround the circos plot
#' @param col Vector of colors of each identity; names of the col vector are supposed to be assigned to indicate each color for each identity
#' @param ident To highlight the interaction between a specific identity class and others; if 'NULL', plot interaction for all identity classes
#' @param name.vert Should the group annotation be vertical to the grid? Defualt is FALSE
#' @return Circos plot showing the ligand-receptor interaction
#' @export
circosPlot <- function(Interact, order=NULL, col=NULL, ident=NULL, name.vert=FALSE){
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
		if (!name.vert){
			chordDiagram(Interact.num.dat, order=all.ident, grid.col=col, annotationTrack=c("name","grid"), transparency=0.1, directional=1, direction.type='arrows', link.arr.type = "big.arrow",preAllocateTracks = list(track.height=mm_h(5)), annotationTrackHeight=convert_height(c(1, 2), "mm"))
		}else{
			chordDiagram(Interact.num.dat, order=all.ident, grid.col=col, annotationTrack=c("grid"), transparency=0.1, directional=1, direction.type='arrows', link.arr.type = "big.arrow", preAllocateTracks = list(track.height=mm_h(5)), annotationTrackHeight=convert_height(c(2, 2), "mm"))
			circos.track(track.index=1, panel.fun=function(x, y){circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)}, bg.border = NA)
		}

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
		
		link.zindex <- 1:nrow(Interact.num.dat)
		ident.line <- (Interact.num.dat$Cell.From==ident) | (Interact.num.dat$Cell.To==ident)
		link.zindex[!ident.line] <- rank(-(Interact.num.dat[!ident.line,'LR.Number']))
		link.zindex[ident.line] <- rank(-(Interact.num.dat[ident.line,'LR.Number']))+length(which(!ident.line))
		if (!name.vert){
			chordDiagram(Interact.num.dat, order=all.ident, grid.col=col, annotationTrack=c("name","grid"), transparency=0.1, directional=1, direction.type='arrows', link.arr.type = "big.arrow", link.zindex=link.zindex, col=line.col, preAllocateTracks = list(track.height=mm_h(5)), annotationTrackHeight=convert_height(c(1, 2), "mm"))
		}else{
			chordDiagram(Interact.num.dat, order=all.ident, grid.col=col, annotationTrack="grid", transparency=0.1, directional=1, direction.type='arrows', link.arr.type = "big.arrow", link.zindex=link.zindex, col=line.col, preAllocateTracks = list(track.height=mm_h(5)), annotationTrackHeight=convert_height(c(2, 2), "mm"))
			circos.track(track.index=1, panel.fun=function(x, y){circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)}, bg.border = NA)
		}
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
	Interact <- list(InteractNumer=Interact.num.dat,InteractGene=Interact.gene.dat, InteractGeneUnfold=lr.unfold.dat, markerL=marker.lig.dat, markerR=marker.rep.dat, logFC.thre=logFC.thre, p.thre=p.thre, species=species)
	return(Interact)
}

#' To present a dot plot for specific ligand-receptor pairs in specific clusters
#' @param Interact Interact list returned by findLRpairs
#' @param ligand.ident Vector containing the ligand ident
#' @param receptor.ident Vector containing the receptor ident
#' @param ident.levels Vector of levels of the identities 
#' @param top.n.inter Show the dotplot for the top n LR pairs with the largest product of Log2FC
#' @param return.data Logical value indicating whether to return the data for the plot or not
#' @return Dotplot showing the ligand-receptor interaction between the selected ligand.ident and receptor.ident
#' @export
dotPlot <- function(Interact, ligand.ident=NULL, receptor.ident=NULL, ident.levels=NULL, top.n.inter=10, return.data=FALSE){
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
		#shape <- 16
	}else{
		inter.gene.dat$Xaxis <- inter.gene.dat$Cell.From
		#shape <- 17
	}

	logFC.thre <- Interact$logFC.thre
	p.thre <- Interact$p.thre

	inter.ident.dat <- inter.gene.dat
	if (!is.null(receptor.ident)){
		if (!all(receptor.ident %in% inter.ident.dat$Cell.To)){
			stop(paste0('select receptor.ident from ',pasteIdent(inter.ident.dat$Cell.To)))
		}
		inter.ident.dat <- inter.ident.dat[inter.ident.dat$Cell.To %in% receptor.ident,]
	}
	if (!is.null(ligand.ident)){
		if (!all(ligand.ident %in% inter.ident.dat$Cell.From)){
			stop(paste0('select ligand.ident from ',pasteIdent(inter.ident.dat$Cell.From)))
		}
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
		#print(each.row)
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

	inter.ident.unfold.dat <- inter.ident.unfold.dat[,c('Xaxis', 'LR.Info', 'Log2FC_LR', 'Log10_P_adj')]

	max.fc <- by(data=inter.ident.unfold.dat$Log2FC_LR, INDICES = inter.ident.unfold.dat$LR.Info, function(x){max(x,na.rm=T)})
	max.LRname.fc <- names(max.fc)[order(max.fc, decreasing=TRUE)]
	
	if (top.n.inter > length(max.LRname.fc)){
		warning(paste0('there is(are) ', length(max.LRname.fc),' significant LR pair(s) for the selected ident, and the input top.n.inter is ', top.n.inter))
		top.n.inter <- length(max.LRname.fc)
	}

	max.LRname.fc <- max.LRname.fc[1:top.n.inter]
	inter.ident.unfold.dat <- subset(inter.ident.unfold.dat, LR.Info %in% max.LRname.fc)

	inter.ident.unfold.dat[which(inter.ident.unfold.dat$Log10_P_adj > 30), 'Log10_P_adj'] <- 30
	plot <- ggplot(inter.ident.unfold.dat,aes(Xaxis,LR.Info)) + 
		geom_point(aes(size=Log2FC_LR,col=Log10_P_adj)) +
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
#' @param category Character to indicate which pathway to investigate; one of "go" (GO terms), "kegg" (for KEGG pathways), 'wiki' (for WikiPathways), and "reactome" (for reactome pathways), or "all" for all pathways
#' @return Interact list containing the ligand-receptor interaction information and the pathways showing overlap with the marker ligand and receptor genes in the dataset
#' @export
findLRpath <- function(Interact, category='all'){
	options(stringsAsFactors=F)
	if (category!='all' & category!='go' & category!='kegg' & category!='wiki' & category!='reactome'){
		stop("wrong category selected. Select one of 'go', 'kegg', 'wiki', and 'reactome', or 'all' for all pathways")
	}

	### select one category
	if (category=='all'){
		path.go.list <- commpathData$DataGOterm[[Interact$species]]
		path.kegg.list <- commpathData$DataKEGGpathway[[Interact$species]]
		path.wiki.list <- commpathData$DataWikiPathway[[Interact$species]]
		path.reac.list <- commpathData$DataReactome[[Interact$species]]
		path.list <- c(path.go.list, path.kegg.list, path.wiki.list, path.reac.list)
		path.list <- path.list[which(!duplicated(names(path.list)))]

	}else if (category=='go'){
		path.list <- commpathData$DataGOterm[[Interact$species]]
	}else if(category=='kegg'){
		path.list <- commpathData$DataKEGGpathway[[Interact$species]]
	}else if(category=='wiki'){
		path.list <- commpathData$DataWikiPathway[[Interact$species]]
	}else if(category=='reactome'){
		path.list <- commpathData$DataReactome[[Interact$species]]
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
			return(paste(ident.lig.vec[ident.lig.vec %in% each.set],collapse=';'))
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

#' To plot a heatmap of those differentially enriched pathways for each cluster
#' @param gsva.mat Matrix containing the pathway enrichment sorces, with rows representing pathways and columns representing cells. Pathway scores are usually computed from gsva, or other methods aiming to measure the pathway enrichment in cells
#' @param ident.label Vector indicating the identity labels of cells, and the order of labels are required to match order of cells (columns) in the gsva.mat
#' @param all.path.dat data.frame of differential enrichment test result from diffAllPath
#' @param top.n.pathway Show the heatmap of top n most significant pathways
#' @param sort Sort criteria used to select the top n pathways, either p.val.adj, p.val, mean.diff, or median.diff
#' @param col Vector of colors used to generate a series of gradient colors to show the enrichment score of pathways; provide a vector containing at least two colors
#' @param show.legend Whether to show the legend or not; default is FALSE
#' @param cell.label.size Text size of the label of cell types
#' @param cell.label.angle Rotation angle of the label of cell types
#' @param pathway.label.size Text size of the label of pathways
#' @param scale Whether to scale the enrichment sorces matrix among cells or not; default is TRUE
#' @param truncation Truncation fold value; scores > (the third quartiles + truncation * interquartile range) and scores < (the first quartiles - truncation * interquartile range) will be adjusted; either a value to indicate the specific truncation value or 'none' to indicate no truncation; default is 1.5
#' @return Heatmap plot showing the top enriched patways in each cluster
#' @export
pathHeatmap <- function(gsva.mat, ident.label, all.path.dat, top.n.pathway=10, sort='p.val.adj', col=NULL, show.legend=FALSE, cell.label.size=NULL, cell.label.angle=0, pathway.label.size=NULL, scale=TRUE, truncation=1.5){
	### select the highly enriched pathways
	if ('mean.diff' %in% colnames(all.path.dat)){
		all.path.dat <- subset(all.path.dat, mean.diff > 0)
	}else if('median.diff' %in% colnames(all.path.dat)){
		all.path.dat <- subset(all.path.dat, median.diff > 0)
	}else{
		stop('please input the intact test result computed from diffAllPath')
	}

	### select top n pathways
	if (sort=='p.val.adj' | sort=='p.val'){
		all.path.dat <- all.path.dat[order(all.path.dat[[sort]], decreasing=FALSE), ]
	}else if (sort=='mean.diff' | sort=='median.diff'){
		all.path.dat <- all.path.dat[order(all.path.dat[[sort]], decreasing=TRUE), ]
	}else{
		stop('select one sort criteria from p.val.adj, p.val, mean.diff, and median.diff')
	}

	if (!is.factor(ident.label)){
		ident.label <- factor(ident.label)
	}
	all.ident <- levels(ident.label)

	path <- sapply(all.ident, function(each.label){
		ident.path.dat <- subset(all.path.dat, cluster==each.label)
		if (top.n.pathway > nrow(ident.path.dat)){
			top.n.pathway <- nrow(ident.path.dat)
			warning(paste0('there is(are) ',nrow(ident.path.dat),' significant pathway(s) for cluster ',each.label,', and the input top.n.pathway is ',top.n.pathway))
		}
		ident.path.dat[1:top.n.pathway, 'description']
	})
	path <- unique(as.character(path))
	
	gsva.mat <- gsva.mat[path,]

	### scale or not
	if (scale){
		gsva.scale.mat <- t(scale(t(gsva.mat)))
	}else{
		gsva.scale.mat <- gsva.mat
	}

	gsva.dat <- melt(gsva.scale.mat, varnames=c('Pathway','Cell'),value.name="Score",  na.rm=TRUE)
	gsva.dat$Pathway <- factor(as.character(gsva.dat$Pathway), levels=rev(path))
	gsva.dat$Cell <- as.character(gsva.dat$Cell)

	cell.anno.dat <- data.frame(Cellname=colnames(gsva.scale.mat),Celltype=ident.label,stringsAsFactors=F)
	gsva.dat$Celltype <- cell.anno.dat[match(gsva.dat$Cell, cell.anno.dat$Cellname),'Celltype']
	
	if (is.numeric(truncation)){
		IQR <- quantile(gsva.dat$Score, 0.75) - quantile(gsva.dat$Score, 0.25)
		max.score <- quantile(gsva.dat$Score, 0.75) + IQR * truncation
		min.score <- quantile(gsva.dat$Score, 0.25) - IQR * truncation
		gsva.dat[gsva.dat$Score > max.score, 'Score'] <- max.score
		gsva.dat[gsva.dat$Score < min.score, 'Score'] <- min.score
	}else if(truncation!='none'){
		stop('either a specific value or "none" is needed for the truncation parameter')
	}
	
	### basic plot
	path.heatmap <- ggplot(data=gsva.dat, aes(x=Cell, y=Pathway, fill=Score)) + 
		geom_tile() +
		facet_grid(facets=~Celltype,scales="free",space="free",switch='x') +
		labs(x='',y='') +
		theme(axis.text.x=element_blank(),
			axis.text.y=element_text(colour='black'),
			axis.ticks=element_blank(),
			axis.line=element_blank(),
			strip.background = element_blank(),
			strip.text.x=element_text(colour="black", angle=cell.label.angle))

	### gradient color
	if (is.null(col)){
		colours <- c("#440154" ,"#21908C", "#FDE725")
	}else if(length(col)==1){
		stop('select at least 2 colors to generate a series of gradient colors')
	}else{
		colours <- col
	}
	path.heatmap <- path.heatmap + scale_fill_gradientn(colours=colours)

	### show legend or not
	if (show.legend){
		path.heatmap <- path.heatmap + theme(legend.title=element_text(colour='black',size=14),
			legend.text=element_text(colour='black',size=14))
	}else{
		path.heatmap <- path.heatmap + theme(legend.position='none')
	}

	### cell and pathway label
	if (!is.null(cell.label.size)){
		path.heatmap <- path.heatmap + theme(strip.text.x=element_text(size=cell.label.size))
	}
	if (!is.null(pathway.label.size)){
		path.heatmap <- path.heatmap + theme(axis.text.y=element_text(size=pathway.label.size))
	}
	return(path.heatmap)
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
findLigand <- function(Interact, select.ident=NULL, select.receptor=NULL){
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

#' To present the interact and variable pathways in a line plot, for a specific ident
#' @param Interact Interact list returned by findLRpairs
#' @param select.ident select.ident
#' @param ident.path.dat ident.path.dat for the selected ident
#' @param top.n.receptor top n receptor with the largest log2FCs to plot
#' @param order order of source clusters
#' @param top.n.path top n pathways with the smallest adjusted p values to plot
#' @param p.thre threshold for adjust p values; Only pathways with a adjust p valua < p.thre would be considered
#' @param dot.ident.col color of the dots representing upstream clusters
#' @param dot.receptor.col color of the dots representing reptors
#' @param bar.pathway.col color of the bars reprsenting pathways
#' @param bar.pathway.width width of the bars reprsenting pathways
#' @param dot.ident.size size of the dots representing source clusters
#' @param dot.receptor.size size of the dots representing receptors
#' @param label.text.size text size in the plot
#' @param label.title.size text size of the title annotation of the plot 
#' @param line.ident.width text size of pathway labels
#' @param line.path.width width of lines connecting receptors and pathways
#' @return Network plot showing receptors in the selected cluster, the upstream clusters which show L-R connections with the selected cluster, and the significant pathways involved in the receptors
#' @export
receptorPathPlot <- function(Interact, select.ident, ident.path.dat, top.n.receptor=5, order=NULL, top.n.path=10, p.thre=0.05, dot.ident.col=NULL, dot.receptor.col=NULL, bar.pathway.col=NULL, bar.pathway.width=10, dot.ident.size=1, dot.receptor.size=1, label.text.size=1, label.title.size=1, line.ident.width=1, line.path.width=1){
	options(stringsAsFactors=F)

	### check the input and select significant pathways
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
	up.ident <- as.vector(unlist(sapply(ident.path.dat$cell.up, function(x){ strsplit(x, split=';') })))
	cur.rep <- as.vector(unlist(sapply(ident.path.dat$receptor.in.path, function(x){ strsplit(x, split=';') })))

	n.elem.each.row <- as.vector(unlist(lapply(sapply(ident.path.dat$receptor.in.path, function(x){ strsplit(x, split=';') }),length)))
	path.all <- rep(ident.path.dat$description, times=n.elem.each.row)
	plot.dat <- data.frame(up.ident=up.ident, cur.rep=cur.rep, path.name=path.all)
	
	### select the highly expressed receptors
	markerR.dat <- subset(Interact$markerR, subset=(cluster==select.ident & gene %in% cur.rep))
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
	### those pathways with a larger position index in the all.sig.path are more significant
	### the following codes is to get the pathways with the top n largest position index
	path.uniq.name <- unique(plot.dat$path.name)
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
	plot.ident.to.receprtor.dat <- unique(plot.dat[,c('up.ident','cur.rep')])
	up.ident <- plot.ident.to.receprtor.dat$up.ident
	cur.rep <- plot.ident.to.receprtor.dat$cur.rep
	### the idents which release ligands
	up.uniq.ident <- unique(up.ident)
	### to check the order parameter
	if (!is.null(order)){
		up.uniq.ident <- orderCheck(up.uniq.ident, order)
	}
	up.uniq.ident <- up.uniq.ident[order(up.uniq.ident, decreasing=TRUE)]
	n.up.ident <- length(up.uniq.ident)
	if (is.null(dot.ident.col)){
		dot.ident.col <- rev(scales::hue_pal(c=100)(n.up.ident))
	}else if(length(dot.ident.col) < n.up.ident){
		stop(paste0('there is(are) ',n.up.ident,' cluster(s), but only ',length(dot.ident.col),' colors provide.'))
	}else{
		dot.ident.col <- rev(dot.ident.col[1:n.up.ident])
	}
	up.ident.coor <- 1:n.up.ident

	### the receptor and their coordinate, and their color
	cur.uniq.rep <- factor(unique(cur.rep))
	cur.uniq.rep <- cur.uniq.rep[order(cur.uniq.rep, decreasing=TRUE)]
	n.rep <- length(cur.uniq.rep)
	rep.coor <- seq(from=1, to=n.up.ident, length.out=n.rep)

	logfc <- top.receptor.dat[match(cur.uniq.rep, top.receptor.dat$gene),'avg_log2FC']
	rep.size <- 5 * dot.receptor.size * logfc
	rep.pct <-  top.receptor.dat[match(cur.uniq.rep, top.receptor.dat$gene),'pct.1']
	rep.pval <-  top.receptor.dat[match(cur.uniq.rep, top.receptor.dat$gene),'p_val_adj']
	rep.pval <- p.remove.inf(rep.pval)
	dot.rep.col <- LRcolor(rep.pval, user.set.col=dot.receptor.col)

	### lines between ligand idents and receptors
	line.ident.to.rep.y.coor <- match(up.ident, up.uniq.ident)
	line.ident.to.rep.yend.coor <- rep.coor[match(cur.rep, cur.uniq.rep)]
	line.col <- dot.ident.col[match(up.ident, up.uniq.ident)]

	### pathways coordinate and color
	plot.ligand.to.path.dat <- unique(subset(plot.dat, subset=!(is.na(path.name)), select=c('cur.rep','path.name')))
	path.name <- plot.ligand.to.path.dat$path.name
	path.uniq.name <- unique(path.name)
	path.uniq.name <- factor(path.uniq.name, levels=top.path)
	path.uniq.name <- path.uniq.name[order(path.uniq.name, decreasing=TRUE)]
	pathway.coor <- seq(from=1, to=n.up.ident, length.out=length(path.uniq.name))
	
	### lines between receptors and pathways
	cur.rep <- plot.ligand.to.path.dat$cur.rep
	line.rep.to.path.y.coor <- rep.coor[match(cur.rep, cur.uniq.rep)]
	line.rep.to.path.yend.coor <- pathway.coor[match(path.name, path.uniq.name)]
	path.line.col <- dot.rep.col[match(cur.rep, cur.uniq.rep)]

	#### rect for pathway
	if ('mean.diff' %in% colnames(ident.path.dat)){
		path.rect.length <- ident.path.dat[match(path.uniq.name, ident.path.dat$description),'mean.diff']
	}else{
		path.rect.length <- ident.path.dat[match(path.uniq.name, ident.path.dat$description),'median.diff']
	}
	path.rect.pval <-  ident.path.dat[match(path.uniq.name,ident.path.dat$description),'p.val.adj']
	path.rect.pval <- p.remove.inf(path.rect.pval)
	bar.path.col <- LRcolor(path.rect.pval, user.set.col=bar.pathway.col)

	### adjust the coordinate of each element
	max.limit <- max(n.up.ident, n.rep, top.n.path)
	# coors of three types of dot
	up.ident.new.coor <- seq(from=1, to=max.limit, length.out=n.up.ident)
	rep.new.coor <- seq(from=1, to=max.limit, length.out=n.rep)
	pathway.new.coor <- seq(from=1, to=max.limit, length.out=top.n.path)
	# coors of lines between up idents and receptors
	line.ident.to.rep.y.coor <- plyr::mapvalues(line.ident.to.rep.y.coor, from=1:n.up.ident, to=up.ident.new.coor, warn_missing=FALSE)
	line.ident.to.rep.yend.coor <- plyr::mapvalues(line.ident.to.rep.yend.coor, from=rep.coor, to=rep.new.coor, warn_missing=FALSE)
	# coors of lines between receptors and pathways
	line.rep.to.path.y.coor <- plyr::mapvalues(line.rep.to.path.y.coor, from=rep.coor, to=rep.new.coor, warn_missing=FALSE)
	line.rep.to.path.yend.coor <- plyr::mapvalues(line.rep.to.path.yend.coor, from=pathway.coor, to=pathway.new.coor, warn_missing=FALSE)
	
	pathplot.theme <- theme(axis.title=element_blank(), axis.text=element_blank(),axis.ticks=element_blank(), axis.line=element_blank(), panel.background=element_rect(fill="white"))
	

	### point
	point_up_ident <- geom_point(aes(x=1,y=up.ident.new.coor),size=8*dot.ident.size,color=dot.ident.col,shape=16) 
	point_receptor <- geom_point(aes(x=2,y=rep.new.coor),size=rep.size,color=dot.rep.col) 

	if (length(pathway.new.coor)==1){
		bar_pathway_width <- 0.3
	}else{
		bar_pathway_width <- 0.3 * (pathway.new.coor[2]-pathway.new.coor[1])
	}
	path.name.max.width <- max(strwidth(path.uniq.name, units="inches", cex=1))
	bar_pathway <- geom_rect(aes(xmin=3, xmax=3+bar.pathway.width*path.name.max.width*path.rect.length, ymin=pathway.new.coor-bar_pathway_width, ymax=pathway.new.coor+bar_pathway_width), size=1, fill=bar.path.col, alpha=1, show.legend=F)
	#print(bar.pathway.width*path.name.max.width*path.rect.length)
	
	### line
	line_ident_to_receptor <- geom_segment(aes(x=1, xend=2, y=line.ident.to.rep.y.coor,yend=line.ident.to.rep.yend.coor), size=line.ident.width, color=line.col, show.legend=F)
	line_receptor_to_pathway <- geom_segment(aes(x=2, xend=3, y=line.rep.to.path.y.coor,yend=line.rep.to.path.yend.coor),size=line.path.width, color=path.line.col, show.legend=F) 	
	### text and label
	label_up_ident <- annotate('text', hjust=0.5, x=1, y=up.ident.new.coor-0.4,label=up.uniq.ident, size=5*label.text.size)
	label_receptor <- annotate('text', hjust=0.5, x=2, y=rep.new.coor-0.4,label=cur.uniq.rep, size=5*label.text.size)
	text_pathway <- geom_text(aes(x=3, y=pathway.new.coor), hjust=0, label=path.uniq.name, size=5*label.text.size)
	label_title <- annotate('text', hjust=c(0.5,0.5,0), x=1:3, y=max(up.ident.new.coor)+0.8,label=c('Upstream', 'Receptor', 'Pathway annotation'), size=5*label.title.size) 
	
	plot <- ggplot() + 
	line_ident_to_receptor +
	line_receptor_to_pathway + # error!!!
	
	label_up_ident +
	label_receptor +
	label_title +

	point_up_ident +
	point_receptor + 
	bar_pathway +
	text_pathway +
	
	scale_x_continuous(limits=c(0,3+max(path.name.max.width))) + 
	pathplot.theme

	return(plot)


}

#' To present the interactions for a selected cluster, including both the upstream and downstream clusters which are connected by specific pathways in the selected cluster
#' @param Interact Interact list returned by findLRpairs
#' @param select.ident select.ident
#' @param ident.path.dat ident.path.dat for the selected ident
#' @param top.n.receptor top n receptor with the largest log2FCs to plot
#' @param order order of source clusters
#' @param top.n.path top n pathways with the smallest adjusted p values to plot
#' @param p.thre threshold for adjust p values; Only pathways with a adjust p valua < p.thre would be considered
#' @param top.n.ligand top n receptor with the largest log2FCs to plot
#' @param dot.ident.col color of the dots representing upstream clusters
#' @param dot.down.ident.col color of the dots representing downstream clusters
#' @param dot.receptor.col color of the dots representing reptors
#' @param bar.pathway.col color of the bars reprsenting pathways
#' @param bar.pathway.width width of the bars reprsenting pathways
#' @param dot.ident.size size of the dots representing source clusters
#' @param dot.receptor.size size of the dots representing receptors
#' @param label.text.size text size in the plot
#' @param label.title.size text size of the title annotation of the plot
#' @param line.ident.width text size of pathway labels
#' @param line.path.width width of lines connecting receptors and pathways
#' @param dot.pathway.size width of lines connecting receptors and pathways
#' @return Network plot showing receptors in the selected cluster, the upstream clusters which show L-R connections with the selected cluster, and the significant pathways involved in the receptors
#' @export
pathInterPlot <- function(Interact, select.ident, ident.path.dat, top.n.receptor=5, order=NULL, top.n.path=10, p.thre = 0.05, top.n.ligand=10, dot.ident.col=NULL, dot.down.ident.col=NULL, dot.receptor.col=NULL, bar.pathway.col=NULL, bar.pathway.width=10, dot.ident.size=1, dot.receptor.size=1, label.text.size=1, label.title.size=1, line.ident.width=1, line.path.width=1, dot.pathway.size=1){
	options(stringsAsFactors=F)
	Interact.num.dat <- subset(Interact$InteractNumer, LR.Number!=0)
	all.ident <- unique(c(Interact.num.dat$Cell.From,Interact.num.dat$Cell.To))
	if (!is.null(order)){
		all.ident <- all.ident[all.ident %in% order]
		all.ident <- factor(all.ident, levels=order)
		all.ident <- all.ident[order(all.ident)]
	}
	all.col <- scales::hue_pal(c=100)(length(all.ident))
	#all.col[1:4] <- c('red','black','green','blue')
	names(all.col) <- as.character(all.ident)

	### check the input and select significant pathways
	if ('t' %in% colnames(ident.path.dat)){
		ident.path.dat <- subset(ident.path.dat, p.val.adj < p.thre & t > 0)
	}else if ('W' %in% colnames(ident.path.dat)){
		ident.path.dat <- subset(ident.path.dat, p.val.adj < p.thre & median.diff > 0)
	}else{
		stop('please input the intact test result computed from diffPath')
	}
	ident.path.dat <- ident.path.dat[order(ident.path.dat$p.val.adj, decreasing=TRUE),]
	all.sig.path <- ident.path.dat$description

	### preprocess and extract useful information
	up.ident <- as.vector(unlist(sapply(ident.path.dat$cell.up, function(x){ strsplit(x, split=';') })))
	cur.rep <- as.vector(unlist(sapply(ident.path.dat$receptor.in.path, function(x){ strsplit(x, split=';') })))

	n.elem.each.row <- as.vector(unlist(lapply(sapply(ident.path.dat$receptor.in.path, function(x){ strsplit(x, split=';') }),length)))
	path.all <- rep(ident.path.dat$description, times=n.elem.each.row)
	plot.dat <- data.frame(up.ident=up.ident, cur.rep=cur.rep, path.name=path.all)
	
	### select the highly expressed receptors
	markerR.dat <- subset(Interact$markerR, subset=(cluster==select.ident & gene %in% cur.rep))
	markerR.dat <- markerR.dat[order(markerR.dat$avg_log2FC, decreasing=TRUE), ]
	if (top.n.receptor > nrow(markerR.dat)){
		warning(paste0('there is(are) ', nrow(markerR.dat),' marker receptor(s) in the selected ident, and the input top.n.receptor is ', top.n.receptor))
		top.n.receptor <- nrow(markerR.dat)
	}
	top.receptor.dat <- markerR.dat[1:top.n.receptor,]
	plot.dat <- subset(plot.dat, cur.rep %in% top.receptor.dat$gene)
	
	### to select the top n pathways associated with the selected receptors
	### those pathways with a larger position index in the all.sig.path are more significant
	### the following codes is to get the pathways with the top n largest position index
	path.uniq.name <- unique(plot.dat$path.name)
	path.position <- match(path.uniq.name, all.sig.path)
	path.position <- path.position[order(path.position, decreasing=TRUE)]
	if (top.n.path > length(path.position)){
		warning(paste0('there is(are) ', length(path.position),' significant pathway(s) for the selected ident, and the input top.n.path is ', top.n.path))
		top.n.path <- length(path.position)
	}
	path.position <- path.position[1:top.n.path]
	top.path <- all.sig.path[path.position]
	plot.dat$path.name[which(!(plot.dat$path.name %in% top.path))] <- NA
	
	### data for plot of upstream idents and receptors
	plot.ident.to.receprtor.dat <- unique(plot.dat[,c('up.ident','cur.rep')])
	up.ident <- plot.ident.to.receprtor.dat$up.ident
	cur.rep <- plot.ident.to.receprtor.dat$cur.rep
	### the idents which release ligands
	up.uniq.ident <- unique(up.ident)
	### to check the order parameter
	if (!is.null(order)){
		up.uniq.ident <- orderCheck(up.uniq.ident, order)
	}
	up.uniq.ident <- up.uniq.ident[order(up.uniq.ident, decreasing=TRUE)]
	n.up.ident <- length(up.uniq.ident)
	if (is.null(dot.ident.col)){
		dot.ident.col <- all.col[as.character(up.uniq.ident)]
	}else if(length(dot.ident.col) < n.up.ident){
		stop(paste0('there is(are) ',n.up.ident,' cluster(s), but only ',length(dot.ident.col),' colors provide.'))
	}else{
		dot.ident.col <- rev(dot.ident.col[1:n.up.ident])
	}
	up.ident.coor <- 1:n.up.ident

	### the receptor and their coordinate, and their color
	cur.uniq.rep <- factor(unique(cur.rep))
	cur.uniq.rep <- cur.uniq.rep[order(cur.uniq.rep, decreasing=TRUE)]
	n.rep <- length(cur.uniq.rep)
	rep.coor <- seq(from=1, to=n.up.ident, length.out=n.rep)

	logfc <- top.receptor.dat[match(cur.uniq.rep, top.receptor.dat$gene),'avg_log2FC']
	rep.size <- 5 * dot.receptor.size * logfc
	rep.pct <-  top.receptor.dat[match(cur.uniq.rep, top.receptor.dat$gene),'pct.1']
	rep.pval <-  top.receptor.dat[match(cur.uniq.rep, top.receptor.dat$gene),'p_val_adj']
	rep.pval <- p.remove.inf(rep.pval)
	dot.rep.col <- LRcolor(rep.pval, user.set.col=dot.receptor.col)

	### lines between ligand idents and receptors
	line.ident.to.rep.y.coor <- match(up.ident, up.uniq.ident)
	line.ident.to.rep.yend.coor <- rep.coor[match(cur.rep, cur.uniq.rep)]
	line.col <- dot.ident.col[match(up.ident, up.uniq.ident)]

	### pathways coordinate
	plot.ligand.to.path.dat <- unique(subset(plot.dat, subset=!(is.na(path.name)), select=c('cur.rep','path.name')))
	path.name <- plot.ligand.to.path.dat$path.name
	path.uniq.name <- unique(path.name)
	path.uniq.name <- factor(path.uniq.name, levels=top.path)
	path.uniq.name <- path.uniq.name[order(path.uniq.name, decreasing=TRUE)]
	pathway.coor <- seq(from=1, to=n.up.ident, length.out=length(path.uniq.name))
	
	### lines between receptors and pathways
	cur.rep <- plot.ligand.to.path.dat$cur.rep
	line.rep.to.path.y.coor <- rep.coor[match(cur.rep, cur.uniq.rep)]
	line.rep.to.path.yend.coor <- pathway.coor[match(path.name, path.uniq.name)]
	path.line.col <- dot.rep.col[match(cur.rep, cur.uniq.rep)]

	#### rect for pathway
	if ('mean.diff' %in% colnames(ident.path.dat)){
		path.rect.length <- ident.path.dat[match(path.uniq.name, ident.path.dat$description),'mean.diff']
	}else{
		path.rect.length <- ident.path.dat[match(path.uniq.name, ident.path.dat$description),'median.diff']
	}
	path.rect.pval <-  ident.path.dat[match(path.uniq.name, ident.path.dat$description),'p.val.adj']
	path.rect.pval <- p.remove.inf(path.rect.pval)
	bar.path.col <- LRcolor(path.rect.pval, user.set.col=bar.pathway.col)

	### data for plot of ligands and downstream clusters
	ligand.in.path <- ident.path.dat[match(path.uniq.name, ident.path.dat$description), 'ligand.in.path']
	if (all(is.na(ligand.in.path))){
		stop('there is no marker ligand in the same pathways with the selected receptors\nselect more receptors or more pathways and try again')
	}
	names(ligand.in.path) <- path.uniq.name
	ligand.in.path <- ligand.in.path[which(!is.na(ligand.in.path))]
	cur.lig <- as.vector(unlist(sapply(ligand.in.path, function(x) {strsplit(x, split=';')})))
	path.times <- sapply(ligand.in.path, function(x) {length(strsplit(x, split=';')[[1]])})
	path.for.lig <- rep(names(path.times), times=path.times)
	plot.pathway.to.ligand.dat <- data.frame(cur.lig=cur.lig, path.for.lig=path.for.lig)

	### select the top n ligand
	cur.uniq.lig <- unique(plot.pathway.to.ligand.dat$cur.lig)
	if (length(cur.uniq.lig) > top.n.ligand){
		markerL.dat <- subset(Interact$markerL, subset=(cluster==select.ident & gene %in% cur.uniq.lig))
		markerL.dat <- markerL.dat[match(cur.uniq.lig, markerL.dat$gene),]
		logfc <- markerL.dat$avg_log2FC
		names(logfc) <- cur.uniq.lig
		top.ligand.name <- names(logfc)[order(logfc, decreasing=TRUE)][1:top.n.ligand]
		plot.pathway.to.ligand.dat <- subset(plot.pathway.to.ligand.dat, cur.lig %in% top.ligand.name)
	}else if (length(cur.uniq.lig) < top.n.ligand){
		warning(paste0('there is(are) ', length(cur.uniq.lig),' marker ligand(s) in the selected ident associated with the selected pathways, and the input top.n.ligand is ', top.n.ligand))
	}

	path.for.lig <- plot.pathway.to.ligand.dat$path.for.lig
	cur.lig <- plot.pathway.to.ligand.dat$cur.lig
	cur.uniq.lig <- unique(cur.lig)
	n.cur.lig <- length(cur.uniq.lig)
	### lines between pathway and ligand
	line.path.to.lig.y.coor <- pathway.coor[match(path.for.lig, path.uniq.name)]
	dot.cur.lig.coor <- seq(from=1, to=n.up.ident, length.out=n.cur.lig)
	line.path.to.lig.yend.coor <- dot.cur.lig.coor[match(cur.lig, cur.uniq.lig)]
	
	### color and size for ligand
	markerL.dat <- subset(Interact$markerL, subset=(cluster==select.ident & gene %in% cur.uniq.lig))
	markerL.dat <- markerL.dat[match(cur.uniq.lig, markerL.dat$gene),]
	logfc <- markerL.dat$avg_log2FC
	lig.size <- 5 * dot.receptor.size * logfc
	lig.pct <-  markerL.dat$pct.1
	lig.pval <-  p.remove.inf(markerL.dat$p_val_adj)
	dot.ligand.col <- LRcolor(lig.pval, user.set.col=dot.receptor.col)
	line.path.to.lig.col <- dot.ligand.col[match(cur.lig, cur.uniq.lig)]

	### lines between ligand and downstream ident
	plot.down.lig.ident.dat <- do.call(rbind,lapply(cur.uniq.lig, function(x){findReceptor(Interact=Interact, select.ident=select.ident, select.ligand=x)}))
	if (nrow(plot.down.lig.ident.dat)==0){
		warning('there is no downstream cluster detected for the selected pathways and the corresponding ligands')
	}

	down.up.ident <- plot.down.lig.ident.dat$Cell.To
	down.uniq.ident <- unique(down.up.ident)
	### to check the order parameter
	if (!is.null(order)){
		down.uniq.ident <- orderCheck(down.uniq.ident, order)
	}
	down.uniq.ident <- down.uniq.ident[order(down.uniq.ident, decreasing=TRUE)]
	### color of downstream ident
	n.down.ident <- length(down.uniq.ident)
	if (is.null(dot.down.ident.col)){
		dot.down.ident.col <- all.col[as.character(down.uniq.ident)]
	}else if(length(dot.down.ident.col) < n.down.ident){
		stop(paste0('there is(are) ',n.down.ident,' cluster(s) in the down stream, but only ',length(dot.down.ident.col),' colors provide.'))
	}else{
		dot.down.ident.col <- rev(dot.down.ident.col[1:n.down.ident])
	}

	down.ident.coor <- seq(from=1, to=n.up.ident, length.out=n.down.ident)
	line.lig.to.ident.yend.coor <- down.ident.coor[match(down.up.ident, down.uniq.ident)]
	line.lig.to.ident.col <- dot.down.ident.col[match(down.up.ident, down.uniq.ident)]
	line.lig.to.ident.y.coor <- dot.cur.lig.coor[match(plot.down.lig.ident.dat$Ligand, cur.uniq.lig)]

	### adjust the coordinate of each element
	max.limit <- max(n.up.ident, n.rep, top.n.path, n.cur.lig, n.down.ident)
	# coors of five types of dot
	up.ident.new.coor <- seq(from=1, to=max.limit, length.out=n.up.ident)
	rep.new.coor <- seq(from=1, to=max.limit, length.out=n.rep)
	pathway.new.coor <- seq(from=1, to=max.limit, length.out=top.n.path)
	dot.cur.lig.new.coor <- seq(from=1, to=max.limit, length.out=n.cur.lig)
	down.ident.new.coor <- seq(from=1, to=max.limit, length.out=n.down.ident)
	# coors of lines between up idents and receptors
	line.ident.to.rep.y.coor <- plyr::mapvalues(line.ident.to.rep.y.coor, from=1:n.up.ident, to=up.ident.new.coor, warn_missing=FALSE)
	line.ident.to.rep.yend.coor <- plyr::mapvalues(line.ident.to.rep.yend.coor, from=rep.coor, to=rep.new.coor, warn_missing=FALSE)
	# coors of lines between receptors and pathways
	line.rep.to.path.y.coor <- plyr::mapvalues(line.rep.to.path.y.coor, from=rep.coor, to=rep.new.coor, warn_missing=FALSE)
	line.rep.to.path.yend.coor <- plyr::mapvalues(line.rep.to.path.yend.coor, from=pathway.coor, to=pathway.new.coor, warn_missing=FALSE)
	# coors of lines between pathways and ligands
	line.path.to.lig.y.coor <- plyr::mapvalues(line.path.to.lig.y.coor, from=pathway.coor, to=pathway.new.coor, warn_missing=FALSE)
	line.path.to.lig.yend.coor <- plyr::mapvalues(line.path.to.lig.yend.coor, from=dot.cur.lig.coor, to=dot.cur.lig.new.coor, warn_missing=FALSE)
	# coors of lines between ligands and idents
	line.lig.to.ident.y.coor <- plyr::mapvalues(line.lig.to.ident.y.coor, from=dot.cur.lig.coor, to=dot.cur.lig.new.coor, warn_missing=FALSE)
	line.lig.to.ident.yend.coor <- plyr::mapvalues(line.lig.to.ident.yend.coor, from=down.ident.coor, to=down.ident.new.coor, warn_missing=FALSE)

	pathplot.theme <- theme(axis.title=element_blank(), axis.text=element_blank(),axis.ticks=element_blank(), axis.line=element_blank(), panel.background=element_rect(fill="white"))
	### point
	point_up_ident <- geom_point(aes(x=1,y=up.ident.new.coor),size=8*dot.ident.size,color=dot.ident.col,shape=16) 
	point_receptor <- geom_point(aes(x=2,y=rep.new.coor),size=rep.size,color=dot.rep.col) 
	point_pathway <- geom_point(aes(x=3,y=pathway.new.coor),size=8*dot.ident.size,color='black',shape=21,fill='white')
	point_ligand <- geom_point(aes(x=4,y=dot.cur.lig.new.coor),size=lig.size,color=dot.ligand.col)
	point_down_ident <- geom_point(aes(x=5,y=down.ident.new.coor),size=8*dot.ident.size,color=dot.down.ident.col) 
	if (length(pathway.new.coor)==1){
		bar_pathway_width <- 0.3
	}else{
		bar_pathway_width <- 0.3 * (pathway.new.coor[2]-pathway.new.coor[1])
	}
	path.name.max.width <- max(strwidth(path.uniq.name, units="inches", cex=1))
	bar_pathway <- geom_rect(aes(xmin=6, xmax=6+bar.pathway.width*path.name.max.width*path.rect.length, ymin=pathway.new.coor-bar_pathway_width, ymax=pathway.new.coor+bar_pathway_width), size=1, fill=bar.path.col, alpha=1, show.legend=F)
	#print(bar.pathway.width*path.name.max.width*path.rect.length)
	
	### line
	line_ident_to_receptor <- geom_segment(aes(x=1, xend=2, y=line.ident.to.rep.y.coor,yend=line.ident.to.rep.yend.coor), size=line.ident.width, color=line.col, show.legend=F)
	line_receptor_to_pathway <- geom_segment(aes(x=2, xend=3, y=line.rep.to.path.y.coor,yend=line.rep.to.path.yend.coor),size=line.path.width, color=path.line.col, show.legend=F) 	
	line_pathway_to_ligand <- geom_segment(aes(x=3, xend=4, y=line.path.to.lig.y.coor,yend=line.path.to.lig.yend.coor), size=line.ident.width, color=line.path.to.lig.col, show.legend=F)
	line_ligand_to_ident <- geom_segment(aes(x=4, xend=5, y=line.lig.to.ident.y.coor,yend=line.lig.to.ident.yend.coor), size=line.ident.width, color=line.lig.to.ident.col, show.legend=F)
	
	### text and label
	label_up_ident <- annotate('text', hjust=0.5, x=1, y=up.ident.new.coor-0.4,label=up.uniq.ident, size=5*label.text.size)
	label_receptor <- annotate('text', hjust=0.5, x=2, y=rep.new.coor-0.4,label=cur.uniq.rep, size=5*label.text.size)
	#text_pathway <- geom_text(aes(x=6.1+path.rect.length, y=pathway.new.coor), hjust=0, label=path.uniq.name, size=5*label.text.size)
	text_pathway <- geom_text(aes(x=6, y=pathway.new.coor), hjust=0, label=path.uniq.name, size=5*label.text.size)
	label_pathway_circle <- annotate('text', hjust=0.5, x=3, y=pathway.new.coor, label=length(path.uniq.name):1, size=5*label.text.size)
	label_pathway_bar <- annotate('text', hjust=1, x=6, y=pathway.new.coor, label=length(path.uniq.name):1, size=5*label.text.size)
	label_ligand <- annotate('text', hjust=0.5, x=4, y=dot.cur.lig.new.coor-0.4,label=cur.uniq.lig, size=5*label.text.size)
	label_down_ident <- annotate('text', hjust=0.5, x=5, y=down.ident.new.coor-0.4,label=down.uniq.ident, size=5*label.text.size)
	label_title <- annotate('text', hjust=c(0.5,0.5,0.5,0.5,0.5,0), x=1:6, y=max(up.ident.new.coor)+0.8,label=c('Upstream', 'Receptor', 'Pathway', 'Ligand', 'Downstream','Pathway annotation'), size=5*label.title.size) 

	# the max limit of x axis is dependent on the longest pathway
	# this limit may need to be adjusted
	plot <- ggplot() + 

	line_pathway_to_ligand +
	line_ident_to_receptor +
	line_receptor_to_pathway + 
	line_ligand_to_ident +
	label_up_ident +
	label_receptor +
	label_ligand +
	label_down_ident +

	label_receptor +
	label_title +
	
	point_up_ident +
	point_receptor + 
	point_pathway +
	point_ligand +
	point_down_ident +
	bar_pathway +
	text_pathway +

	label_pathway_circle +
	label_pathway_bar +
	scale_x_continuous(limits=c(0,6+max(path.name.max.width))) + 

	pathplot.theme

	return(plot)
}

