#' To present a circos plot
#' @param object CommPath object
#' @param plot To present a circos plot for LR count ("count") or overall interaction intensity ("intensity") among cell clusters
#' @param ident.col Vector of colors of each cluster; names of the ident.col vector are supposed to be assigned to indicate each color for each cluster
#' @param filter Logical value indicating to present the circos plot for the filtered LR interactions or not; default is TRUE
#' @param select.ident To highlight the interaction between a specific cluster class and others; if 'NULL', plot interaction for all clusters
#' @param name.vert Should the group annotation be vertical to the grid? Defualt is FALSE
#' @importFrom circlize circos.clear circos.par chordDiagram circos.track circos.text
#' @return Circos plot showing the ligand-receptor interaction
#' @export
circosPlot <- function(object, plot='count', ident.col=NULL, filter=TRUE, select.ident=NULL, name.vert=FALSE){
	options(stringsAsFactors=F)
	Cluster <- object@cell.info$Cluster
	if (!is.factor(Cluster)){ Cluster <- factor(Cluster) }
	ident.label <- levels(Cluster)

	### select LR count or cluster intensity to plot
	if(filter){
		Interact.num.dat <- object@interact.filter$InteractNumber
		if(is.null(Interact.num.dat)){
			stop('No filtered LR interaction detected, set the parameter "filter" as FALSE and try again')
		}
	}else{
		Interact.num.dat <- object@interact$InteractNumber
	}

	if (plot=='count'){
		Interact.num.dat$Plot <- Interact.num.dat$LR.count
	}else if(plot=='intensity'){
		Interact.num.dat$Plot <- Interact.num.dat$Intensity
	}else{
		stop('Select "count" or "intensity" for the parameter "plot"!')
	}
	Interact.num.dat <- Interact.num.dat[, c('cell.from','cell.to','Plot')]
	
	Interact.num.dat = Interact.num.dat[Interact.num.dat$Plot!=0,]
	all.ident <- unique(c(Interact.num.dat$cell.from,Interact.num.dat$cell.to))
	all.ident <- orderCheck(all.ident, ident.label)
	all.ident <- all.ident[order(all.ident)]
	
	### to check the ident.col parameter
	if (is.null(ident.col)){ 
		ident.col <- scales::hue_pal(c=100)(length(all.ident))
		names(ident.col) <- all.ident
	}else{
		if (is.null(names(ident.col))){
			stop(paste0('Wrong "col" parameter! The colors should be named according to the clusters'))
		}else{
			colo.name <- names(ident.col)
			ident.missed <- all.ident[which(!(all.ident %in% colo.name))]
			ident.missed <- pasteIdent(ident.missed)
			stop(paste0('The color of ',ident.missed,' may be missed in the input color'))
		}
	}

	circos.clear()
	circos.par(start.degree=90, clock.wise=F)

	### plot interaction for all identity classes
	if (is.null(select.ident)){
		if (!name.vert){
			chordDiagram(Interact.num.dat, order=all.ident, grid.col=ident.col, annotationTrack=c("name","grid"), transparency=0.1, directional=1, direction.type='arrows', link.arr.type = "big.arrow",preAllocateTracks = list(track.height=mm_h(5)), annotationTrackHeight=convert_height(c(1, 2), "mm"))
		}else{
			chordDiagram(Interact.num.dat, order=all.ident, grid.col=ident.col, annotationTrack=c("grid"), transparency=0.1, directional=1, direction.type='arrows', link.arr.type = "big.arrow", preAllocateTracks = list(track.height=mm_h(5)), annotationTrackHeight=convert_height(c(2, 2), "mm"))
			circos.track(track.index=1, panel.fun=function(x, y){circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)}, bg.border = NA)
		}

	}else{
		### highlight the interaction for a specific identity class
		### check the select.ident parameter
		if (!(all(select.ident %in% all.ident)) | length(select.ident)>1){
			stop(paste0('Select one identity class from ', pasteIdent(all.ident)))
		}

		line.col <- rep('gray',nrow(Interact.num.dat))
		line.col[(Interact.num.dat$cell.from==select.ident)] <- ident.col[as.character(select.ident)]
		ident.lig = Interact.num.dat[Interact.num.dat$cell.to==select.ident,'cell.from']
		line.col[(Interact.num.dat$cell.to==select.ident)] <- ident.col[as.character(ident.lig)]
		
		link.zindex <- 1:nrow(Interact.num.dat)
		ident.line <- (Interact.num.dat$cell.from==select.ident) | (Interact.num.dat$cell.to==select.ident)
		link.zindex[!ident.line] <- rank(-(Interact.num.dat[!ident.line,'Plot']))
		link.zindex[ident.line] <- rank(-(Interact.num.dat[ident.line,'Plot']))+length(which(!ident.line))
		if (!name.vert){
			chordDiagram(Interact.num.dat, order=all.ident, grid.col=ident.col, annotationTrack=c("name","grid"), transparency=0.1, directional=1, direction.type='arrows', link.arr.type = "big.arrow", link.zindex=link.zindex, col=line.col, preAllocateTracks = list(track.height=mm_h(5)), annotationTrackHeight=convert_height(c(1, 2), "mm"))
		}else{
			chordDiagram(Interact.num.dat, order=all.ident, grid.col=ident.col, annotationTrack="grid", transparency=0.1, directional=1, direction.type='arrows', link.arr.type = "big.arrow", link.zindex=link.zindex, col=line.col, preAllocateTracks = list(track.height=mm_h(5)), annotationTrackHeight=convert_height(c(2, 2), "mm"))
			circos.track(track.index=1, panel.fun=function(x, y){circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)}, bg.border = NA)
		}
	}
}

#' To present a dot plot for specific ligand-receptor pairs in specific clusters
#' @param object CommPath object
#' @param ligand.ident Vector containing the ligand ident
#' @param receptor.ident Vector containing the receptor ident
#' @param ident.levels Vector of levels of the identities 
#' @param top.n.inter Show the dotplot for the top n LR pairs with the largest product of Log2FC
#' @param filter Logical value indicating to present the circos plot for the filtered LR interacyion or not; default is TRUE
#' @param return.data Logical value indicating whether to return the data for the plot or not
#' @importFrom ggplot2 ggplot geom_point scale_color_manual labs theme element_text element_rect element_line aes scale_color_gradientn
#' @return Dotplot showing the ligand-receptor interaction between the selected ligand.ident and receptor.ident
#' @export
dotPlot.LR <- function(object, ligand.ident=NULL, receptor.ident=NULL, ident.levels=NULL, top.n.inter=10, filter=TRUE, return.data=FALSE){
	options(stringsAsFactors=F)
	if (is.null(ligand.ident) & is.null(receptor.ident)){
		stop("Either ligand.ident or ligand.ident need to be asigned")
	}
	if (length(ligand.ident)>1 & length(receptor.ident)>1){
		stop("Specify one cluster for ligand or receptor analysis")
	}

	# get the InteractGene dataframe
	if(filter){
		inter.gene.dat <- object@interact.filter$InteractGene
		if(is.null(inter.gene.dat)){
			stop('No filtered LR interaction detected, set the parameter "filter" as FALSE and try again')
		}
	}else{
		inter.gene.dat <- object@interact$InteractGene		
	}
	
	if (!is.null(receptor.ident)){
		which.rep.not.in <- which(!(receptor.ident %in% inter.gene.dat$cell.to))
		which.rep.in <- which(receptor.ident %in% inter.gene.dat$cell.to)

		if (length(which.rep.not.in) > 0){
			if (length(which.rep.in) > 0){
				warning(paste0('No ligand and upstream cluster detected for cluster ', pasteIdent(unique(receptor.ident[which.rep.not.in]))))
			}else{
				stop(paste0('No ligand and upstream cluster detected for cluster ', pasteIdent(unique(receptor.ident[which.rep.not.in]))))
			}
		}
		receptor.ident <- receptor.ident[which.rep.in]
		inter.gene.dat <- inter.gene.dat[inter.gene.dat$cell.to %in% receptor.ident,]
	}
	if (!is.null(ligand.ident)){
		which.lig.not.in <- which(!(ligand.ident %in% inter.gene.dat$cell.from))
		which.lig.in <- which(ligand.ident %in% inter.gene.dat$cell.from)

		if (length(which.lig.not.in) > 0){
			if (length(which.lig.in) > 0){
				warning(paste0('No ligand and upstream cluster detected for cluster ', pasteIdent(unique(ligand.ident[which.lig.not.in]))))
			}else{
				stop(paste0('No ligand and upstream cluster detected for cluster ', pasteIdent(unique(ligand.ident[which.lig.not.in]))))
			}
		}
		ligand.ident <- ligand.ident[which.lig.in]
		inter.gene.dat <- inter.gene.dat[which(inter.gene.dat$cell.from %in% ligand.ident),]
	}

	inter.gene.dat$LR.info <- paste(inter.gene.dat$ligand, inter.gene.dat$receptor, sep='/')

	### to find those LR pairs with largest log2FC.LR
	max.fc <- by(data=inter.gene.dat$log2FC.LR, INDICES=inter.gene.dat$LR.info, function(x){max(x,na.rm=T)})
	max.LRname.fc <- names(max.fc)[order(max.fc, decreasing=TRUE)]
	
	if (top.n.inter > length(max.LRname.fc)){
		warning(paste0('There is(are) ', length(max.LRname.fc),' significant LR pair(s) for the selected ident, and the input top.n.inter is ', top.n.inter))
		top.n.inter <- length(max.LRname.fc)
	}

	max.LRname.fc <- max.LRname.fc[1:top.n.inter]
	inter.gene.dat <- subset(inter.gene.dat, LR.info %in% max.LRname.fc)

	if (return.data){
		return(inter.gene.dat[,c('cell.from', 'cell.to', 'LR.info', 'log2FC.LR', 'P.val.LR', 'P.val.adj.LR')])
	}
	
	inter.gene.dat$Log10.P.adj <- -log10(inter.gene.dat$P.val.adj.LR)
	inter.gene.dat$Log10.P.adj <- p.remove.inf(inter.gene.dat$Log10.P.adj)
	if(all(inter.gene.dat$Log10.P.adj==0)){inter.gene.dat$Log10.P.adj='Infinite'}

	#inter.gene.dat[which(inter.gene.dat$Log10.P.adj > 30), 'Log10.P.adj'] <- 30
	
	### adjust the x label and title
	if (length(ligand.ident)==1){
		inter.gene.dat$Xaxis <- inter.gene.dat$cell.to
		x.title <- paste0('Downstream clusters for cluster ',ligand.ident)
	}else{
		inter.gene.dat$Xaxis <- inter.gene.dat$cell.from
		x.title <- paste0('Upstream clusters for cluster ',receptor.ident)
	}
	### you may want to adjust the order of the clusters in the x axis
	inter.gene.dat$Xaxis <- as.character(inter.gene.dat$Xaxis)
	if (is.null(ident.levels)){
		all.ident <- object@cell.info$Cluster
		if (!is.factor(all.ident)){ all.ident <- factor(all.ident) }
		ident.levels <- levels(all.ident)
	}
	x.levels <- ident.levels[ident.levels %in% inter.gene.dat$Xaxis]
	if (all(inter.gene.dat$Xaxis %in% x.levels)){
		inter.gene.dat$Xaxis <- factor(inter.gene.dat$Xaxis, levels=x.levels)
	}else{
		ident.missed <- unique(inter.gene.dat$Xaxis[!(inter.gene.dat$Xaxis %in% x.levels)])
		ident.missed <- pasteIdent(ident.missed)
		stop(paste0('The ident class ',ident.missed,' may be missed in the input ident.levels'))
	}

	if (all(inter.gene.dat$Log10.P.adj=='Infinite')){
		plot <- ggplot(data=inter.gene.dat,aes(x=Xaxis,y=LR.info)) + 
			geom_point(aes(size=log2FC.LR,col='Infinite')) +
			scale_color_manual(values='red') + 
			labs(color='-log10(p.adj)',size='log2FC.LR',x=x.title,y="") + 
			theme(axis.text=element_text(colour='black',size=12),
			axis.title=element_text(colour='black',size=12),
			panel.background=element_rect(fill="white",color="black"),
			panel.grid=element_line(size=0.5,colour='gray')
			)
		warning('Adjusted p values for all LR pairs are 0')
	}else{
		plot <- ggplot(data=inter.gene.dat,aes(x=Xaxis,y=LR.info)) + 
			geom_point(aes(size=log2FC.LR,col=Log10.P.adj)) +
			scale_color_gradientn(values = seq(0,1,0.2),colours=c("#0570b0", "grey", "#d7301f")) + 
			labs(color='-log10(p.adj)',size='log2FC.LR',x=x.title,y="") + 
			theme(axis.text=element_text(colour='black',size=12),
			axis.title=element_text(colour='black',size=12),
			panel.background=element_rect(fill="white",color="black"),
			panel.grid=element_line(size=0.5,colour='gray')
			)
	}
	return(plot)
}

#' To present a dot plot for top ligand-receptor pairs involved in the specific pathways in the selected clusters
#' @param object CommPath object
#' @param acti.path.filtered.dat Data frame of differential activation test result of filtered from diffAllPath
#' @param pathway To present a dot plot for ligand-receptor pairs involved in which pathway
#' @param ligand.ident Vector containing the ligand ident
#' @param receptor.ident Vector containing the receptor ident
#' @param ident.levels Vector of levels of the identities 
#' @param top.n.inter Show the dotplot for the top n LR pairs with the largest product of Log2FC
#' @param return.data Logical value indicating whether to return the data for the plot or not
#' @importFrom ggplot2 ggplot geom_point scale_color_manual labs theme element_text element_rect element_line aes scale_color_gradientn
#' @return Dotplot showing the ligand-receptor interaction involved in the specific pathways in the selected clusters
#' @export
dotPlot.pathway <- function(object, acti.path.filtered.dat, pathway, ligand.ident=NULL, receptor.ident=NULL, ident.levels=NULL, top.n.inter=10, return.data=FALSE){
	options(stringsAsFactors=F)
	if (is.null(ligand.ident) & is.null(receptor.ident)){
		stop("Either ligand.ident or ligand.ident need to be asigned")
	}
	if (length(ligand.ident)>1 & length(receptor.ident)>1){
		stop("Specify one cluster for ligand or receptor analysis")
	}

	if('t' %in% colnames(acti.path.filtered.dat)){
		acti.path.filtered.dat <- acti.path.filtered.dat[which(acti.path.filtered.dat$mean.diff > 0 & acti.path.filtered.dat$P.val.adj < 0.05), c('description','ligand.in.path','receptor.in.path','cluster')]
	}else{
		acti.path.filtered.dat <- acti.path.filtered.dat[which(acti.path.filtered.dat$median.diff > 0 & acti.path.filtered.dat$P.val.adj < 0.05), c('description','ligand.in.path','receptor.in.path','cluster')]
	}

	if (!is.null(receptor.ident)){
		cur.path.dat <- acti.path.filtered.dat[which(acti.path.filtered.dat$cluster==receptor.ident), ]
		if(nrow(cur.path.dat)==0){
			stop(paste0('There is no significantly up-regulatged ligand-containing pathways for cluster ', receptor.ident))
		}

		cur.rep.char <- unique(cur.path.dat[which(cur.path.dat$description==pathway), 'receptor.in.path'])
		if(length(cur.rep.char)==0){
			stop(paste0('The selected pathway is not up-regulatged in the cluster', receptor.ident))
		}else if(length(cur.rep.char)==1 & is.na(cur.rep.char)){
			stop('The selected pathway does not contain any marker receptor')
		}

		cur.rep.vec <- strsplit(cur.rep.char, split=';')[[1]]
		up.ligand.dat <- findLigand(object, select.ident=receptor.ident, select.receptor=cur.rep.vec, filter=TRUE)
		up.ligand.dat$Xaxis <- up.ligand.dat$cell.from
		up.ligand.dat$Yaxis <- paste(up.ligand.dat$ligand, up.ligand.dat$receptor, sep='/')
		plot.lr.dat <- up.ligand.dat
		x.title <- paste0('Upstream clusters for cluster ',receptor.ident)

	}else{
		cur.path.dat <- acti.path.filtered.dat[which(acti.path.filtered.dat$cluster==ligand.ident), ]
		if(nrow(cur.path.dat)==0){
			stop(paste0('There is no significantly up-regulatged ligand-containing pathways for cluster ', ligand.ident))
		}

		lig.rep.char <- unique(cur.path.dat[which(cur.path.dat$description==pathway), 'ligand.in.path'])
		if(length(lig.rep.char)==0){
			stop(paste0('The selected pathway is not up-regulatged in the cluster', ligand.ident))
		}else if(length(lig.rep.char)==1 & is.na(lig.rep.char)){
			stop('The selected pathway does not contain any marker ligand')
		}
		
		cur.lig.vec <- strsplit(lig.rep.char, split=';')[[1]]
		down.ligand.dat <- findReceptor(object, select.ident=ligand.ident, select.ligand=cur.lig.vec, filter=TRUE)
		down.ligand.dat$Xaxis <- down.ligand.dat$cell.to
		down.ligand.dat$Yaxis <- paste(down.ligand.dat$ligand, down.ligand.dat$receptor, sep='/')
		plot.lr.dat <- down.ligand.dat
		x.title <- paste0('Downstream clusters for cluster ',ligand.ident)
	}

	if(top.n.inter <= length(unique(plot.lr.dat$Yaxis))){
		maxfc.lr <- by(data=plot.lr.dat$log2FC.LR, INDICES=plot.lr.dat$Yaxis, FUN=function(x){max(x)})
		maxfc.lr.name <- names(maxfc.lr[order(maxfc.lr, decreasing=TRUE)][1:top.n.inter])
		plot.lr.dat <- plot.lr.dat[which(plot.lr.dat$Yaxis %in% maxfc.lr.name), ]
	}else{
		warning(paste0('There is(are) ', length(unique(plot.lr.dat$Yaxis)),' significant LR pair(s) involved in the selected pathway for the selected ident, and the input top.n.inter is ', top.n.inter))
	}

	if (return.data){
		return.dat <- plot.lr.dat[,c('cell.from', 'cell.to', 'Yaxis', 'log2FC.LR', 'P.val.LR', 'P.val.adj.LR')]
		colnames(return.dat) <- c('cell.from', 'cell.to', 'LR.info', 'log2FC.LR', 'P.val.LR', 'P.val.adj.LR')
		return(return.dat)
	}

	plot.lr.dat$Log10.P.adj <- -log10(plot.lr.dat$P.val.adj.LR)
	plot.lr.dat$Log10.P.adj <- p.remove.inf(plot.lr.dat$Log10.P.adj)
	if(all(plot.lr.dat$Log10.P.adj==0)){ plot.lr.dat$Log10.P.adj='Infinite' }

	### you may want to adjust the order of the clusters in the x axis
	plot.lr.dat$Xaxis <- as.character(plot.lr.dat$Xaxis)
	if (is.null(ident.levels)){
		all.ident <- object@cell.info$Cluster
		if (!is.factor(all.ident)){ all.ident <- factor(all.ident) }
		ident.levels <- levels(all.ident)
	}
	x.levels <- ident.levels[ident.levels %in% plot.lr.dat$Xaxis]
	if (all(plot.lr.dat$Xaxis %in% x.levels)){
		plot.lr.dat$Xaxis <- factor(plot.lr.dat$Xaxis, levels=x.levels)
	}else{
		ident.missed <- unique(plot.lr.dat$Xaxis[!(plot.lr.dat$Xaxis %in% x.levels)])
		ident.missed <- pasteIdent(ident.missed)
		stop(paste0('The ident class ',ident.missed,' may be missed in the input ident.levels'))
	}

	if (all(plot.lr.dat$Log10.P.adj=='Infinite')){
		plot <- ggplot(data=plot.lr.dat,aes(x=Xaxis,y=Yaxis)) + 
			geom_point(aes(size=log2FC.LR,col='Infinite')) +
			scale_color_manual(values='red') + 
			labs(color='-log10(p.adj)',size='log2FC.LR',x=x.title,y="") + 
			theme(axis.text=element_text(colour='black',size=12),
			axis.title=element_text(colour='black',size=12),
			panel.background=element_rect(fill="white",color="black"),
			panel.grid=element_line(size=0.5,colour='gray')
			)
		warning('Adjusted p values for all LR pairs are 0')
	}else{
		plot <- ggplot(data=plot.lr.dat,aes(x=Xaxis,y=Yaxis)) + 
			geom_point(aes(size=log2FC.LR,col=Log10.P.adj)) +
			scale_color_gradientn(values = seq(0,1,0.2),colours=c("#0570b0", "grey", "#d7301f")) + 
			labs(color='-log10(p.adj)',size='log2FC.LR',x=x.title,y="") + 
			theme(axis.text=element_text(colour='black',size=12),
			axis.title=element_text(colour='black',size=12),
			panel.background=element_rect(fill="white",color="black"),
			panel.grid=element_line(size=0.5,colour='gray')
			)
	}
	return(plot)
}


#' To plot a heatmap of those differentially enriched pathways for each cluster
#' @param object CommPath object
#' @param acti.path.dat Data frame of differential enrichment test result from diffAllPath; if NULL, diffAllPath would be run to get the acti.path.dat
#' @param top.n.pathway Show the heatmap of top n most significant pathways
#' @param path.order Sort criteria used to select the top n pathways, either 'P.val' or 'P.val.adj', which represent the original and adjusted p values, or 'diff' which represents the mean (in t test) or median (in wilcox test) difference
#' @param col Vector of colors used to generate a series of gradient colors to show the enrichment score of pathways; provide a vector containing at least two colors
#' @param cell.aver Whether to display averaged pathway enrichment scores among cells from the same clusters or to display scores for all cells; default is FALSE, which means to display scores for all cells
#' @param cell.label.size Text size of the label of cell types
#' @param cell.label.angle Rotation angle of the label of cell types
#' @param pathway.label.size Text size of the label of pathways
#' @param scale Whether to scale the enrichment sorces matrix among cells or not; default is TRUE
#' @param truncation Truncation fold value; scores > (the third quartiles + truncation * interquartile range) and scores < (the first quartiles - truncation * interquartile range) will be adjusted; either a value to indicate the specific truncation value or 'none' to indicate no truncation; default is 1.5
#' @param show.legend Whether to show the legend or not; default is FALSE
#' @importFrom ggplot2 geom_tile facet_grid element_blank
#' @importFrom stats quantile
#' @return Heatmap plot showing the top enriched patways in each cluster
#' @export
pathHeatmap <- function(object, acti.path.dat=NULL, top.n.pathway=10, path.order='P.val.adj', col=NULL, cell.aver=FALSE, cell.label.size=NULL, cell.label.angle=45, pathway.label.size=NULL, scale=TRUE, truncation=1.5, show.legend=TRUE){
	if (is.null(acti.path.dat)){ acti.path.dat <- diffAllPath(object) }
	
	acti.path.dat <- subset(acti.path.dat, P.val.adj < 0.05)
	### select the highly enriched pathways
	if ('mean.diff' %in% colnames(acti.path.dat)){
		acti.path.dat <- subset(acti.path.dat, mean.diff > 0)
		acti.path.dat$stat <- acti.path.dat$t
		acti.path.dat$diff <- acti.path.dat$mean.diff
	}else if('median.diff' %in% colnames(acti.path.dat)){
		acti.path.dat <- subset(acti.path.dat, median.diff > 0)
		acti.path.dat$stat <- acti.path.dat$W
		acti.path.dat$diff <- acti.path.dat$median.diff
	}else{
		stop('Please input the intact test result computed from diffAllPath')
	}

	### select top n pathways
	if (path.order=='P.val.adj' | path.order=='P.val'){
		acti.path.dat <- acti.path.dat[order(acti.path.dat[[path.order]], -acti.path.dat[['stat']]), ]
	}else if (path.order=='diff'){
		acti.path.dat <- acti.path.dat[order(-acti.path.dat[['diff']], acti.path.dat[['P.val.adj']],), ]
	}else{
		stop('Select one order criteria from P.val.adj, P.val, and diff')
	}

	ident.label <- object@cell.info$Cluster
	gsva.mat <- object@pathway$acti.score

	if (!is.factor(ident.label)){
		ident.label <- factor(ident.label)
	}
	all.ident <- levels(ident.label)

	path <- sapply(all.ident, function(each.label){
		ident.path.dat <- subset(acti.path.dat, cluster==each.label)
		if (nrow(ident.path.dat)==0){ 
			warning(paste0('There is no significantly up-regulatged pathway for cluster ',each.label))
			return(FALSE) 
		}
		if (top.n.pathway > nrow(ident.path.dat)){
			warning(paste0('There is(are) ',nrow(ident.path.dat),' significantly up-regulatged pathway(s) for cluster ',each.label,', and the input top.n.pathway is ',top.n.pathway))
			top.n.pathway <- nrow(ident.path.dat)
		}
		return(ident.path.dat[1:top.n.pathway, 'description'])
	})
	path <- unique(as.character(unlist(path)))
	path <- path[which(path!='FALSE')]
	
	if (length(path)==0){
		stop('There is no significantly up-regulatged pathway(s) for all clusters!')
	}
	
	gsva.mat <- gsva.mat[path,]
	
	if (cell.aver){
		gsva.cell.mat <- t(gsva.mat)
		path.aver <- by(data=gsva.cell.mat, INDICES=ident.label, FUN=colMeans)
		path.mat <- matrix(unlist(path.aver), nrow=ncol(gsva.cell.mat))
		rownames(path.mat) <- colnames(gsva.cell.mat)
		colnames(path.mat) <- names(path.aver)

		### scale or not
		if (scale){
			gsva.scale.mat <- t(scale(t(path.mat)))
		}else{
			gsva.scale.mat <- path.mat
		}

		gsva.dat <- melt(gsva.scale.mat, varnames=c('Pathway','Cell'),value.name="Score",  na.rm=TRUE)
		gsva.dat$Pathway <- factor(as.character(gsva.dat$Pathway), levels=rev(path))
		### Cell is the Celltype
		gsva.dat$Cell <- factor(as.character(gsva.dat$Cell), levels=levels(ident.label))
		
		if (is.numeric(truncation)){
			IQR <- quantile(gsva.dat$Score, 0.75) - quantile(gsva.dat$Score, 0.25)
			max.score <- quantile(gsva.dat$Score, 0.75) + IQR * truncation
			min.score <- quantile(gsva.dat$Score, 0.25) - IQR * truncation
			gsva.dat[gsva.dat$Score > max.score, 'Score'] <- max.score
			gsva.dat[gsva.dat$Score < min.score, 'Score'] <- min.score
		}else if(truncation!='none'){
			stop('Either a specific value or "none" is needed for the truncation parameter')
		}
		
		### basic plot
		path.heatmap <- ggplot(data=gsva.dat, aes(x=Cell, y=Pathway, fill=Score)) + 
			geom_tile() +
			labs(x='',y='') +
			theme(axis.text.y=element_text(colour='black'),
				axis.text.x=element_text(colour='black'),
				axis.ticks=element_blank(),
				axis.line=element_blank(),
				panel.background=element_rect(fill="white"))
		path.heatmap <- path.heatmap + theme(axis.text.x=rotated.axis.element.text(cell.label.angle,'x'))

		### gradient color
		if (is.null(col)){
			colours <- c("#440154" ,"#21908C", "#FDE725")
		}else if(length(col)==1){
			stop('Select at least 2 colors to generate a series of gradient colors')
		}else{
			colours <- col
		}
		path.heatmap <- path.heatmap + scale_fill_gradientn(colours=colours)

		### show legend or not
		if (!show.legend){
			path.heatmap <- path.heatmap + theme(legend.position='none')
		}

		### cell and pathway label
		if (!is.null(cell.label.size)){
			path.heatmap <- path.heatmap + theme(axis.text.x=element_text(size=cell.label.size))
		}
		if (!is.null(pathway.label.size)){
			path.heatmap <- path.heatmap + theme(axis.text.y=element_text(size=pathway.label.size))
		}
		return(path.heatmap)

	}else{
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
			stop('Either a specific value or "none" is needed for the truncation parameter')
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
				strip.text.x=element_text(colour="black", angle=cell.label.angle),
				panel.background=element_rect(fill="white"))

		### gradient color
		if (is.null(col)){
			colours <- c("#440154" ,"#21908C", "#FDE725")
		}else if(length(col)==1){
			stop('Select at least 2 colors to generate a series of gradient colors')
		}else{
			colours <- col
		}
		path.heatmap <- path.heatmap + scale_fill_gradientn(colours=colours)

		### show legend or not
		if (!show.legend){
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
	
}

#' To present the interactions for a selected cluster, including the upstream clusters and the activated pathways in the selected cluster
#' @param object CommPath object
#' @param select.ident Plot the activated pathways for which cluster or cell type?
#' @param up.ident Upstream clusters in the communication chain; default is all clusters
#' @param acti.path.dat Data frame of differential activation test result from diffAllPath
#' @param top.n.path Top n pathways with the smallest adjusted p values to plot
#' @param path.order Sort criteria used to select the top n pathways, either 'P.val' or 'P.val.adj', which represent the original and adjusted p values, or 'diff' which represents the mean (in t test) or median (in wilcox test) difference
#' @param select.path Vector of pathways for which users want to present a chain plot; default is all, which means to plot all pathways. The parameter top.n.path would be masked if this one is set with specific pathways.
#' @param p.thre Threshold for adjust p values; Only pathways with a adjust p valua < p.thre would be considered
#' @param top.n.receptor Top n receptor with the largest log2FCs to plot
#' @param dot.ident.col Color of the dots representing clusters
#' @param dot.ident.size Size of the dots representing clusters
#' @param dot.gene.col Color of the dots representing receptors
#' @param dot.gene.size Size of the dots representing receptors
#' @param bar.pathway.col Color of the bars reprsenting pathways
#' @param label.text.size Text size in the plot
#' @param label.title.size Text size of the title annotation of the plot 
#' @importFrom ggplot2 geom_rect geom_segment geom_text annotate scale_x_continuous
#' @importFrom graphics strwidth
#' @return Network plot showing the significantly activated pathways in the select cluster compared to all other clusters, receptors involved in the pathways, and the upstream clusters which show LR connections with the selected cluster
#' @export
pathPlot <- function(object, select.ident, up.ident=NULL, acti.path.dat=NULL, top.n.path=5, path.order='P.val.adj', select.path='all', p.thre=0.05, top.n.receptor=10, dot.ident.col=NULL, dot.ident.size=1, dot.gene.col=NULL, dot.gene.size=1, bar.pathway.col=NULL, label.text.size=1, label.title.size=1){
	options(stringsAsFactors=F)
	all.ident <- object@cell.info$Cluster
	if (!is.factor(all.ident)){ all.ident <- factor(all.ident) }
	ident.label <- levels(all.ident)
	
	# check the input up.ident
	if (!is.null(up.ident)){
		up.in.ident <- up.ident[which(up.ident %in% ident.label)]
		if (length(up.in.ident)==0){
			stop('Wrong up.ident parameter! All selected up.idents are not present in the dataset!')
		}
		up.missed.ident <- up.ident[which(!(up.ident %in% ident.label))]
		if (length(up.missed.ident) > 0){ warning(paste0('These selected up.idents are not present in the dataset: ', pasteIdent(up.missed.ident), '.')) }
	}else{
		up.in.ident <- ident.label
	}

	### check the input and select significant pathways
	if (is.null(acti.path.dat)) {
		ident.path.dat <- diffPath(object, select.ident.1=select.ident)
	}else{
		ident.path.dat <- subset(acti.path.dat, cluster==select.ident)
	}
	if ('t' %in% colnames(ident.path.dat)){
		ident.path.dat <- subset(ident.path.dat, P.val.adj < p.thre & t > 0)
		ident.path.dat$stat <- ident.path.dat$t
		ident.path.dat$diff <- ident.path.dat$mean.diff
	}else if ('W' %in% colnames(ident.path.dat)){
		ident.path.dat <- subset(ident.path.dat, P.val.adj < p.thre & median.diff > 0)
		ident.path.dat$stat <- ident.path.dat$W
		ident.path.dat$diff <- ident.path.dat$median.diff
	}else{
		stop('Please input the integrate ident.path.dat computed from diffPath')
	}

	if (nrow(ident.path.dat)==0){ stop(paste0('There is no significantly up-regulatged pathways for cluster ', select.ident )) }

	### remove pathways with NA receptor
	ident.path.dat <- subset(ident.path.dat, !is.na(receptor.in.path))
	if (nrow(ident.path.dat)==0){ stop(paste0('There is no marker receptor in the significantly up-regulatged pathways for cluster ', select.ident )) }

	if (path.order=='P.val.adj' | path.order=='P.val'){
		ident.path.dat <- ident.path.dat[order(ident.path.dat[,path.order], -ident.path.dat[,'stat']),]
	}else if(path.order=='diff'){
		ident.path.dat <- ident.path.dat[order(-ident.path.dat[,'diff'], ident.path.dat[,'P.val.adj']),]
	}

	### select user-defined pathways, with the select.path parameter
	if (length(select.path) > 1 | (length(select.path)==1 & select.path[1]!='all')){
		path.num.plot <- which(ident.path.dat$description %in% select.path)
		if (length(path.num.plot) > 0){
			ident.path.dat <- ident.path.dat[path.num.plot, ]
		}else{
			stop('All pathways selected are not up-regulatged in the selected cluster. Run diffPath to retrieve all up-regulatged pathways.')
		}
		top.n.path <- length(path.num.plot)
		path.missed <- select.path[which(!(select.path %in% ident.path.dat$description))]
		if (length(path.missed) > 0){
			warning(paste0('These pathways are not up-regulatged in the selected cluster: ', pasteIdent(path.missed), '. Run diffPath to retrieve all up-regulatged pathways.'))
		}
	}else{
		### select top n pathways, if the select.path parameter is not set
		if (top.n.path > nrow(ident.path.dat)){
			warning(paste0('There is(are) ', nrow(ident.path.dat),' significant pathway(s) for the selected ident, and the input top.n.path is ', top.n.path))
			top.n.path <- nrow(ident.path.dat)
		}
		ident.path.dat <- ident.path.dat[1:top.n.path, ]
	}
	all.sig.path <- as.vector(ident.path.dat$description)

	### data for plot of receptors and upstream clusters
	receptor.in.path <- as.character(ident.path.dat[, 'receptor.in.path'])
	names(receptor.in.path) <- all.sig.path
	cur.rep <- as.vector(unlist(sapply(receptor.in.path, function(x) {strsplit(x, split=';')})))
	path.times <- sapply(receptor.in.path, function(x) {length(strsplit(x, split=';')[[1]])})
	cur.path <- rep(names(path.times), times=path.times)
	plot.receptor.to.pathway.dat <- data.frame(cur.rep=cur.rep, cur.path=cur.path)

	### select the top n receptors
	cur.uniq.rep <- as.vector(unique(plot.receptor.to.pathway.dat$cur.rep))
	cur.rep.LR.inten <- cluster.lr.inten(cur.uniq.rep, object, select.ident, ident.label, find='ligand')
	cur.rep.LR.inten <- subsetMatrix(cur.rep.LR.inten, row.select=up.in.ident)
	keep.rep <- names(which(colSums(is.na(cur.rep.LR.inten)) < nrow(cur.rep.LR.inten)))
	cur.rep.LR.inten <- subsetMatrix(cur.rep.LR.inten, col.select=keep.rep)

	cur.rep.max.LR.inten <- apply(cur.rep.LR.inten, 2, function(x){max(x, na.rm=TRUE)})
	cur.rep.max.LR.inten <- cur.rep.max.LR.inten[order(cur.rep.max.LR.inten, decreasing=TRUE)]
	
	if (length(cur.rep.max.LR.inten) < top.n.receptor){
		warning(paste0('There is(are) ', length(cur.rep.max.LR.inten),' marker receptors(s) in the selected ident associated with the selected pathways, and the input top.n.receptor is ', top.n.receptor))
		top.n.receptor <- length(cur.rep.max.LR.inten)
	}
	
	top.rep.name <- names(cur.rep.max.LR.inten)[1:top.n.receptor]
	plot.receptor.to.pathway.dat <- subset(plot.receptor.to.pathway.dat, cur.rep %in% top.rep.name)
	top.rep.LR.inten <- subsetMatrix(cur.rep.LR.inten, col.select=top.rep.name)
	plot.ident.to.receprtor.dat <- melt(top.rep.LR.inten, varnames=c('cell.from','receptor'),value.name="LR.inten",  na.rm=TRUE)
	plot.ident.to.receprtor.dat <- factor.to.character(plot.ident.to.receprtor.dat)

	### pathways and their color
	path.uniq.name <- factor(all.sig.path, levels=all.sig.path)
	path.uniq.name <- path.uniq.name[order(path.uniq.name, decreasing=TRUE)]
	
	#### rect for pathway
	path.rect.length <- ident.path.dat[match(path.uniq.name, ident.path.dat$description),'diff']
	path.rect.stat <- ident.path.dat[match(path.uniq.name, ident.path.dat$description),'stat']
	path.rect.pval <-  ident.path.dat[match(path.uniq.name, ident.path.dat$description),'P.val.adj']
	if (all(path.rect.pval==0)){
		path.rect.pval <- path.rect.stat
		warning('All adjusted p values for the selected pathways are 0. The colors of bars representing pathways are adjusted to indicate the t values (in t test) or W values (in wilcox test) from difference tests')
	}else{
		path.rect.pval <- -log10(path.rect.pval)
		path.rect.pval <- p.remove.inf(path.rect.pval)
	}
	bar.path.col <- LRcolor(path.rect.pval, user.set.col=bar.pathway.col)

	
	### the receptor and their coordinate, and their color
	cur.uniq.rep <- factor(top.rep.name, levels=top.rep.name)
	cur.uniq.rep <- cur.uniq.rep[order(cur.uniq.rep, decreasing=TRUE)]
	n.rep <- length(cur.uniq.rep)
	
	top.markerR.dat <- subset(object@interact.filter$markerR, subset=(cluster==select.ident & gene %in% top.rep.name))

	fc.rep <- top.markerR.dat[match(cur.uniq.rep, top.markerR.dat$gene),'avg_log2FC']
	rep.size <- 10 * dot.gene.size * fc.rep
	rep.pct <-  top.markerR.dat[match(cur.uniq.rep, top.markerR.dat$gene),'pct.1']
	rep.pval <-  top.markerR.dat[match(cur.uniq.rep, top.markerR.dat$gene),'p_val_adj']
	if (all(rep.pval==0)){
		rep.pval <- fc.rep
		warning('All adjusted p values for the top.n.receptors are 0. The colors of dots representing receptors are adjusted to indicate the average log2FC of receptors')
	}else{
		rep.pval <- -log10(rep.pval)
		rep.pval <- p.remove.inf(rep.pval)
	}
	dot.rep.col <- LRcolor(rep.pval, user.set.col=dot.gene.col)


	### the upstream ident and their coordinate, and their color
	up.uniq.ident <- factor(up.in.ident, levels=up.in.ident)
	up.uniq.ident <- up.uniq.ident[order(up.uniq.ident, decreasing=TRUE)]
	n.up.ident <- length(up.uniq.ident)
	if (is.null(dot.ident.col)){
		dot.ident.col <- rev(scales::hue_pal(c=100)(length(ident.label)))
		names(dot.ident.col) <- ident.label
	}else{
		dot.ident.names.col <- names(dot.ident.col)
		if (is.null(dot.ident.names.col)){
			stop(paste0('Wrong dot.ident.col parameter! The colors should be named according to their clusters'))
		}else{
			all.missed.ident <- up.in.ident[which(!(up.in.ident %in% dot.ident.names.col))]
			if (length(all.missed.ident) > 0){ stop(paste0('Wrong dot.ident.col parameter! Please add colors for these clusters: ', pasteIdent(all.missed.ident), '.')) }
		}
	}
	dot.up.ident.col <- rev(dot.ident.col[up.in.ident])
	
	### adjust the coordinate of each element
	max.limit <- max(n.up.ident, n.rep, top.n.path)
	up.ident.new.coor <- seq(from=1, to=max.limit, length.out=n.up.ident)
	rep.new.coor <- seq(from=1, to=max.limit, length.out=n.rep)
	pathway.new.coor <- seq(from=1, to=max.limit, length.out=top.n.path)

	if (n.up.ident==1){ up.ident.new.coor <- median(1:max.limit) }
	if (n.rep==1){ rep.new.coor <- median(1:max.limit) }
	if (top.n.path==1){ pathway.new.coor <- median(1:max.limit) }

	if (length(pathway.new.coor)==1){
		bar_pathway_width <- 0.3
	}else{
		bar_pathway_width <- 0.3 * (pathway.new.coor[2]-pathway.new.coor[1]) * top.n.path / max.limit
	}
	path.name.max.width <- max(strwidth(path.uniq.name, units="inches", cex=1))
	#bar.pathway.length <- bar.pathway.width*path.name.max.width*path.rect.length
	bar.pathway.length <- 2 * path.rect.length / max(path.rect.length)
	bar_pathway <- geom_rect(aes(xmin=3, xmax=3+bar.pathway.length, ymin=pathway.new.coor-bar_pathway_width, ymax=pathway.new.coor+bar_pathway_width), size=1, fill=bar.path.col, alpha=1, show.legend=F)
	point_receptor <- geom_point(aes(x=2,y=rep.new.coor),size=rep.size,color=dot.rep.col) 
	point_up_ident <- geom_point(aes(x=1,y=up.ident.new.coor),size=8*dot.ident.size,color=dot.up.ident.col,shape=16) 

	### lines between elements
	### lines between up idents and receptors
	line.ident.to.rep.y.coor <- up.ident.new.coor[match(plot.ident.to.receprtor.dat$cell.from, up.uniq.ident)]
	line.ident.to.rep.yend.coor <- rep.new.coor[match(plot.ident.to.receprtor.dat$receptor, cur.uniq.rep)]
	line.col <- dot.up.ident.col[match(plot.ident.to.receprtor.dat$cell.from, up.uniq.ident)]
	line.ident.width <- LRinten.to.width(plot.ident.to.receprtor.dat$LR.inten)

	line_ident_to_receptor <- geom_segment(aes(x=1, xend=2, y=line.ident.to.rep.y.coor,yend=line.ident.to.rep.yend.coor), size=line.ident.width, color=line.col, show.legend=F, alpha=0.6)
	### lines between receptors and pathways
	line.rep.to.path.y.coor <- rep.new.coor[match(plot.receptor.to.pathway.dat$cur.rep, cur.uniq.rep)]
	line.rep.to.path.yend.coor <- pathway.new.coor[match(plot.receptor.to.pathway.dat$cur.path, path.uniq.name)]
	path.line.col <- dot.rep.col[match(plot.receptor.to.pathway.dat$cur.rep, cur.uniq.rep)]

	line_receptor_to_pathway <- geom_segment(aes(x=2, xend=3, y=line.rep.to.path.y.coor,yend=line.rep.to.path.yend.coor),size=1, color=path.line.col, show.legend=F)
	
	### text and label
	label_up_ident <- annotate('text', hjust=0.5, x=1, y=up.ident.new.coor-0.4,label=up.uniq.ident, size=5*label.text.size)
	label_receptor <- annotate('text', hjust=0.5, x=2, y=rep.new.coor-0.4,label=cur.uniq.rep, size=5*label.text.size)
	text_pathway <- geom_text(aes(x=3, y=pathway.new.coor), hjust=0, label=path.uniq.name, size=5*label.text.size)
	label_title <- annotate('text', hjust=c(0.5,0.5,0), x=1:3, y=max.limit+0.8,label=c('Upstream', 'Receptor', 'Pathway annotation'), size=5*label.title.size) 
	
	
	pathplot.theme <- theme(axis.title=element_blank(), axis.text=element_blank(),axis.ticks=element_blank(), axis.line=element_blank(), panel.background=element_rect(fill="white"))
	plot <- ggplot() + 
	line_ident_to_receptor +
	line_receptor_to_pathway +
	point_up_ident +
	point_receptor + 
	bar_pathway +
	text_pathway +
	label_up_ident +
	label_receptor +
	label_title +
	scale_x_continuous(limits=c(0,3+max(path.name.max.width,bar.pathway.length))) + 
	pathplot.theme

	return(plot)
}


#' To compare differentially activated pathways and the involved receptors between the selected clusters of two CommPath object
#' @param object.1 CommPath object 1
#' @param object.2 CommPath object 2 for comparison
#' @param select.ident Plot the activated pathways for which cluster or cell type?
#' @param up.ident Upstream clusters in the communication chain; default is all clusters
#' @param diff.path.dat Data frame of defferential activation test result from comparePath; if NULL, comparePath would be run to get diff.path.dat
#' @param top.n.path Top n pathways with the smallest adjusted p values to plot
#' @param path.order Sort criteria used to select the top n pathways, either 'P.val' or 'P.val.adj', which represent the original and adjusted p values, or 'diff' which represents the mean (in t test) or median (in wilcox test) difference
#' @param select.path Vector of pathways for which users want to present a chain plot; default is all, which means to plot all pathways. The parameter top.n.path would be masked if this one is set with specific pathways.
#' @param p.thre Threshold for adjust p values; Only pathways with a adjust p valua < p.thre would be considered
#' @param top.n.receptor Top n receptor with the largest log2FCs to plot
#' @param dot.ident.col Color of the dots representing clusters
#' @param dot.ident.size Size of the dots representing clusters
#' @param dot.gene.col Color of the dots representing receptors
#' @param dot.gene.size Size of the dots representing receptors
#' @param bar.pathway.col Color of the bars reprsenting pathways
#' @param label.text.size Text size in the plot
#' @param label.title.size Text size of the title annotation of the plot
#' @return Network plot showing the significantly differentially activated pathways between the selected clusters of two CommPath object, receptors involved in the pathways, and the upstream clusters which show LR connections with the selected cluster
#' @export
pathPlot.compare <- function(object.1, object.2, select.ident, up.ident=NULL, diff.path.dat=NULL, top.n.path=5, path.order='P.val.adj', select.path='all', p.thre=0.05, top.n.receptor=10, dot.ident.col=NULL, dot.ident.size=1, dot.gene.col=NULL, dot.gene.size=1, bar.pathway.col=NULL, label.text.size=1, label.title.size=1){
	options(stringsAsFactors=F)
	if (is.null(diff.path.dat)){
		message(paste0('Identifying differentially activated pathways between cluster ',select.ident,' in object 1 and object 2'))
		diff.path.dat <- comparePath(object.1, object.2, select.ident, method='t.test', p.adjust='BH', min.size=10, only.posi=FALSE, only.sig=TRUE)
	}
	all.ident <- object.1@cell.info$Cluster
	if (!is.factor(all.ident)){ all.ident <- factor(all.ident) }
	ident.label <- levels(all.ident)
	
	# check the input up.ident
	if (!is.null(up.ident)){
		up.in.ident <- up.ident[which(up.ident %in% ident.label)]
		if (length(up.in.ident)==0){
			stop('Wrong up.ident parameter! All selected up.idents are not present in the dataset!')
		}
		up.missed.ident <- up.ident[which(!(up.ident %in% ident.label))]
		if (length(up.missed.ident) > 0){ warning(paste0('These selected up.idents are not present in the dataset: ', pasteIdent(up.missed.ident), '.')) }
	}else{
		up.in.ident <- ident.label
	}

	### check the input and select significant pathways
	ident.path.dat <- diff.path.dat
	if ('t' %in% colnames(ident.path.dat)){
		ident.path.dat <- subset(ident.path.dat, P.val.adj < p.thre & t > 0)
		ident.path.dat$stat <- ident.path.dat$t
		ident.path.dat$diff <- ident.path.dat$mean.diff
	}else if ('W' %in% colnames(ident.path.dat)){
		ident.path.dat <- subset(ident.path.dat, P.val.adj < p.thre & median.diff > 0)
		ident.path.dat$stat <- ident.path.dat$W
		ident.path.dat$diff <- ident.path.dat$median.diff
	}else{
		stop('Please input the integrate ident.path.dat computed from diffPath')
	}

	if (nrow(ident.path.dat)==0){ stop(paste0('There is no significantly up-regulatged pathways for cluster ', select.ident )) }

	### remove pathways with NA receptor
	ident.path.dat <- subset(ident.path.dat, !is.na(receptor.in.path))
	if (nrow(ident.path.dat)==0){ stop(paste0('There is no marker receptor in the significantly up-regulatged pathways for cluster ', select.ident )) }

	if (path.order=='P.val.adj' | path.order=='P.val'){
		ident.path.dat <- ident.path.dat[order(ident.path.dat[,path.order], -ident.path.dat[,'stat']),]
	}else if(path.order=='diff'){
		ident.path.dat <- ident.path.dat[order(-ident.path.dat[,'diff'], ident.path.dat[,'P.val.adj']),]
	}

	### select user-defined pathways, with the select.path parameter
	if (length(select.path) > 1 | (length(select.path)==1 & select.path[1]!='all')){
		path.num.plot <- which(ident.path.dat$description %in% select.path)
		if (length(path.num.plot) > 0){
			ident.path.dat <- ident.path.dat[path.num.plot, ]
		}else{
			stop('All pathways selected are not up-regulatged in the selected cluster. Run comparePath to retrieve all up-regulatged pathways.')
		}
		top.n.path <- length(path.num.plot)
		path.missed <- select.path[which(!(select.path %in% ident.path.dat$description))]
		if (length(path.missed) > 0){
			warning(paste0('These pathways are not up-regulatged in the selected cluster: ', pasteIdent(path.missed), '. Run comparePath to retrieve all up-regulatged pathways.'))
		}
	}else{
		### select top n pathways, if the select.path parameter is not set
		if (top.n.path > nrow(ident.path.dat)){
			warning(paste0('There is(are) ', nrow(ident.path.dat),' significant pathway(s) for the selected ident, and the input top.n.path is ', top.n.path))
			top.n.path <- nrow(ident.path.dat)
		}
		ident.path.dat <- ident.path.dat[1:top.n.path, ]
	}
	all.sig.path <- as.vector(ident.path.dat$description)

	### data for plot of receptors and upstream clusters
	receptor.in.path <- as.character(ident.path.dat[, 'receptor.in.path'])
	names(receptor.in.path) <- all.sig.path
	# receptor.in.path <- receptor.in.path[which(!is.na(receptor.in.path))]
	cur.rep <- as.vector(unlist(sapply(receptor.in.path, function(x) {strsplit(x, split=';')})))
	path.times <- sapply(receptor.in.path, function(x) {length(strsplit(x, split=';')[[1]])})
	cur.path <- rep(names(path.times), times=path.times)
	plot.receptor.to.pathway.dat <- data.frame(cur.rep=cur.rep, cur.path=cur.path)

	### select the top n receptors
	cur.uniq.rep <- as.vector(unique(plot.receptor.to.pathway.dat$cur.rep))
	cur.rep.LR.inten <- cluster.lr.inten(cur.uniq.rep, object.1, select.ident, ident.label, find='ligand')
	cur.rep.LR.inten <- subsetMatrix(cur.rep.LR.inten, row.select=up.in.ident)
	keep.rep <- names(which(colSums(is.na(cur.rep.LR.inten)) < nrow(cur.rep.LR.inten)))
	cur.rep.LR.inten <- subsetMatrix(cur.rep.LR.inten, col.select=keep.rep)

	cur.rep.max.LR.inten <- apply(cur.rep.LR.inten, 2, function(x){max(x, na.rm=TRUE)})
	cur.rep.max.LR.inten <- cur.rep.max.LR.inten[order(cur.rep.max.LR.inten, decreasing=TRUE)]
	
	if (length(cur.rep.max.LR.inten) < top.n.receptor){
		warning(paste0('There is(are) ', length(cur.rep.max.LR.inten),' marker receptors(s) in the selected ident associated with the selected pathways, and the input top.n.receptor is ', top.n.receptor))
		top.n.receptor <- length(cur.rep.max.LR.inten)
	}
	
	top.rep.name <- names(cur.rep.max.LR.inten)[1:top.n.receptor]
	plot.receptor.to.pathway.dat <- subset(plot.receptor.to.pathway.dat, cur.rep %in% top.rep.name)
	top.rep.LR.inten <- subsetMatrix(cur.rep.LR.inten, col.select=top.rep.name)
	plot.ident.to.receprtor.obj1.dat <- melt(top.rep.LR.inten, varnames=c('cell.from','receptor'),value.name="LR.inten",  na.rm=FALSE)
	plot.ident.to.receprtor.obj1.dat <- factor.to.character(plot.ident.to.receprtor.obj1.dat)

	### LR.inten for obj.2
	cur.rep.LR.obj2.inten <- cluster.lr.inten(top.rep.name, object.2, select.ident, ident.label, find='ligand')
	cur.rep.LR.obj2.inten <- subsetMatrix(cur.rep.LR.obj2.inten, row.select=up.in.ident, col.select=top.rep.name)

	plot.ident.to.receprtor.obj2.dat <- melt(cur.rep.LR.obj2.inten, varnames=c('cell.from','receptor'),value.name="LR.inten",  na.rm=FALSE)
	plot.ident.to.receprtor.obj2.dat <- factor.to.character(plot.ident.to.receprtor.obj2.dat)

	LR.1.Dens <- plot.ident.to.receprtor.obj1.dat$LR.inten
	LR.2.Dens <- plot.ident.to.receprtor.obj2.dat$LR.inten

	plot.ident.to.receprtor.obj1.dat$LR.inten.diff <- as.numeric(LR.1.Dens > LR.2.Dens)
	plot.ident.to.receprtor.obj1.dat[which((is.na(LR.1.Dens)) & (!is.na(LR.2.Dens))), 'LR.inten.diff'] <- 0
	plot.ident.to.receprtor.obj1.dat[which((!is.na(LR.1.Dens)) & (is.na(LR.2.Dens))), 'LR.inten.diff'] <- 1
	plot.ident.to.receprtor.obj1.dat <- subset(plot.ident.to.receprtor.obj1.dat, !is.na(LR.inten))
	plot.ident.to.receprtor.dat <- plot.ident.to.receprtor.obj1.dat
	
	### pathways and their color
	path.uniq.name <- factor(all.sig.path, levels=all.sig.path)
	path.uniq.name <- path.uniq.name[order(path.uniq.name, decreasing=TRUE)]
	
	#### rect for pathway
	path.rect.length <- ident.path.dat[match(path.uniq.name, ident.path.dat$description),'diff']
	path.rect.stat <- ident.path.dat[match(path.uniq.name, ident.path.dat$description),'stat']
	path.rect.pval <-  ident.path.dat[match(path.uniq.name, ident.path.dat$description),'P.val.adj']

	if (all(path.rect.pval==0)){
		path.rect.pval <- path.rect.stat
		warning('All adjusted p values for the selected pathways are 0. The colors of bars representing pathways are adjusted to indicate the t values (in t test) or W values (in wilcox test) from difference tests')
	}else{
		path.rect.pval <- -log10(path.rect.pval)
		path.rect.pval <- p.remove.inf(path.rect.pval)
	}
	bar.path.col <- LRcolor(path.rect.pval, user.set.col=bar.pathway.col)

	
	### the receptor and their coordinate, and their color
	cur.uniq.rep <- factor(top.rep.name, levels=top.rep.name)
	cur.uniq.rep <- cur.uniq.rep[order(cur.uniq.rep, decreasing=TRUE)]
	n.rep <- length(cur.uniq.rep)
	
	top.markerR.dat <- subset(object.1@interact$markerR, subset=(cluster==select.ident & gene %in% top.rep.name))
	fc.rep <- top.markerR.dat[match(cur.uniq.rep, top.markerR.dat$gene),'avg_log2FC']
	rep.size <- 10 * dot.gene.size * fc.rep
	rep.pct <-  top.markerR.dat[match(cur.uniq.rep, top.markerR.dat$gene),'pct.1']
	
	rep.pval <-  top.markerR.dat[match(cur.uniq.rep, top.markerR.dat$gene),'p_val_adj']
	if (all(rep.pval==0)){
		rep.pval <- fc.rep
		warning('All adjusted p values for the top.n.receptors are 0. The colors of dots representing receptors are adjusted to indicate the average log2FC of receptors')
	}else{
		rep.pval <- -log10(rep.pval)
		rep.pval <- p.remove.inf(rep.pval)
	}
	dot.rep.col <- LRcolor(rep.pval, user.set.col=dot.gene.col)

	
	### the upstream ident and their coordinate, and their color
	up.uniq.ident <- factor(up.in.ident, levels=up.in.ident)
	up.uniq.ident <- up.uniq.ident[order(up.uniq.ident, decreasing=TRUE)]
	n.up.ident <- length(up.uniq.ident)
	### color for dots representing ident
	if (is.null(dot.ident.col)){
		dot.ident.col <- rev(scales::hue_pal(c=100)(length(ident.label)))
		names(dot.ident.col) <- ident.label
	}else{
		dot.ident.names.col <- names(dot.ident.col)
		if (is.null(dot.ident.names.col)){
			stop(paste0('Wrong dot.ident.col parameter! The colors should be named according to their clusters'))
		}else{
			all.missed.ident <- up.in.ident[which(!(up.in.ident %in% dot.ident.names.col))]
			if (length(all.missed.ident) > 0){ stop(paste0('Wrong dot.ident.col parameter! Please add colors for these clusters: ', pasteIdent(all.missed.ident), '.')) }
		}

	}
	dot.up.ident.col <- rev(dot.ident.col[up.in.ident])

	### adjust the coordinate of each element
	max.limit <- max(n.up.ident, n.rep, top.n.path)
	up.ident.new.coor <- seq(from=1, to=max.limit, length.out=n.up.ident)
	rep.new.coor <- seq(from=1, to=max.limit, length.out=n.rep)
	pathway.new.coor <- seq(from=1, to=max.limit, length.out=top.n.path)

	if (n.up.ident==1){ up.ident.new.coor <- median(1:max.limit) }
	if (n.rep==1){ rep.new.coor <- median(1:max.limit) }
	if (top.n.path==1){ pathway.new.coor <- median(1:max.limit) }

	if (length(pathway.new.coor)==1){
		bar_pathway_width <- 0.3
	}else{
		bar_pathway_width <- 0.3 * (pathway.new.coor[2]-pathway.new.coor[1]) * top.n.path / max.limit
	}
	path.name.max.width <- max(strwidth(path.uniq.name, units="inches", cex=1))
	#bar.pathway.length <- bar.pathway.width*path.name.max.width*path.rect.length
	bar.pathway.length <- 2 * path.rect.length / max(path.rect.length)
	bar_pathway <- geom_rect(aes(xmin=3, xmax=3+bar.pathway.length, ymin=pathway.new.coor-bar_pathway_width, ymax=pathway.new.coor+bar_pathway_width), size=1, fill=bar.path.col, alpha=1, show.legend=F)
	point_receptor <- geom_point(aes(x=2,y=rep.new.coor),size=rep.size,color=dot.rep.col) 
	point_up_ident <- geom_point(aes(x=1,y=up.ident.new.coor),size=8*dot.ident.size,color=dot.up.ident.col,shape=16) 

	### lines between elements
	### lines between up idents and receptors
	line.ident.to.rep.y.coor <- up.ident.new.coor[match(plot.ident.to.receprtor.dat$cell.from, up.uniq.ident)]
	line.ident.to.rep.yend.coor <- rep.new.coor[match(plot.ident.to.receprtor.dat$receptor, cur.uniq.rep)]
	#line.col <- dot.ident.col[match(plot.ident.to.receprtor.dat$cell.from, up.uniq.ident)]
	line.ident.to.rep.col <- ifelse(plot.ident.to.receprtor.dat$LR.inten.diff > 0, 'red', 'blue')
	line.ident.to.rep.width <- LRinten.to.width(plot.ident.to.receprtor.dat$LR.inten)

	line_ident_to_receptor <- geom_segment(aes(x=1, xend=2, y=line.ident.to.rep.y.coor,yend=line.ident.to.rep.yend.coor), size=line.ident.to.rep.width, color=line.ident.to.rep.col, show.legend=F, alpha=0.6)
	### lines between receptors and pathways
	line.rep.to.path.y.coor <- rep.new.coor[match(plot.receptor.to.pathway.dat$cur.rep, cur.uniq.rep)]
	line.rep.to.path.yend.coor <- pathway.new.coor[match(plot.receptor.to.pathway.dat$cur.path, path.uniq.name)]
	path.line.col <- dot.rep.col[match(plot.receptor.to.pathway.dat$cur.rep, cur.uniq.rep)]

	line_receptor_to_pathway <- geom_segment(aes(x=2, xend=3, y=line.rep.to.path.y.coor,yend=line.rep.to.path.yend.coor),size=1, color=path.line.col, show.legend=F)
	
	### text and label
	label_up_ident <- annotate('text', hjust=0.5, x=1, y=up.ident.new.coor-0.4,label=up.uniq.ident, size=5*label.text.size)
	label_receptor <- annotate('text', hjust=0.5, x=2, y=rep.new.coor-0.4,label=cur.uniq.rep, size=5*label.text.size)
	text_pathway <- geom_text(aes(x=3, y=pathway.new.coor), hjust=0, label=path.uniq.name, size=5*label.text.size)
	label_title <- annotate('text', hjust=c(0.5,0.5,0), x=1:3, y=max.limit+0.8,label=c('Upstream', 'Receptor', 'Pathway annotation'), size=5*label.title.size) 

	pathplot.theme <- theme(axis.title=element_blank(), axis.text=element_blank(),axis.ticks=element_blank(), axis.line=element_blank(), panel.background=element_rect(fill="white"))
	plot <- ggplot() + 
		line_ident_to_receptor + line_receptor_to_pathway +
		point_up_ident + point_receptor + 
		bar_pathway + text_pathway +
		label_up_ident + label_receptor + label_title +
		scale_x_continuous(limits=c(0,3+max(path.name.max.width,bar.pathway.length))) + 
		pathplot.theme

	return(plot)
}


#' To present the interactions for a selected cluster, including both the upstream and downstream clusters which are connected by specific pathways in the selected cluster
#' @param object CommPath object
#' @param select.ident Plot the activated pathways for which cluster or cell type?
#' @param up.ident Upstream clusters in the communication chain; default is all clusters
#' @param down.ident down.ident clusters in the communication chain; default is all clusters
#' @param acti.path.dat Data frame of differential activation test result from diffAllPath
#' @param top.n.path Top n pathways with the smallest adjusted p values to plot
#' @param path.order Sort criteria used to select the top n pathways, either 'P.val' or 'P.val.adj', which represent the original and adjusted p values, or 'diff' which represents the mean (in t test) or median (in wilcox test) difference
#' @param select.path Vector of pathways for which users want to present a chain plot; default is all, which means to plot all pathways. The parameter top.n.path would be masked if this one is set with specific pathways.
#' @param p.thre Threshold for adjust p values; Only pathways with a adjust p valua < p.thre would be considered
#' @param top.n.receptor Top n receptor with the largest log2FCs to plot
#' @param top.n.ligand Top n ligand with the largest log2FCs to plot
#' @param dot.ident.col Color of the dots representing clusters
#' @param dot.ident.size Size of the dots representing clusters
#' @param dot.gene.col Color of the dots representing receptors and ligands
#' @param dot.gene.size Size of the dots representing receptors and ligands
#' @param bar.pathway.col Color of the bars reprsenting pathways
#' @param label.text.size Text size in the plot
#' @param label.title.size Text size of the title annotation of the plot
#' @return Network plot showing receptors in the selected cluster, the upstream clusters which show L-R connections with the selected cluster, and the significant pathways involved in the receptors
#' @export
pathChainPlot <- function(object, select.ident, up.ident=NULL, down.ident=NULL, acti.path.dat=NULL, top.n.path=5, path.order='P.val.adj', select.path='all', p.thre=0.05, top.n.receptor=10, top.n.ligand=10, dot.ident.col=NULL, dot.ident.size=1, dot.gene.col=NULL, dot.gene.size=1, bar.pathway.col=NULL, label.text.size=1, label.title.size=1){
	options(stringsAsFactors=F)
	all.ident <- object@cell.info$Cluster
	if (!is.factor(all.ident)){ all.ident <- factor(all.ident) }
	ident.label <- levels(all.ident)
	
	# check the input up.ident
	if (!is.null(up.ident)){
		up.in.ident <- up.ident[which(up.ident %in% ident.label)]
		if (length(up.in.ident)==0){
			stop('Wrong up.ident parameter! All selected up.idents are not present in the dataset!')
		}
		up.missed.ident <- up.ident[which(!(up.ident %in% ident.label))]
		if (length(up.missed.ident) > 0){ warning(paste0('These selected up.idents are not present in the dataset: ', pasteIdent(up.missed.ident), '.')) }
	}else{
		up.in.ident <- ident.label
	}

	# check the input down.ident
	if (!is.null(down.ident)){
		down.in.ident <- down.ident[which(down.ident %in% ident.label)]
		if (length(down.in.ident)==0){
			stop('Wrong down.ident parameter! All selected down.idents are not present in the dataset!')
		}
		down.missed.ident <- down.ident[which(!(down.ident %in% ident.label))]
		if (length(down.missed.ident) > 0){ warning(paste0('These selected down.ident are not present in the dataset: ', pasteIdent(down.missed.ident), '.')) }
	}else{
		down.in.ident <- ident.label
	}

	### check the input and select significant pathways
	if (is.null(acti.path.dat)) {
		ident.path.dat <- diffPath(object, select.ident.1=select.ident)
	}else{
		ident.path.dat <- subset(acti.path.dat, cluster==select.ident)
	}

	if ('t' %in% colnames(ident.path.dat)){
		ident.path.dat <- subset(ident.path.dat, P.val.adj < p.thre & t > 0)
		ident.path.dat$stat <- ident.path.dat$t
		ident.path.dat$diff <- ident.path.dat$mean.diff
	}else if ('W' %in% colnames(ident.path.dat)){
		ident.path.dat <- subset(ident.path.dat, P.val.adj < p.thre & median.diff > 0)
		ident.path.dat$stat <- ident.path.dat$W
		ident.path.dat$diff <- ident.path.dat$median.diff
	}else{
		stop('Please input the integrate acti.path.dat computed from diffPath')
	}

	if (nrow(ident.path.dat)==0){ stop(paste0('There is no significantly up-regulatged pathways for cluster ', select.ident )) }

	if (path.order=='P.val.adj' | path.order=='P.val'){
		ident.path.dat <- ident.path.dat[order(ident.path.dat[,path.order], -ident.path.dat[,'stat']),]
	}else if(path.order=='diff'){
		ident.path.dat <- ident.path.dat[order(-ident.path.dat[,'diff'], ident.path.dat[,'P.val.adj']),]
	}

	### select user-defined pathways, with the select.path parameter
	if (length(select.path) > 1 | (length(select.path)==1 & select.path[1]!='all')){
		path.num.plot <- which(ident.path.dat$description %in% select.path)
		if (length(path.num.plot) > 0){
			ident.path.dat <- ident.path.dat[path.num.plot, ]
		}else{
			stop('All pathways selected are not up-regulatged in the selected cluster. Run diffPath to retrieve all up-regulatged pathways.')
		}
		top.n.path <- length(path.num.plot)
		path.missed <- select.path[which(!(select.path %in% ident.path.dat$description))]
		if (length(path.missed) > 0){
			warning(paste0('These pathways are not up-regulatged in the selected cluster: ', pasteIdent(path.missed), '. Run diffPath to retrieve all up-regulatged pathways.'))
		}
	}else{
		### select top n pathways, if the select.path parameter is not set
		if (top.n.path > nrow(ident.path.dat)){
			warning(paste0('There is(are) ', nrow(ident.path.dat),' significant pathway(s) for the selected ident, and the input top.n.path is ', top.n.path))
			top.n.path <- nrow(ident.path.dat)
		}
		ident.path.dat <- ident.path.dat[1:top.n.path, ]
	}
	all.sig.path <- as.vector(ident.path.dat$description)

	### data for plot of receptors and upstream clusters
	receptor.in.path <- as.character(ident.path.dat[, 'receptor.in.path'])
	if (all(is.na(receptor.in.path))){
		stop('There is no marker receptor in the selected pathways\nselect more pathways and try again')
	}
	names(receptor.in.path) <- all.sig.path
	receptor.in.path <- receptor.in.path[which(!is.na(receptor.in.path))]
	cur.rep <- as.vector(unlist(sapply(receptor.in.path, function(x) {strsplit(x, split=';')})))
	path.rep.times <- sapply(receptor.in.path, function(x) {length(strsplit(x, split=';')[[1]])})
	cur.rep.path <- rep(names(path.rep.times), times=path.rep.times)
	plot.receptor.to.pathway.dat <- data.frame(cur.rep=cur.rep, cur.path=cur.rep.path)
	### select the top n receptors
	cur.uniq.rep <- as.vector(unique(plot.receptor.to.pathway.dat$cur.rep))
	cur.rep.LR.inten <- cluster.lr.inten(cur.uniq.rep, object, select.ident, ident.label, find='ligand')
	cur.rep.LR.inten <- subsetMatrix(cur.rep.LR.inten, row.select=up.in.ident)
	keep.rep <- names(which(colSums(is.na(cur.rep.LR.inten)) < nrow(cur.rep.LR.inten)))
	cur.rep.LR.inten <- subsetMatrix(cur.rep.LR.inten, col.select=keep.rep)

	cur.rep.max.LR.inten <- apply(cur.rep.LR.inten, 2, function(x){max(x, na.rm=TRUE)})
	cur.rep.max.LR.inten <- cur.rep.max.LR.inten[order(cur.rep.max.LR.inten, decreasing=TRUE)]
	
	if (length(cur.rep.max.LR.inten) < top.n.receptor){
		warning(paste0('There is(are) ', length(cur.rep.max.LR.inten),' marker receptor(s) in the selected ident associated with the selected pathways, and the input top.n.receptor is ', top.n.receptor))
		top.n.receptor <- length(cur.rep.max.LR.inten)
	}
	
	top.rep.name <- names(cur.rep.max.LR.inten)[1:top.n.receptor]
	plot.receptor.to.pathway.dat <- subset(plot.receptor.to.pathway.dat, cur.rep %in% top.rep.name)
	top.rep.LR.inten <- subsetMatrix(cur.rep.LR.inten, col.select=top.rep.name)

	plot.ident.to.receprtor.dat <- melt(top.rep.LR.inten, varnames=c('cell.from','receptor'),value.name="LR.inten",  na.rm=TRUE)
	plot.ident.to.receprtor.dat <- factor.to.character(plot.ident.to.receprtor.dat)
	
	### data for plot of ligands and downstream clusters
	ligand.in.path <- as.character(ident.path.dat[, 'ligand.in.path'])
	if (all(is.na(ligand.in.path))){
		stop('There is no marker ligand in the selected pathways\nselect more pathways and try again\nor you may want to try pathPlot instead')
	}
	names(ligand.in.path) <- all.sig.path
	ligand.in.path <- ligand.in.path[which(!is.na(ligand.in.path))]
	cur.lig <- as.vector(unlist(sapply(ligand.in.path, function(x) {strsplit(x, split=';')})))
	path.lig.times <- sapply(ligand.in.path, function(x) {length(strsplit(x, split=';')[[1]])})
	cur.lig.path <- rep(names(path.lig.times), times=path.lig.times)
	plot.pathway.to.ligand.dat <- data.frame(cur.lig=cur.lig, cur.path=cur.lig.path)
	### select the top n ligands
	cur.uniq.lig <- as.vector(unique(plot.pathway.to.ligand.dat$cur.lig))
	cur.lig.LR.inten <- cluster.lr.inten(cur.uniq.lig, object, select.ident, ident.label, find='receptor')

	cur.lig.LR.inten <- subsetMatrix(cur.lig.LR.inten, row.select=down.in.ident)
	keep.lig <- names(which(colSums(is.na(cur.lig.LR.inten)) < nrow(cur.lig.LR.inten)))
	cur.lig.LR.inten <- subsetMatrix(cur.lig.LR.inten, col.select=keep.lig)

	cur.lig.max.LR.inten <- apply(cur.lig.LR.inten, 2, function(x){max(x, na.rm=TRUE)})
	cur.lig.max.LR.inten <- cur.lig.max.LR.inten[order(cur.lig.max.LR.inten, decreasing=TRUE)]
	
	if (length(cur.lig.max.LR.inten) < top.n.ligand){
		warning(paste0('There is(are) ', length(cur.lig.max.LR.inten),' marker ligand(s) in the selected ident associated with the selected pathways, and the input top.n.ligand is ', top.n.ligand))
		top.n.ligand <- length(cur.lig.max.LR.inten)
	}
	
	top.lig.name <- names(cur.lig.max.LR.inten)[1:top.n.ligand]
	plot.pathway.to.ligand.dat <- subset(plot.pathway.to.ligand.dat, cur.lig %in% top.lig.name)
	top.lig.LR.inten <- subsetMatrix(cur.lig.LR.inten, col.select=top.lig.name)

	plot.ligand.to.ident.dat <- melt(top.lig.LR.inten, varnames=c('cell.to','ligand'),value.name="LR.inten",  na.rm=TRUE)
	plot.ligand.to.ident.dat <- factor.to.character(plot.ligand.to.ident.dat)

	### pathways and their color
	path.uniq.name <- factor(all.sig.path, levels=all.sig.path)
	path.uniq.name <- path.uniq.name[order(path.uniq.name, decreasing=TRUE)]
	
	#### rect for pathway
	path.rect.length <- ident.path.dat[match(path.uniq.name, ident.path.dat$description),'diff']
	path.rect.stat <- ident.path.dat[match(path.uniq.name, ident.path.dat$description),'stat']
	path.rect.pval <-  ident.path.dat[match(path.uniq.name, ident.path.dat$description),'P.val.adj']

	if (all(path.rect.pval==0)){
		path.rect.pval <- path.rect.stat
		warning('All adjusted p values for the selected pathways are 0. The colors of bars representing pathways are adjusted to indicate the t values (in t test) or W values (in wilcox test) from difference tests')
	}else{
		path.rect.pval <- -log10(path.rect.pval)
		path.rect.pval <- p.remove.inf(path.rect.pval)
	}
	bar.path.col <- LRcolor(path.rect.pval, user.set.col=bar.pathway.col)

	
	### the receptor and their coordinate, and their color
	cur.uniq.rep <- factor(top.rep.name, levels=top.rep.name)
	cur.uniq.rep <- cur.uniq.rep[order(cur.uniq.rep, decreasing=TRUE)]
	n.rep <- length(cur.uniq.rep)
	
	top.markerR.dat <- subset(object@interact.filter$markerR, subset=(cluster==select.ident & gene %in% top.rep.name))

	fc.rep <- top.markerR.dat[match(cur.uniq.rep, top.markerR.dat$gene),'avg_log2FC']
	rep.size <- 10 * dot.gene.size * fc.rep
	rep.pct <-  top.markerR.dat[match(cur.uniq.rep, top.markerR.dat$gene),'pct.1']
	rep.pval <-  top.markerR.dat[match(cur.uniq.rep, top.markerR.dat$gene),'p_val_adj']

	if (all(rep.pval==0)){
		rep.pval <- fc.rep
		warning('All adjusted p values for the top.n.receptors are 0. The colors of dots representing receptors are adjusted to indicate the average log2FC of receptors')
	}else{
		rep.pval <- -log10(rep.pval)
		rep.pval <- p.remove.inf(rep.pval)
	}
	dot.rep.col <- LRcolor(rep.pval, user.set.col=dot.gene.col)


	### the ligand and their coordinate, and their color
	cur.uniq.lig <- factor(top.lig.name, levels=top.lig.name)
	cur.uniq.lig <- cur.uniq.lig[order(cur.uniq.lig, decreasing=TRUE)]
	n.lig <- length(cur.uniq.lig)
	
	top.markerL.dat <- subset(object@interact.filter$markerL, subset=(cluster==select.ident & gene %in% top.lig.name))

	fc.lig <- top.markerL.dat[match(cur.uniq.lig, top.markerL.dat$gene),'avg_log2FC']
	lig.size <- 10 * dot.gene.size * fc.lig
	lig.pct <-  top.markerL.dat[match(cur.uniq.lig, top.markerL.dat$gene),'pct.1']
	lig.pval <-  top.markerL.dat[match(cur.uniq.lig, top.markerL.dat$gene),'p_val_adj']

	if (all(lig.pval==0)){
		lig.pval <- fc.lig
		warning('All adjusted p values for the top.n.ligands are 0. The colors of dots representing ligands are adjusted to indicate the average log2FC of ligands')
	}else{
		lig.pval <- -log10(lig.pval)
		lig.pval <- p.remove.inf(lig.pval)
	}
	dot.lig.col <- LRcolor(lig.pval, user.set.col=dot.gene.col)

	### the upstream ident and their coordinate, and their color
	up.uniq.ident <- factor(up.in.ident, levels=up.in.ident)
	up.uniq.ident <- up.uniq.ident[order(up.uniq.ident, decreasing=TRUE)]
	### n.up.ident is also the n.down.ident
	n.up.ident <- length(up.uniq.ident)
	if (is.null(dot.ident.col)){
		dot.ident.col <- rev(scales::hue_pal(c=100)(length(ident.label)))
		names(dot.ident.col) <- ident.label
	}else{
		dot.ident.names.col <- names(dot.ident.col)
		if (is.null(dot.ident.names.col)){
			stop(paste0('Wrong dot.ident.col parameter! The colors should be named according to their clusters'))
		}else{
			all.in.ident <- unique(c(up.in.ident, down.in.ident))
			all.missed.ident <- all.in.ident[which(!(all.in.ident %in% dot.ident.names.col))]
			if (length(all.missed.ident) > 0){ stop(paste0('Wrong dot.ident.col parameter! Please add colors for these clusters: ', pasteIdent(all.missed.ident), '.')) }
		}

	}
	dot.up.ident.col <- rev(dot.ident.col[up.in.ident])
	### the downstream ident and their coordinate, and their color
	down.uniq.ident <- factor(down.in.ident, levels=down.in.ident)
	down.uniq.ident <- down.uniq.ident[order(down.uniq.ident, decreasing=TRUE)]
	### n.up.ident is also the n.down.ident
	n.down.ident <- length(down.uniq.ident)
	dot.down.ident.col <- rev(dot.ident.col[down.in.ident])

	### adjust the coordinate of each element
	max.limit <- max(n.up.ident, n.rep, top.n.path, n.lig, n.down.ident)
	up.ident.new.coor <- seq(from=1, to=max.limit, length.out=n.up.ident)
	rep.new.coor <- seq(from=1, to=max.limit, length.out=n.rep)
	pathway.new.coor <- seq(from=1, to=max.limit, length.out=top.n.path)
	lig.new.coor <- seq(from=1, to=max.limit, length.out=n.lig)
	down.ident.new.coor <- seq(from=1, to=max.limit, length.out=n.down.ident)

	if (n.up.ident==1){ up.ident.new.coor <- median(1:max.limit) }
	if (n.rep==1){ rep.new.coor <- median(1:max.limit) }
	if (top.n.path==1){ pathway.new.coor <- median(1:max.limit) }
	if (n.lig==1){ lig.new.coor <- median(1:max.limit) }
	if (n.down.ident==1){ down.ident.new.coor <- median(1:max.limit) }
	
	if (length(pathway.new.coor)==1){
		bar_pathway_width <- 0.3
	}else{
		bar_pathway_width <- 0.3 * (pathway.new.coor[2]-pathway.new.coor[1]) * top.n.path / max.limit
	}
	path.name.max.width <- max(strwidth(path.uniq.name, units="inches", cex=1))
	#bar.pathway.length <- bar.pathway.width*path.name.max.width*path.rect.length
	bar.pathway.length <- 2 * path.rect.length / max(path.rect.length)
	bar_pathway <- geom_rect(aes(xmin=6, xmax=6+bar.pathway.length, ymin=pathway.new.coor-bar_pathway_width, ymax=pathway.new.coor+bar_pathway_width), size=1, fill=bar.path.col, alpha=1, show.legend=F)
	point_up_ident <- geom_point(aes(x=1,y=up.ident.new.coor),size=8*dot.ident.size,color=dot.up.ident.col,shape=16) 
	point_receptor <- geom_point(aes(x=2,y=rep.new.coor),size=rep.size,color=dot.rep.col)
	point_pathway <- geom_point(aes(x=3,y=pathway.new.coor),size=8*dot.ident.size,color='black',shape=21,fill='white')
	point_ligand <- geom_point(aes(x=4,y=lig.new.coor),size=lig.size,color=dot.lig.col) 
	point_down_ident <- geom_point(aes(x=5,y=down.ident.new.coor),size=8*dot.ident.size,color=dot.down.ident.col,shape=16) 

	### lines between elements
	### lines between up idents and receptors
	line.ident.to.rep.y.coor <- up.ident.new.coor[match(plot.ident.to.receprtor.dat$cell.from, up.uniq.ident)]
	line.ident.to.rep.yend.coor <- rep.new.coor[match(plot.ident.to.receprtor.dat$receptor, cur.uniq.rep)]
	line.ident.to.rep.col <- dot.up.ident.col[match(plot.ident.to.receprtor.dat$cell.from, up.uniq.ident)]
	line.ident.to.rep.width <- LRinten.to.width(plot.ident.to.receprtor.dat$LR.inten)
	line_ident_to_receptor <- geom_segment(aes(x=1, xend=2, y=line.ident.to.rep.y.coor,yend=line.ident.to.rep.yend.coor), size=line.ident.to.rep.width, color=line.ident.to.rep.col, show.legend=F, alpha=0.6)
	### lines between receptors and pathways
	line.rep.to.path.y.coor <- rep.new.coor[match(plot.receptor.to.pathway.dat$cur.rep, cur.uniq.rep)]
	line.rep.to.path.yend.coor <- pathway.new.coor[match(plot.receptor.to.pathway.dat$cur.path, path.uniq.name)]
	path.rep.line.col <- dot.rep.col[match(plot.receptor.to.pathway.dat$cur.rep, cur.uniq.rep)]
	line_receptor_to_pathway <- geom_segment(aes(x=2, xend=3, y=line.rep.to.path.y.coor,yend=line.rep.to.path.yend.coor),size=1, color=path.rep.line.col, show.legend=F)
	### lines between pathways and ligands
	line.path.to.lig.y.coor <- pathway.new.coor[match(plot.pathway.to.ligand.dat$cur.path, path.uniq.name)]
	line.path.to.lig.yend.coor <- lig.new.coor[match(plot.pathway.to.ligand.dat$cur.lig, cur.uniq.lig)]
	path.lig.line.col <- dot.lig.col[match(plot.pathway.to.ligand.dat$cur.lig, cur.uniq.lig)]
	line_pathway_to_ligand <- geom_segment(aes(x=3, xend=4, y=line.path.to.lig.y.coor,yend=line.path.to.lig.yend.coor),size=1, color=path.lig.line.col, show.legend=F)
	### lines between ligands and down idents
	line.lig.to.ident.y.coor <- lig.new.coor[match(plot.ligand.to.ident.dat$ligand, cur.uniq.lig)]
	line.lig.to.ident.yend.coor <- down.ident.new.coor[match(plot.ligand.to.ident.dat$cell.to, down.uniq.ident)]
	line.lig.to.ident.col <- dot.down.ident.col[match(plot.ligand.to.ident.dat$cell.to, down.uniq.ident)]
	line.lig.to.ident.width <- LRinten.to.width(plot.ligand.to.ident.dat$LR.inten)
	line_ligand_to_ident <- geom_segment(aes(x=4, xend=5, y=line.lig.to.ident.y.coor,yend=line.lig.to.ident.yend.coor), size=line.lig.to.ident.width, color=line.lig.to.ident.col, show.legend=F, alpha=0.6)

	
	### text and label
	label_up_ident <- annotate('text', hjust=0.5, x=1, y=up.ident.new.coor-0.4,label=up.uniq.ident, size=5*label.text.size)
	label_receptor <- annotate('text', hjust=0.5, x=2, y=rep.new.coor-0.4,label=cur.uniq.rep, size=5*label.text.size)
	#text_pathway <- geom_text(aes(x=6.1+path.rect.length, y=pathway.new.coor), hjust=0, label=path.uniq.name, size=5*label.text.size)
	text_pathway <- geom_text(aes(x=6, y=pathway.new.coor), hjust=0, label=path.uniq.name, size=5*label.text.size)
	label_pathway_circle <- annotate('text', hjust=0.5, x=3, y=pathway.new.coor, label=length(path.uniq.name):1, size=5*label.text.size)
	label_pathway_bar <- annotate('text', hjust=1, x=6, y=pathway.new.coor, label=length(path.uniq.name):1, size=5*label.text.size)
	label_ligand <- annotate('text', hjust=0.5, x=4, y=lig.new.coor-0.4,label=cur.uniq.lig, size=5*label.text.size)
	label_down_ident <- annotate('text', hjust=0.5, x=5, y=down.ident.new.coor-0.4,label=down.uniq.ident, size=5*label.text.size)
	label_title <- annotate('text', hjust=c(0.5,0.5,0.5,0.5,0.5,0), x=1:6, y=max.limit+0.8,label=c('Upstream', 'Receptor', 'Pathway', 'Ligand', 'Downstream','Pathway annotation'), size=5*label.title.size)

	pathplot.theme <- theme(axis.title=element_blank(), axis.text=element_blank(),axis.ticks=element_blank(), axis.line=element_blank(), panel.background=element_rect(fill="white"))

	# the max limit of x axis is dependent on the longest pathway
	# this limit may need to be adjusted
	plot <- ggplot() +
	line_ident_to_receptor +
	line_receptor_to_pathway + 
	line_pathway_to_ligand +
	line_ligand_to_ident +
	
	point_up_ident +
	point_receptor +
	point_pathway +
	point_ligand +
	point_down_ident +
	bar_pathway +
	text_pathway +
	
	label_up_ident +
	label_receptor +
	label_ligand +
	label_down_ident +
	label_title +
	label_pathway_circle +
	label_pathway_bar +

	scale_x_continuous(limits=c(0,6+max(path.name.max.width,bar.pathway.length))) +
	pathplot.theme

	return(plot)
}


#' To compare the pathway mediated cell-cell communication flow for a specific cluster between two CommPath object
#' @param object.1 CommPath object 1
#' @param object.2 CommPath object 2 for comparison
#' @param select.ident Plot the activated pathways for which cluster or cell type?
#' @param up.ident Upstream clusters in the communication chain; default is all clusters
#' @param down.ident down.ident clusters in the communication chain; default is all clusters
#' @param diff.path.dat Data frame of defferential activation test result from comparePath; if NULL, comparePath would be run to get diff.path.dat
#' @param top.n.path Top n pathways with the smallest adjusted p values to plot
#' @param path.order Sort criteria used to select the top n pathways, either 'P.val' or 'P.val.adj', which represent the original and adjusted p values, or 'diff' which represents the mean (in t test) or median (in wilcox test) difference
#' @param select.path Vector of pathways for which users want to present a chain plot; default is all, which means to plot all pathways. The parameter top.n.path would be masked if this one is set with specific pathways.
#' @param p.thre Threshold for adjust p values; Only pathways with a adjust p valua < p.thre would be considered
#' @param top.n.receptor Top n receptor with the largest log2FCs to plot
#' @param top.n.ligand Top n ligand with the largest log2FCs to plot
#' @param dot.ident.col Color of the dots representing clusters
#' @param dot.ident.size Size of the dots representing clusters
#' @param dot.gene.col Color of the dots representing receptors
#' @param dot.gene.size Size of the dots representing receptors
#' @param bar.pathway.col Color of the bars reprsenting pathways
#' @param label.text.size Text size in the plot
#' @param label.title.size Text size of the title annotation of the plot
#' @return Network plot showing the significantly differentially activated pathways between the selected clusters of two CommPath object, receptors involved in the pathways, and the upstream clusters which show LR connections with the selected cluster
#' @export
pathChainPlot.compare <- function(object.1, object.2, select.ident, up.ident=NULL, down.ident=NULL, diff.path.dat=NULL, top.n.path=5, path.order='P.val.adj', select.path='all', p.thre=0.05, top.n.receptor=10, top.n.ligand=10, dot.ident.col=NULL, dot.ident.size=1, dot.gene.col=NULL, dot.gene.size=1, bar.pathway.col=NULL, label.text.size=1, label.title.size=1){
	options(stringsAsFactors=F)
	if (is.null(diff.path.dat)){
		message(paste0('Identifying differentially activated pathways between cluster ',select.ident,' in object 1 and object 2'))
		diff.path.dat <- comparePath(object.1, object.2, select.ident, method='t.test', p.adjust='BH', min.size=10, only.posi=FALSE, only.sig=TRUE)
	}
	
	all.ident <- object.1@cell.info$Cluster
	if (!is.factor(all.ident)){ all.ident <- factor(all.ident) }
	ident.label <- levels(all.ident)
	
	# check the input up.ident
	if (!is.null(up.ident)){
		up.in.ident <- up.ident[which(up.ident %in% ident.label)]
		if (length(up.in.ident)==0){
			stop('Wrong up.ident parameter! All selected up.idents are not present in the dataset!')
		}
		up.missed.ident <- up.ident[which(!(up.ident %in% ident.label))]
		if (length(up.missed.ident) > 0){ warning(paste0('These selected up.idents are not present in the dataset: ', pasteIdent(up.missed.ident), '.')) }
	}else{
		up.in.ident <- ident.label
	}
	
	# check the input down.ident
	if (!is.null(down.ident)){
		down.in.ident <- down.ident[which(down.ident %in% ident.label)]
		if (length(down.in.ident)==0){
			stop('Wrong down.ident parameter! All selected down.idents are not present in the dataset!')
		}
		down.missed.ident <- down.ident[which(!(down.ident %in% ident.label))]
		if (length(down.missed.ident) > 0){ warning(paste0('These selected down.ident are not present in the dataset: ', pasteIdent(down.missed.ident), '.')) }
	}else{
		down.in.ident <- ident.label
	}

	### check the input and select significant pathways
	ident.path.dat <- diff.path.dat
	if ('t' %in% colnames(ident.path.dat)){
		ident.path.dat <- subset(ident.path.dat, P.val.adj < p.thre & t > 0)
		ident.path.dat$stat <- ident.path.dat$t
		ident.path.dat$diff <- ident.path.dat$mean.diff
	}else if ('W' %in% colnames(ident.path.dat)){
		ident.path.dat <- subset(ident.path.dat, P.val.adj < p.thre & median.diff > 0)
		ident.path.dat$stat <- ident.path.dat$W
		ident.path.dat$diff <- ident.path.dat$median.diff
	}else{
		stop('Please input the integrate acti.path.dat computed from comparePath')
	}

	if (nrow(ident.path.dat)==0){ stop(paste0('There is no significantly up-regulatged pathways for cluster ', select.ident )) }

	if (path.order=='P.val.adj' | path.order=='P.val'){
		ident.path.dat <- ident.path.dat[order(ident.path.dat[,path.order], -ident.path.dat[,'stat']),]
	}else if(path.order=='diff'){
		ident.path.dat <- ident.path.dat[order(-ident.path.dat[,'diff'], ident.path.dat[,'P.val.adj']),]
	}

	### select user-defined pathways, with the select.path parameter
	if (length(select.path) > 1 | (length(select.path)==1 & select.path[1]!='all')){
		path.num.plot <- which(ident.path.dat$description %in% select.path)
		if (length(path.num.plot) > 0){
			ident.path.dat <- ident.path.dat[path.num.plot, ]
		}else{
			stop('All pathways selected are not up-regulatged in the selected cluster. Run comparePath to retrieve all up-regulatged pathways.')
		}
		top.n.path <- length(path.num.plot)
		path.missed <- select.path[which(!(select.path %in% ident.path.dat$description))]
		if (length(path.missed) > 0){
			warning(paste0('These pathways are not up-regulatged in the selected cluster: ', pasteIdent(path.missed), '. Run comparePath to retrieve all up-regulatged pathways.'))
		}
	}else{
		### select top n pathways, if the select.path parameter is not set
		if (top.n.path > nrow(ident.path.dat)){
			warning(paste0('There is(are) ', nrow(ident.path.dat),' significant pathway(s) for the selected ident, and the input top.n.path is ', top.n.path))
			top.n.path <- nrow(ident.path.dat)
		}
		ident.path.dat <- ident.path.dat[1:top.n.path, ]
	}
	all.sig.path <- as.vector(ident.path.dat$description)

	### data for plot of receptors and upstream clusters
	receptor.in.path <- as.character(ident.path.dat[, 'receptor.in.path'])
	if (all(is.na(receptor.in.path))){
		stop('There is no marker receptor in the selected pathways\nselect more pathways and try again\nor you may want to try pathPlot.compare instead')
	}
	names(receptor.in.path) <- all.sig.path
	receptor.in.path <- receptor.in.path[which(!is.na(receptor.in.path))]
	cur.rep <- as.vector(unlist(sapply(receptor.in.path, function(x) {strsplit(x, split=';')})))
	path.rep.times <- sapply(receptor.in.path, function(x) {length(strsplit(x, split=';')[[1]])})
	cur.rep.path <- rep(names(path.rep.times), times=path.rep.times)
	plot.receptor.to.pathway.dat <- data.frame(cur.rep=cur.rep, cur.path=cur.rep.path)
	### select the top n receptors
	cur.uniq.rep <- as.vector(unique(plot.receptor.to.pathway.dat$cur.rep))
	cur.rep.LR.inten <- cluster.lr.inten(cur.uniq.rep, object.1, select.ident, ident.label, find='ligand')

	cur.rep.LR.inten <- subsetMatrix(cur.rep.LR.inten, row.select=up.in.ident)
	keep.rep <- names(which(colSums(is.na(cur.rep.LR.inten)) < nrow(cur.rep.LR.inten)))
	cur.rep.LR.inten <- subsetMatrix(cur.rep.LR.inten, col.select=keep.rep)

	cur.rep.max.LR.inten <- apply(cur.rep.LR.inten, 2, function(x){max(x, na.rm=TRUE)})
	cur.rep.max.LR.inten <- cur.rep.max.LR.inten[order(cur.rep.max.LR.inten, decreasing=TRUE)]
	
	if (length(cur.rep.max.LR.inten) < top.n.receptor){
		warning(paste0('There is(are) ', length(cur.rep.max.LR.inten),' marker receptor(s) in the selected ident associated with the selected pathways, and the input top.n.receptor is ', top.n.receptor))
		top.n.receptor <- length(cur.rep.max.LR.inten)
	}
	
	top.rep.name <- names(cur.rep.max.LR.inten)[1:top.n.receptor]
	plot.receptor.to.pathway.dat <- subset(plot.receptor.to.pathway.dat, cur.rep %in% top.rep.name)

	top.rep.LR.inten <- subsetMatrix(cur.rep.LR.inten, col.select=top.rep.name)
	plot.ident.to.receprtor.obj1.dat <- melt(top.rep.LR.inten, varnames=c('cell.from','receptor'),value.name="LR.inten",  na.rm=FALSE)
	plot.ident.to.receprtor.obj1.dat <- factor.to.character(plot.ident.to.receprtor.obj1.dat)
	
	### LR.inten for obj.2 for receptor
	cur.rep.LR.obj2.inten <- cluster.lr.inten(top.rep.name, object.2, select.ident, ident.label, find='ligand')
	cur.rep.LR.obj2.inten <- subsetMatrix(cur.rep.LR.obj2.inten, row.select=up.in.ident, col.select=top.rep.name)

	plot.ident.to.receprtor.obj2.dat <- melt(cur.rep.LR.obj2.inten, varnames=c('cell.from','receptor'),value.name="LR.inten",  na.rm=FALSE)
	plot.ident.to.receprtor.obj2.dat <- factor.to.character(plot.ident.to.receprtor.obj2.dat)

	LR.1.Dens.rep <- plot.ident.to.receprtor.obj1.dat$LR.inten
	LR.2.Dens.rep <- plot.ident.to.receprtor.obj2.dat$LR.inten

	plot.ident.to.receprtor.obj1.dat$LR.inten.diff <- as.numeric(LR.1.Dens.rep > LR.2.Dens.rep)
	plot.ident.to.receprtor.obj1.dat[which((is.na(LR.1.Dens.rep)) & (!is.na(LR.2.Dens.rep))), 'LR.inten.diff'] <- 0
	plot.ident.to.receprtor.obj1.dat[which((!is.na(LR.1.Dens.rep)) & (is.na(LR.2.Dens.rep))), 'LR.inten.diff'] <- 1
	plot.ident.to.receprtor.obj1.dat <- subset(plot.ident.to.receprtor.obj1.dat, !is.na(LR.inten))
	plot.ident.to.receprtor.dat <- plot.ident.to.receprtor.obj1.dat
	
	### data for plot of ligands and downstream clusters
	ligand.in.path <- as.character(ident.path.dat[, 'ligand.in.path'])
	if (all(is.na(ligand.in.path))){
		stop('There is no marker ligand in the selected pathways\nselect more pathways and try again')
	}
	names(ligand.in.path) <- all.sig.path
	ligand.in.path <- ligand.in.path[which(!is.na(ligand.in.path))]
	cur.lig <- as.vector(unlist(sapply(ligand.in.path, function(x) {strsplit(x, split=';')})))
	path.lig.times <- sapply(ligand.in.path, function(x) {length(strsplit(x, split=';')[[1]])})
	cur.lig.path <- rep(names(path.lig.times), times=path.lig.times)
	plot.pathway.to.ligand.dat <- data.frame(cur.lig=cur.lig, cur.path=cur.lig.path)
	### select the top n ligands
	cur.uniq.lig <- as.vector(unique(plot.pathway.to.ligand.dat$cur.lig))
	cur.lig.LR.inten <- cluster.lr.inten(cur.uniq.lig, object.1, select.ident, ident.label, find='receptor')

	cur.lig.LR.inten <- subsetMatrix(cur.lig.LR.inten, row.select=down.in.ident)
	keep.lig <- names(which(colSums(is.na(cur.lig.LR.inten)) < nrow(cur.lig.LR.inten)))
	cur.lig.LR.inten <- subsetMatrix(cur.lig.LR.inten, col.select=keep.lig)

	cur.lig.max.LR.inten <- apply(cur.lig.LR.inten, 2, function(x){max(x, na.rm=TRUE)})
	cur.lig.max.LR.inten <- cur.lig.max.LR.inten[order(cur.lig.max.LR.inten, decreasing=TRUE)]
	
	if (length(cur.lig.max.LR.inten) < top.n.ligand){
		warning(paste0('There is(are) ', length(cur.lig.max.LR.inten),' marker ligand(s) in the selected ident associated with the selected pathways, and the input top.n.ligand is ', top.n.ligand))
		top.n.ligand <- length(cur.lig.max.LR.inten)
	}
	
	top.lig.name <- names(cur.lig.max.LR.inten)[1:top.n.ligand]
	plot.pathway.to.ligand.dat <- subset(plot.pathway.to.ligand.dat, cur.lig %in% top.lig.name)
	
	top.lig.LR.inten <- subsetMatrix(cur.lig.LR.inten, col.select=top.lig.name)
	plot.ligand.to.ident.obj1.dat <- melt(top.lig.LR.inten, varnames=c('cell.to','ligand'),value.name="LR.inten",  na.rm=FALSE)
	plot.ligand.to.ident.obj1.dat <- factor.to.character(plot.ligand.to.ident.obj1.dat)

	### LR.inten for obj.2 for ligand
	cur.lig.LR.obj2.inten <- cluster.lr.inten(top.lig.name, object.2, select.ident, ident.label, find='receptor')
	cur.lig.LR.obj2.inten <- subsetMatrix(cur.lig.LR.obj2.inten, row.select=down.in.ident, col.select=top.lig.name)

	plot.ligand.to.ident.obj2.dat <- melt(cur.lig.LR.obj2.inten, varnames=c('cell.to','ligand'),value.name="LR.inten",  na.rm=FALSE)
	plot.ligand.to.ident.obj2.dat <- factor.to.character(plot.ligand.to.ident.obj2.dat)

	LR.1.Dens.lig <- plot.ligand.to.ident.obj1.dat$LR.inten
	LR.2.Dens.lig <- plot.ligand.to.ident.obj2.dat$LR.inten

	plot.ligand.to.ident.obj1.dat$LR.inten.diff <- as.numeric(LR.1.Dens.lig > LR.2.Dens.lig)
	plot.ligand.to.ident.obj1.dat[which((is.na(LR.1.Dens.lig)) & (!is.na(LR.2.Dens.lig))), 'LR.inten.diff'] <- 0
	plot.ligand.to.ident.obj1.dat[which((!is.na(LR.1.Dens.lig)) & (is.na(LR.2.Dens.lig))), 'LR.inten.diff'] <- 1
	plot.ligand.to.ident.obj1.dat <- subset(plot.ligand.to.ident.obj1.dat, !is.na(LR.inten))
	plot.ligand.to.ident.dat <- plot.ligand.to.ident.obj1.dat
	
	### pathways and their color
	path.uniq.name <- factor(all.sig.path, levels=all.sig.path)
	path.uniq.name <- path.uniq.name[order(path.uniq.name, decreasing=TRUE)]
	
	#### rect for pathway
	path.rect.length <- ident.path.dat[match(path.uniq.name, ident.path.dat$description),'diff']
	path.rect.stat <- ident.path.dat[match(path.uniq.name, ident.path.dat$description),'stat']
	path.rect.pval <-  ident.path.dat[match(path.uniq.name, ident.path.dat$description),'P.val.adj']

	if (all(path.rect.pval==0)){
		path.rect.pval <- path.rect.stat
		warning('All adjusted p values for the selected pathways are 0. The colors of bars representing pathways are adjusted to indicate the t values (in t test) or W values (in wilcox test) from difference tests')
	}else{
		path.rect.pval <- -log10(path.rect.pval)
		path.rect.pval <- p.remove.inf(path.rect.pval)
	}
	bar.path.col <- LRcolor(path.rect.pval, user.set.col=bar.pathway.col)

	
	### the receptor and their coordinate, and their color
	cur.uniq.rep <- factor(top.rep.name, levels=top.rep.name)
	cur.uniq.rep <- cur.uniq.rep[order(cur.uniq.rep, decreasing=TRUE)]
	n.rep <- length(cur.uniq.rep)
	
	top.markerR.dat <- subset(object.1@interact$markerR, subset=(cluster==select.ident & gene %in% top.rep.name))
	fc.rep <- top.markerR.dat[match(cur.uniq.rep, top.markerR.dat$gene),'avg_log2FC']
	rep.size <- 10 * dot.gene.size * fc.rep
	rep.pct <-  top.markerR.dat[match(cur.uniq.rep, top.markerR.dat$gene),'pct.1']
	rep.pval <-  top.markerR.dat[match(cur.uniq.rep, top.markerR.dat$gene),'p_val_adj']

	if (all(rep.pval==0)){
		rep.pval <- fc.rep
		warning('All adjusted p values for the top.n.receptors are 0. The colors of dots representing receptors are adjusted to indicate the average log2FC of receptors')
	}else{
		rep.pval <- -log10(rep.pval)
		rep.pval <- p.remove.inf(rep.pval)
	}
	dot.rep.col <- LRcolor(rep.pval, user.set.col=dot.gene.col)

	### the ligand and their coordinate, and their color
	cur.uniq.lig <- factor(top.lig.name, levels=top.lig.name)
	cur.uniq.lig <- cur.uniq.lig[order(cur.uniq.lig, decreasing=TRUE)]
	n.lig <- length(cur.uniq.lig)
	
	top.markerL.dat <- subset(object.1@interact$markerL, subset=(cluster==select.ident & gene %in% top.lig.name))
	fc.lig <- top.markerL.dat[match(cur.uniq.lig, top.markerL.dat$gene),'avg_log2FC']
	lig.size <- 10 * dot.gene.size * fc.lig
	lig.pct <- top.markerL.dat[match(cur.uniq.lig, top.markerL.dat$gene),'pct.1']
	lig.pval <- top.markerL.dat[match(cur.uniq.lig, top.markerL.dat$gene),'p_val_adj']

	if (all(lig.pval==0)){
		lig.pval <- fc.lig
		warning('All adjusted p values for the top.n.ligands are 0. The colors of dots representing ligands are adjusted to indicate the average log2FC of ligands')
	}else{
		lig.pval <- -log10(lig.pval)
		lig.pval <- p.remove.inf(lig.pval)
	}
	dot.lig.col <- LRcolor(lig.pval, user.set.col=dot.gene.col)

	### the upstream and downstream ident and their coordinate, and their color
	up.uniq.ident <- factor(up.in.ident, levels=up.in.ident)
	up.uniq.ident <- up.uniq.ident[order(up.uniq.ident, decreasing=TRUE)]
	n.up.ident <- length(up.uniq.ident)
	### color for dots representing ident
	if (is.null(dot.ident.col)){
		dot.ident.col <- rev(scales::hue_pal(c=100)(length(ident.label)))
		names(dot.ident.col) <- ident.label
	}else{
		dot.ident.names.col <- names(dot.ident.col)
		if (is.null(dot.ident.names.col)){
			stop(paste0('Wrong dot.ident.col parameter! The colors should be named according to their clusters'))
		}else{
			all.in.ident <- unique(c(up.in.ident, down.in.ident))
			all.missed.ident <- all.in.ident[which(!(all.in.ident %in% dot.ident.names.col))]
			if (length(all.missed.ident) > 0){ stop(paste0('Wrong dot.ident.col parameter! Please add colors for these clusters: ', pasteIdent(all.missed.ident), '.')) }
		}

	}
	dot.up.ident.col <- rev(dot.ident.col[up.in.ident])
	### the downstream ident and their coordinate, and their color
	down.uniq.ident <- factor(down.in.ident, levels=down.in.ident)
	down.uniq.ident <- down.uniq.ident[order(down.uniq.ident, decreasing=TRUE)]
	n.down.ident <- length(down.uniq.ident)
	dot.down.ident.col <- rev(dot.ident.col[down.in.ident])

	### adjust the coordinate of each element
	max.limit <- max(n.up.ident, n.rep, top.n.path, n.lig, n.down.ident)
	up.ident.new.coor <- seq(from=1, to=max.limit, length.out=n.up.ident)
	rep.new.coor <- seq(from=1, to=max.limit, length.out=n.rep)
	pathway.new.coor <- seq(from=1, to=max.limit, length.out=top.n.path)
	lig.new.coor <- seq(from=1, to=max.limit, length.out=n.lig)
	down.ident.new.coor <- seq(from=1, to=max.limit, length.out=n.down.ident)

	if (n.up.ident==1){ up.ident.new.coor <- median(1:max.limit) }
	if (n.rep==1){ rep.new.coor <- median(1:max.limit) }
	if (top.n.path==1){ pathway.new.coor <- median(1:max.limit) }
	if (n.lig==1){ lig.new.coor <- median(1:max.limit) }
	if (n.down.ident==1){ down.ident.new.coor <- median(1:max.limit) }

	if (length(pathway.new.coor)==1){
		bar_pathway_width <- 0.3
	}else{
		bar_pathway_width <- 0.3 * (pathway.new.coor[2]-pathway.new.coor[1]) * top.n.path / max.limit
	}
	path.name.max.width <- max(strwidth(path.uniq.name, units="inches", cex=1))
	bar.pathway.length <- 2 * path.rect.length / max(path.rect.length)
	bar_pathway <- geom_rect(aes(xmin=6, xmax=6+bar.pathway.length, ymin=pathway.new.coor-bar_pathway_width, ymax=pathway.new.coor+bar_pathway_width), size=1, fill=bar.path.col, alpha=1, show.legend=F)
	point_up_ident <- geom_point(aes(x=1,y=up.ident.new.coor),size=8*dot.ident.size,color=dot.up.ident.col,shape=16) 
	point_receptor <- geom_point(aes(x=2,y=rep.new.coor),size=rep.size,color=dot.rep.col)
	point_pathway <- geom_point(aes(x=3,y=pathway.new.coor),size=8*dot.ident.size,color='black',shape=21,fill='white')
	point_ligand <- geom_point(aes(x=4,y=lig.new.coor),size=lig.size,color=dot.lig.col) 
	point_down_ident <- geom_point(aes(x=5,y=down.ident.new.coor),size=8*dot.ident.size,color=dot.down.ident.col,shape=16) 

	### lines between elements
	### lines between up idents and receptors
	line.ident.to.rep.y.coor <- up.ident.new.coor[match(plot.ident.to.receprtor.dat$cell.from, up.uniq.ident)]
	line.ident.to.rep.yend.coor <- rep.new.coor[match(plot.ident.to.receprtor.dat$receptor, cur.uniq.rep)]
	line.ident.to.rep.col <- ifelse(plot.ident.to.receprtor.dat$LR.inten.diff > 0, 'red', 'blue')
	line.ident.to.lig.width <- LRinten.to.width(plot.ident.to.receprtor.dat$LR.inten)
	line_ident_to_receptor <- geom_segment(aes(x=1, xend=2, y=line.ident.to.rep.y.coor,yend=line.ident.to.rep.yend.coor), size=line.ident.to.lig.width, color=line.ident.to.rep.col, show.legend=F, alpha=0.6)
	### lines between receptors and pathways
	line.rep.to.path.y.coor <- rep.new.coor[match(plot.receptor.to.pathway.dat$cur.rep, cur.uniq.rep)]
	line.rep.to.path.yend.coor <- pathway.new.coor[match(plot.receptor.to.pathway.dat$cur.path, path.uniq.name)]
	path.rep.line.col <- dot.rep.col[match(plot.receptor.to.pathway.dat$cur.rep, cur.uniq.rep)]
	line_receptor_to_pathway <- geom_segment(aes(x=2, xend=3, y=line.rep.to.path.y.coor,yend=line.rep.to.path.yend.coor),size=1, color=path.rep.line.col, show.legend=F)
	### lines between pathways and ligands
	line.path.to.lig.y.coor <- pathway.new.coor[match(plot.pathway.to.ligand.dat$cur.path, path.uniq.name)]
	line.path.to.lig.yend.coor <- lig.new.coor[match(plot.pathway.to.ligand.dat$cur.lig, cur.uniq.lig)]
	path.lig.line.col <- dot.lig.col[match(plot.pathway.to.ligand.dat$cur.lig, cur.uniq.lig)]
	line_pathway_to_ligand <- geom_segment(aes(x=3, xend=4, y=line.path.to.lig.y.coor,yend=line.path.to.lig.yend.coor),size=1, color=path.lig.line.col, show.legend=F)
	### lines between ligands and down idents
	line.lig.to.ident.y.coor <- lig.new.coor[match(plot.ligand.to.ident.dat$ligand, cur.uniq.lig)]
	line.lig.to.ident.yend.coor <- down.ident.new.coor[match(plot.ligand.to.ident.dat$cell.to, down.uniq.ident)]
	line.lig.to.ident.col <- ifelse(plot.ligand.to.ident.dat$LR.inten.diff > 0, 'red', 'blue')
	line.lig.to.ident.width <- LRinten.to.width(plot.ligand.to.ident.dat$LR.inten)
	line_ligand_to_ident <- geom_segment(aes(x=4, xend=5, y=line.lig.to.ident.y.coor,yend=line.lig.to.ident.yend.coor), size=line.lig.to.ident.width, color=line.lig.to.ident.col, show.legend=F, alpha=0.6)
	
	### text and label
	label_up_ident <- annotate('text', hjust=0.5, x=1, y=up.ident.new.coor-0.4,label=up.uniq.ident, size=5*label.text.size)
	label_receptor <- annotate('text', hjust=0.5, x=2, y=rep.new.coor-0.4,label=cur.uniq.rep, size=5*label.text.size)
	#text_pathway <- geom_text(aes(x=6.1+path.rect.length, y=pathway.new.coor), hjust=0, label=path.uniq.name, size=5*label.text.size)
	text_pathway <- geom_text(aes(x=6, y=pathway.new.coor), hjust=0, label=path.uniq.name, size=5*label.text.size)
	label_pathway_circle <- annotate('text', hjust=0.5, x=3, y=pathway.new.coor, label=length(path.uniq.name):1, size=5*label.text.size)
	label_pathway_bar <- annotate('text', hjust=1, x=6, y=pathway.new.coor, label=length(path.uniq.name):1, size=5*label.text.size)
	label_ligand <- annotate('text', hjust=0.5, x=4, y=lig.new.coor-0.4,label=cur.uniq.lig, size=5*label.text.size)
	label_down_ident <- annotate('text', hjust=0.5, x=5, y=down.ident.new.coor-0.4,label=down.uniq.ident, size=5*label.text.size)
	label_title <- annotate('text', hjust=c(0.5,0.5,0.5,0.5,0.5,0), x=1:6, y=max.limit+0.8,label=c('Upstream', 'Receptor', 'Pathway', 'Ligand', 'Downstream','Pathway annotation'), size=5*label.title.size) 

	pathplot.theme <- theme(axis.title=element_blank(), axis.text=element_blank(),axis.ticks=element_blank(), axis.line=element_blank(), panel.background=element_rect(fill="white"))

	# the max limit of x axis is dependent on the longest pathway
	# this limit may need to be adjusted
	plot <- ggplot() + 
	line_pathway_to_ligand +
	line_ident_to_receptor +
	line_receptor_to_pathway + 
	line_ligand_to_ident +
	
	point_up_ident +
	point_receptor + 
	point_pathway +
	point_ligand +
	point_down_ident +
	bar_pathway +
	text_pathway +
	
	label_up_ident +
	label_receptor +
	label_ligand +
	label_down_ident +
	label_title +
	label_pathway_circle +
	label_pathway_bar +

	scale_x_continuous(limits=c(0,6+max(path.name.max.width,bar.pathway.length))) +
	pathplot.theme
	return(plot)
}

#' To present a net plot for pathways and upstream or downstream associated LR interactions
#' @param object CommPath object
#' @param select.ident Cluster of interest
#' @param plot To show the network of activated pathways and the upstream or downstream LR interactions; default is 'upstream'
#' @param ident.col Vector of colors of each cluster; names of the col vector are supposed to be assigned to indicate each color for each cluster
#' @param vert.size.attr Which statistical measures should be mapped to the size of vertex; Type "getNetAttr(object)" to retrieve the available statistical measures for pathways
#' @param vert.size.LR Size of node representing LR interactions
#' @param vert.size.path.adj pseudocount to adjust the size of node representing pathways
#' @param top.n.path Top n activated pathways in the selected cluster of interest; if NULL, slelect all activated pathways
#' @param path.order Sort criteria used to select the top n pathways, either "P.val" or "P.val.adj", which represent the original and adjusted p values, or "diff" which represents the mean (in t test) or median (in wilcox test) difference; default is "P.val.adj". This parameter would be masked if top.n.path is set as NULL.
#' @param select.path Vector of pathways for which users want to present a chain plot; default is all, which means to plot all pathways. The parameters top.n.path and path.order would be masked if this one is set with specific pathways.
#' @param layout Network layout; defualt is Fruchterman-Reingold ('layout.fruchterman.reingold'). Refer to the igraph package for more information. Note: an error would be raised if the selected layout do not apply to the graph.
#' @param LR.label To show which label on the node representing LR interactions. Available options: LR (to show both the ligand and receptor label), L (to show only the ligand label), R (to show only the receptor label), or NA (no label); default is NA
#' @param pathway.label Logical value indicating to display the label of pathways or not; default is TRUE
#' @param edge.arrow.size Size of arrow
#' @param edge.color Color of edge
#' @param vertex.label.cex.LR Font size of LR label
#' @param vertex.label.cex.path Font size of pathway label
#' @param vertex.label.color color of vertex label
#' @param vertex.label.family Font family of the label
#' @param vertex.frame.color Color of Node border
#' @param node.pie Logical value indicating to show the node representing pathways in pie charts or not; default is TRUE
#' @param return.net Logical value indicating to return the igraph object or not; default is FALSE
#' @param return.data Logical value indicating to return the link dataframe or not; default is FALSE. This parameter allows users to export the link data in a dataframe, which could be analyzed in other network visualization tools, for example Cytoscape.
#' @importFrom igraph graph_from_data_frame V E degree layout_with_fr
#' @importFrom plyr rbind.fill mapvalues
#' @return Net plot showing the activated pathways and associated LR interactions
#' @export
pathNetPlot <- function(object, select.ident, plot='upstream', ident.col=NULL, vert.size.attr="degree", vert.size.LR=0.5, vert.size.path.adj=5, top.n.path=NULL, path.order="P.val.adj", select.path='all', layout="layout.auto", LR.label='LR', pathway.label=TRUE, edge.arrow.size=0.2, edge.color='gray75', vertex.label.cex.LR=0.5,vertex.label.cex.path=0.5, vertex.label.color="black", vertex.label.family="Helvetica", vertex.frame.color="#ffffff", node.pie=TRUE, return.net=FALSE, return.data=FALSE){
	if(plot=='upstream'){
		pathNetPlot.upstream(object=object, select.ident=select.ident, ident.col=ident.col, vert.size.attr=vert.size.attr, vert.size.LR=vert.size.LR, vert.size.path.adj=vert.size.path.adj, top.n.path=top.n.path, select.path=select.path, path.order=path.order, layout=layout, LR.label=LR.label, pathway.label=pathway.label, edge.arrow.size=edge.arrow.size, edge.color=edge.color, vertex.label.cex.LR=vertex.label.cex.LR, vertex.label.cex.path=vertex.label.cex.path, vertex.label.color=vertex.label.color, vertex.label.family=vertex.label.family, vertex.frame.color=vertex.frame.color, node.pie=node.pie, return.net=return.net, return.data=return.data)
	}else if(plot=='downstream'){
		pathNetPlot.downstream(object=object, select.ident=select.ident, ident.col=ident.col, vert.size.attr=vert.size.attr, vert.size.LR=vert.size.LR, vert.size.path.adj=vert.size.path.adj, top.n.path=top.n.path, select.path=select.path, path.order=path.order, layout=layout, LR.label=LR.label, pathway.label=pathway.label, edge.arrow.size=edge.arrow.size, edge.color=edge.color, vertex.label.cex.LR=vertex.label.cex.LR, vertex.label.cex.path=vertex.label.cex.path, vertex.label.color=vertex.label.color, vertex.label.family=vertex.label.family, vertex.frame.color=vertex.frame.color, node.pie=node.pie, return.net=return.net, return.data=return.data)
	}else{
		stop('Wrong "plot" parameter! Available options: upstream (to show the network of activated pathways and the upstream LR interactions) or downstream (to show the network of activated pathways and the downstream LR interactions)')
	}
}

#' To present a net plot for pathways and upstream associated LR interactions
#' @param object CommPath object
#' @param select.ident Cluster of interest
#' @param ident.col Vector of colors of each cluster; names of the col vector are supposed to be assigned to indicate each color for each cluster
#' @param vert.size.attr Which statistical measures should be mapped to the size of vertex; Type "getNetAttr(object)" to retrieve the available statistical measures for pathways
#' @param vert.size.LR Size of node representing LR interactions
#' @param vert.size.path.adj pseudocount to adjust the size of node representing pathways
#' @param top.n.path Top n activated pathways in the selected cluster of interest; if NULL, slelect all activated pathways
#' @param path.order Sort criteria used to select the top n pathways, either "P.val" or "P.val.adj", which represent the original and adjusted p values, or "diff" which represents the mean (in t test) or median (in wilcox test) difference; default is "P.val.adj". This parameter would be masked if top.n.path is set as NULL.
#' @param select.path Vector of pathways for which users want to present a chain plot; default is all, which means to plot all pathways. The parameters top.n.path and path.order would be masked if this one is set with specific pathways.
#' @param layout Network layout; defualt is Fruchterman-Reingold ('layout.fruchterman.reingold'). Refer to the igraph package for more information. Note: an error would be raised if the selected layout do not apply to the graph.
#' @param LR.label To show which label on the node representing LR interactions. Available options: LR (to show both the ligand and receptor label), L (to show only the ligand label), R (to show only the receptor label), or NA (no label); default is NA
#' @param pathway.label Logical value indicating to display the label of pathways or not; default is TRUE
#' @param edge.arrow.size Size of arrow
#' @param edge.color Color of edge
#' @param vertex.label.cex.LR Font size of LR label
#' @param vertex.label.cex.path Font size of pathway label
#' @param vertex.label.color color of vertex label
#' @param vertex.label.family Font family of the label
#' @param vertex.frame.color Color of Node border
#' @param node.pie Logical value indicating to show the node representing pathways in pie charts or not; default is TRUE
#' @param return.net Logical value indicating to return the igraph object or not; default is FALSE
#' @param return.data Logical value indicating to return the link dataframe or not; default is FALSE. This parameter allows users to export the link data in a dataframe, which could be analyzed in other network visualization tools, for example Cytoscape.
#' @importFrom igraph graph_from_data_frame V E degree layout_with_fr
#' @importFrom plyr rbind.fill mapvalues
#' @return Net plot showing the activated pathways and associated LR interactions
#' @export
pathNetPlot.upstream <- function(object, select.ident, ident.col=NULL, vert.size.attr="degree", vert.size.LR=0.5, vert.size.path.adj=5, top.n.path=NULL, path.order="P.val.adj", select.path='all', layout="layout.auto", LR.label='LR', pathway.label=TRUE, edge.arrow.size=0.2,  edge.color='gray75', vertex.label.cex.LR=0.5, vertex.label.cex.path=0.5,vertex.label.color="black", vertex.label.family="Helvetica", vertex.frame.color="#ffffff", node.pie=TRUE, return.net=FALSE, return.data=FALSE){
	options(stringsAsFactors=F)

	path.net.dat <- object@pathway.net$upstream
	if (!(select.ident %in% path.net.dat$cell.to)){
		stop(paste0('There is no significantly up-regulatged pathways for cluster ', select.ident, ' or the up-regulatged pathways contain no marker receptor'))
	}
	path.curcell.dat <- path.net.dat[which(path.net.dat$cell.to==select.ident), ]

	### node of pathway
	if('t.path' %in% colnames(path.curcell.dat)){
		node.path <- unique(path.curcell.dat[, c('path.contain.rep.unfold','description.path','mean.diff.path', 't.path', 'P.val.path', 'P.val.adj.path')])
		colnames(node.path) <- c('ID','description.path','mean.diff.path', 't.path', 'P.val.path', 'P.val.adj.path')
		node.path$stat <- node.path$t.path
		node.path$diff <- node.path$mean.diff
	}else{
		node.path <- unique(path.curcell.dat[, c('path.contain.rep.unfold','description.path','mean.diff.path', 'W.path', 'P.val.path', 'P.val.adj.path')])
		colnames(node.path) <- c('ID','description.path','mean.diff.path', 'W.path', 'P.val.path', 'P.val.adj.path')
		node.path$stat <- node.path$W.path
		node.path$diff <- node.path$median.diff
	}
	node.path$type <- 'Pathway'

	### select the user-defined pathways, with the select.path parameter
	if (length(select.path) > 1 | (length(select.path)==1 & select.path[1]!='all')){
		path.num.plot <- which(node.path$ID %in% select.path)
		if (length(path.num.plot) > 0){
			node.path <- node.path[path.num.plot, ]
		}else{
			stop('All pathways selected are not up-regulatged in the selected cluster. Set return.data as TRUE to retrieve all up-regulatged pathways.')
		}

		ident.missed <- select.path[which(!(select.path %in% node.path$ID))]
		if (length(ident.missed) > 0){
			warning(paste0('These pathways are not up-regulatged in the selected cluster: ', pasteIdent(ident.missed), '. Set return.data as TRUE to retrieve all up-regulatged pathways.'))
		}
	}else{
		### select the top n pathways
		if (is.null(top.n.path)){
			top.n.path <- nrow(node.path)
		}else if (is.numeric(top.n.path) & (top.n.path%%1==0)){
			if (path.order=='P.val.adj'){
				node.path <- node.path[order(node.path[,'P.val.adj.path'], -node.path[,'stat']),]
			}else if(path.order=='P.val'){
				node.path <- node.path[order(node.path[,'P.val.path'], -node.path[,'stat']),]
			}else if(path.order=='diff'){
				node.path <- node.path[order(-node.path[,'diff'], node.path[,'P.val.adj.path']),]
			}

			if (top.n.path > nrow(node.path)){
				warning(paste0('There is(are) ', nrow(node.path),' significant pathway(s) for the selected ident, and the input top.n.path is ', top.n.path))
				top.n.path <- nrow(node.path)
			}
		}else{
			stop('Wrong "top.n.path" parameter! Input an positive integer to select the corresponding number of up-regulatged pathway(s)')
		}
		node.path <- node.path[1:top.n.path, ]
	}
	node.path$stat <- NULL
	path.curcell.dat <- path.curcell.dat[which(path.curcell.dat$path.contain.rep.unfold %in% node.path$ID), ]

	### node of LR
	node.LR <- unique(path.curcell.dat[, c('LR.type','LR.pair','cell.from','cell.to','ligand','receptor','log2FC.LR', 'P.val.LR', 'P.val.adj.LR')])	
	colnames(node.LR) <- c('ID','LR.pair','cell.from','cell.to', 'ligand', 'receptor', 'log2FC.LR', 'P.val.LR', 'P.val.adj.LR')
	node.LR$type <- 'LR'
	node.LR <- node.LR[which(!duplicated(node.LR$ID)),]

	### get the node of all vertex and construct the net
	node.dat <- rbind.fill(node.LR, node.path)
	link.dat <- path.curcell.dat[, c('LR.type', 'path.contain.rep.unfold')]
	net <- graph_from_data_frame(d=link.dat, vertices=node.dat, directed=T)
	### return the igraph object if ordered by user
	if(return.net){
		return(net)
	}

	if(return.data){
		if('t.path' %in% colnames(path.curcell.dat)){
			return(path.curcell.dat[,c('cell.from', 'cell.to', 'ligand', 'receptor', 'log2FC.LR', 'P.val.LR', 'P.val.adj.LR','description.path','mean.diff.path', 'mean.1.path', 'mean.2.path', 't.path', 'df.path', 'P.val.path', 'P.val.adj.path')])
		}else{
			return(path.curcell.dat[,c('cell.from', 'cell.to', 'ligand', 'receptor', 'log2FC.LR', 'P.val.LR', 'P.val.adj.LR','description.path','median.diff.path', 'median.1.path', 'median.2.path', 'W.path', 'P.val.path', 'P.val.adj.path')])
		}
	}

	### network plot
	vert.attr.list <- get.vertex.attribute(net)
	node.dat.id <- vert.attr.list$name
	### vertex label
	if(is.na(LR.label)){
		vert.label <- rep(NA,length(node.dat.id))
	}else if(LR.label=='L'){
		vert.label <- vert.attr.list$ligand
	}else if(LR.label=='R'){
		vert.label <- vert.attr.list$receptor
	}else if(LR.label=='LR'){
		vert.label <- vert.attr.list$LR.pair
	}else{
		stop('Wrong LR.label parameter! Available options: LR (to show both the ligand and receptor label), L (to show only the ligand label), R (to show only the receptor label), or NA (no label) ')
	}

	if(pathway.label){
		if(all(is.na(vert.label))){
			vert.label <- vert.attr.list$description.path
		}else{
			vert.label[which(is.na(vert.label))] <- node.dat.id[which(is.na(vert.label))]
		}
	}

	### vertex size
	if(vert.size.attr=='degree'){
		vert.size <- igraph::degree(net, mode="in")
		vert.size[which(vert.size==0)] <- NA
	}else{
		if(vert.size.attr=='logP'){
			vert.size <- vert.attr.list[['P.val.path']]
			vert.size <- -log10(vert.size)
		}else if(vert.size.attr=='logP.adj'){
				vert.size <- vert.attr.list[['P.val.adj.path']]
				vert.size <- -log10(vert.size)
		}else{
			vert.size <- vert.attr.list[[paste0(vert.size.attr,'.path')]]
			if(is.null(vert.size)){
				stop(paste0('Wrong "vert.size.attr" parameter! There is no valid attribute named as ',vert.size.attr,'\nType "getPathAttr(object)" for available attribute'))
			}
		}
	}
	vert.size <- scale_1(vert.size)
	vert.size <- vert.size * vert.size.path.adj
	vert.size[which(is.na(vert.size))] <- vert.size.LR
	
	### vertex color
	all.ident <- object@cell.info$Cluster
	if (!is.factor(all.ident)){ all.ident <- factor(all.ident) }
	ident.label <- levels(all.ident)
	vert.cluster.ident <- ident.label[which(ident.label %in% as.character(vert.attr.list$cell.from))]
	if(is.null(ident.col)){
		scales.color <- scales::hue_pal(c=100)(length(ident.label))
		vert.cluster.color <- scales.color[which(ident.label %in% as.character(vert.attr.list$cell.from))]
	}else{
		col.name <- names(ident.col)
		if(is.null(col.name)){
			stop(paste0('Wrong "ident.col" parameter! The colors should be named according to the clusters'))
		}else{
			ident.missed <- vert.cluster.ident[which(!(vert.cluster.ident %in% col.name))]
			if(length(ident.missed) > 0){
				stop(paste0('The color of ',pasteIdent(ident.missed),' may be missed in the input color'))
			}
			vert.cluster.color <- as.vector(ident.col[vert.cluster.ident])
		}
	}
	vertex.pie.color <- list(vert.cluster.color)
	### pie border
	vertex.frame.na.color <- mapvalues(vert.attr.list$type, from=c('LR','Pathway'), to=c(NA,vertex.frame.color), warn_missing=FALSE) 
	
	### pie chart composition
	pie.comp.path <- as.data.frame.array(table(path.curcell.dat[,c('cell.from','path.contain.rep.unfold')]))
	pie.comp.path$rownames <- rownames(pie.comp.path)
	pie.comp.path <- pie.comp.path[vert.cluster.ident, ]
	pie.comp.path$rownames <- NULL
	pie.comp.node <- as.data.frame.array(table(node.dat[,c('cell.from','ID')]))
	pie.comp.node <- pie.comp.node[vert.cluster.ident, node.dat$ID]

	for(each.path in colnames(pie.comp.path)){
		pie.comp.node[ ,each.path] <- pie.comp.path[ ,each.path] 
	}
	pie.comp.node <- lapply(pie.comp.node, function(x){x})

	if(!node.pie){
		scales.color <- scales::hue_pal(c=100)(length(ident.label))
		vert.cluster.color <- scales.color[which(ident.label %in% as.character(vert.attr.list$cell.from))]
		vert.color <- mapvalues(vert.attr.list$cell.from, from=c(vert.cluster.ident,NA), to=c(vert.cluster.color,'tomato'), warn_missing=FALSE) 
	}

	### layout
	l <- do.call(layout, list(net))
	vertex.label.cex <- as.numeric(mapvalues(vert.attr.list$type, from=c('LR','Pathway'), to=c(vertex.label.cex.LR, vertex.label.cex.path), warn_missing=FALSE))
	if(node.pie){
		plot(net, vertex.shape="pie", vertex.pie=pie.comp.node, edge.arrow.size=edge.arrow.size, edge.color=edge.color, layout=l, vertex.size=vert.size, vertex.label=vert.label, vertex.label.color=vertex.label.color, vertex.label.cex=vertex.label.cex, vertex.label.family=vertex.label.family,vertex.frame.color=vertex.frame.na.color, vertex.pie.color=vertex.pie.color)
		legend(x=1.1, y=1, legend=vert.cluster.ident, pch=21, col=vertex.frame.color, pt.bg=vert.cluster.color, pt.cex=1, cex=.5, bty="n", ncol=1)
	}else{
		plot(net, edge.arrow.size=edge.arrow.size, edge.color=edge.color, layout=l, vertex.size=vert.size, vertex.label=vert.label, vertex.label.color=vertex.label.color, vertex.label.cex=vertex.label.cex, vertex.label.family=vertex.label.family,vertex.frame.color=vertex.frame.color, vertex.color=vert.color)
		legend(x=1.1, y=1, legend=vert.cluster.ident, pch=21, col=vertex.frame.color, pt.bg=vert.cluster.color, pt.cex=1, cex=.5, bty="n", ncol=1)
	}
}

#' To present a net plot for pathways and downstream associated LR interactions
#' @param object CommPath object
#' @param select.ident Cluster of interest
#' @param ident.col Vector of colors of each cluster; names of the col vector are supposed to be assigned to indicate each color for each cluster
#' @param vert.size.attr Which statistical measures should be mapped to the size of vertex; Type "getNetAttr(object)" to retrieve the available statistical measures for pathways
#' @param vert.size.LR Size of node representing LR interactions
#' @param vert.size.path.adj pseudocount to adjust the size of node representing pathways
#' @param top.n.path Top n activated pathways in the selected cluster of interest; if NULL, slelect all activated pathways
#' @param path.order Sort criteria used to select the top n pathways, either "P.val" or "P.val.adj", which represent the original and adjusted p values, or "diff" which represents the mean (in t test) or median (in wilcox test) difference; default is "P.val.adj". This parameter would be masked if top.n.path is set as NULL.
#' @param select.path Vector of pathways for which users want to present a chain plot; default is all, which means to plot all pathways. The parameters top.n.path and path.order would be masked if this one is set with specific pathways.
#' @param layout Network layout; defualt is Fruchterman-Reingold ('layout.fruchterman.reingold'). Refer to the igraph package for more information. Note: an error would be raised if the selected layout do not apply to the graph.
#' @param LR.label To show which label on the node representing LR interactions. Available options: LR (to show both the ligand and receptor label), L (to show only the ligand label), R (to show only the receptor label), or NA (no label); default is NA
#' @param pathway.label Logical value indicating to display the label of pathways or not; default is TRUE
#' @param edge.arrow.size Size of arrow
#' @param edge.color Color of edge
#' @param vertex.label.cex.LR Font size of LR label
#' @param vertex.label.cex.path Font size of pathway label
#' @param vertex.label.color color of vertex label
#' @param vertex.label.family Font family of the label
#' @param vertex.frame.color Color of Node border
#' @param node.pie Logical value indicating to show the node representing pathways in pie charts or not; default is TRUE
#' @param return.net Logical value indicating to return the igraph object or not; default is FALSE
#' @param return.data Logical value indicating to return the link dataframe or not; default is FALSE. This parameter allows users to export the link data in a dataframe, which could be analyzed in other network visualization tools, for example Cytoscape.
#' @importFrom igraph graph_from_data_frame V E degree layout_with_fr
#' @importFrom plyr rbind.fill mapvalues
#' @return Net plot showing the activated pathways and associated LR interactions
#' @export
pathNetPlot.downstream <- function(object, select.ident, ident.col=NULL, vert.size.attr="degree", vert.size.LR=0.5, vert.size.path.adj=5, top.n.path=NULL, path.order="P.val.adj", select.path='all', layout="layout.auto", LR.label='LR', pathway.label=TRUE, edge.arrow.size=0.2, edge.color='gray75', vertex.label.cex.LR=0.5, vertex.label.cex.path=0.5, vertex.label.color="black", vertex.label.family="Helvetica", vertex.frame.color="#ffffff", node.pie=TRUE, return.net=FALSE, return.data=FALSE){
	options(stringsAsFactors=F)

	path.net.dat <- object@pathway.net$downstream
	if (!(select.ident %in% path.net.dat$cell.from)){
		stop(paste0('There is no significantly up-regulatged pathways for cluster ', select.ident, ' or the up-regulatged pathways contain no marker ligand'))
	}
	path.curcell.dat <- path.net.dat[which(path.net.dat$cell.from==select.ident), ]

	### node of pathway
	if('t.path' %in% colnames(path.curcell.dat)){
		node.path <- unique(path.curcell.dat[, c('path.contain.lig.unfold','description.path','mean.diff.path', 't.path', 'P.val.path', 'P.val.adj.path')])
		colnames(node.path) <- c('ID','description.path','mean.diff.path', 't.path', 'P.val.path', 'P.val.adj.path')
		node.path$stat <- node.path$t.path
		node.path$diff <- node.path$mean.diff
	}else{
		node.path <- unique(path.curcell.dat[, c('path.contain.lig.unfold','description.path','mean.diff.path', 'W.path', 'P.val.path', 'P.val.adj.path')])
		colnames(node.path) <- c('ID','description.path','mean.diff.path', 'W.path', 'P.val.path', 'P.val.adj.path')
		node.path$stat <- node.path$W.path
		node.path$diff <- node.path$median.diff
	}
	node.path$type <- 'Pathway'

	### select the user-defined pathways, with the select.path parameter
	if (length(select.path) > 1 | (length(select.path)==1 & select.path[1]!='all')){
		path.num.plot <- which(node.path$ID %in% select.path)
		if (length(path.num.plot) > 0){
			node.path <- node.path[path.num.plot, ]
		}else{
			stop('All pathways selected are not up-regulatged in the selected cluster. Set return.data as TRUE to retrieve all up-regulatged pathways.')
		}

		ident.missed <- select.path[which(!(select.path %in% node.path$ID))]
		if (length(ident.missed) > 0){
			warning(paste0('These pathways are not up-regulatged in the selected cluster: ', pasteIdent(ident.missed), '. Set return.data as TRUE to retrieve all up-regulatged pathways.'))
		}
	}else{
		### select the top n pathways
		if (is.null(top.n.path)){
			top.n.path <- nrow(node.path)
		}else if (is.numeric(top.n.path) & (top.n.path%%1==0)){
			if (path.order=='P.val.adj'){
				node.path <- node.path[order(node.path[,'P.val.adj.path'], -node.path[,'stat']),]
			}else if(path.order=='P.val'){
				node.path <- node.path[order(node.path[,'P.val.path'], -node.path[,'stat']),]
			}else if(path.order=='diff'){
				node.path <- node.path[order(-node.path[,'diff'], node.path[,'P.val.adj.path']),]
			}

			if (top.n.path > nrow(node.path)){
				warning(paste0('There is(are) ', nrow(node.path),' significant pathway(s) for the selected ident, and the input top.n.path is ', top.n.path))
				top.n.path <- nrow(node.path)
			}
		}else{
			stop('Wrong "top.n.path" parameter! Input an positive integer to select the corresponding number of up-regulatged pathway(s)')
		}	
		node.path <- node.path[1:top.n.path, ]
	}
	node.path$stat <- NULL
	path.curcell.dat <- path.curcell.dat[which(path.curcell.dat$path.contain.lig.unfold %in% node.path$ID), ]

	### node of LR
	node.LR <- unique(path.curcell.dat[, c('LR.type','LR.pair','cell.from','cell.to','ligand','receptor','log2FC.LR', 'P.val.LR', 'P.val.adj.LR')])	
	colnames(node.LR) <- c('ID','LR.pair','cell.from','cell.to','ligand', 'receptor', 'log2FC.LR', 'P.val.LR', 'P.val.adj.LR')
	node.LR$type <- 'LR'
	node.LR <- node.LR[which(!duplicated(node.LR$ID)),]

	### get the node of all vertex and construct the net
	node.dat <- rbind.fill(node.LR, node.path)
	link.dat <- path.curcell.dat[, c('path.contain.lig.unfold','LR.type')]
	net <- graph_from_data_frame(d=link.dat, vertices=node.dat, directed=T)
	### return the igraph object if ordered by user
	if(return.net){
		return(net)
	}

	if(return.data){
		if('t.path' %in% colnames(path.curcell.dat)){
			return(path.curcell.dat[,c('cell.from', 'cell.to', 'ligand', 'receptor', 'log2FC.LR', 'P.val.LR', 'P.val.adj.LR','description.path','mean.diff.path', 'mean.1.path', 'mean.2.path', 't.path', 'df.path', 'P.val.path', 'P.val.adj.path')])
		}else{
			return(path.curcell.dat[,c('cell.from', 'cell.to', 'ligand', 'receptor', 'log2FC.LR', 'P.val.LR', 'P.val.adj.LR','description.path','median.diff.path', 'median.1.path', 'median.2.path', 'W.path', 'P.val.path', 'P.val.adj.path')])
		}
	}

	### network plot
	vert.attr.list <- get.vertex.attribute(net)
	node.dat.id <- vert.attr.list$name
	### vertex label
	if(is.na(LR.label)){
		vert.label <- rep(NA,length(node.dat.id))
	}else if(LR.label=='L'){
		vert.label <- vert.attr.list$ligand
	}else if(LR.label=='R'){
		vert.label <- vert.attr.list$receptor
	}else if(LR.label=='LR'){
		vert.label <- vert.attr.list$LR.pair
	}else{
		stop('Wrong LR.label parameter! Available options: LR (to show both the ligand and receptor label), L (to show only the ligand label), R (to show only the receptor label), or NA (no label) ')
	}

	if(pathway.label){
		if(all(is.na(vert.label))){
			vert.label <- vert.attr.list$description.path
		}else{
			vert.label[which(is.na(vert.label))] <- node.dat.id[which(is.na(vert.label))]
		}
	}

	### vertex size
	if(vert.size.attr=='degree'){
		vert.size <- igraph::degree(net, mode="out")
		vert.size[which(vert.size==0)] <- NA
	}else{
		if(vert.size.attr=='logP'){
			vert.size <- vert.attr.list[['P.val.path']]
			vert.size <- -log10(vert.size)
		}else if(vert.size.attr=='logP.adj'){
				vert.size <- vert.attr.list[['P.val.adj.path']]
				vert.size <- -log10(vert.size)
		}else{
			vert.size <- vert.attr.list[[paste0(vert.size.attr,'.path')]]
			if(is.null(vert.size)){
				stop(paste0('Wrong "vert.size.attr" parameter! There is no valid attribute named as ',vert.size.attr,'\nType "getPathAttr(object)" for available attribute'))
			}
		}
	}
	vert.size <- scale_1(vert.size)
	vert.size <- vert.size * vert.size.path.adj
	vert.size[which(is.na(vert.size))] <- vert.size.LR
	
	### vertex color
	all.ident <- object@cell.info$Cluster
	if (!is.factor(all.ident)){ all.ident <- factor(all.ident) }
	ident.label <- levels(all.ident)
	vert.cluster.ident <- ident.label[which(ident.label %in% as.character(vert.attr.list$cell.to))]
	if(is.null(ident.col)){
		scales.color <- scales::hue_pal(c=100)(length(ident.label))
		vert.cluster.color <- scales.color[which(ident.label %in% as.character(vert.attr.list$cell.to))]
	}else{
		col.name <- names(ident.col)
		if(is.null(col.name)){
			stop(paste0('Wrong "ident.col" parameter! The colors should be named according to the clusters'))
		}else{
			ident.missed <- vert.cluster.ident[which(!(vert.cluster.ident %in% col.name))]
			if(length(ident.missed) > 0){
				stop(paste0('The color of ',pasteIdent(ident.missed),' may be missed in the input color'))
			}
			vert.cluster.color <- as.vector(ident.col[vert.cluster.ident])
		}
	}
	vertex.pie.color <- list(vert.cluster.color)
	### pie border
	vertex.frame.na.color <- mapvalues(vert.attr.list$type, from=c('LR','Pathway'), to=c(NA,vertex.frame.color), warn_missing=FALSE) 
	
	### pie chart composition
	pie.comp.path <- as.data.frame.array(table(path.curcell.dat[,c('cell.to','path.contain.lig.unfold')]))
	pie.comp.path$rownames <- rownames(pie.comp.path)
	pie.comp.path <- pie.comp.path[vert.cluster.ident, ]
	pie.comp.path$rownames <- NULL
	pie.comp.node <- as.data.frame.array(table(node.dat[,c('cell.to','ID')]))
	pie.comp.node <- pie.comp.node[vert.cluster.ident, node.dat$ID]

	for(each.path in colnames(pie.comp.path)){
		pie.comp.node[ ,each.path] <- pie.comp.path[ ,each.path] 
	}
	pie.comp.node <- lapply(pie.comp.node, function(x){x})

	if(!node.pie){
		scales.color <- scales::hue_pal(c=100)(length(ident.label))
		vert.cluster.color <- scales.color[which(ident.label %in% as.character(vert.attr.list$cell.to))]
		vert.color <- mapvalues(vert.attr.list$cell.to, from=c(vert.cluster.ident,NA), to=c(vert.cluster.color,'tomato'), warn_missing=FALSE) 
	}

	### layout
	l <- do.call(layout, list(net))
	vertex.label.cex <- as.numeric(mapvalues(vert.attr.list$type, from=c('LR','Pathway'), to=c(vertex.label.cex.LR, vertex.label.cex.path), warn_missing=FALSE))
	if(node.pie){
		plot(net, vertex.shape="pie", vertex.pie=pie.comp.node, edge.arrow.size=edge.arrow.size, edge.color=edge.color, layout=l, vertex.size=vert.size, vertex.label=vert.label, vertex.label.color=vertex.label.color, vertex.label.cex=vertex.label.cex, vertex.label.family=vertex.label.family,vertex.frame.color=vertex.frame.na.color, vertex.pie.color=vertex.pie.color)
		legend(x=1.1, y=1, legend=vert.cluster.ident, pch=21, col=vertex.frame.color, pt.bg=vert.cluster.color, pt.cex=1, cex=.5, bty="n", ncol=1)
	}else{
		plot(net, edge.arrow.size=edge.arrow.size, edge.color=edge.color, layout=l, vertex.size=vert.size, vertex.label=vert.label, vertex.label.color=vertex.label.color, vertex.label.cex=vertex.label.cex, vertex.label.family=vertex.label.family,vertex.frame.color=vertex.frame.color, vertex.color=vert.color)
		legend(x=1.1, y=1, legend=vert.cluster.ident, pch=21, col=vertex.frame.color, pt.bg=vert.cluster.color, pt.cex=1, cex=.5, bty="n", ncol=1)
	}
}
