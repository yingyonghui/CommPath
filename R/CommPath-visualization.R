#' To present a circos plot
#' @param object CommPath object
#' @param plot To present a circos plot for LR count ("count") or overall interaction intensity ("intensity") among cell clusters
#' @param col Vector of colors of each identity; names of the col vector are supposed to be assigned to indicate each color for each identity
#' @param select.ident To highlight the interaction between a specific identity class and others; if 'NULL', plot interaction for all identity classes
#' @param name.vert Should the group annotation be vertical to the grid? Defualt is FALSE
#' @importFrom circlize circos.clear circos.par chordDiagram circos.track circos.text
#' @return Circos plot showing the ligand-receptor interaction
#' @export
circosPlot <- function(object, plot='count', col=NULL, select.ident=NULL, name.vert=FALSE){
	options(stringsAsFactors=F)
	Cluster <- object@cell.info$Cluster
	if (!is.factor(Cluster)){ Cluster <- factor(Cluster) }
	ident.label <- levels(Cluster)

	### select LR count or cluster intensity to plot
	Interact.num.dat <- object@interact$InteractNumer
	if (plot=='count'){
		Interact.num.dat$Plot <- Interact.num.dat$LR.Count
	}else if(plot=='intensity'){
		Interact.num.dat$Plot <- Interact.num.dat$Intensity
	}else{
		stop('Select "count" or "intensity" for the parameter "plot"!')
	}
	Interact.num.dat <- Interact.num.dat[, c('Cell.From','Cell.To','Plot')]
	
	Interact.num.dat = Interact.num.dat[Interact.num.dat$Plot!=0,]
	all.ident <- unique(c(Interact.num.dat$Cell.From,Interact.num.dat$Cell.To))
	all.ident <- orderCheck(all.ident, ident.label)
	all.ident <- all.ident[order(all.ident)]
	
	### to check the col parameter
	if (is.null(col)){ 
		col <- scales::hue_pal(c=100)(length(all.ident))
		names(col) <- all.ident
	}else{
		if (is.null(names(col))){
			stop('The cols should be named according to the identity class')
		}else{
			colo.name <- names(col)
			ident.missed <- all.ident[!(all.ident %in% colo.name)]
			ident.missed <- pasteIdent(ident.missed)
			stop(paste0('The ident class ',ident.missed,' may be missed in the input color'))
		}
	}

	circos.clear()
	circos.par(start.degree=90, clock.wise=F)

	### plot interaction for all identity classes
	if (is.null(select.ident)){
		if (!name.vert){
			chordDiagram(Interact.num.dat, order=all.ident, grid.col=col, annotationTrack=c("name","grid"), transparency=0.1, directional=1, direction.type='arrows', link.arr.type = "big.arrow",preAllocateTracks = list(track.height=mm_h(5)), annotationTrackHeight=convert_height(c(1, 2), "mm"))
		}else{
			chordDiagram(Interact.num.dat, order=all.ident, grid.col=col, annotationTrack=c("grid"), transparency=0.1, directional=1, direction.type='arrows', link.arr.type = "big.arrow", preAllocateTracks = list(track.height=mm_h(5)), annotationTrackHeight=convert_height(c(2, 2), "mm"))
			circos.track(track.index=1, panel.fun=function(x, y){circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)}, bg.border = NA)
		}

	}else{
		### highlight the interaction for a specific identity class
		### check the select.ident parameter
		if (!(all(select.ident %in% all.ident)) | length(select.ident)>1){
			stop(paste0('Select one identity class from ', pasteIdent(all.ident)))
		}

		line.col <- rep('gray',nrow(Interact.num.dat))
		line.col[(Interact.num.dat$Cell.From==select.ident)] <- col[as.character(select.ident)]
		ident.lig = Interact.num.dat[Interact.num.dat$Cell.To==select.ident,'Cell.From']
		line.col[(Interact.num.dat$Cell.To==select.ident)] <- col[as.character(ident.lig)]
		
		link.zindex <- 1:nrow(Interact.num.dat)
		ident.line <- (Interact.num.dat$Cell.From==select.ident) | (Interact.num.dat$Cell.To==select.ident)
		link.zindex[!ident.line] <- rank(-(Interact.num.dat[!ident.line,'Plot']))
		link.zindex[ident.line] <- rank(-(Interact.num.dat[ident.line,'Plot']))+length(which(!ident.line))
		if (!name.vert){
			chordDiagram(Interact.num.dat, order=all.ident, grid.col=col, annotationTrack=c("name","grid"), transparency=0.1, directional=1, direction.type='arrows', link.arr.type = "big.arrow", link.zindex=link.zindex, col=line.col, preAllocateTracks = list(track.height=mm_h(5)), annotationTrackHeight=convert_height(c(1, 2), "mm"))
		}else{
			chordDiagram(Interact.num.dat, order=all.ident, grid.col=col, annotationTrack="grid", transparency=0.1, directional=1, direction.type='arrows', link.arr.type = "big.arrow", link.zindex=link.zindex, col=line.col, preAllocateTracks = list(track.height=mm_h(5)), annotationTrackHeight=convert_height(c(2, 2), "mm"))
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
#' @param return.data Logical value indicating whether to return the data for the plot or not
#' @importFrom ggplot2 ggplot geom_point scale_color_manual labs theme element_text element_rect element_line aes scale_color_gradientn
#' @return Dotplot showing the ligand-receptor interaction between the selected ligand.ident and receptor.ident
#' @export
dotPlot <- function(object, ligand.ident=NULL, receptor.ident=NULL, ident.levels=NULL, top.n.inter=10, return.data=FALSE){
	options(stringsAsFactors=F)
	if (is.null(ligand.ident) & is.null(receptor.ident)){
		stop("Either ligand.ident or ligand.ident need to be asigned")
	}
	if (length(ligand.ident)>1 & length(receptor.ident)>1){
		stop("Specify one cluster for ligand or receptor analysis")
	}

	# get the InteractGene dataframe
	inter.gene.dat <- object@interact$InteractGene
	
	if (!is.null(receptor.ident)){
		which.rep.not.in <- which(!(receptor.ident %in% inter.gene.dat$Cell.To))
		which.rep.in <- which(receptor.ident %in% inter.gene.dat$Cell.To)

		if (length(which.rep.not.in) > 0){
			if (length(which.rep.in) > 0){
				warning(paste0('No ligand and upstream cluster detected for cluster ', pasteIdent(unique(receptor.ident[which.rep.not.in]))))
			}else{
				stop(paste0('No ligand and upstream cluster detected for cluster ', pasteIdent(unique(receptor.ident[which.rep.not.in]))))
			}
		}
		receptor.ident <- receptor.ident[which.rep.in]
		inter.gene.dat <- inter.gene.dat[inter.gene.dat$Cell.To %in% receptor.ident,]
	}
	if (!is.null(ligand.ident)){
		which.lig.not.in <- which(!(ligand.ident %in% inter.gene.dat$Cell.From))
		which.lig.in <- which(ligand.ident %in% inter.gene.dat$Cell.From)

		if (length(which.lig.not.in) > 0){
			if (length(which.lig.in) > 0){
				warning(paste0('No ligand and upstream cluster detected for cluster ', pasteIdent(unique(ligand.ident[which.lig.not.in]))))
			}else{
				stop(paste0('No ligand and upstream cluster detected for cluster ', pasteIdent(unique(ligand.ident[which.lig.not.in]))))
			}
		}
		ligand.ident <- ligand.ident[which.lig.in]
		inter.gene.dat <- inter.gene.dat[inter.gene.dat$Cell.From %in% ligand.ident,]
	}

	inter.gene.dat$LR.Info <- paste(inter.gene.dat$Ligand, ' --> ', inter.gene.dat$Receptor)

	### to find those LR pairs with largest Log2FC.LR
	max.fc <- by(data=inter.gene.dat$Log2FC.LR, INDICES=inter.gene.dat$LR.Info, function(x){max(x,na.rm=T)})
	max.LRname.fc <- names(max.fc)[order(max.fc, decreasing=TRUE)]
	
	if (top.n.inter > length(max.LRname.fc)){
		warning(paste0('There is(are) ', length(max.LRname.fc),' significant LR pair(s) for the selected ident, and the input top.n.inter is ', top.n.inter))
		top.n.inter <- length(max.LRname.fc)
	}

	max.LRname.fc <- max.LRname.fc[1:top.n.inter]
	inter.gene.dat <- subset(inter.gene.dat, LR.Info %in% max.LRname.fc)

	if (return.data){
		return(inter.gene.dat[,c('Cell.From', 'Cell.To', 'LR.Info', 'Log2FC.LR', 'P.val.LR', 'P.val.adj.LR')])
	}
	
	inter.gene.dat$Log10.P.adj <- -log10(inter.gene.dat$P.val.adj.LR)
	inter.gene.dat$Log10.P.adj <- p.remove.inf(inter.gene.dat$Log10.P.adj)
	if(all(inter.gene.dat$Log10.P.adj==0)){inter.gene.dat$Log10.P.adj='Infinite'}

	#inter.gene.dat[which(inter.gene.dat$Log10.P.adj > 30), 'Log10.P.adj'] <- 30
	
	### adjust the x label and title
	if (length(ligand.ident)==1){
		inter.gene.dat$Xaxis <- inter.gene.dat$Cell.To
		x.title <- paste0('Receptor clusters for cluster ',ligand.ident)
	}else{
		inter.gene.dat$Xaxis <- inter.gene.dat$Cell.From
		x.title <- paste0('Ligand clusters to cluster ',receptor.ident)
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
		plot <- ggplot(data=inter.gene.dat,aes(x=Xaxis,y=LR.Info)) + 
			geom_point(aes(size=Log2FC.LR,col='Infinite')) +
			scale_color_manual(values='red') + 
			labs(color='-log10(p.adj)',size='Log2FC.LR',x=x.title,y="") + 
			theme(axis.text=element_text(colour='black',size=12),
			axis.title=element_text(colour='black',size=12),
			panel.background=element_rect(fill="white",color="black"),
			panel.grid=element_line(size=0.5,colour='gray')
			)
		warning('Adjusted p values for all LR pairs are 0')
	}else{
		plot <- ggplot(data=inter.gene.dat,aes(x=Xaxis,y=LR.Info)) + 
			geom_point(aes(size=Log2FC.LR,col=Log10.P.adj)) +
			scale_color_gradientn(values = seq(0,1,0.2),colours=c("#0570b0", "grey", "#d7301f")) + 
			labs(color='-log10(p.adj)',size='Log2FC.LR',x=x.title,y="") + 
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
#' @param path.order Sort criteria used to select the top n pathways, either 'p.val' or 'p.val.adj', which represent the original and adjusted p values, or 'diff' which represents the mean (in t test) or median (in wilcox test) difference
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
pathHeatmap <- function(object, acti.path.dat=NULL, top.n.pathway=10, path.order='p.val.adj', col=NULL, cell.aver=FALSE, cell.label.size=NULL, cell.label.angle=45, pathway.label.size=NULL, scale=TRUE, truncation=1.5, show.legend=TRUE){
	if (is.null(acti.path.dat)){
		acti.path.dat <- diffAllPath(object)
	}
	
	acti.path.dat <- subset(acti.path.dat, p.val.adj < 0.05)
	### select the highly enriched pathways
	if ('mean.diff' %in% colnames(acti.path.dat)){
		acti.path.dat <- subset(acti.path.dat, mean.diff > 0)
		acti.path.dat$diff <- acti.path.dat$mean.diff
	}else if('median.diff' %in% colnames(acti.path.dat)){
		acti.path.dat <- subset(acti.path.dat, median.diff > 0)
		acti.path.dat$diff <- acti.path.dat$median.diff
	}else{
		stop('Please input the intact test result computed from diffAllPath')
	}

	### select top n pathways
	if (path.order=='p.val.adj' | path.order=='p.val'){
		acti.path.dat <- acti.path.dat[order(acti.path.dat[[path.order]], -acti.path.dat[['diff']]), ]
	}else if (path.order=='diff'){
		acti.path.dat <- acti.path.dat[order(-acti.path.dat[['diff']], acti.path.dat[['p.val.adj']],), ]
	}else{
		stop('Select one order criteria from p.val.adj, p.val, and diff')
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
#' @param acti.path.dat Data frame of differential activation test result from diffAllPath
#' @param top.n.path Top n pathways with the smallest adjusted p values to plot
#' @param path.order Sort criteria used to select the top n pathways, either 'p.val' or 'p.val.adj', which represent the original and adjusted p values, or 'diff' which represents the mean (in t test) or median (in wilcox test) difference
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
pathPlot <- function(object, select.ident, acti.path.dat=NULL, top.n.path=5, path.order='p.val.adj', p.thre=0.05, top.n.receptor=10, dot.ident.col=NULL, dot.ident.size=1, dot.gene.col=NULL, dot.gene.size=1, bar.pathway.col=NULL, label.text.size=1, label.title.size=1){
	options(stringsAsFactors=F)
	all.ident <- object@cell.info$Cluster
	if (!is.factor(all.ident)){ all.ident <- factor(all.ident) }
	ident.label <- levels(all.ident)

	### check the input and select significant pathways
	if (is.null(acti.path.dat)) {
		ident.path.dat <- diffPath(object, select.ident.1=select.ident)
	}else{
		ident.path.dat <- subset(acti.path.dat, cluster==select.ident)
	}
	if ('t' %in% colnames(ident.path.dat)){
		ident.path.dat <- subset(ident.path.dat, p.val.adj < p.thre & t > 0)
		ident.path.dat$diff <- ident.path.dat$mean.diff
	}else if ('W' %in% colnames(ident.path.dat)){
		ident.path.dat <- subset(ident.path.dat, p.val.adj < p.thre & median.diff > 0)
		ident.path.dat$diff <- ident.path.dat$median.diff
	}else{
		stop('Please input the integrate ident.path.dat computed from diffPath')
	}

	if (nrow(ident.path.dat)==0){ paste0(stop('There is no significantly up-regulatged pathways for cluster ', select.ident )) }

	### remove pathways with NA receptor
	ident.path.dat <- subset(ident.path.dat, !is.na(receptor.in.path))
	if (nrow(ident.path.dat)==0){ paste0(stop('There is no marker receptor in the significantly up-regulatged pathways for cluster ', select.ident )) }

	if (path.order=='p.val.adj' | path.order=='p.val'){
		ident.path.dat <- ident.path.dat[order(ident.path.dat[,path.order], -ident.path.dat[,'diff']),]
	}else if(path.order=='diff'){
		ident.path.dat <- ident.path.dat[order(-ident.path.dat[,'diff'], ident.path.dat[,'p.val.adj']),]
	}

	if (top.n.path > nrow(ident.path.dat)){
		warning(paste0('There is(are) ', nrow(ident.path.dat),' significant pathway(s) for the selected ident, and the input top.n.path is ', top.n.path))
		top.n.path <- nrow(ident.path.dat)
	}
	ident.path.dat <- ident.path.dat[1:top.n.path, ]
	all.sig.path <- as.vector(ident.path.dat$description)

	### data for plot of receptors and upstream clusters
	receptor.in.path <- ident.path.dat[, 'receptor.in.path']
	names(receptor.in.path) <- all.sig.path
	# receptor.in.path <- receptor.in.path[which(!is.na(receptor.in.path))]
	cur.rep <- as.vector(unlist(sapply(receptor.in.path, function(x) {strsplit(x, split=';')})))
	path.times <- sapply(receptor.in.path, function(x) {length(strsplit(x, split=';')[[1]])})
	cur.path <- rep(names(path.times), times=path.times)
	plot.receptor.to.pathway.dat <- data.frame(cur.rep=cur.rep, cur.path=cur.path)

	### select the top n receptors
	cur.uniq.rep <- as.vector(unique(plot.receptor.to.pathway.dat$cur.rep))
	cur.rep.LR.inten <- cluster.lr.inten(cur.uniq.rep, object, select.ident, ident.label, find='ligand')
	cur.rep.max.LR.inten <- apply(cur.rep.LR.inten, 2, function(x){max(x, na.rm=TRUE)})
	cur.rep.max.LR.inten <- cur.rep.max.LR.inten[order(cur.rep.max.LR.inten, decreasing=TRUE)]
	
	if (length(cur.uniq.rep) < top.n.receptor){
		warning(paste0('There is(are) ', length(cur.uniq.rep),' marker receptors(s) in the selected ident associated with the selected pathways, and the input top.n.receptor is ', top.n.receptor))
		top.n.receptor <- length(cur.uniq.rep)
	}
	
	top.rep.name <- names(cur.rep.max.LR.inten)[1:top.n.receptor]
	top.rep.LR.inten <- cur.rep.LR.inten[, top.rep.name]

	### if only one top.rep.name
	if (is.null(dim(top.rep.LR.inten))){
		top.rep.LR.inten <- as.matrix(top.rep.LR.inten)
		colnames(top.rep.LR.inten) <- top.rep.name
	}
	plot.receptor.to.pathway.dat <- subset(plot.receptor.to.pathway.dat, cur.rep %in% top.rep.name)
	plot.ident.to.receprtor.dat <- melt(top.rep.LR.inten, varnames=c('Cell.From','Receptor'),value.name="LR.inten",  na.rm=TRUE)
	plot.ident.to.receprtor.dat <- factor.to.character(plot.ident.to.receprtor.dat)

	### pathways and their color
	path.uniq.name <- factor(all.sig.path, levels=all.sig.path)
	path.uniq.name <- path.uniq.name[order(path.uniq.name, decreasing=TRUE)]
	
	#### rect for pathway
	path.rect.length <- ident.path.dat[match(path.uniq.name, ident.path.dat$description),'diff']
	path.rect.pval <-  ident.path.dat[match(path.uniq.name, ident.path.dat$description),'p.val.adj']
	if (all(path.rect.pval==0)){
		path.rect.pval <- path.rect.length
		warning('All adjusted p values for the selected pathways are 0. The colors of bars representing pathways are adjusted to indicate the mean (in t test) or median (in wilcox test) difference of pathways')
	}else{
		path.rect.pval <- -log10(path.rect.pval)
		path.rect.pval <- p.remove.inf(path.rect.pval)
	}
	bar.path.col <- LRcolor(path.rect.pval, user.set.col=bar.pathway.col)

	
	### the receptor and their coordinate, and their color
	cur.uniq.rep <- factor(top.rep.name, levels=top.rep.name)
	cur.uniq.rep <- cur.uniq.rep[order(cur.uniq.rep, decreasing=TRUE)]
	n.rep <- length(cur.uniq.rep)
	
	top.markerR.dat <- subset(object@interact$markerR, subset=(cluster==select.ident & gene %in% top.rep.name))

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
	up.uniq.ident <- factor(ident.label, levels=ident.label)
	up.uniq.ident <- up.uniq.ident[order(up.uniq.ident, decreasing=TRUE)]
	n.up.ident <- length(up.uniq.ident)
	if (is.null(dot.ident.col)){
		dot.ident.col <- rev(scales::hue_pal(c=100)(n.up.ident))
	}else if(length(dot.ident.col) < n.up.ident){
		stop(paste0('There is(are) ',n.up.ident,' cluster(s), but only ',length(dot.ident.col),' colors provided.'))
	}else{
		dot.ident.col <- rev(dot.ident.col[1:n.up.ident])
	}
	
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
	point_up_ident <- geom_point(aes(x=1,y=up.ident.new.coor),size=8*dot.ident.size,color=dot.ident.col,shape=16) 

	### lines between elements
	### lines between up idents and receptors
	line.ident.to.rep.y.coor <- up.ident.new.coor[match(plot.ident.to.receprtor.dat$Cell.From, up.uniq.ident)]
	line.ident.to.rep.yend.coor <- rep.new.coor[match(plot.ident.to.receprtor.dat$Receptor, cur.uniq.rep)]
	line.col <- dot.ident.col[match(plot.ident.to.receprtor.dat$Cell.From, up.uniq.ident)]
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
#' @param diff.marker.dat Data frame of defferential expression test result from compareMarker; if NULL, compareMarker would be run to get diff.marker.dat
#' @param diff.path.dat Data frame of defferential activation test result from comparePath; if NULL, comparePath would be run to get diff.path.dat
#' @param top.n.path Top n pathways with the smallest adjusted p values to plot
#' @param path.order Sort criteria used to select the top n pathways, either 'p.val' or 'p.val.adj', which represent the original and adjusted p values, or 'diff' which represents the mean (in t test) or median (in wilcox test) difference
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
pathPlot.compare <- function(object.1, object.2, select.ident, diff.marker.dat=NULL, diff.path.dat=NULL, top.n.path=5, path.order='p.val.adj', p.thre=0.05, top.n.receptor=10, dot.ident.col=NULL, dot.ident.size=1, dot.gene.col=NULL, dot.gene.size=1, bar.pathway.col=NULL, label.text.size=1, label.title.size=1){
	options(stringsAsFactors=F)
	# if(is.null(diff.marker.dat)){
	# 	message(paste0('Identifying differentially expressed ligands(receptors) between cluster ',select.ident,' in object 1 and object 2'))
	# 	diff.marker.dat <- compareMarker(object.1, object.2, select.ident, method='wilcox.test', p.adjust='BH', only.posi=FALSE, only.sig=TRUE)
	# }
	if (is.null(diff.path.dat)){
		message(paste0('Identifying differentially activated pathways between cluster ',select.ident,' in object 1 and object 2'))
		diff.path.dat <- comparePath(object.1, object.2, select.ident, method='wilcox.test', p.adjust='BH', min.size=10, only.posi=FALSE, only.sig=TRUE)
	}
	all.ident <- object.1@cell.info$Cluster
	if (!is.factor(all.ident)){ all.ident <- factor(all.ident) }
	ident.label <- levels(all.ident)

	### check the input and select significant pathways
	ident.path.dat <- diff.path.dat
	if ('t' %in% colnames(ident.path.dat)){
		ident.path.dat <- subset(ident.path.dat, p.val.adj < p.thre & t > 0)
		ident.path.dat$diff <- ident.path.dat$mean.diff
	}else if ('W' %in% colnames(ident.path.dat)){
		ident.path.dat <- subset(ident.path.dat, p.val.adj < p.thre & median.diff > 0)
		ident.path.dat$diff <- ident.path.dat$median.diff
	}else{
		stop('Please input the integrate ident.path.dat computed from diffPath')
	}

	if (nrow(ident.path.dat)==0){ paste0(stop('There is no significantly up-regulatged pathways for cluster ', select.ident )) }

	### remove pathways with NA receptor
	ident.path.dat <- subset(ident.path.dat, !is.na(receptor.in.path))
	if (nrow(ident.path.dat)==0){ paste0(stop('There is no marker receptor in the significantly up-regulatged pathways for cluster ', select.ident )) }

	if (path.order=='p.val.adj' | path.order=='p.val'){
		ident.path.dat <- ident.path.dat[order(ident.path.dat[,path.order], -ident.path.dat[,'diff']),]
	}else if(path.order=='diff'){
		ident.path.dat <- ident.path.dat[order(-ident.path.dat[,'diff'], ident.path.dat[,'p.val.adj']),]
	}

	if (top.n.path > nrow(ident.path.dat)){
		warning(paste0('There is(are) ', nrow(ident.path.dat),' significant pathway(s) for the selected ident, and the input top.n.path is ', top.n.path))
		top.n.path <- nrow(ident.path.dat)
	}
	ident.path.dat <- ident.path.dat[1:top.n.path, ]
	all.sig.path <- as.vector(ident.path.dat$description)

	### data for plot of receptors and upstream clusters
	receptor.in.path <- ident.path.dat[, 'receptor.in.path']
	names(receptor.in.path) <- all.sig.path
	# receptor.in.path <- receptor.in.path[which(!is.na(receptor.in.path))]
	cur.rep <- as.vector(unlist(sapply(receptor.in.path, function(x) {strsplit(x, split=';')})))
	path.times <- sapply(receptor.in.path, function(x) {length(strsplit(x, split=';')[[1]])})
	cur.path <- rep(names(path.times), times=path.times)
	plot.receptor.to.pathway.dat <- data.frame(cur.rep=cur.rep, cur.path=cur.path)

	### select the top n receptors
	cur.uniq.rep <- as.vector(unique(plot.receptor.to.pathway.dat$cur.rep))
	cur.rep.LR.inten <- cluster.lr.inten(cur.uniq.rep, object.1, select.ident, ident.label, find='ligand')
	cur.rep.max.LR.inten <- apply(cur.rep.LR.inten, 2, function(x){max(x, na.rm=TRUE)})
	cur.rep.max.LR.inten <- cur.rep.max.LR.inten[order(cur.rep.max.LR.inten, decreasing=TRUE)]
	
	if (length(cur.uniq.rep) < top.n.receptor){
		warning(paste0('There is(are) ', length(cur.uniq.rep),' marker receptors(s) in the selected ident associated with the selected pathways, and the input top.n.receptor is ', top.n.receptor))
		top.n.receptor <- length(cur.uniq.rep)
	}
	
	top.rep.name <- names(cur.rep.max.LR.inten)[1:top.n.receptor]
	plot.receptor.to.pathway.dat <- subset(plot.receptor.to.pathway.dat, cur.rep %in% top.rep.name)

	top.rep.LR.inten <- cur.rep.LR.inten[, top.rep.name]
	if (is.null(dim(top.rep.LR.inten))){
		top.rep.LR.inten <- as.matrix(top.rep.LR.inten)
		colnames(top.rep.LR.inten) <- top.rep.name
	}
	plot.ident.to.receprtor.obj1.dat <- melt(top.rep.LR.inten, varnames=c('Cell.From','Receptor'),value.name="LR.inten",  na.rm=FALSE)
	plot.ident.to.receprtor.obj1.dat <- factor.to.character(plot.ident.to.receprtor.obj1.dat)

	### LR.inten for obj.2
	cur.rep.LR.obj2.inten <- cluster.lr.inten(top.rep.name, object.2, select.ident, ident.label, find='ligand')
	plot.ident.to.receprtor.obj2.dat <- melt(cur.rep.LR.obj2.inten, varnames=c('Cell.From','Receptor'),value.name="LR.inten",  na.rm=FALSE)
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
	path.rect.pval <-  ident.path.dat[match(path.uniq.name, ident.path.dat$description),'p.val.adj']

	if (all(path.rect.pval==0)){
		path.rect.pval <- path.rect.length
		warning('All adjusted p values for the selected pathways are 0. The colors of bars representing pathways are adjusted to indicate the mean (in t test) or median (in wilcox test) difference of pathways')
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
	up.uniq.ident <- factor(ident.label, levels=ident.label)
	up.uniq.ident <- up.uniq.ident[order(up.uniq.ident, decreasing=TRUE)]
	
	n.up.ident <- length(up.uniq.ident)
	if (is.null(dot.ident.col)){
		dot.ident.col <- rev(scales::hue_pal(c=100)(n.up.ident))
	}else if(length(dot.ident.col) < n.up.ident){
		stop(paste0('There is(are) ',n.up.ident,' cluster(s), but only ',length(dot.ident.col),' colors provided.'))
	}else{
		dot.ident.col <- rev(dot.ident.col[1:n.up.ident])
	}
	
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
	point_up_ident <- geom_point(aes(x=1,y=up.ident.new.coor),size=8*dot.ident.size,color=dot.ident.col,shape=16) 

	### lines between elements
	### lines between up idents and receptors
	line.ident.to.rep.y.coor <- up.ident.new.coor[match(plot.ident.to.receprtor.dat$Cell.From, up.uniq.ident)]
	line.ident.to.rep.yend.coor <- rep.new.coor[match(plot.ident.to.receprtor.dat$Receptor, cur.uniq.rep)]
	#line.col <- dot.ident.col[match(plot.ident.to.receprtor.dat$Cell.From, up.uniq.ident)]
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


#' To present the interactions for a selected cluster, including both the upstream and downstream clusters which are connected by specific pathways in the selected cluster
#' @param object CommPath object
#' @param select.ident Plot the activated pathways for which cluster or cell type?
#' @param acti.path.dat Data frame of differential activation test result from diffAllPath
#' @param top.n.path Top n pathways with the smallest adjusted p values to plot
#' @param path.order Sort criteria used to select the top n pathways, either 'p.val' or 'p.val.adj', which represent the original and adjusted p values, or 'diff' which represents the mean (in t test) or median (in wilcox test) difference
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
pathInterPlot <- function(object, select.ident, acti.path.dat=NULL, top.n.path=5, path.order='p.val.adj', p.thre=0.05, top.n.receptor=10, top.n.ligand=10, dot.ident.col=NULL, dot.ident.size=1, dot.gene.col=NULL, dot.gene.size=1, bar.pathway.col=NULL, label.text.size=1, label.title.size=1){
	options(stringsAsFactors=F)
	all.ident <- object@cell.info$Cluster
	if (!is.factor(all.ident)){ all.ident <- factor(all.ident) }
	ident.label <- levels(all.ident)

	### check the input and select significant pathways
	if (is.null(acti.path.dat)) {
		ident.path.dat <- diffPath(object, select.ident.1=select.ident)
	}else{
		ident.path.dat <- subset(acti.path.dat, cluster==select.ident)
	}

	if ('t' %in% colnames(ident.path.dat)){
		ident.path.dat <- subset(ident.path.dat, p.val.adj < p.thre & t > 0)
		ident.path.dat$diff <- ident.path.dat$mean.diff
	}else if ('W' %in% colnames(ident.path.dat)){
		ident.path.dat <- subset(ident.path.dat, p.val.adj < p.thre & median.diff > 0)
		ident.path.dat$diff <- ident.path.dat$median.diff
	}else{
		stop('Please input the integrate acti.path.dat computed from diffPath')
	}

	if (nrow(ident.path.dat)==0){ paste0(stop('There is no significantly up-regulatged pathways for cluster ', select.ident )) }

	if (path.order=='p.val.adj' | path.order=='p.val'){
		ident.path.dat <- ident.path.dat[order(ident.path.dat[,path.order], -ident.path.dat[,'diff']),]
	}else if(path.order=='diff'){
		ident.path.dat <- ident.path.dat[order(-ident.path.dat[,'diff'], ident.path.dat[,'p.val.adj']),]
	}
	if (top.n.path > nrow(ident.path.dat)){
		warning(paste0('There is(are) ', nrow(ident.path.dat),' significant pathway(s) for the selected ident, and the input top.n.path is ', top.n.path))
		top.n.path <- nrow(ident.path.dat)
	}
	ident.path.dat <- ident.path.dat[1:top.n.path, ]
	all.sig.path <- as.vector(ident.path.dat$description)

	### data for plot of receptors and upstream clusters
	receptor.in.path <- ident.path.dat[, 'receptor.in.path']
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
	cur.rep.max.LR.inten <- apply(cur.rep.LR.inten, 2, function(x){max(x, na.rm=TRUE)})
	cur.rep.max.LR.inten <- cur.rep.max.LR.inten[order(cur.rep.max.LR.inten, decreasing=TRUE)]
	
	if (length(cur.uniq.rep) < top.n.receptor){
		warning(paste0('There is(are) ', length(cur.uniq.rep),' marker receptor(s) in the selected ident associated with the selected pathways, and the input top.n.receptor is ', top.n.receptor))
		top.n.receptor <- length(cur.uniq.rep)
	}
	
	top.rep.name <- names(cur.rep.max.LR.inten)[1:top.n.receptor]
	plot.receptor.to.pathway.dat <- subset(plot.receptor.to.pathway.dat, cur.rep %in% top.rep.name)

	top.rep.LR.inten <- cur.rep.LR.inten[, top.rep.name]
	if (is.null(dim(top.rep.LR.inten))){
		top.rep.LR.inten <- as.matrix(top.rep.LR.inten)
		colnames(top.rep.LR.inten) <- top.rep.name
	}
	plot.ident.to.receprtor.dat <- melt(top.rep.LR.inten, varnames=c('Cell.From','Receptor'),value.name="LR.inten",  na.rm=TRUE)
	plot.ident.to.receprtor.dat <- factor.to.character(plot.ident.to.receprtor.dat)
	
	### data for plot of ligands and downstream clusters
	ligand.in.path <- ident.path.dat[, 'ligand.in.path']
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
	cur.lig.LR.inten <- cluster.lr.inten(cur.uniq.lig, object, select.ident, ident.label, find='receptor')
	cur.lig.max.LR.inten <- apply(cur.lig.LR.inten, 2, function(x){max(x, na.rm=TRUE)})
	cur.lig.max.LR.inten <- cur.lig.max.LR.inten[order(cur.lig.max.LR.inten, decreasing=TRUE)]
	
	if (length(cur.uniq.lig) < top.n.ligand){
		warning(paste0('There is(are) ', length(cur.uniq.lig),' marker ligand(s) in the selected ident associated with the selected pathways, and the input top.n.ligand is ', top.n.ligand))
		top.n.ligand <- length(cur.uniq.lig)
	}
	
	top.lig.name <- names(cur.lig.max.LR.inten)[1:top.n.ligand]
	plot.pathway.to.ligand.dat <- subset(plot.pathway.to.ligand.dat, cur.lig %in% top.lig.name)

	top.lig.LR.inten <- cur.lig.LR.inten[, top.lig.name]
	if (is.null(dim(top.lig.LR.inten))){
		top.lig.LR.inten <- as.matrix(top.lig.LR.inten)
		colnames(top.lig.LR.inten) <- top.rep.name
	}
	plot.ligand.to.ident.dat <- melt(top.lig.LR.inten, varnames=c('Cell.To','Ligand'),value.name="LR.inten",  na.rm=TRUE)
	plot.ligand.to.ident.dat <- factor.to.character(plot.ligand.to.ident.dat)

	### pathways and their color
	path.uniq.name <- factor(all.sig.path, levels=all.sig.path)
	path.uniq.name <- path.uniq.name[order(path.uniq.name, decreasing=TRUE)]
	
	#### rect for pathway
	path.rect.length <- ident.path.dat[match(path.uniq.name, ident.path.dat$description),'diff']
	path.rect.pval <-  ident.path.dat[match(path.uniq.name, ident.path.dat$description),'p.val.adj']

	if (all(path.rect.pval==0)){
		path.rect.pval <- path.rect.length
		warning('All adjusted p values for the selected pathways are 0. The colors of bars representing pathways are adjusted to indicate the mean (in t test) or median (in wilcox test) difference of pathways')
	}else{
		path.rect.pval <- -log10(path.rect.pval)
		path.rect.pval <- p.remove.inf(path.rect.pval)
	}
	bar.path.col <- LRcolor(path.rect.pval, user.set.col=bar.pathway.col)

	
	### the receptor and their coordinate, and their color
	cur.uniq.rep <- factor(top.rep.name, levels=top.rep.name)
	cur.uniq.rep <- cur.uniq.rep[order(cur.uniq.rep, decreasing=TRUE)]
	n.rep <- length(cur.uniq.rep)
	
	top.markerR.dat <- subset(object@interact$markerR, subset=(cluster==select.ident & gene %in% top.rep.name))

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
	
	top.markerL.dat <- subset(object@interact$markerL, subset=(cluster==select.ident & gene %in% top.lig.name))

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


	### the upstream and downstream ident and their coordinate, and their color
	up.uniq.ident <- factor(ident.label, levels=ident.label)
	up.uniq.ident <- up.uniq.ident[order(up.uniq.ident, decreasing=TRUE)]
	### n.up.ident is also the n.down.ident
	n.up.ident <- length(up.uniq.ident)
	if (is.null(dot.ident.col)){
		dot.ident.col <- rev(scales::hue_pal(c=100)(n.up.ident))
	}else if(length(dot.ident.col) < n.up.ident){
		stop(paste0('There is(are) ',n.up.ident,' cluster(s), but only ',length(dot.ident.col),' colors provided.'))
	}else{
		dot.ident.col <- rev(dot.ident.col[1:n.up.ident])
	}
	
	### adjust the coordinate of each element
	max.limit <- max(n.up.ident, n.rep, top.n.path, n.lig)
	up.ident.new.coor <- seq(from=1, to=max.limit, length.out=n.up.ident)
	rep.new.coor <- seq(from=1, to=max.limit, length.out=n.rep)
	pathway.new.coor <- seq(from=1, to=max.limit, length.out=top.n.path)
	lig.new.coor <- seq(from=1, to=max.limit, length.out=n.lig)
	down.ident.new.coor <- seq(from=1, to=max.limit, length.out=n.up.ident)

	if (n.up.ident==1){ 
		up.ident.new.coor <- median(1:max.limit) 
		down.ident.new.coor <- median(1:max.limit) 
	}
	if (n.rep==1){ rep.new.coor <- median(1:max.limit) }
	if (top.n.path==1){ pathway.new.coor <- median(1:max.limit) }
	if (n.lig==1){ lig.new.coor <- median(1:max.limit) }
	
	if (length(pathway.new.coor)==1){
		bar_pathway_width <- 0.3
	}else{
		bar_pathway_width <- 0.3 * (pathway.new.coor[2]-pathway.new.coor[1]) * top.n.path / max.limit
	}
	path.name.max.width <- max(strwidth(path.uniq.name, units="inches", cex=1))
	#bar.pathway.length <- bar.pathway.width*path.name.max.width*path.rect.length
	bar.pathway.length <- 2 * path.rect.length / max(path.rect.length)
	bar_pathway <- geom_rect(aes(xmin=6, xmax=6+bar.pathway.length, ymin=pathway.new.coor-bar_pathway_width, ymax=pathway.new.coor+bar_pathway_width), size=1, fill=bar.path.col, alpha=1, show.legend=F)
	point_up_ident <- geom_point(aes(x=1,y=up.ident.new.coor),size=8*dot.ident.size,color=dot.ident.col,shape=16) 
	point_receptor <- geom_point(aes(x=2,y=rep.new.coor),size=rep.size,color=dot.rep.col)
	point_pathway <- geom_point(aes(x=3,y=pathway.new.coor),size=8*dot.ident.size,color='black',shape=21,fill='white')
	point_ligand <- geom_point(aes(x=4,y=lig.new.coor),size=lig.size,color=dot.lig.col) 
	point_down_ident <- geom_point(aes(x=5,y=down.ident.new.coor),size=8*dot.ident.size,color=dot.ident.col,shape=16) 

	### lines between elements
	### lines between up idents and receptors
	line.ident.to.rep.y.coor <- up.ident.new.coor[match(plot.ident.to.receprtor.dat$Cell.From, up.uniq.ident)]
	line.ident.to.rep.yend.coor <- rep.new.coor[match(plot.ident.to.receprtor.dat$Receptor, cur.uniq.rep)]
	line.ident.to.rep.col <- dot.ident.col[match(plot.ident.to.receprtor.dat$Cell.From, up.uniq.ident)]
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
	line.lig.to.ident.y.coor <- lig.new.coor[match(plot.ligand.to.ident.dat$Ligand, cur.uniq.lig)]
	line.lig.to.ident.yend.coor <- up.ident.new.coor[match(plot.ligand.to.ident.dat$Cell.To, up.uniq.ident)]
	line.lig.to.ident.col <- dot.ident.col[match(plot.ligand.to.ident.dat$Cell.To, up.uniq.ident)]
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
	label_down_ident <- annotate('text', hjust=0.5, x=5, y=down.ident.new.coor-0.4,label=up.uniq.ident, size=5*label.text.size)
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


#' To compare the pathway mediated cell-cell communication flow for a specific cluster between two CommPath object
#' @param object.1 CommPath object 1
#' @param object.2 CommPath object 2 for comparison
#' @param select.ident Plot the activated pathways for which cluster or cell type?
#' @param diff.marker.dat Data frame of defferential expression test result from compareMarker; if NULL, compareMarker would be run to get diff.marker.dat
#' @param diff.path.dat Data frame of defferential activation test result from comparePath; if NULL, comparePath would be run to get diff.path.dat
#' @param top.n.path Top n pathways with the smallest adjusted p values to plot
#' @param path.order Sort criteria used to select the top n pathways, either 'p.val' or 'p.val.adj', which represent the original and adjusted p values, or 'diff' which represents the mean (in t test) or median (in wilcox test) difference
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
pathInterPlot.compare <- function(object.1, object.2, select.ident, diff.marker.dat=NULL, diff.path.dat=NULL, top.n.path=5, path.order='p.val.adj', p.thre=0.05, top.n.receptor=10, top.n.ligand=10, dot.ident.col=NULL, dot.ident.size=1, dot.gene.col=NULL, dot.gene.size=1, bar.pathway.col=NULL, label.text.size=1, label.title.size=1){
	options(stringsAsFactors=F)
	# if(is.null(diff.marker.dat)){
	# 	message(paste0('Identifying differentially expressed ligands(receptors) between cluster ',select.ident,' in object 1 and object 2'))
	# 	diff.marker.dat <- compareMarker(object.1, object.2, select.ident, method='wilcox.test', p.adjust='BH', only.posi=FALSE, only.sig=TRUE)
	# }
	if (is.null(diff.path.dat)){
		message(paste0('Identifying differentially activated pathways between cluster ',select.ident,' in object 1 and object 2'))
		diff.path.dat <- comparePath(object.1, object.2, select.ident, method='wilcox.test', p.adjust='BH', min.size=10, only.posi=FALSE, only.sig=TRUE)
	}
	
	all.ident <- object.1@cell.info$Cluster
	if (!is.factor(all.ident)){ all.ident <- factor(all.ident) }
	ident.label <- levels(all.ident)

	### check the input and select significant pathways
	ident.path.dat <- diff.path.dat
	if ('t' %in% colnames(ident.path.dat)){
		ident.path.dat <- subset(ident.path.dat, p.val.adj < p.thre & t > 0)
		ident.path.dat$diff <- ident.path.dat$mean.diff
	}else if ('W' %in% colnames(ident.path.dat)){
		ident.path.dat <- subset(ident.path.dat, p.val.adj < p.thre & median.diff > 0)
		ident.path.dat$diff <- ident.path.dat$median.diff
	}else{
		stop('Please input the integrate acti.path.dat computed from comparePath')
	}

	if (nrow(ident.path.dat)==0){ paste0(stop('There is no significantly up-regulatged pathways for cluster ', select.ident )) }

	if (path.order=='p.val.adj' | path.order=='p.val'){
		ident.path.dat <- ident.path.dat[order(ident.path.dat[,path.order], -ident.path.dat[,'diff']),]
	}else if(path.order=='diff'){
		ident.path.dat <- ident.path.dat[order(-ident.path.dat[,'diff'], ident.path.dat[,'p.val.adj']),]
	}

	if (top.n.path > nrow(ident.path.dat)){
		warning(paste0('There is(are) ', nrow(ident.path.dat),' significant pathway(s) for the selected ident, and the input top.n.path is ', top.n.path))
		top.n.path <- nrow(ident.path.dat)
	}
	ident.path.dat <- ident.path.dat[1:top.n.path, ]
	all.sig.path <- as.vector(ident.path.dat$description)

	### data for plot of receptors and upstream clusters
	receptor.in.path <- ident.path.dat[, 'receptor.in.path']
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
	cur.rep.LR.inten <- cluster.lr.inten(cur.uniq.rep, object.1, select.ident, ident.label, find='ligand')
	cur.rep.max.LR.inten <- apply(cur.rep.LR.inten, 2, function(x){max(x, na.rm=TRUE)})
	cur.rep.max.LR.inten <- cur.rep.max.LR.inten[order(cur.rep.max.LR.inten, decreasing=TRUE)]
	
	if (length(cur.uniq.rep) < top.n.receptor){
		warning(paste0('There is(are) ', length(cur.uniq.rep),' marker receptor(s) in the selected ident associated with the selected pathways, and the input top.n.receptor is ', top.n.receptor))
		top.n.receptor <- length(cur.uniq.rep)
	}
	
	top.rep.name <- names(cur.rep.max.LR.inten)[1:top.n.receptor]
	plot.receptor.to.pathway.dat <- subset(plot.receptor.to.pathway.dat, cur.rep %in% top.rep.name)

	top.rep.LR.inten <- cur.rep.LR.inten[, top.rep.name]
	if (is.null(dim(top.rep.LR.inten))){
		top.rep.LR.inten <- as.matrix(top.rep.LR.inten)
		colnames(top.rep.LR.inten) <- top.rep.name
	}
	plot.ident.to.receprtor.obj1.dat <- melt(top.rep.LR.inten, varnames=c('Cell.From','Receptor'),value.name="LR.inten",  na.rm=FALSE)
	plot.ident.to.receprtor.obj1.dat <- factor.to.character(plot.ident.to.receprtor.obj1.dat)
	
	### LR.inten for obj.2 for receptor
	cur.rep.LR.obj2.inten <- cluster.lr.inten(top.rep.name, object.2, select.ident, ident.label, find='ligand')
	plot.ident.to.receprtor.obj2.dat <- melt(cur.rep.LR.obj2.inten, varnames=c('Cell.From','Receptor'),value.name="LR.inten",  na.rm=FALSE)
	plot.ident.to.receprtor.obj2.dat <- factor.to.character(plot.ident.to.receprtor.obj2.dat)

	LR.1.Dens.rep <- plot.ident.to.receprtor.obj1.dat$LR.inten
	LR.2.Dens.rep <- plot.ident.to.receprtor.obj2.dat$LR.inten

	plot.ident.to.receprtor.obj1.dat$LR.inten.diff <- as.numeric(LR.1.Dens.rep > LR.2.Dens.rep)
	plot.ident.to.receprtor.obj1.dat[which((is.na(LR.1.Dens.rep)) & (!is.na(LR.2.Dens.rep))), 'LR.inten.diff'] <- 0
	plot.ident.to.receprtor.obj1.dat[which((!is.na(LR.1.Dens.rep)) & (is.na(LR.2.Dens.rep))), 'LR.inten.diff'] <- 1
	plot.ident.to.receprtor.obj1.dat <- subset(plot.ident.to.receprtor.obj1.dat, !is.na(LR.inten))
	plot.ident.to.receprtor.dat <- plot.ident.to.receprtor.obj1.dat
	
	### data for plot of ligands and downstream clusters
	ligand.in.path <- ident.path.dat[, 'ligand.in.path']
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
	cur.lig.max.LR.inten <- apply(cur.lig.LR.inten, 2, function(x){max(x, na.rm=TRUE)})
	cur.lig.max.LR.inten <- cur.lig.max.LR.inten[order(cur.lig.max.LR.inten, decreasing=TRUE)]
	
	if (length(cur.uniq.lig) < top.n.ligand){
		warning(paste0('There is(are) ', length(cur.uniq.lig),' marker ligand(s) in the selected ident associated with the selected pathways, and the input top.n.ligand is ', top.n.ligand))
		top.n.ligand <- length(cur.uniq.lig)
	}
	
	top.lig.name <- names(cur.lig.max.LR.inten)[1:top.n.ligand]
	plot.pathway.to.ligand.dat <- subset(plot.pathway.to.ligand.dat, cur.lig %in% top.lig.name)

	top.lig.LR.inten <- cur.lig.LR.inten[, top.lig.name]
	if (is.null(dim(top.lig.LR.inten))){
		top.lig.LR.inten <- as.matrix(top.lig.LR.inten)
		colnames(top.lig.LR.inten) <- top.rep.name
	}
	plot.ligand.to.ident.obj1.dat <- melt(top.lig.LR.inten, varnames=c('Cell.To','Ligand'),value.name="LR.inten",  na.rm=FALSE)
	plot.ligand.to.ident.obj1.dat <- factor.to.character(plot.ligand.to.ident.obj1.dat)

	### LR.inten for obj.2 for ligand
	cur.lig.LR.obj2.inten <- cluster.lr.inten(top.lig.name, object.2, select.ident, ident.label, find='receptor')
	plot.ligand.to.ident.obj2.dat <- melt(cur.lig.LR.obj2.inten, varnames=c('Cell.To','Ligand'),value.name="LR.inten",  na.rm=FALSE)
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
	path.rect.pval <-  ident.path.dat[match(path.uniq.name, ident.path.dat$description),'p.val.adj']

	if (all(path.rect.pval==0)){
		path.rect.pval <- path.rect.length
		warning('All adjusted p values for the selected pathways are 0. The colors of bars representing pathways are adjusted to indicate the mean (in t test) or median (in wilcox test) difference of pathways')
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


	### the upstream and downstream ident and their coordinate, and their color
	up.uniq.ident <- factor(ident.label, levels=ident.label)
	up.uniq.ident <- up.uniq.ident[order(up.uniq.ident, decreasing=TRUE)]
	### n.up.ident is also the n.down.ident
	n.up.ident <- length(up.uniq.ident)
	if (is.null(dot.ident.col)){
		dot.ident.col <- rev(scales::hue_pal(c=100)(n.up.ident))
	}else if(length(dot.ident.col) < n.up.ident){
		stop(paste0('There is(are) ',n.up.ident,' cluster(s), but only ',length(dot.ident.col),' colors provided.'))
	}else{
		dot.ident.col <- rev(dot.ident.col[1:n.up.ident])
	}
	
	### adjust the coordinate of each element
	max.limit <- max(n.up.ident, n.rep, top.n.path, n.lig)
	up.ident.new.coor <- seq(from=1, to=max.limit, length.out=n.up.ident)
	rep.new.coor <- seq(from=1, to=max.limit, length.out=n.rep)
	pathway.new.coor <- seq(from=1, to=max.limit, length.out=top.n.path)
	lig.new.coor <- seq(from=1, to=max.limit, length.out=n.lig)
	down.ident.new.coor <- seq(from=1, to=max.limit, length.out=n.up.ident)


	if (n.up.ident==1){ 
		up.ident.new.coor <- median(1:max.limit) 
		down.ident.new.coor <- median(1:max.limit) 
	}
	if (n.rep==1){ rep.new.coor <- median(1:max.limit) }
	if (top.n.path==1){ pathway.new.coor <- median(1:max.limit) }
	if (n.lig==1){ lig.new.coor <- median(1:max.limit) }
	
	if (length(pathway.new.coor)==1){
		bar_pathway_width <- 0.3
	}else{
		bar_pathway_width <- 0.3 * (pathway.new.coor[2]-pathway.new.coor[1]) * top.n.path / max.limit
	}
	path.name.max.width <- max(strwidth(path.uniq.name, units="inches", cex=1))
	#bar.pathway.length <- bar.pathway.width*path.name.max.width*path.rect.length
	bar.pathway.length <- 2 * path.rect.length / max(path.rect.length)
	bar_pathway <- geom_rect(aes(xmin=6, xmax=6+bar.pathway.length, ymin=pathway.new.coor-bar_pathway_width, ymax=pathway.new.coor+bar_pathway_width), size=1, fill=bar.path.col, alpha=1, show.legend=F)
	point_up_ident <- geom_point(aes(x=1,y=up.ident.new.coor),size=8*dot.ident.size,color=dot.ident.col,shape=16) 
	point_receptor <- geom_point(aes(x=2,y=rep.new.coor),size=rep.size,color=dot.rep.col)
	point_pathway <- geom_point(aes(x=3,y=pathway.new.coor),size=8*dot.ident.size,color='black',shape=21,fill='white')
	point_ligand <- geom_point(aes(x=4,y=lig.new.coor),size=lig.size,color=dot.lig.col) 
	point_down_ident <- geom_point(aes(x=5,y=down.ident.new.coor),size=8*dot.ident.size,color=dot.ident.col,shape=16) 

	### lines between elements
	### lines between up idents and receptors
	line.ident.to.rep.y.coor <- up.ident.new.coor[match(plot.ident.to.receprtor.dat$Cell.From, up.uniq.ident)]
	line.ident.to.rep.yend.coor <- rep.new.coor[match(plot.ident.to.receprtor.dat$Receptor, cur.uniq.rep)]
	#line.ident.to.rep.col <- dot.ident.col[match(plot.ident.to.receprtor.dat$Cell.From, up.uniq.ident)]
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
	line.lig.to.ident.y.coor <- lig.new.coor[match(plot.ligand.to.ident.dat$Ligand, cur.uniq.lig)]
	line.lig.to.ident.yend.coor <- up.ident.new.coor[match(plot.ligand.to.ident.dat$Cell.To, up.uniq.ident)]
	#line.lig.to.ident.col <- dot.ident.col[match(plot.ligand.to.ident.dat$Cell.To, up.uniq.ident)]
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
	label_down_ident <- annotate('text', hjust=0.5, x=5, y=down.ident.new.coor-0.4,label=up.uniq.ident, size=5*label.text.size)
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
