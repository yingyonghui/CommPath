#' To present a circos plot
#' @param object Commpath object
#' @param order Vector of orders of the identities arround the circos plot
#' @param col Vector of colors of each identity; names of the col vector are supposed to be assigned to indicate each color for each identity
#' @param ident To highlight the interaction between a specific identity class and others; if 'NULL', plot interaction for all identity classes
#' @param name.vert Should the group annotation be vertical to the grid? Defualt is FALSE
#' @return Circos plot showing the ligand-receptor interaction
#' @export
circosPlot <- function(object, order=NULL, col=NULL, ident=NULL, name.vert=FALSE){
	options(stringsAsFactors=F)

	Interact.num.dat <- object@interact$InteractNumer
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

#' To present a dot plot for specific ligand-receptor pairs in specific clusters
#' @param object Commpath object
#' @param ligand.ident Vector containing the ligand ident
#' @param receptor.ident Vector containing the receptor ident
#' @param ident.levels Vector of levels of the identities 
#' @param top.n.inter Show the dotplot for the top n LR pairs with the largest product of Log2FC
#' @param return.data Logical value indicating whether to return the data for the plot or not
#' @return Dotplot showing the ligand-receptor interaction between the selected ligand.ident and receptor.ident
#' @export
dotPlot <- function(object, ligand.ident=NULL, receptor.ident=NULL, ident.levels=NULL, top.n.inter=10, return.data=FALSE){
	options(stringsAsFactors=F)
	if (is.null(ligand.ident) & is.null(receptor.ident)){
		stop("either ligand.ident or ligand.ident need to be asigned")
	}
	if (length(ligand.ident)>1 & length(receptor.ident)>1){
		stop("specify one cluster for ligand or receptor analysis")
	}

	# get the InteractGene dataframe
	inter.gene.dat <- object@interact$InteractGeneUnfold
	
	if (!is.null(receptor.ident)){
		if (!all(receptor.ident %in% inter.gene.dat$Cell.To)){
			stop(paste0('select receptor.ident from ',pasteIdent(inter.gene.dat$Cell.To)))
		}
		inter.gene.dat <- inter.gene.dat[inter.gene.dat$Cell.To %in% receptor.ident,]
	}
	if (!is.null(ligand.ident)){
		if (!all(ligand.ident %in% inter.gene.dat$Cell.From)){
			stop(paste0('select ligand.ident from ',pasteIdent(inter.gene.dat$Cell.From)))
		}
		inter.gene.dat <- inter.gene.dat[inter.gene.dat$Cell.From %in% ligand.ident,]
	}

	inter.gene.dat$LR.Info <- paste(inter.gene.dat$Ligand, ' --> ', inter.gene.dat$Receptor)

	### to find those LR pairs with largest Log2FC.LR
	max.fc <- by(data=inter.gene.dat$Log2FC.LR, INDICES=inter.gene.dat$LR.Info, function(x){max(x,na.rm=T)})
	max.LRname.fc <- names(max.fc)[order(max.fc, decreasing=TRUE)]
	
	if (top.n.inter > length(max.LRname.fc)){
		warning(paste0('there is(are) ', length(max.LRname.fc),' significant LR pair(s) for the selected ident, and the input top.n.inter is ', top.n.inter))
		top.n.inter <- length(max.LRname.fc)
	}

	max.LRname.fc <- max.LRname.fc[1:top.n.inter]
	inter.gene.dat <- subset(inter.gene.dat, LR.Info %in% max.LRname.fc)

	if (return.data){
		return(inter.gene.dat[,c('Cell.From', 'Cell.To', 'LR.Info', 'Log2FC.LR', 'P.val.LR', 'P.val.adj.LR')])
	}
	
	inter.gene.dat$Log10.P.adj <- p.remove.inf(inter.gene.dat$P.val.adj.LR)
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
		all.ident <- object@meta.info$cell.info$Cluster
		if (!is.factor(all.ident)){ all.ident <- factor(all.ident) }
		ident.levels <- levels(all.ident)
	}
	x.levels <- ident.levels[ident.levels %in% inter.gene.dat$Xaxis]
	if (all(inter.gene.dat$Xaxis %in% x.levels)){
		inter.gene.dat$Xaxis <- factor(inter.gene.dat$Xaxis, levels=x.levels)
	}else{
		ident.missed <- unique(inter.gene.dat$Xaxis[!(inter.gene.dat$Xaxis %in% x.levels)])
		ident.missed <- pasteIdent(ident.missed)
		stop(paste0('the ident class ',ident.missed,' may be missed in the input ident.levels'))
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
		warning('adjusted p values for all LR pairs are 0')
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
#' @param object Commpath object
#' @param all.path.dat data.frame of differential enrichment test result from diffAllPath; if NULL, diffAllPath would be run to get the all.path.dat
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
pathHeatmap <- function(object, all.path.dat=NULL, top.n.pathway=10, sort='p.val.adj', col=NULL, show.legend=FALSE, cell.label.size=NULL, cell.label.angle=0, pathway.label.size=NULL, scale=TRUE, truncation=1.5){
	if (is.null(all.path.dat)){
		all.path.dat <- diffAllPath(object)
	}

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

	ident.label <- object@meta.info$cell.info$Cluster
	gsva.mat <- object@pathway$acti.score

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



#' To present the interact and variable pathways in a line plot, for a specific ident
#' @param object Commpath object
#' @param select.ident select.ident
#' @param ident.path.dat ident.path.dat for the selected ident; if NULL, diffPath would be run to get the ident.path.dat
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
# path.order='p.val.adj';top.n.receptor=10; order=NULL; top.n.path=5; p.thre=0.05; dot.ident.col=NULL; dot.receptor.col=NULL; bar.pathway.col=NULL; bar.pathway.width=10; dot.ident.size=1; dot.receptor.size=1; label.text.size=1; label.title.size=1; line.ident.width=1; line.path.width=1

receptorPathPlot <- function(object, select.ident, ident.path.dat=NULL, path.order='p.val.adj',top.n.receptor=5, order=NULL, top.n.path=10, p.thre=0.05, dot.ident.col=NULL, dot.receptor.col=NULL, bar.pathway.col=NULL, bar.pathway.width=10, dot.ident.size=1, dot.receptor.size=1, label.text.size=1, label.title.size=1, line.ident.width=1, line.path.width=1){
	options(stringsAsFactors=F)
	all.ident <- object@meta.info$cell.info$Cluster
	if (!is.factor(all.ident)){ all.ident <- factor(all.ident) }
	ident.label <- levels(all.ident)

	### check the input and select significant pathways
	if (is.null(ident.path.dat)){
		ident.path.dat <- diffPath(object, select.ident.1=select.ident)
	}
	if ('t' %in% colnames(ident.path.dat)){
		ident.path.dat <- subset(ident.path.dat, p.val.adj < p.thre & t > 0)
	}else if ('W' %in% colnames(ident.path.dat)){
		ident.path.dat <- subset(ident.path.dat, p.val.adj < p.thre & median.diff > 0)
	}else{
		stop('please input the integrate ident.path.dat computed from diffPath')
	}

	if (nrow(ident.path.dat)==0){ paste0(stop('there is no significantly up-regulatged pathways for cluster ', select.ident )) }

	### remove pathways with NA receptor
	ident.path.dat <- subset(ident.path.dat, !is.na(receptor.in.path))
	if (nrow(ident.path.dat)==0){ paste0(stop('there is no marker receptor in the significantly up-regulatged pathways for cluster ', select.ident )) }

	if (path.order=='p.val.adj'){
		ident.path.dat <- ident.path.dat[order(ident.path.dat[,'p.val.adj'], decreasing=FALSE),]
	}

	if (top.n.path > nrow(ident.path.dat)){
		warning(paste0('there is(are) ', nrow(ident.path.dat),' significant pathway(s) for the selected ident, and the input top.n.path is ', top.n.path))
		top.n.path <- nrow(ident.path.dat)
	}
	ident.path.dat <- ident.path.dat[1:top.n.path, ]
	all.sig.path <- as.vector(ident.path.dat$description)

	### data for plot of receptors and upstream clusters
	receptor.in.path <- ident.path.dat[, 'receptor.in.path']
	# if (all(is.na(receptor.in.path))){
	# 	stop('there is no marker receptor in the selected pathways\nselect more pathways and try again')
	# }
	names(receptor.in.path) <- all.sig.path
	# receptor.in.path <- receptor.in.path[which(!is.na(receptor.in.path))]
	cur.rep <- as.vector(unlist(sapply(receptor.in.path, function(x) {strsplit(x, split=';')})))
	path.times <- sapply(receptor.in.path, function(x) {length(strsplit(x, split=';')[[1]])})
	cur.path <- rep(names(path.times), times=path.times)
	plot.receptor.to.pathway.dat <- data.frame(cur.rep=cur.rep, cur.path=cur.path)

	### select the top n receptors
	cur.uniq.rep <- unique(plot.receptor.to.pathway.dat$cur.rep)
	cur.rep.LR.dens <- cluster.lr.dens(cur.uniq.rep, object, select.ident, ident.label, find='ligand')
	cur.rep.max.LR.dens <- apply(cur.rep.LR.dens, 2, function(x){max(x, na.rm=TRUE)})
	cur.rep.max.LR.dens <- cur.rep.max.LR.dens[order(cur.rep.max.LR.dens, decreasing=TRUE)]
	
	if (length(cur.uniq.rep) < top.n.receptor){
		warning(paste0('there is(are) ', length(cur.uniq.lig),' marker ligand(s) in the selected ident associated with the selected pathways, and the input top.n.ligand is ', top.n.ligand))
		top.n.receptor <- length(cur.uniq.rep)
	}
	
	top.rep.name <- names(cur.rep.max.LR.dens)[1:top.n.receptor]
	top.rep.LR.dens <- cur.rep.LR.dens[, top.rep.name]
	plot.receptor.to.pathway.dat <- subset(plot.receptor.to.pathway.dat, cur.rep %in% top.rep.name)
	plot.ident.to.receprtor.dat <- melt(top.rep.LR.dens, varnames=c('Cell.From','Receptor'),value.name="LR.Dens",  na.rm=TRUE)
	plot.ident.to.receprtor.dat <- factor.to.character(plot.ident.to.receprtor.dat)



	### preprocess and extract useful information
	plot.dat <- extract.info(ident.path.dat)
	up.ident <- plot.dat$up.ident
	cur.rep <- plot.dat$cur.rep
	### select the highly expressed receptors
	markerR.dat <- subset(object@interact$markerR, subset=(cluster==select.ident & gene %in% cur.rep))
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
#' @param object Commpath object
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
pathInterPlot <- function(object, select.ident, ident.path.dat=NULL, top.n.receptor=5, order=NULL, top.n.path=10, p.thre = 0.05, top.n.ligand=10, dot.ident.col=NULL, dot.down.ident.col=NULL, dot.receptor.col=NULL, bar.pathway.col=NULL, bar.pathway.width=10, dot.ident.size=1, dot.receptor.size=1, label.text.size=1, label.title.size=1, line.ident.width=1, line.path.width=1, dot.pathway.size=1){
	options(stringsAsFactors=F)
	Interact.num.dat <- subset(object@interact$InteractNumer, LR.Number!=0)
	all.ident <- unique(c(Interact.num.dat$Cell.From,Interact.num.dat$Cell.To))
	if (!is.null(order)){
		all.ident <- all.ident[all.ident %in% order]
		all.ident <- factor(all.ident, levels=order)
		all.ident <- all.ident[order(all.ident)]
	}
	all.col <- scales::hue_pal(c=100)(length(all.ident))
	names(all.col) <- as.character(all.ident)

	### check the input and select significant pathways
	if (is.null(ident.path.dat)){
		ident.path.dat <- diffPath(object, select.ident.1=select.ident)
	}
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
	plot.dat <- extract.info(ident.path.dat)
	up.ident <- plot.dat$up.ident
	cur.rep <- plot.dat$cur.rep
	### select the highly expressed receptors
	markerR.dat <- subset(object@interact$markerR, subset=(cluster==select.ident & gene %in% cur.rep))
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
		markerL.dat <- subset(object@interact$markerL, subset=(cluster==select.ident & gene %in% cur.uniq.lig))
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
	markerL.dat <- subset(object@interact$markerL, subset=(cluster==select.ident & gene %in% cur.uniq.lig))
	markerL.dat <- markerL.dat[match(cur.uniq.lig, markerL.dat$gene),]
	logfc <- markerL.dat$avg_log2FC
	lig.size <- 5 * dot.receptor.size * logfc
	lig.pct <-  markerL.dat$pct.1
	lig.pval <-  p.remove.inf(markerL.dat$p_val_adj)
	dot.ligand.col <- LRcolor(lig.pval, user.set.col=dot.receptor.col)
	line.path.to.lig.col <- dot.ligand.col[match(cur.lig, cur.uniq.lig)]

	### lines between ligand and downstream ident
	plot.down.lig.ident.dat <- do.call(rbind,lapply(cur.uniq.lig, function(x){findReceptor(object=object, select.ident=select.ident, select.ligand=x)}))
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

