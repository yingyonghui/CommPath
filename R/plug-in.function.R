#' To paste a vector of idents into a ',' separated string
#' @param ident.missed Vector of idents
#' @return String with idents pasted
pasteIdent <- function(ident.missed){
	if (length(ident.missed) > 1){ 
		ident.missed <- paste0(paste(ident.missed[1:(length(ident.missed)-1)], collapse=', '), ' and ', ident.missed[length(ident.missed)])
	}
	return(ident.missed)
}

#' To adjust those p values==0
#' @param pvalues pvalues
#' @param scale scale
#' @return p values with inf adjusted
p.remove.inf <- function(pvalues, scale=0.1){
	if (all(is.infinite(pvalues))){ 
		pvalues <- rep(0, length(pvalues))
	}else if (is.infinite(max(pvalues))){
		secend.max <- max(pvalues[which(is.finite(pvalues))])
		pvalues[which(is.infinite(pvalues))] <- secend.max + min(pvalues) * scale
	}
	return(pvalues)
}

#' To get the color for the dots representing receptors and ligands
#' @param LRpvalue pvalues of ligands or receptors
#' @param user.set.col colors set by user
#' @return Colors matching each p value
LRcolor <- function(LRpvalue, user.set.col){
	if (length(LRpvalue) > 1){
		if (is.null(user.set.col)){
			dot.col <- circlize::colorRamp2(breaks=seq(min(LRpvalue),max(LRpvalue), length=2),colors=c("#feb24c", "#bd0026"))(LRpvalue)
		}else if(length(user.set.col)==1){
			dot.col <- user.set.col
		}else{
			dot.col <- circlize::colorRamp2(breaks=seq(min(LRpvalue),max(LRpvalue), length=length(user.set.col)),colors=user.set.col)(LRpvalue)
		}
	}else if (length(LRpvalue) == 1){
		dot.col <- "#bd0026"

	}else{
		dot.col <- c()
	}
	return(dot.col)
}

#' To check whether the specified order of idents by users contain all idents that present in the dataset
#' @param all.ident Vector of all idents that present in the dataset
#' @param order Vector of the specified order of idents by users
#' @return if all idents are not contained, stop the procedure
orderCheck <- function(all.ident, order){
	if (all(all.ident %in% order)){
		order.overlap <- order[which(order %in% all.ident)]
		all.ident <- factor(all.ident, levels=order.overlap)
		return(all.ident)
	}else{
		ident.missed <- all.ident[which(!(all.ident %in% order))]
		ident.missed <- pasteIdent(ident.missed)
		stop(paste0('the ident class ',ident.missed,' may be missed in the input order'))
	}
}

#' To convert the variables in a data.frame from factor to character
#' @param x data.frame with columns named as Cell.From and Cell.To
#' @param column which column
#' @return data.frame with variables converted to character
factor.to.character <- function(x, column=c(1, 2)){
	for (each.col in column){
		x[, each.col] <- as.character(x[, each.col])
	}
	return(x)
}

#' To subset a CommPath object
#' @param object CommPath object
#' @param ident.keep idents to keep
#' @return a subset CommPath object
subsetCommPath <- function(object, ident.keep){
	cell.info <- object@cell.info
	cell.info <- subset(cell.info, Cluster %in% ident.keep)
	object@cell.info <- cell.info
	object@data <- object@data[, rownames(cell.info)]
	object@data <- object@data[, rownames(cell.info)]
	object@LR.marker <- subset(object@LR.marker, cluster %in% ident.keep)
	object@pathway$acti.score <- object@pathway$acti.score[, rownames(cell.info)]
	return(object)
}

#' To extract information from ident.path.dat
#' @param ident.path.dat ident.path.dat
#' @return data.frame
extract.info <- function(ident.path.dat){
	up.ident <- as.vector(unlist(sapply(ident.path.dat$cell.up, function(x){ strsplit(x, split=';') })))
	cur.rep <- as.vector(unlist(sapply(ident.path.dat$receptor.in.path, function(x){ strsplit(x, split=';') })))

	n.elem.each.row <- as.vector(unlist(lapply(sapply(ident.path.dat$receptor.in.path, function(x){ strsplit(x, split=';') }),length)))
	path.all <- rep(ident.path.dat$description, times=n.elem.each.row)
	plot.dat <- data.frame(up.ident=up.ident, cur.rep=cur.rep, path.name=path.all)
	return(plot.dat)
}

#' To extract top n element
#' @param x vector of elements with names
#' @param n top a
#' @return top n elements in x
order.and.top <- function(x, n){
	x <- x[order(x, decreasing=TRUE)]
	return(x[1:n])
}

#' To get interaction intensity between clusters and receptors
#' @param top.rep.name top.rep.name
#' @param object CommPath object
#' @param select.ident select.ident
#' @param ident.label ident.label
#' @param find find
#' @return matrix of clusters * genes, elements are interaction intensity between genes and clusters
cluster.lr.inten <- function(top.rep.name, object, select.ident, ident.label, find) {
	if(find=='ligand'){
		top.rep.LR.inten <- sapply(top.rep.name, function(eachrep){
			up.ligand.dat <- findLigand(object, select.ident=select.ident, select.receptor=eachrep, filter=TRUE)
			max.cluster.sum.LR <- by(data=up.ligand.dat$log2FC.LR, INDICES=up.ligand.dat$cell.from, sum)
			each.cluster.inten <- max.cluster.sum.LR[ident.label]
			names(each.cluster.inten) <- ident.label
			return(each.cluster.inten)
		})
	}else if(find=='receptor'){
		top.rep.LR.inten <- sapply(top.rep.name, function(eachrep){
			up.ligand.dat <- findReceptor(object, select.ident=select.ident, select.ligand=eachrep, filter=TRUE)
			max.cluster.sum.LR <- by(data=up.ligand.dat$log2FC.LR, INDICES=up.ligand.dat$cell.to, sum)
			each.cluster.inten <- max.cluster.sum.LR[ident.label]
			names(each.cluster.inten) <- ident.label
			return(each.cluster.inten)
		})
	}
	
	return(top.rep.LR.inten)
}

#' To transform LR intensity to line width
#' @param x LR intensity
#' @return Line width
LRinten.to.width <- function(x){ return(x/max(x) + 1) }

#' To set up the hjust and vjust of text on axises
#' @param angle rotation angle
#' @param position 'x' or 'y' axis
#' @importFrom ggplot2 element_text
#' @return Theme element_text
rotated.axis.element.text <- function(angle,position='x'){
	angle <- angle[1]
	position <- position[1]
	positions <- list(x=0,y=90,top=180,right=270)
	if(!is.numeric(angle)){
		stop("cell.label.angle must be numeric",call.=FALSE)
	}
	rads  <- (-angle - positions[[ position ]])*pi/180
	hjust <- 0.5*(1 - sin(rads))
	vjust <- 0.5*(1 + cos(rads))
	return(element_text(angle=angle,vjust=vjust,hjust=hjust))
}

#' To scale x with minimum equal to 1 and maximum equal to 2
#' @param x numeric vector
#' @return scaled x with minimum equal to 1 and maximum equal to 2
scale_1 <- function(x){
	if (length(which(!is.na(x)))==1) {
		x.scale <- rep(NA, length(x))
		x.scale[which(!is.na(x))] <- 1
	}else{
		x.min <- min(x, na.rm=TRUE)
		x.max <- max(x, na.rm=TRUE)
		x.scale <- (x - x.min) / (x.max - x.min) + 1
	}
	return(x.scale)
}

#' To retrieve the available statistical measures for pathways
#' @param object CommPath object
#' @return print available statistical measures for pathways
getPathAttr <- function(object){
	if(is.null(object@pathway.net)){
		stop('Please run "pathNet" before run getPathAttr')
	}

	if('t.path' %in% colnames(object@pathway.net)){
		print(c('mean.diff','mean','t','P.val','P.val.adj'))
	}else{
		print(c('median.diff','median','W','P.val','P.val.adj'))
	}

}

#' To subset a matrix while keeping the dimension
#' @param mat matrix
#' @param row.select rows slected
#' @param col.select columns selected
#' @return subset of the matrix
subsetMatrix <- function(mat, row.select=NULL, col.select=NULL){
	if (is.null(row.select)){ row.select <- rownames(mat) }
	if (is.null(col.select)){ col.select <- colnames(mat) }
	new.mat <- as.matrix(mat[row.select, col.select])
	if (nrow(mat)==1 | length(row.select)==1){ new.mat <- t(new.mat) }
	rownames(new.mat) <- row.select
	colnames(new.mat) <- col.select
	return(new.mat)
}