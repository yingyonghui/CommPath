#' This is a plug-in function, aimming to paste a vector of idents into a ',' separated string
#' @param ident.missed Vector of idents
#' @return String with idents pasted
pasteIdent <- function(ident.missed){
	if (length(ident.missed) > 1){ 
		ident.missed <- paste0(paste(ident.missed[1:(length(ident.missed)-1)], collapse=', '), ' and ', ident.missed[length(ident.missed)])
	}
	return(ident.missed)
}

#' this is a plug-in function, aiming to adjust those p values==0
#' @param pvalues pvalues
#' @param scale scale
#' @return p values with inf adjusted
p.remove.inf <- function(pvalues, scale=0.1){
	pvalues <- -log10(pvalues)
	if (all(is.infinite(pvalues))){ 
		pvalues <- rep(0, length(pvalues))
	}else if (is.infinite(max(pvalues))){
		secend.max <- max(pvalues[which(is.finite(pvalues))])
		pvalues[which(is.infinite(pvalues))] <- secend.max + min(pvalues) * scale
	}
	return(pvalues)
}


#' this is a plug-in function, aiming to get the color for the dots representing receptors and ligands
#' @param LRpvalue pvalues of ligands or receptors
#' @param user.set.col colors set by user
#' @return Colors matching each p value
LRcolor <- function(LRpvalue, user.set.col){
	if (is.null(user.set.col)){
		dot.col <- circlize::colorRamp2(breaks=seq(min(LRpvalue),max(LRpvalue), length=2),colors=c("#feb24c", "#bd0026"))(LRpvalue)
	}else if(length(user.set.col)==1){
		dot.col <- user.set.col
	}else{
		dot.col <- circlize::colorRamp2(breaks=seq(min(LRpvalue),max(LRpvalue), length=length(user.set.col)),colors=user.set.col)(LRpvalue)
	}
	return(dot.col)
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

#' This is a plug-in function, aimming to convert the variables in a data.frame from factor to character
#' @param x data.frame with columns named as Cell.From and Cell.To
#' @return data.frame with variables converted to character
factor.to.character <- function(x){
	x[['Cell.From']] <- as.character(x[['Cell.From']])
	x[['Cell.To']] <- as.character(x[['Cell.To']])
	return(x)
}
