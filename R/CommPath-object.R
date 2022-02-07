#setClassUnion(name='exprMatrix', members=c("matrix", "dgCMatrix"))

#' To define the CommPath object and the key slots
#' @slot data normalized expression matrix
#' @slot cell.info vector of cell labels
#' @slot meta.info list of important parameters used during the analysis
#' @slot LR.marker data frame of the differential expression test result of ligands and receptors
#' @slot interact list containing information of LR interaction among clusters
#' @slot pathway list containing information of pathways associated with ligands and receptors
#' @exportClass CommPath
CommPath <- methods::setClass("CommPath",
	slots = c(data = 'ANY',
	cell.info = 'data.frame',
	meta.info = 'list',
	LR.marker = 'data.frame',
	interact = 'list',
	pathway = 'list')
)

#' To create a CommPath object
#' @param expr.mat Matrix or data frame of expression matrix, with genes in rows and cells in columns
#' @param cell.info Vector of lables indicating identity classes of cells in the expression matrix, and the order of lables should match the order of cells in the expression matrix; or a data frame containing the meta infomation of cells with the row names matching the cells in the expression matrix and a column named as "Cluster" must be included to indicate identity classes of cells
#' @param species Species
#' @return CommPath object
#' @export
createCommPath <- function(expr.mat, cell.info, species){
	if ((length(species) > 1) | (!(species %in% c('hsapiens', 'mmusculus', 'rnorvegicus')))){
		stop("Select one species from 'hsapiens', 'mmusculus', and 'rnorvegicus'")
	}
	if(is.vector(cell.info) | is.factor(cell.info)){
		if(length(cell.info)!=ncol(expr.mat)){
			stop('The input cell.info is vector, and the length of cell.info should match the ncol of expr.mat!')
		}
		cell.info <- data.frame(Cluster=cell.info)
		rownames(cell.info) <- colnames(expr.mat)
	}else if(is.data.frame(cell.info)){
		if(nrow(cell.info)!=ncol(expr.mat)){
			stop('The input cell.info is a data frame, and the row names of cell.info should match the column names of expr.mat!')
		}
		if(any(rownames(cell.info) !=  colnames(expr.mat))){
			stop('The input cell.info is a data frame, and the row names of cell.info should match the column names of expr.mat!')
		}
	}else{
		stop('Either input a vector or a data frame for cell.info')
	}
	object <- methods::new(Class="CommPath",
		data = as(expr.mat, "dgCMatrix"),
		cell.info = cell.info,
		meta.info = list(species=species, logFC.thre=NULL, p.thre=NULL),
		LR.marker = data.frame(),
		interact = list(),
		pathway = list()
	)

	return(object)
}


#' Show method for CommPath
#' @param object A CommPath object
#' @docType methods
setMethod(f="show", signature="CommPath", definition=function(object) {
	cat("An object of class", class(object), "with\n",
		nrow(x = object@data), "genes *", ncol(x = object@data), "cells.\n"
	)
	invisible(x = NULL)
})
