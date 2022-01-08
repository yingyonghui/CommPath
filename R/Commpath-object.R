#' To define the Commpath object and the key slots
#' @slot data normalized expression matrix
#' @slot meta.info vector of cell labels
#' @slot LR.marker data.frame of the differential expression test result of ligands and receptors
#' @slot interact list containing information of LR interaction among clusters
#' @slot pathway list containing information of pathways associated with ligands and receptors
#' @exportClass Commpath
Commpath <- methods::setClass("Commpath",
	slots = c(data = 'matrix',
	meta.info = 'list',
	LR.marker = 'data.frame',
	interact = 'list',
	pathway = 'list')
)

#' To create a Commpath object
#' @param expr.mat Matrix or data frame of expression matrix, with genes in rows and cells in columns
#' @param cell.info data.frame of mata information of cells, and the row names should match the column names of expr.mat
#' @param species Species
#' @return Commpath object
#' @export
createCommpath <- function(expr.mat, cell.info, species){
	### if cell.info is a vector of labels
	cell.info <- data.frame(Cluster=cell.info)
	rownames(cell.info) <- colnames(expr.mat)

	# if (length(species)>1){
	# 	stop("select one species once")
	# }
	# if (!(species %in% c('hsapiens', 'mmusculus', 'rnorvegicus'))){
	# 	stop("select one species from 'hsapiens', 'mmusculus', and 'rnorvegicus'")
	# }

	object <- methods::new(Class="Commpath",
		data = expr.mat,
		meta.info = list(cell.info=cell.info, species=species, logFC.thre=NULL, p.thre=NULL),
		LR.marker = data.frame(),
		interact = list(),
		pathway = list()
	)
	return(object)
}
