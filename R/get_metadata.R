#' @export
#' @docType methods
#' @rdname get_metadata-methods
setGeneric(name = "get_metadata",
           def = function(chroms) {
             standardGeneric("get_metadata")
           }
)

#' @title  Method get_metadata
#' @rdname get_metadata-methods
#' @aliases get_metadata,GCxGC-method
#' @description `get_metadata` retrieves the metadata contained in a
#' joined_chrom object.
#' @details This function acceses to the \emph{groups} slot created by the
#'  joined_chrom function. The \emph{Names} are the names of the
#'  chromatograms.
#'  
#' @param chroms a joined_chrom object created by
#'  joined_chrom function.
#' 
#' @examples
#' data(Myrothecium)
#' myr_data <- get_metadata(Myrothecium)
#' myr_data
setMethod(f = "get_metadata", signature = "joined_chrom",
          definition = function(chroms) return(chroms@groups) )
