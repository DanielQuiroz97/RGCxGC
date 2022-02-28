#' @export
#' @docType methods
#' @rdname print_mpca-methods
setGeneric(name = "print_mpca",
           def = function(Object)
             standardGeneric("print_mpca"))
#' @title  Print MPCA summary
#' @rdname print_mpca-methods
#' @aliases print_mpca,MPCA-method
#' @description `print_mpca` call the MPCA object to print the summary of this
#'  analysis.
#' 
#' @details  The plot function employs the built-in print function and a
#'  precomputed MPCA summary to display the explained and cumulative variance
#'  for each principal component.
#' 
#' @param Object a MPCA object
#' @exportMethod print_mpca
#' @examples 
#' 
#' data(MTBLS579)
#' MTBLS_mpca <- m_prcomp(MTBLS579, center = TRUE)
#' print_mpca(MTBLS_mpca)
setMethod(f = "print_mpca", signature = "MPCA",
          definition = function(Object){
            Object@summary$summary
          })