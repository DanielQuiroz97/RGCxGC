#' @export
#' @docType methods
#' @rdname print-methods
setGeneric(name = "print",
           def = function(Object)
             standardGeneric("print"))
#' @title  Print MPCA summary
#' @rdname print-methods
#' @aliases print,MPCA-method
#' @description `print` call the MPCA object to print the summary of this
#'  analysis.
#' 
#' @details  This plot function employs the built-in print function and a
#'  precomputed MPCA summary.
#' 
#' @param Object a MPCA object
#' @exportMethod print
#' @examples 
#' 
#' data(MTBLS579)
#' MTBLS_mpca <- m_prcomp(MTBLS579)
#' print(MTBLS_mpca)
setMethod(f = "print", signature = "MPCA",
          definition = function(Object){
            Object@summary$summary
          })