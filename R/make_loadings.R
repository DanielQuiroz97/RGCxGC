setGeneric(name = "base_makelds",
           def = function(Object, D1_cols, D2_rows){
             standardGeneric("base_makelds")
           }
)

setMethod(f = "base_makelds",
          signature = "matrix",
          definition = function(Object, D1_cols, D2_rows){
            m_loadings <- lapply(data.frame(Object),
                                 function(x, D1_cols, D2_rows){
                                   m_lds <- matrix(x, nrow = D2_rows,
                                                   ncol = D1_cols, byrow = T)
                                   m_lds
                                 }, D1_cols = D1_cols, D2_rows = D2_rows)
            m_loadings
          }
)

#' @title Import foreign model loadings
#' @description `make_loading` method takes the loading matrix obtained
#'  by a foreing package and fold them into two-dimensional matrix
#'  
#' @details We strongly recommend to use the plsda function in the mixOmics
#'  package to perform partial least squares-discriminant analysis. The result
#'  of this model is a list containing a loading matrix.
#'  The method start the matrix A with \emph{m} and \emph{n} dimensions.
#'  where \emph{m} is the eigenvalues and \emph{n}
#'  is the number of loadings which the model returns.
#'  
#' @param floadings a numeric matrix with loadings in columns and eigenvalues
#'  in rows
#' @param time a vector of length two with the tota time range of the
#'  chromatographic run
#' @param mod_time the modulation time
#' @param acq_rate the acquisition rate of the mass analyzer. If it is not
#'  known at the moment, it is printed when read_chrom function is performed
#' @export


make_loadings <- function(floadings, time, mod_time, acq_rate){
  D1_cols <- mod_time * acq_rate
  D2_rows <- length(x) / D1_cols
  if ( !(D1_cols * D2_rows == nrow(floadings)) )
    stop("The length of the loadings does not match with the two-dimensional
        chromatogram dimensions")
  m_loadings <-  base_makelds(floadings, D1_cols = D1_cols,
                              D2_rows = D2_rows)
  lds <- new("foreign", loadings = m_loadings,
             time = time, mod_time = mod_time)
  return(lds)
}


