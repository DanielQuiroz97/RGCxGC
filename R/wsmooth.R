setGeneric(name = "base_smooth",
           def = function(Object, penalty, lambda){
             standardGeneric("base_smooth")
           })

setMethod(f = "base_smooth",
          signature = "raw_GCxGC",
          definition = function(Object, penalty, lambda){
            if (penalty > 2 | penalty < 1)
              stop("Only 1 or 2 penalty are allowed")
            if (penalty == 1){
              processed_GCxGC <- apply(Object@chromatogram, MARGIN = 2,
                                       ptw::whit1, lambda = lambda)
              return(processed_GCxGC)
            } else{
              processed_GCxGC <- apply(Object@chromatogram, MARGIN = 2,
                                       ptw::whit2, lambda = lambda)
              return(processed_GCxGC)
            }
          })

setGeneric(name = "method_smooth",
           def = function(Object, penalty, lambda){
             standardGeneric("method_smooth")
           })

setMethod(f = "method_smooth",
          definition = function(Object, penalty, lambda){
            if (!is(Object, "raw_GCxGC")) 
              stop("The provided object is not a raw_GCxGC")
            Object@chromatogram <- base_smooth(Object, penalty, lambda)
            return(Object)
          })
#' @title  Two-dimensional smoothing
#' 
#' @description  `wsmooth` weighted whittaker  smoothing.
#' 
#' @details This function takes a raw two-dimensional chromatogram and performs
#'  the weighted wittaker smoothing. It smooths the signal with linear or
#'  quadratic penalty, depending on the provided penalty,
#'  along side the first dimension, based on Whittaker smoother
#'  \insertCite{Eilers2003}{RGCxGC}.
#' 
#' @param chromatogram \emph{raw_GCxGC} or \emph{preproc_GCxGC} object with
#'    \emph{name} and \emph{mod_time} slots.
#' @param  penalty an integer of the order of the penalty. Only penalty of
#'   first (penalty = 1) and second  (penalty = 2) order are allowed. By
#'   default, the smooth function is performed with first penalty order.
#' @param lambda smoothing parameter: larger values lead to more smoothing.
#' @param min_int minimum intensity value. The smoothing routine usually
#'   creates low intensity artifacts which can obscure other compounds signals.
#'   The min intensity value replace signals bellow the given value with 0.
#'   For quadrupole mass detector this artifacts may range from 0-100, while
#'   for TOF mass analyzers it can be 0-1e3.
#' @importFrom methods new is
#' @export
#' @examples 
#' 
#' chrom_name <- system.file("extdata", "08GB.cdf", package = "RGCxGC")
#' chrom_2D <- read_chrom(chrom_name, 5L)
#' chrom_smooth <- wsmooth(chrom_2D, penalty = 1, lambda = 1e1)
#' plot(chrom_smooth, nlevels = 150,
#'            color.palette = colorRamps::matlab.like,
#'            main = expression(paste(lambda, "= 10, penalty = 1")) )
#' # Remove intensities bellow 1.75e5 (too high)    
#' chrom_smooth2 <- wsmooth(chrom_2D, penalty = 1,
#'                          lambda = 1e1, min_int = 1.75e5)
#' plot(chrom_smooth2, nlevels = 150,
#'            color.palette = colorRamps::matlab.like,
#'            main = expression(paste(lambda,
#'             "= 10, penalty = 1, min_int = 1.75e5")) )
#' @references
#'     \insertAllCited{}
wsmooth <- function(chromatogram, penalty = 1, lambda = 1, min_int = 0) {
  preproc_chrom <- new("preproc_GCxGC", name = chromatogram@name,
                       mod_time = chromatogram@mod_time,
                       time = chromatogram@time)
  smoothed_chrom <- method_smooth(Object = chromatogram,
                                  penalty = penalty,
                                  lambda = lambda)
  if (min_int < 0) {  
    stop("Please, provide a positive minimun intensity")
  } else if (min_int > 0){
    smoothed_chrom@chromatogram[smoothed_chrom@chromatogram < min_int] <- 0
    preproc_chrom@chromatogram <- smoothed_chrom@chromatogram
  } else {
    preproc_chrom@chromatogram <- smoothed_chrom@chromatogram
  }
  
  return(preproc_chrom)
}