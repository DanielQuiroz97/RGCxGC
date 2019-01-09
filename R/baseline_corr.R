setGeneric(name = "base_baselineCorr",
           def = function(Object){
             standardGeneric("base_baselineCorr")
           })

setMethod(f = "base_baselineCorr",
          signature = c("raw_GCxGC"),
          definition = function(Object){
            chrom_bc <- apply(Object@chromatogram, MARGIN = 2, baseline.corr)
            return(chrom_bc)
          })

setGeneric(name = "method_baselineCorr",
           def = function(Object){
             standardGeneric("method_baselineCorr")
           })

setMethod(f = "method_baselineCorr",
          signature = c("raw_GCxGC"),
          definition = function(Object){
            Object@chromatogram <- base_baselineCorr(Object)
            return(Object)
          })
#' @title  Bidimensional baseline correction
#' 
#' @description  `baseline_corr` provides a bidimensional baseline correction
#' using asymetric least squares.
#' 
#' @details This function takes a raw bidimensional chromatogram and performs
#'  the baseline correction  with implemented function in
#'  \code{\link[ptw]{baseline.corr}}  \insertCite{Eilers2004}{RGCxGC}.
#' 
#' @param chromatogram a \emph{raw_GCxGC} object like with provided
#'  \emph{name} and \emph{mod_time} slots.
#'  
#' @importFrom ptw baseline.corr
#' @importFrom methods new
#' @importFrom Rdpack reprompt
#' @export
#' @examples 
#' library(colorRamps)
#' chrom_name <- system.file("extdata", "08GB.cdf", package = "RGCxGC")
#' chrom_2D <- read_chrom(chrom_name, 5L)
#' chrom_bsline <- baseline_corr(chrom_2D)
#' plot(chrom_bsline, nlevels = 150,
#'            color.palette = matlab.like)
#' 
#' @references 
#'     \insertAllCited{}
baseline_corr <- function(chromatogram) {
  preproc_chrom <- new('preproc_GCxGC', name = chromatogram@name,
                       mod_time = chromatogram@mod_time,
                       time = chromatogram@time)
  blcorr_chrom <- method_baselineCorr(chromatogram)
  preproc_chrom@chromatogram <- blcorr_chrom@chromatogram
  return(preproc_chrom)
}