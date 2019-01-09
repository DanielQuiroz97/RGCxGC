#' @export
#' @docType methods
#' @rdname plot-methods
setGeneric(name = "plot",
           def = function(Object, ...) {
             standardGeneric("plot")
           }
)
#' @title  Method plot
#' @rdname plot-methods
#' @aliases plot,GCxGC-method
#' @description `plot` plot the bidimensional chromatogram as a
#'  filled contour plot
#' 
#' @details  This plot function employs the built-in countour function. As
#'  mentioned in \insertCite{Reichenbach2004;textual}{RGCxGC}, interpolation
#'  is usedto display non-native GCxGC data.
#' 
#' @param Object a GCxGC chromatogram, it could be a raw, or preprocessed
#'   chromatogram
#' @param ... Other parameters passes to \code{\link[graphics]{filled.contour}}
#'  function.
#' @importFrom colorRamps matlab.like2
#' @importFrom graphics filled.contour axis
#' @exportMethod plot
#' @examples 
#' 
#' library(colorRamps)
#' chrom_name <-  system.file("extdata", "08GB.cdf", package = "RGCxGC")
#' chrom_2D <- read_chrom(chrom_name, 5L)
#' plot(chrom_2D, nlevels = 150, color.palette = matlab.like)
#' 
#' @references 
#'     \insertAllCited{}
setMethod(f = 'plot', signature = 'GCxGC',
          definition = function(Object, ...){
            labx <- round(seq(Object@time[1], Object@time[2],
                              length.out = 5), 2)
            laby <- round(seq(0, Object@mod_time, length.out = 5), 2)
            graphics::filled.contour(t(Object@chromatogram),
                           plot.axes = {
                             axis(1, at = seq(0, 1, length.out = 5),
                                  labels = labx)
                             axis(2, at = seq(0, 1, length.out = 5),
                                  labels = laby)
                           }, xlab = "1D min", ylab = "2D sec", ...)
          })