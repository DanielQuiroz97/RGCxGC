#' @export
#' @docType methods
#' @rdname plot-methods
setGeneric(name = "plot",
           def = function(Object, type = "f", ...) {
             standardGeneric("plot")
           }
)
#' @title  Method plot
#' @rdname plot-methods
#' @aliases plot,GCxGC-method
#' @description `plot` plot the two-dimensional chromatogram as a contour plot.
#' 
#' @details  This plot function employs the built-in countour function. As
#'  mentioned in \insertCite{Reichenbach2004;textual}{RGCxGC}, interpolation
#'  is usedto display non-native GCxGC data.
#' 
#' @param Object a GCxGC chromatogram, it can be a raw, or preprocessed
#'   chromatogram.
#' @param type a character indicating the type of chromatogram representation.
#'  By default, type = "f" for \code{\link[graphics]{filled.contour}} function,
#'  if type = "c" only contours or isolines will be displayed by using the
#'   \code{\link[graphics]{contour}} function.
#' @param ... Other parameters passed to \code{\link[graphics]{filled.contour}}
#'  or \code{\link[graphics]{filled.contour}} function, it depends on the value
#'  of the \emph{type} argument.
#' @importFrom colorRamps matlab.like2
#' @importFrom graphics axis filled.contour contour
#' @exportMethod plot
#' @examples 
#' 
#' library(colorRamps)
#' chrom_name <-  system.file("extdata", "08GB.cdf", package = "RGCxGC")
#' chrom_2D <- read_chrom(chrom_name, 5L)
#' plot(chrom_2D, nlevels = 150, color.palette = matlab.like)
#' plot(chrom_2D, type = "c", nlevels = 50, col = matlab.like(30))
#' @references 
#'     \insertAllCited{}
setMethod(f = 'plot', signature = 'GCxGC',
          definition = function(Object, type = "f", ...){
            if ( !(type %in% c("f", "c")) )
              stop("Only f (filled contour) and c (contour) types are allowed")
            labx <- round(seq(Object@time[1], Object@time[2],
                              length.out = 5), 2)
            laby <- round(seq(Object@mod_time[1], Object@mod_time[2],
                              length.out = 5), 2)
            if (type %in% "f"){
            graphics::filled.contour(t(Object@chromatogram),
                           plot.axes = {
                             axis(1, at = seq(0, 1, length.out = 5),
                                  labels = labx)
                             axis(2, at = seq(0, 1, length.out = 5),
                                  labels = laby)
                           }, xlab = "1D min", ylab = "2D sec",
                           ... = ...)
            } else {
              graphics::contour(t(Object@chromatogram),
                                drawlabels = FALSE,
                                frame.plot = TRUE,
                                axes = FALSE,
                                xlab = "1D min", ylab = "2D sec",
                                ... = ...)
              axis(1, at = seq(0, 1, length.out = 5),
                   labels = labx)
              axis(2, at = seq(0, 1, length.out = 5),
                   labels = laby)
            }
          })