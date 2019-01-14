#' @export
#' @docType methods
#' @rdname plot_loading-methods
setGeneric(name = "plot_loading",
           def = function(Object, type = "n", pc = 1, ...) {
             standardGeneric("plot_loading")
           }
)
#' @title  Plot two dimensional loadings from MPCA
#' @rdname plot_loading-methods
#' @aliases plot_loading,MPCA-method
#' @description `plot_loading` plot the loadings of the previously MPCA performed.
#' 
#' @details  This function takes the loadings of MPCA and eval if a certain
#'   variable was removed previous compute de MPCA and fill the removed
#'   variables with cero. Then, it plots the loadings of a single
#'   principal component in two dimensions.
#' 
#' @param Object a MPCA object
#' @param type the value type of loadings, \emph{p} for positive, 
#'  \emph{n} for negative, and \emph{b} for negative and positive 
#'  loading values.
#' @param pc the number of the principal component to plot
#' @param ... Other parameters passes to \code{\link[graphics]{filled.contour}}
#'  function.
#' @importFrom colorRamps matlab.like matlab.like2
#' @importFrom graphics filled.contour
#' @exportMethod plot_loading
#' @examples 
#' 
#' library(colorRamps)
#' data(MTBLS579)
#' # MPCA with mean-centered and scaled data
#' MTBLS579_mpca <- m_prcomp(MTBLS579)
#' # Negative loadings of the first principal component
#' \donttest{
#' plot_loading(MTBLS579_mpca, type = "n", pc = 1,
#'              color.palette = matlab.like)
#' # Positive loadings of the first principal component
#' plot_loading(MTBLS579_mpca, type = "p", pc = 1,
#'              color.palette = matlab.like)
#' }
setMethod(f = "plot_loading", signature = "MPCA",
          definition = function(Object, type = "b", pc = 1, ...){
            if (length(pc) > 1)
              stop("Only on principal component has to be provided")
            npcs <- length(Object@loadings$loadings)
            if (pc > npcs)
              stop(paste("There are just", npcs, "PCs"))
            if (!(type %in% c("n", "p", "b")))
              stop("Only negative (n), positive (p), or both (b) loadings
                   types are available")
            loading <- Object@loadings$loadings[[pc]]
            if (type %in% "n"){
              loading[loading > 0] <- 0
            } else if(type %in% "p"){
              loading[loading < 0] <- 0
            }
            D1 <- Object@loadings$dimension[1]
            D2 <- Object@loadings$dimension[2]
            if (length(!Object@loadings$var_col) > 1){
              filled_loading <- vector(length = D1 * D2, "numeric")
              load_it <- 1
              for (i in seq_along(filled_loading)) {
                if (all(i != Object@loadings$var_col)) {
                  filled_loading[i] <- loading[load_it]
                  load_it <- load_it + 1
                }
                mloading <- matrix(filled_loading, byrow = T,
                                   nrow = D1, ncol = D2)
              }
            } else {
              mloading <- matrix(loading, nrow = D1,
                                 ncol = D2, byrow = T)
            }
            labx <- round(seq(Object@time[1], Object@time[2],
                              length.out = 5), 2)
            laby <- round(seq(0, Object@mod_time, length.out = 5), 2)
            filled.contour(t(mloading), 
                           plot.axes = {
                             axis(1, at = seq(0, 1, length.out = 5),
                                  labels = labx)
                             axis(2, at = seq(0, 1, length.out = 5),
                                  labels = laby)
                           }, xlab = "1D min", ylab = "2D sec",...)
          } )