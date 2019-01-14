setGeneric(name = "base_MPCA",
           def = function(chrom, center, scale, npcs, ...){
             standardGeneric("base_MPCA")
           })

setMethod(f = "base_MPCA",
          signature = "joined_chrom",
          definition = function(chrom, center, scale, npcs, ...) {
            if (length(chrom@chromatograms) < 2)
              stop("More than one chromatogram is needed")
            chrom_dim <- t(sapply(chrom@chromatograms, dim))
            D1  <- unique(chrom_dim[, 1])
            D2 <- unique(chrom_dim[, 2])
            if (length(D1) > 1 | length(D2) > 1)
              stop('All chromatograms muts have the same dimesions')
            raw_signal <- sapply(chrom@chromatograms,
                                 function(x) as.vector(t(x)) )
            raw_signal <-  t(raw_signal)
            col_0var <- apply(raw_signal, 2, var) != 0
            col_removed <- which(!col_0var)
            raw_signal <- raw_signal[, col_0var]
            pca <- prcomp(raw_signal, center = center, scale. = scale, ...)
            sum_pca <- list(summary = summary(pca))
            lds <-  lapply(seq(npcs), function(i) pca$rotation[, i])
            names(lds) <- paste0("PC", seq(npcs))
            loadings <- list("loadings" = lds, "var_col" = col_removed,
                             "dimension" = c(nrow = D1, ncol = D2))
            mpca <- list(scores = pca$x, loadings = loadings,
                         summ = sum_pca)
            return(mpca)
          })

setGeneric(name = "method_MPCA",
           def = function(chrom, center, scale, npcs, ...){
             standardGeneric("method_MPCA")
           })

setMethod(f = "method_MPCA",
          signature = "joined_chrom",
          definition = function(chrom, center, scale, npcs, ...){
            MPCA <- new("MPCA")
            MPCA_raw <- base_MPCA(chrom = chrom, center = center,
                              scale = scale, npcs = npcs, ...)
            MPCA@scores <- MPCA_raw$scores
            MPCA@loadings <- MPCA_raw$loadings
            MPCA@summary <- MPCA_raw$summ
            if (!is.null(chrom@groups))
              MPCA@groups <- chrom@groups
            MPCA@time <- chrom@time
            MPCA@mod_time <- chrom@mod_time
            return(MPCA)
          })
#' Multiway Principal Component Analysis
#'
#' `MPCA` Performs a multiway principal components analysis on the given
#' bidimensional chromatograms and returns the rusults  as object of class
#' MPCA.
#"
#' Before to perform the calculation, each given chromatogramas are unfolded
#' to a single dimension. All chromatograms are merged and principal component 
#' analysis is performed with the built-in \code{\link[stats]{prcomp}} function.
#' The print method for these objects prints the summary of the analysis.
#' This algorithm was first presented by \insertCite{Wold1987}{RGCxGC}.
#"
#' @param chrom Multiple chromatograms readed or batch aligned
#' @param center A logical value indicating whether the variables should be
#'   shifted to be zero centered. True is set by default and is strongly
#'   seggested not to change to False.   
#' @param scale a logical value indicating whether the variables should
#'  be scaled to have unit variance before the analysis takes place. The
#'  default is True to give the same variable importance in chemometrics.
#' @param npcs an integer indicating how many principals components are
#'  desired to mantain. The default is 3 principal components.
#' @param ... Other argments passed to prcomp function.
#' 
#' @return MPCA returns a list whit class "MPCA" containing the summary of the
#'   analysis, the scores matrix, and unfolded loadings, and the metadata if it
#'   was providen when chromatograms were joined.
#'  
#' @importFrom stats prcomp var
#' @importFrom methods new
#' @export
#' @examples
#' 
#' data(MTBLS579)
#' \donttest{
#' MTBLS579_mpca <- m_prcomp(MTBLS579)
#' print(MTBLS579_mpca)
#' scores(MTBLS579_mpca)
#' plot_loading(MTBLS579_mpca)
#' }
#' 
#' @references 
#'     \insertAllCited{}
m_prcomp <- function(chrom, center = TRUE, scale = TRUE, npcs = 3, ...) {
  multiway_pca <- method_MPCA(chrom = chrom, center = center,
                              scale = scale, npcs = npcs, ...)
  return(multiway_pca)
}