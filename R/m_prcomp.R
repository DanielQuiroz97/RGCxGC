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
            if (length(col_removed) != 0){
              raw_signal <- raw_signal[, -col_removed]
              warning(paste0(length(col_removed), "variables with 0 variance were
                             removed"))
            }
            
            if ( sum( is.na(col_0var) ) ){
              warning("NAs wer found and replaced with 0")
              raw_signal[is.na(raw_signal)] <- 0
            }
            pca <- prcomp(raw_signal, center = center, scale. = scale)
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
#' `m_prcomp` Performs a multiway principal components analysis on a given
#' two-dimensional chromatograms and returns the results as object of class
#' MPCA.
#"
#' Before to perform the calculation, each given chromatogramas are unfolded
#' to a single dimension. All chromatograms are merged and principal component 
#' analysis is performed with the built-in \code{\link[stats]{prcomp}} function.
#' The print method for these objects prints the summary of the analysis.
#' This algorithm was first presented by \insertCite{Wold1987}{RGCxGC}.
#"
#' @param chrom Multiple chromatograms read or batch aligned
#' @param center A logical value indicating whether the variables should be
#'   shifted to be zero centered. FALSE is set by default.
#' @param scale a logical value indicating whether the variables should
#'  be scaled to have unit variance before the analysis takes place. The
#'  default is FALSE.
#' @param npcs an integer indicating how many principals components are
#'  desired to maintain. The default is 3 principal components.
#' @param ... Other arguments passed to \code{\link[stats]{prcomp}}
#'  function.
#' 
#' @return MPCA returns a list whit class "MPCA" containing the summary of the
#'   analysis, the scores matrix, unfolded loadings, and the metadata if it
#'   was provided when chromatograms were joined.
#'  
#' @importFrom stats prcomp var
#' @importFrom methods new
#' @export
#' @examples
#' 
#' data(MTBLS579)
#' \donttest{
#' # Perform multiway principal component analysis
#' MTBLS579_mpca <- m_prcomp(MTBLS579, center = TRUE)
#' # Print MPCA summary
#' print(MTBLS579_mpca)
#' # Retrieve MPCA scores
#' scores(MTBLS579_mpca)
#' # Plot bidimensional scores
#' plot_loading(MTBLS579_mpca)
#' }
#' 
#' @references 
#'     \insertAllCited{}
m_prcomp <- function(chrom, center = FALSE, scale = FALSE, npcs = 3, ...) {
  multiway_pca <- method_MPCA(chrom = chrom, center = center,
                              scale = scale, npcs = npcs, ...)
  return(multiway_pca)
}