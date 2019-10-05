setGeneric(name = "base_reference",
           def = function(joined_list){
  standardGeneric("base_reference")
})

setMethod(f = "base_reference",
          signature = "joined_chrom",
          definition = function(joined_list){
            # Unfold chrom
            chroms <- joined_list@chromatograms
            unfolded <- sapply(chroms,
                                FUN = function(x) as.vector(t(x)) )
            # Get the mean for each pixel
            mean_chrom <- apply(unfolded, 1, mean)
            # Get chromaogram dimensions
            n_col <- ncol( chroms[[1]] )
            n_row <- nrow( chroms[[1]] )
            # Fold the census chromatogram
            census_chrom <- matrix(mean_chrom, ncol = n_col,
                                   nrow = n_row, byrow = T)
            census_chrom
          } )

setGeneric(name = "method_reference",
           def = function(chroms){
             standardGeneric("method_reference")
           })

setMethod(f = "method_reference",
          signature = "joined_chrom",
          definition = function(chroms){
            census_chrom <- base_reference(chroms)
            ref_chrom <- new("preproc_GCxGC")
            ref_chrom@chromatogram <- census_chrom
            ref_chrom@time <- chroms@time
            ref_chrom@mod_time <- chroms@mod_time
            ref_chrom
          })

#' Make reference chromatogram
#' 
#' `reference_chrom` makes a reference chromatogram by calculating the mean
#' of multiple chromatograms.
#' 
#' The aim of this function is to create a consensus chromatogram to be used
#' as a reference one in the peak alignment process. In other words, the new
#' chromatogram will be used as a template and other chromatograms will be
#' aligned against it. This function overlap the provided chromatograms and
#' compute the mean for each pixel in order to create a reference 
#' two-dimensional chromatogram.
#' 
#' @param chromatograms a joined_chrom object.
#' @export
#' @examples 
#' 
#' # Read chromatogram 1
#' GB08_fl <- system.file("extdata", "08GB.cdf", package = "RGCxGC")
#' MTBLS08 <- read_chrom(GB08_fl, mod_time = 5)
#' 
#' # Read chromatogram 2
#' GB09_fl <- system.file("extdata", "09GB.cdf", package = "RGCxGC")
#' MTBLS09 <- read_chrom(GB09_fl, mod_time = 5)
#' 
#' # Join chromatograms
#' joined <-  join_chromatograms(MTBLS08, MTBLS09)
#' reference <- reference_chrom(joined)
#' plot(reference, main = "Reference chromaogram")

reference_chrom <- function(chromatograms) {
  ref_crom <- method_reference(chromatograms)
  ref_crom
}