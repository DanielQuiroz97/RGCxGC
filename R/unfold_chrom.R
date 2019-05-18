setGeneric(name = "base_unfold",
           def = function(x)
             standardGeneric("base_unfold"))

setMethod(f = "base_unfold",
          signature = "list",
          definition = function(x){
            if (any(names(x) %in% "")) 
              stop("A named chromatogram must be provided")
            raw_signal <- sapply(x,
                                 function(x) as.vector(t(x)) )
            raw_signal <-  t(raw_signal)
            rownames(raw_signal) <- names(x)
            raw_signal
          })

#' Unfold two-dimensional chromatograms
#' 
#' `unfold` converts two-dimensional chromatograms into a one dimensional
#' vector. Then, all chromatograms are joined into a matrix.
#' 
#' This function takes a single argument, batch_2DCOW or joined_chrom objects
#' and extracts each chromatogram (matix) and it is unfolded to a one dimensional
#' vector. Then, each one dimensional vector is joined in a single matrix. The
#' column matrix represents a given retention time, while row represents samples.
#' Also, in order to keep information about chromatographic runs, the retention
#' times for both dimensions are also exported.
#' 
#' @param Object a batch_2DCOW or joined_chrom objects.
#' @importFrom methods is new
#' @export
#' @examples  
#' 
#' GB08_fl <- system.file("extdata", "08GB.cdf", package = "RGCxGC")
#' GB09_fl <- system.file("extdata", "09GB.cdf", package = "RGCxGC")
#' GB08 <- read_chrom(GB08_fl, 5L)
#' GB09 <- read_chrom(GB09_fl, 5L)
#' join_gc <- join_chromatograms(GB08, GB09)
#' chromatograms_1D <- unfold_chrom(join_gc)
#' chrom_mt <- chromatograms_1D$chromatogram
#' chromatograms_1D$mod_time
#' chromatograms_1D$$time

unfold_chrom <- function(Object) {
  if (!(is(Object, "batch_2DCOW") || is(Object, "joined_chrom")) )
    stop('Only batch aligned and joined chrom objects are allowed')
  if (is(Object, 'batch_2DCOW')){
    unfolded_mt <- base_unfold(Object@Batch_2DCOW)
  } else {
    unfolded_mt <- base_unfold(Object@chromatograms)
  }
  unfolded <- list(chromatogram = unfolded_mt,
                   time = Object@time, mod_time = Object@mod_time)
  unfolded
}