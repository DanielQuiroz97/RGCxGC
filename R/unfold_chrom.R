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
#' `unfold` converts the two-dimensional chromatograms into a one dimensional
#' vector. Then, all chromatograms are joined into a matrix.
#' 
#' This function takes a single argument, batch_2DCOW or joined_chrom objects
#' and extracts each chromatogram and then it is unfolded into a one-dimensional
#' vector. Then, each one dimensional vector is joined in a single matrix, where
#' each row represent an observation or a chromatogram and each column
#' represent a variable, in our case, each retention time.
#' Also, in order to keep information about chromatographic runs, the retention
#' times for both dimensions are also exported.
#' 
#' @param Object a batch_2DCOW or joined_chrom objects.
#' @importFrom methods is new
#' @export
#' @examples  
#' 
#' data(Myrothecium)
#' # Unfold 2D chromatogram
#' chrom_1D <- unfold_chrom(Myrothecium)
#' # Retrieve retention time for the first dimension
#' time_1D <- chrom_1D$time
#' # Retrieve the modulation time
#' modulation <- chrom_1D$mod_time

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