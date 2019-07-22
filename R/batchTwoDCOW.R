setGeneric(name = "base_batch_2DCOW",
           def = function(reference, sample_chroms,
                          segments, max_warp){
             standardGeneric("base_batch_2DCOW")
           })

setMethod(f = "base_batch_2DCOW",
          signature = c("GCxGC", "list"),
          definition = function(reference, sample_chroms,
                                segments, max_warp){
            multiple_alig <- lapply(sample_chroms, method_TwoDCOW,
                                    ref_chrom = reference,
                                    segments = segments,
                                    max_warp = max_warp)
            all_chromatograms <- c(reference, multiple_alig)
            return(all_chromatograms)
          })

setGeneric(name = "method_batch_2DCOW",
           def = function(reference, sample_chroms,
                          segments, max_warp, ref_name){
             standardGeneric("method_batch_2DCOW")
           })
setMethod(f = "method_batch_2DCOW",
          signature = c("GCxGC", "list"),
          definition = function(reference, sample_chroms,
                                segments, max_warp, ref_name){
            lst_aligned <- base_batch_2DCOW(reference = reference,
                                            sample_chroms = sample_chroms,
                                            segments = segments,
                                            max_warp = max_warp)
            names(lst_aligned) <- c(ref_name, names(lst_aligned)[-1])
            chrom_2DCOW <- new("batch_2DCOW")
            chrom_2DCOW@name <- "batch_2DCOW"
            chrom_2DCOW@mod_time <- reference@mod_time
            chrom_2DCOW@time <- reference@time
            chrom_2DCOW@chromatogram <- reference@chromatogram
            chrom_2DCOW@Batch_2DCOW <- lapply(lst_aligned,
                                              function(x) x@chromatogram)
            return(chrom_2DCOW)
          })

#' Two Dimensional COW in batch.
#'
#' `batch_2DCOW` returns the aligned chromatogram in a named list slot.
#' The first chromatogram is considered as the reference. 
#'
#' This is an adaptation of bidimesional COW alignment, first implementated 
#' in MATLAB. This function takes a sample chromatogram to be aligned 
#' to the reference. The argument [segment] will be used to split the whole
#' chromatogram in n parts in the first and the second dimension respectevily.
#' The [max_warp] argument provides de maximum tolerace of the signal
#' transformation as well to the first and the second dimension.
#'
#' @param reference A GCxGC chromatogram wich will be taken as the reference 
#'   chromatogram
#' @param sample_chroms A named list with the sample chromatograms which will
#'  be aligned to the reference chromatogram
#' @param segments A two integer vector with number of segments
#'  which the first and second dimension will be subdivided, respectively.
#' @param max_warp A two intger vector with the maximum warping parameter.
#'  \emph{Name} to merge with the chromatograms.
#' @importFrom methods new is
#' @export
#' @examples
#' 
#' GB08_fl <- system.file("extdata", "08GB.cdf", package = "RGCxGC")
#' GB09_fl <- system.file("extdata", "09GB.cdf", package = "RGCxGC")
#' chrom_nm <- c(GB08_fl, GB09_fl)
#' \donttest{
#' batch_alignment <- batch_2DCOW(chrom_nm, 5L, c(10, 40), c(1, 10))
#' }
#' 
batch_2DCOW <- function(reference, sample_chroms, segments, max_warp) {
  if (length(sample_chroms) < 2)
    stop("At least two sample chromatograms are needed")
  if (length(names(sample_chroms)) != length(sample_chroms))
    stop('A named list must be provided')
  if (!is(sample_chroms, "list"))
    stop("A list must be provided")
  ref_name <- deparse(substitute(reference))
  batch_aligned <- method_batch_2DCOW(reference = reference,
                                      sample_chroms = sample_chroms,
                                      segments = segments,
                                      max_warp = max_warp,
                                      ref_name = ref_name)
  return(batch_aligned)
}