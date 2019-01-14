setGeneric(name = "base_batch_2DCOW",
           def = function(chrom_names, mod_time,
                          segments, max_warp){
             standardGeneric("base_batch_2DCOW")
           })

setMethod(f = "base_batch_2DCOW",
          signature = c("character", "integer"),
          definition = function(chrom_names, mod_time,
                                segments, max_warp){
            reference_chrom <- read_chrom(chrom_names[1],
                                          mod_time = mod_time)
            sample_chrom <- lapply(chrom_names[-1], read_chrom,
                                   mod_time = mod_time)
            multiple_alig <- lapply(sample_chrom, method_TwoDCOW,
                                    ref_chrom = reference_chrom,
                                    segments = segments,
                                    max_warp = max_warp)
            file_names  <- sapply(basename(chrom_names), 
                                  function(x) {strsplit(x,
                                                        split = "[.]")[[1]][1]})
            names(multiple_alig) <- file_names[-1]
            multiple_alig[[length(chrom_names)]] <- reference_chrom
            names(multiple_alig)[length(chrom_names)] <- file_names[1]
            all_chromatograms <- list(reference_chrom, multiple_alig)
            names(all_chromatograms) <- c(file_names[1], "Aligned")
            all_chromatograms$time <- reference_chrom@time
            return(all_chromatograms)
          })

setGeneric(name = "method_batch_2DCOW",
           def = function(chrom_names, mod_time,
                          segments, max_warp){
             standardGeneric("method_batch_2DCOW")
           })
setMethod(f = "method_batch_2DCOW",
          signature = c("character", "integer"),
          definition = function(chrom_names, mod_time,
                                segments, max_warp){
            lst_aligned <- base_batch_2DCOW(chrom_names = chrom_names,
                                            mod_time = mod_time,
                                            segments = segments,
                                            max_warp = max_warp)
            chrom_2DCOW <- new("batch_2DCOW")
            chrom_2DCOW@name <- "batch_2DCOW"
            chrom_2DCOW@mod_time <- mod_time
            chrom_2DCOW@time <- lst_aligned$time
            lst_aligned <- lst_aligned[- length(lst_aligned)]
            chrom_2DCOW@chromatogram <- lst_aligned[[1]]@chromatogram
            chrom_2DCOW@Batch_2DCOW <- lapply(lst_aligned$Aligned,
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
#' @param chrom_names The names of the chromatograms to be aligned,
#'  the first chromatogram name will be considered as the 
#'  reference chromatogram.
#' @param mod_time The modulation time of the experiment.
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
batch_2DCOW <- function(chrom_names, mod_time, segments, max_warp) {
  if (length(chrom_names) < 2)
    stop("At least two chromatograms are needed")
  batch_aligned <- method_batch_2DCOW(chrom_names = chrom_names,
                                      mod_time = mod_time,
                                      segments = segments,
                                      max_warp = max_warp)
  return(batch_aligned)
}