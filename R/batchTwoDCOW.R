setGeneric(name = "base_batch_2DCOW",
           def = function(reference, sample_chroms,
                          segments, max_warp, add_ref){
             standardGeneric("base_batch_2DCOW")
           })

setMethod(f = "base_batch_2DCOW",
          signature = c("GCxGC", "list"),
          definition = function(reference, sample_chroms,
                                segments, max_warp, add_ref){
            multiple_alig <- lapply(sample_chroms, method_TwoDCOW,
                                    ref_chrom = reference,
                                    segments = segments,
                                    max_warp = max_warp)
            if (add_ref){
              all_chromatograms <- c(reference, multiple_alig)  
            } else all_chromatograms <- multiple_alig
            
            return(all_chromatograms)
          })

setGeneric(name = "method_batch_2DCOW",
           def = function(reference, sample_chroms,
                          segments, max_warp, ref_name, add_ref){
             standardGeneric("method_batch_2DCOW")
           })
setMethod(f = "method_batch_2DCOW",
          signature = c("GCxGC", "list"),
          definition = function(reference, sample_chroms,
                                segments, max_warp, ref_name, add_ref){
            lst_aligned <- base_batch_2DCOW(reference = reference,
                                            sample_chroms = sample_chroms,
                                            segments = segments,
                                            max_warp = max_warp,
                                            add_ref = add_ref)
            if (add_ref){
              names(lst_aligned) <- c(ref_name, names(lst_aligned)[-1])  
            }
            chrom_2DCOW <- new("batch_2DCOW")
            chrom_2DCOW@name <- "batch_2DCOW"
            chrom_2DCOW@mod_time <- reference@mod_time
            chrom_2DCOW@time <- reference@time
            chrom_2DCOW@chromatogram <- reference@chromatogram
            chrom_2DCOW@Batch_2DCOW <- lapply(lst_aligned,
                                              function(x) x@chromatogram)
            return(chrom_2DCOW)
          })

#' Two-dimensional COW in batch.
#'
#' `batch_2DCOW` perform two-dicmensional correlation optimized warping
#' alignment in batch.
#' 
#' The first argument is the reference chromatogram which other chromatograms
#' will aligned against. Then, a named list is needed for the sample_chroms
#' argumnet. The chromatogram in this list will be aligned using the reference
#' chromatogram. By default, the ference chromatogram will be not included in
#' the subsequent analysis, such as MPCA. If you would like to add the reference
#' chromatogram, then add_ref = T.
#' 
#' This is an adaptation of two-dimesional COW alignment, firstly implemented
#' in MATLAB. This function takes a set of samples chromatogram to be aligned 
#' against to the reference. The argument [segment] will be used 
#' to split the whole chromatogram in \emph{n} and \emph{m} parts in the first
#' and the second dimension, respectively. The [max_warp] argument provides
#' the maximum tolerance of the signal transformation for the first and the
#' second dimension, respectively.
#'
#' @param reference a GCxGC chromatogram wich will be taken as the reference 
#'   chromatogram.
#' @param sample_chroms a named list with the sample chromatograms which will
#'  be aligned against to the reference chromatogram.
#' @param segments a two integer vector with the number of segments
#'  which the first and second dimension will be divided, respectively.
#' @param max_warp a two intger vector with the maximum warping parameter for
#'  the first and second dimension
#' @param add_ref a logical indicating if the reference chromatogram will
#'  be joined together with the sample chromatograms. By the fault add_ref = F.
#'  If add_ref is set T, the provide reference chromatogram will be included as
#'  another sample chromatogram in the downstream analysis.
#' @importFrom methods new is
#' @export
#' @examples
#' 
#' # Read Sample chromatogram
#' GB08_fl <- system.file("extdata", "08GB.cdf", package = "RGCxGC")
#' MTBLS08 <- read_chrom(GB08_fl, mod_time = 5)
#' 
#' # Read reference chromatogram
#' GB09_fl <- system.file("extdata", "09GB.cdf", package = "RGCxGC")
#' MTBLS09 <- read_chrom(GB09_fl, mod_time = 5)
#' 
#' # Create a named list
#' # MTBLS08 will be repeated for exemplification
#' # Considerer that chromatograms are renamed considering the list names
#' batch_samples <- list(Chrom1 = MTBLS08, Chrom2 = MTBLS08)
#' 
#' \donttest{
#' # Perform batch 2DCOW alignment
#' # Add the reference chromatogram as another sample
#' batch_alignment <- batch_2DCOW(MTBLS09, batch_samples,
#'                                c(10, 40), c(1, 10), add_ref = T)
#' # Exclude the reference chromatogram in the sample chromatogram set
#' batch_alignment <- batch_2DCOW(MTBLS09, batch_samples, c(10, 40), c(1, 10))
#' }
#' 
batch_2DCOW <- function(reference, sample_chroms, segments, max_warp,
                        add_ref = F) {
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
                                      ref_name = ref_name,
                                      add_ref = add_ref)
  return(batch_aligned)
}