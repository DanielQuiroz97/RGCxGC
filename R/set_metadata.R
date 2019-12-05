#' @export
#' @docType methods
#' @rdname set_metadata-methods
setGeneric(name = "set_metadata",
           def = function(Object, metadata)
             standardGeneric("set_metadata"))

#' @title  Set the metadata for a joined_chrom
#' @rdname set_metadata-methods
#' @description `set_metadata` fill metadata slot of a joined chrom object.
#' @aliases set_metadata,GCxGC-method
#' @param Object a joined_chrom object
#' @param metadata a data.frame containing the metadata. It must have a column
#'  named as \emph{Names} to merge with the chromatograms.
#' @examples 
#' 
#' GB08_fl <- system.file("extdata", "08GB.cdf", package = "RGCxGC")
#' GB09_fl <- system.file("extdata", "09GB.cdf", package = "RGCxGC")
#' GB08 <- read_chrom(GB08_fl, 5L)
#' GB09 <- read_chrom(GB09_fl, 5L)
#' extra_info <- data.frame(Names = c("GB08", "GB09"),
#'                        Type = c("Control", "Treatment"))
#' join_chrom <- join_chromatograms(GB08, GB09)
#' join_metadata <- set_metadata(join_chrom, metadata = extra_info)

setMethod(f = 'set_metadata',
          signature = c("joined_chrom", "data.frame"),
          definition = function(Object, metadata){
            Object@groups <- metadata
            return(Object)
          })