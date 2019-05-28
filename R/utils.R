
setGeneric(name = "get_metadata",
           def = function(Object)
             standardGeneric("get_metadata"))
#' Extract metadata from a joinded_chrom
#' 
#' `get_metadata` if a metadata was provided when join_chromatogram
#' funcion was used, the matadata will be retrieve from de object.
#' 
#' @param Object a joined_chrom object
#' @export
#' @examples 
#' GB08_fl <- system.file("extdata", "08GB.cdf", package = "RGCxGC")
#' GB09_fl <- system.file("extdata", "09GB.cdf", package = "RGCxGC")
#' GB08 <- read_chrom(GB08_fl, 5L)
#' GB09 <- read_chrom(GB09_fl, 5L)
#' metadata <- data.frame(Names = c("GB08", "GB09"),
#'                        Type = c("Control", "Treatment"))
#' join_metadata <- join_chromatograms(GB08, GB09, groups = metadata)
#' get_metadata(join_metadata)
setMethod(f = "get_metadata",
          signature = "joined_chrom",
          definition = function(Object){
            return(Object@groups)
          })

setGeneric(name = "set_metadata",
           def = function(Object, metadata)
             standardGeneric("set_metadata"))

#' Set the metadata for a joined_chrom
#' 
#' `set_metadata` fill metadata slot of a joined chrom.
#' 
#' @param Object a joined_chrom object
#' @param metadata a data.frame containing the metadata. It must have a column
#' @export
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


setGeneric(name = "dephase_chrom",
           def = function(Object, rel_dephase)
             standardGeneric("dephase_chrom"))

#' Dephase the 2D retention time
#' 
#' `dephase_chrom`  move the retention time in the second dimension of a
#'  bidimensional chromatogram. This procedure is ussually applied in cases
#'  when peaks are near at the final of the second dimension of the 
#'  chromatogram. The dephasing is performing by splitting the chromagram
#'  with the relative value provided.
#'  
#'  @param Object a GCxGC class object
#'  @param rel_dephase a numeric value from 0 to 100 with the relative dephasing
#'    position.
#'  @export
#'  @examples 
#'  
#'  GB08_fl <- system.file("extdata", "08GB.cdf", package = "RGCxGC")
#'  GB08 <- read_chrom(GB08_fl, 5L)
#'  plot(GB08, nlevels = 150, color.palette = matlab.like,
#'       main = "No dephased chromatogram")
#'  GB08_d25 <- dephase_chrom(GB08, 25)
#'  plot(GB08_d25, nlevels = 150, color.palette = matlab.like,
#'       main = "25% dephased chromatogram")
setMethod(f = "dephase_chrom",
          signature = c("GCxGC"),
          definition =  function(Object, rel_dephase) {
            if (rel_dephase < 0  | rel_dephase > 100)
              stop("A relative value from 0 to 100 must be provided")
            chrom_rows <- nrow(Object@chromatogram)
            dephase_index <- seq(1, chrom_rows * rel_dephase / 100)
            chrom1 <- Object@chromatogram[dephase_index, ]
            chrom2 <- Object@chromatogram[-dephase_index, ]
            chrom <- rbind(chrom2, chrom1)
            Object@chromatogram <- chrom
            return(Object)
          })