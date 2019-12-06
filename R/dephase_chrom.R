#' @export
#' @docType methods
#' @rdname dephase_chrom-methods
setGeneric(name = "dephase_chrom",
           def = function(Object, rel_dephase)
             standardGeneric("dephase_chrom"))

#' @title  Method dephase_chrom
#' @rdname dephase_chrom-methods
#' @aliases dephase_chrom,GCxGC-method
#' @description `dephase_chrom` shifts the retention time in the second
#'  dimension of the two-dimensional chromatogram. This procedure is usually
#'  applied in cases when part of peaks is splited in at the final and beginning
#'  of the second dimension. The dephasing procedure is performing by splitting
#'  the chromagram with the relative value provided.
#'  
#' @param Object a GCxGC class object
#' @param rel_dephase a numeric value from 0 to 100 with the relative dephasing
#'    position.
#' @export
#' @examples 
#'  library(colorRamps)
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
