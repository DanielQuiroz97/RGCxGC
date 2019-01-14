setGeneric(name = "base_GCxGC",
           def = function(Object, sam_rate, per_eval) {
             standardGeneric("base_GCxGC")
           }
)

setMethod(f = "base_GCxGC",
          signature = "GCxGC",
          definition = function(Object, sam_rate, per_eval) {
            
            raw_1d_chrom <- RNetCDF::open.nc(Object@name)
            scan_time <- RNetCDF::var.get.nc(raw_1d_chrom,
                                             "scan_acquisition_time")
            scan_first <- scan_time[ seq(length(scan_time) * per_eval) ]
            scan_floor <- floor(scan_first / 60)
            scan_minut <- table(scan_floor)
            scan_minut <- scan_minut[-c(1, length(scan_minut))]
            homogeneous <- all(scan_minut[1] == scan_minut)
            if (!homogeneous) stop("Sampling rate is not homogeneuous")
            
            tic <- RNetCDF::var.get.nc(raw_1d_chrom,
                                       "total_intensity")
            if (missing(sam_rate)){
              sam_rate <- var.get.nc(raw_1d_chrom,
                                     "scan_acquisition_time")[1:2]
              sam_rate <- 1 / abs(diff(sam_rate)) 
            }
            tic_length <- length(tic)
            time_rn <- range(scan_time) / 60
            len_1d <- floor(sam_rate * Object@mod_time)
            len_2d <- floor(tic_length / len_1d)
            if (len_2d * len_1d < tic_length)
              warning(paste('The last', tic_length - len_2d * len_1d,
                            'signals will be omitted'))
            bidim_chrom <- list(chromatogram = matrix(tic, nrow = len_1d,
                                                      ncol = len_2d),
                                time = time_rn)
            return(bidim_chrom)
            
          }
)

setGeneric(name = "Mread_GCxGC",
           def = function(Object, sam_rate, per_eval) {
             standardGeneric("Mread_GCxGC")
           }
)

setMethod(f = "Mread_GCxGC",
          signature = "raw_GCxGC",
          definition = function(Object, sam_rate, per_eval) {
            readed_gc <- base_GCxGC(Object = Object,
                                    sam_rate = sam_rate,
                                    per_eval = per_eval)
            Object@chromatogram <- readed_gc$chromatogram
            Object@time <- readed_gc$time
            return(Object)
          }
)
#' Read bidimensional total ion current chromatogram.
#' 
#' `read_GCxGC` returns a \emph{raw_GCxGC} with the sample name, the modulation
#' time and the bidimensional chromatogram.
#' 
#' This function reads the netCDF file and retrieve the values in the
#' \emph{total_intensity} variable. Then, with the provided sampling rate and
#' modulation time, it is folded into a numerical matrix (bidimensional
#' chromatogram). This function is an adaptation of the presented routine
#' from \insertCite{Skov2008;textual}{RGCxGC}.
#' 
#' For unusual retention times, more than 60 seconds, the chemical equipment
#' converts the measured points per minute in function of sampling rate.
#' 
#' @param name A name of the netCDF file to which the data will be retrieved.
#' @param mod_time The modulation time of the chromatographic run.
#' @param per_eval An integer with the percentage of the run time to be
#'  evaluate, if the sampling rate is homogeneous.
#' @param sam_rate the sampling rate of the equipment. If sam_rate is missing,
#'  the sampling rate is calculated by the dividing one by the
#'  diference of two adjacent scan time.
#' @importFrom RNetCDF open.nc var.get.nc
#' @importFrom methods is
#' @export 
#' @examples
#'  
#' GB08_fl <- system.file("extdata", "08GB.cdf", package = "RGCxGC")
#' GB08 <- read_chrom(GB08_fl, 5L)
#' 
#' @references
#'     \insertAllCited{}
read_chrom <- function(name, mod_time, sam_rate, per_eval = .10){
  chromatogram <- new('raw_GCxGC', name = name, mod_time = mod_time)
  chromatogram <- Mread_GCxGC(chromatogram, sam_rate, per_eval)
  return(chromatogram)
}
