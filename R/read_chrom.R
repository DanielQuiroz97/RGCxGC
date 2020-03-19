setGeneric(name = "base_GCxGC",
           def = function(Object, mod_time, sam_rate, per_eval,
                          x_cut, y_cut, verbose) {
             standardGeneric("base_GCxGC")
           }
)

setMethod(f = "base_GCxGC",
          signature = "GCxGC",
          definition = function(Object, mod_time, sam_rate, per_eval,
                                x_cut, y_cut, verbose) {
            
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
            #close.nc(raw_1d_chrom)
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
            # Fold into two-dimensional chromatogram
            chromatogram = matrix(tic, nrow = len_1d, ncol = len_2d)
            # Cut chromatogram
            if (!is.null(x_cut) || !is.null(y_cut)){
              if (!is.null(x_cut)){
                raw_time <- seq(time_rn[1], time_rn[2], length.out = len_2d)
                cut_time <- raw_time[raw_time > x_cut[1] &
                                     raw_time < x_cut[2]]
                if(length(cut_time) == 0)
                  stop('Please provide a congruet time range to cut
                       chromatogram')
                cut_index <- raw_time %in% cut_time
                chromatogram <- chromatogram[, cut_index]
                time_rn <- x_cut
              }
              if (!is.null(y_cut)) {
                raw_time <- seq(0, mod_time, length.out = len_1d)
                cut_time <- raw_time[raw_time > y_cut[1] &
                                     raw_time < y_cut[2]]
                if(length(cut_time) == 0)
                  stop('Please provide a congruent time range to cut
                       chromatogram')
                cut_index <- raw_time %in% cut_time
                chromatogram <- chromatogram[cut_index, ]
                n_mod_time <- y_cut
              }
            }
            bidim_chrom <- list(chromatogram = chromatogram, time = time_rn)
            if (!is.null(y_cut))
              bidim_chrom <- c(bidim_chrom, mod_time = n_mod_time)
            else n_mod_time <- c(0, mod_time)
            if (missing(verbose)) verbose <- TRUE
            if (verbose){
              cat("Retention time ranges:\n")
              cat(paste("1D (min):", round(time_rn[1], 2),
                        round(time_rn[2], 2),"\n"))
              cat(paste("2D (sec):", n_mod_time[1], n_mod_time[2], "\n"))
              cat(paste("Acquisition rate:", round(sam_rate, 0), "\n"))
            }
            return(bidim_chrom)
          }
)

setGeneric(name = "Mread_GCxGC",
           def = function(Object, mod_time, sam_rate, per_eval,
                          x_cut, y_cut, verbose) {
             standardGeneric("Mread_GCxGC")
           }
)

setMethod(f = "Mread_GCxGC",
          signature = "raw_GCxGC",
          definition = function(Object, mod_time, sam_rate, per_eval,
                                x_cut, y_cut, verbose) {
            readed_gc <- base_GCxGC(Object = Object,
                                    mod_time = mod_time,
                                    sam_rate = sam_rate,
                                    per_eval = per_eval,
                                    x_cut = x_cut,
                                    y_cut = y_cut,
                                    verbose = verbose)
            Object@chromatogram <- readed_gc$chromatogram
            Object@time <- readed_gc$time
            if (length(readed_gc) > 2)
              Object@mod_time <- c(readed_gc$mod_time1, readed_gc$mod_time2)
            else Object@mod_time <- c(0, mod_time)
            return(Object)
          }
)
#' Read two-dimensional chromatogram
#' 
#' `read_GCxGC` returns a \emph{raw_GCxGC} with the sample name, the modulation
#' time, the chromatographic time range and the two-dimensional chromatogram.
#' 
#' This function reads the NetCDF file and retrieves values in the
#' \emph{total_intensity} array. Then, with the provided sampling rate and
#' modulation time, the chromatogram is folded into a numerical matrix,
#' representing the two-dimensional chromatogram. This function is an
#' adaptation of the presented routine by \insertCite{Skov2008;textual}{RGCxGC}.
#' 
#' @param name a name of the NetCDF file where the data is alocated.
#' @param mod_time a integer, the modulation time of the chromatographic run.
#' @param sam_rate a integer, the sampling rate of the equipment.
#'  If sam_rate is missing, the sampling rate is calculated by the dividing 1 by
#'  the difference of two adjacent scan time.
#' @param per_eval a double between 0 and 1, with the percentage of the run time
#'  records to be evaluated if the sampling rate is homogeneous.
#' @param x_cut a vector with two elements representing the retention time range
#'  to be mantained in the first dimension.
#' @param y_cut a vector with two elements representing the retention time range
#'  to be mantained in the second dimension.
#' @param verbose a logical indicating if the information of chromatogram is
#'  printted in the console.
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
read_chrom <- function(name, mod_time, sam_rate, per_eval = .10,
                       x_cut = NULL, y_cut = NULL, verbose = TRUE){
  chromatogram <- new('raw_GCxGC', name = name, mod_time = mod_time)
  chromatogram <- Mread_GCxGC(chromatogram, mod_time, sam_rate, per_eval,
                              x_cut, y_cut, verbose = verbose)
  return(chromatogram)
}
