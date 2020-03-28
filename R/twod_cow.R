setGeneric(name = "minMax",
           def = function(data, tresh1 = 1, tresh2){
             standardGeneric("minMax")
           })
setMethod(f = "minMax", signature = "numeric",
          definition = function(data, tresh1 = 1, tresh2){
            data[data < tresh1] <- tresh1
            data[data > tresh2] <- tresh2
            data
          })

setGeneric(name = "base_TwoDCOW",
           def = function(sample_chrom, ref_chrom, segments, max_warp){
             standardGeneric('base_TwoDCOW')
           })

setMethod(f = "base_TwoDCOW", signature = c("matrix", 'matrix'),
          definition = function(sample_chrom, ref_chrom, segments, max_warp) {
            dim_x <- dim(sample_chrom)
            dim_ref <- dim(ref_chrom)
            if ( (dim_x[1] != dim_ref[1]) | (dim_x[2] != dim_ref[2])) {
              print("Two images has different dimensions\n")
              print(paste("Reference image:", dim_x[1], "x", dim_x[2], "\n"))
              print(paste("Work image:", dim_ref[1], "x", dim_ref[2], "\n"))
              row_x <- min(dim_x[1], dim_ref[1])
              col_x <- min(dim_x[2], dim_ref[2])
              sample_chrom <- sample_chrom[seq(1, row_x), seq(1, col_x)]
              ref_chrom <- ref_chrom[seq(1, row_x), seq(1, col_x)]
            } else{
              row_x <- dim_ref[1]
              col_x <- dim_ref[2]
            }
            n_row_segs <- floor( row_x / segments[1] )
            if (n_row_segs == 0) {
              segments[1] <- floor( n_row_segs / 2)
              n_row_segs <- 2
              print("The segment length is forced to be smaller
                    than warping-vector length\n")
              print(paste0("segments[1] = ", segments[1]))
            }
            if (max_warp[1] >= (segments[1] * 0.8)) {
              max_warp <- round(segments[1] * 0.3)
              print(paste("The slack for the 1D was changed to:", max_warp[1]))
            }
            ind_d1 <- seq(1, floor(row_x / segments[1]) *
                            segments[1], by = segments[1] )
            len_d1 <- length(ind_d1)
            if (ind_d1[len_d1] != row_x ){
              ind_d1[len_d1 + 1] <- row_x
              len_d1 <- length(ind_d1)
            }
            n_col_segs <- floor(col_x / segments[2])
            if (n_col_segs == 0){
              segments[2] <- floor(n_col_segs / 2)
              n_col_segs <- 2
              print("The segment lenght is forced to be smaller
                    than warping-vector lenght\n")
              print(paste("segments[2] =", segments[2]))
            }
            if (max_warp[2] >= (segments[2] * 0.8)){
              max_warp[2] <- round(segments[2] * 0.3)
              print(paste("The slack for the 2D was changed to:", max_warp[2]))
            }
            ind_d2 <- seq(1, floor(col_x / segments[2]) *
                            segments[2], by = segments[2])
            len_d2 <- length(ind_d2)
            if (ind_d2[len_d2] != col_x) {
              ind_d2[len_d2 + 1] <- col_x
              len_d2 <- length(ind_d2)
            }
            w_dim <- c(len_d1, len_d2)
            ret_warp_x <- matrix(nrow = w_dim[1], ncol = w_dim[2])
            ret_warp_y <- NULL
            Ind_2D <- list(D1 = ind_d1, D2 = ind_d2)
            for (i in c(1, 2)) {
              tmp_seg <- if (i == 1) rep(ind_d2, 2) else rep(ind_d1, 2)
              tmp_seg <- matrix(tmp_seg, nrow = 2, byrow = T)
              jlen <- ifelse(i == 1, len_d1, len_d2)
              for (j in  seq(jlen)) {
                time_step <-  Ind_2D[[i]][j]
                tmp_ref <- ProfileImage(ref_chrom, i, time_step, segments[i])
                tmp_x <- ProfileImage(sample_chrom, i, time_step, segments[i])
                if ( i == 1){
                  tmp_xw <- cow(tmp_ref, tmp_x, tmp_seg, max_warp[2])
                  ret_warp_x[j, ] <- tmp_xw$Warping
                } else{
                  tmp_xw <- cow(tmp_ref, tmp_x, tmp_seg, max_warp[1])
                  if (j  == 1 ) {
                    ret_warp_y <- tmp_xw$Warping[1, ]
                  } else {
                    index <- seq(w_dim[1] * j + 1, w_dim[1] * (j + 1) )
                    ret_warp_y[index] <- tmp_xw$Warping[1, ]
                  }
                }
              }
            }
            ret_warp_y <- na.omit(ret_warp_y)
            ret_warp_y <- matrix(ret_warp_y, nrow = w_dim[1], ncol = w_dim[2])
            x_warp <- matrix(nrow = row_x, ncol = col_x)
            for (i in seq(1, len_d1 - 1) ) {
              for (j in seq(ind_d1[i], ind_d1[i + 1] - 1) ) {
                tmp_xw <- ret_warp_x[i, ] + (ret_warp_x[i + 1, ] -
                                               ret_warp_x[i, ]) *
                  (j - ind_d1[i]) / (ind_d1[i + 1] - ind_d1[i])
                for (k in seq(1, len_d2 - 1)) {
                  tmp_idxr <- seq(ind_d2[k], ind_d2[k + 1] - 1)
                  tmp_slen <- ind_d2[k + 1] - ind_d2[k]
                  x_warp[j, tmp_idxr] <- tmp_xw[k] + seq(0, tmp_slen - 1) *
                    (tmp_xw[k + 1] - tmp_xw[k] - 1) / (tmp_slen - 1)
                }
                x_warp[j, col_x] <- tmp_xw[n_col_segs + 1]
              }
            }
            tmp_xw <- ret_warp_x[n_row_segs + 1, ]
            for (i in seq(1, len_d2 - 1)) {
              tmp_idxr <- seq(ind_d2[i], ind_d2[i + 1] - 1)
              tmp_slen <- ind_d2[i + 1] - ind_d2[i]
              x_warp[row_x, tmp_idxr] <- tmp_xw[i] + seq(0, tmp_slen - 1 ) *
                (tmp_xw[i + 1] - tmp_xw[i] - 1) / (tmp_slen - 1)
            }
            x_warp[row_x, col_x] <- tmp_xw[n_col_segs + 1]
            y_warp <- matrix(nrow = row_x, ncol = col_x)
            for (i in seq(1, len_d2 - 1)) {
              for (j in seq(ind_d2[i], ind_d2[i + 1] - 1)) {
                tmp_xw <- ret_warp_y[, i] + (ret_warp_y[, i + 1] -
                                               ret_warp_y[, i]) *
                  (j - ind_d2[i]) / (ind_d2[i + 1] - ind_d2[i])
                for (k in seq(1, len_d1 - 1)) {
                  tmp_idxr <- seq(ind_d1[k], ind_d1[k + 1] - 1)
                  tmp_slen <- ind_d1[k + 1] - ind_d1[k]
                  y_warp[tmp_idxr, j] <- tmp_xw[k] + seq(0, tmp_slen - 1) *
                    (tmp_xw[k + 1] - tmp_xw[k] - 1) / (tmp_slen - 1)
                }
                y_warp[row_x, j] <- tmp_xw[n_row_segs + 1]
              }
            }
            tmp_xw <- ret_warp_y[, n_col_segs + 1]
            for (i  in seq(1, len_d1 - 1)) {
              tmp_idrx <- seq(ind_d1[i], ind_d1[i + 1] - 1)
              tmp_slen <- ind_d1[i + 1] - ind_d1[i]
              y_warp[tmp_idrx, col_x] <- tmp_xw[i] + seq(0, tmp_slen - 1) *
                (tmp_xw[i + 1] - tmp_xw[i] - 1) / (tmp_slen - 1)
            }
            y_warp[row_x, col_x] <- tmp_xw[n_row_segs + 1]
            aligned <- matrix(NA_integer_, nrow = row_x, ncol = col_x)
            for (i in seq(1, row_x)) {
              for (j in seq(1, col_x)) {
                r <- y_warp[i, j]
                c <- x_warp[i, j]
                rl <- minMax(floor(r), tresh2 = row_x)
                ru <- minMax(ceiling(r), tresh2 = row_x)
                cl <- minMax(floor(c), tresh2 = col_x)
                cu <- minMax(ceiling(c), tresh2 = col_x)
                x <- c - cl
                y <- r - rl
                aligned[i, j] <- sample_chrom[rl, cl] * (1 - y)  * (1 - x) +
                  sample_chrom[rl, cu] * x * (1 - y) + sample_chrom[ru, cl] *
                  (1 - x) * y + sample_chrom[ru, cu] * x * y
              }
            }
            return(aligned)
          })

setGeneric(name = "method_TwoDCOW",
           def = function(sample_chrom, ref_chrom, segments, max_warp){
             standardGeneric("method_TwoDCOW")
           })
setMethod(f = "method_TwoDCOW", signature = c("raw_GCxGC"),
          definition = function(sample_chrom, ref_chrom, segments, max_warp){
            if (all(sample_chrom@mod_time != ref_chrom@mod_time))
              stop('The modulation time of chromatograms are not the same')
            al_chrom <- new("aligned_GCxGC")
            al_chrom@mod_time <- sample_chrom@mod_time
            al_chrom@name <- sample_chrom@name
            al_chrom@time <- sample_chrom@time
            al_chrom@chromatogram <- base_TwoDCOW(sample_chrom = sample_chrom@chromatogram,
                                                  ref_chrom = ref_chrom@chromatogram,
                                                  segments = segments,
                                                  max_warp = max_warp)
            return(al_chrom)
          })
#' Two-dimensional correlation optimized warping alignment
#'
#' This is an adaptation of two-dimesional COW alignment, first implemented 
#' in MATLAB \insertCite{Tomasi2004}{RGCxGC}. 
#' This functions takes a sample chromatogram to be aligned 
#' against a reference. The argument [segment] will be used to split the whole
#' chromatogram in \emph{n} and \emph{m} parts the first and the second
#' dimension, respectively. The [max_warp] argument provides de maximum
#' tolerance of the signal transformation for the first and the second dimension
#' \insertCite{DabaoZhang2008}{RGCxGC}.
#'
#' @param sample_chrom A GCxGC class chromatogram imported by read_chrom 
#'  function or a preprocessed chromatogram.
#' @param ref_chrom A representative GCxGC chromatogram chosen to be the 
#'   template which sample_chrom will be aligned.
#' @param segments A two integer vector with number of segments
#'  which the first and second dimension will be divided, respectively.
#' @param max_warp A two integer vector with the maximum warping parameter.
#' 
#' @importFrom stats na.omit approx
#' @export
#' @examples
#' 
#' GB08_fl <- system.file("extdata", "08GB.cdf", package = "RGCxGC")
#' GB09_fl <- system.file("extdata", "09GB.cdf", package = "RGCxGC")
#' GB08 <- read_chrom(GB08_fl, 5L)
#' GB09 <- read_chrom(GB09_fl, 5L)
#' \donttest{
#' GB09_al <- twod_cow(sample_chrom = GB09, ref_chrom = GB08,
#'                     segments = c(20, 40), max_warp = c(2, 8))
#' }
#' 
#' @references
#'     \insertAllCited{}
twod_cow <- function(sample_chrom, ref_chrom, segments, max_warp) {
  aligned_chrom <- method_TwoDCOW(sample_chrom = sample_chrom,
                                  ref_chrom = ref_chrom,
                                  segments = segments,
                                  max_warp = max_warp)
  return(aligned_chrom)
}