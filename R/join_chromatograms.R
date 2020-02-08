setGeneric(name = "base_joinChrom",
           def = function(x)
             standardGeneric("base_joinChrom"))

setMethod(f = "base_joinChrom",
          signature = "raw_GCxGC",
          definition = function(x){
            if (!is(x, "batch_2DCOW")) { 
              chrom_name <- deparse(substitute(x))
              joinedChrom <- list(x@chromatogram)
              names(joinedChrom) <- chrom_name
            } else {
              chrom_names <- names(x@Batch_2DCOW)
              joinedChrom <- x@Batch_2DCOW
              names(joinedChrom) <- chrom_names
            }
            return(joinedChrom)
          })

setGeneric(name = "method_joinChrom",
           def = function(x, y, groups, xy_names, ...){
             standardGeneric("method_joinChrom")}
           )

setMethod(f = "method_joinChrom",
          signature = "raw_GCxGC",
          definition = function(x, y, groups, xy_names, ...){
            all_chrom <- base_joinChrom(x)
            all_chrom <- c(all_chrom, base_joinChrom(y))
            if (!is(x, "batch_2DCOW") & !is(y, "batch_2DCOW"))
              names(all_chrom) <- xy_names
            others_chrom <- list(...)
            test_chrom <- sapply(others_chrom, is, "GCxGC")
            if (!all(test_chrom))
              stop('Only GCxGC objects are acepted')
            if (length(others_chrom)  > 0) {
              complex_condi <- sapply(others_chrom, is, "batch_2DCOW")
              # Single chromatograms only
              if (sum(!complex_condi) == length(complex_condi)){
                others_nm <- names(others_chrom)
                others_joined <- lapply(others_chrom, base_joinChrom)
                others_joined <- lapply(rapply(others_joined, enquote,
                                               how = "unlist"), eval)
                names(others_joined) <- others_nm
                # Only batch_2DCOW chromatograms
              } else if (sum(complex_condi) == length(complex_condi)){
                others_joined <- lapply(others_chrom, base_joinChrom)
                others_joined <- lapply(rapply(others_joined, enquote,
                                               how = "unlist"), eval)
                # Both, single and batch chromatogramas are provided
              } else {
                others_chrom_batch <- others_chrom[complex_condi]
                others_chrom_basic <- others_chrom[!complex_condi]
                others_joined_batch <- lapply(others_chrom_batch,
                                              base_joinChrom)
                others_batch_nm <- unlist(lapply(others_joined_batch, names))
                others_basic_nm <- names(others_chrom_basic)
                other_joined_basic <- lapply(others_chrom_basic, base_joinChrom)
                others_joined <- c(others_joined_batch, other_joined_basic)
                others_joined <- lapply(rapply(others_joined, enquote,
                                               how = "unlist"), eval)
                names(others_joined) <- c(others_batch_nm, others_basic_nm)
              }
              # Missing names asignation for single chroms
              all_chrom <- c(all_chrom, others_joined)
            }
            joined_chrom <- new("joined_chrom")
            joined_chrom@chromatograms <- all_chrom
            if (!missing(groups)) {
              if (!inherits(groups, "data.frame"))
                stop("A data frame containing groups are requiered")
              chrom_names <- names(all_chrom)
              if(length(chrom_names) != nrow(groups))
                stop("provided chromatograms and groups are
                     not the same length")
              metadata <- data.frame(Names = chrom_names)
              metadata <- merge(metadata, groups, by = "Names")
              joined_chrom@groups <- metadata
            }
            joined_chrom@time <- x@time
            joined_chrom@mod_time <- x@mod_time
            return(joined_chrom)
          })
#' @title  Join multiple two-dimensional chromatograms into a single R object
#' 
#' @description `join_chromatograms` saves chromatograms in a
#'  named list slot. Also, it saves information like metadata and 
#'  retention times.
#' 
#' @param x,y a GCxGC object, either single or batch chromatograms.
#' @param groups a data.frame containing the metadata. It must have a column
#'  named as \emph{Names} to merge with the imported chromatograms.
#' @param ... other GCxGC objects to be merged
#' @importFrom methods is new
#' @export
#' @examples 
#' 
#' GB08_fl <- system.file("extdata", "08GB.cdf", package = "RGCxGC")
#' GB09_fl <- system.file("extdata", "09GB.cdf", package = "RGCxGC")
#' GB08 <- read_chrom(GB08_fl, 5L)
#' GB09 <- read_chrom(GB09_fl, 5L)
#' join_gc <- join_chromatograms(GB08, GB09)
#' metadata <- data.frame(Names = c("GB08", "GB09"),
#'                        Type = c("Control", "Treatment"))
#' join_metadata <- join_chromatograms(GB08, GB09, groups = metadata)
join_chromatograms <- function(x, y, groups, ...) {
  
  if (is(x, "batch_2DCOW") & is(y, "batch_2DCOW")){
    xy_names <- NULL
  } else { 
    xy_names <- c(deparse(substitute(x)), deparse(substitute(y))) 
  }
  joinedChrom <- method_joinChrom(x = x, y = y, groups = groups,
                                  xy_names = xy_names, ... = ...)
  return(joinedChrom)
}