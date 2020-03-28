#' @export
#' @docType methods
#' @rdname import_mixOmics-methods
setGeneric(name = "import_mixOmics",
           def = function(chromatogram, model, mod_time, 
                          time_range, sampling_rt){
             standardGeneric("import_mixOmics")
           })

#' @title Import mixOmics discriminant analysis
#' @rdname import_mixOmics-methods
#' @aliases import_mixOmics,PLSDA-method
#' @description `import_mixOmics` transform a mixOmics discriminant analysis
#'   eigenvalues/eigenvectors into a familiar structure to be handled in the
#'   RGCxGC package.
#' @details This function takes a model built through
#'   \code{\link[mixOmics]{plsda}} and \code{\link[mixOmics]{splsda}}, then,
#'   access to loading values and transform each dimension loading into
#'   a two-dimensional matrix. This matrix represents a typical GCxGC
#'   chromatogram. By default, user can provide a chromatogram where the
#'   required information will be retrieved from. On the other hand,
#'   user can also provide all the needed information
#'   (\emph{mod_time}, \emph{time_range}, \emph{sampling_rt})
#'   to fold eigenvectors into a typical GCxGC chromatogram.
#'   
#' @param chromatogram a typicial GCxGC imported or preprocessed chromatogram.
#' @param model a partial least square discriminant analysis based model, built
#'   by mixOmics package.
#' @param mod_time modulation time of the second dimension
#' @param time_range an atomic vector of lenght two with the time range
#'   of chromatographic  run.
#' @param sampling_rt the sampling rate of the equipment.
#' @importFrom mixOmics splsda tune.splsda
#' @exportMethod import_mixOmics
#' @examples 
#' \donttest{
#' 
#' #### Preparing data ####
#' # Load libraries
#' library(mixOmics)
#' library(caret)
#' # Load chromatograms
#' data(Myrothecium)
#' 
#' # Unfold chromatograms
#' list_chrom <- unfold_chrom(Myrothecium)
#' unfolded_chrom <- list_chrom$chromatogram
#' colnames(unfolded_chrom) <- paste0("RT", seq(dim(unfolded_chrom)[2]))
#' metadata <- get_metadata(Myrothecium)
#' 
#' index <- get_metadata(Myrothecium)
#' # Create a response variable
#' Y <- factor(index$Type)
#' 
#' #### Build the PLS-DA model ####
#' # For reprosucibility
#' set.seed(10)
#  #Tune pls-da
#' chrom_dim <- dim(unfolded_chrom)[2]
#' list.keepX <- seq(chrom_dim/3, chrom_dim, by = 5000)
#' tune.splsda <- tune.splsda(unfolded_chrom, Y, ncomp = 2, validation = 'loo',
#'                           progressBar = TRUE, dist = 'max.dist',
#'                           cpus = 1, # Set cpus according with your pc
#'                           test.keepX = list.keepX)
#' # Number of variables per component
#' tune.splsda$choice.keepX
#' # Remove zero variance predictor variables
#' zero_var <- caret::nearZeroVar(unfolded_chrom)
#' unfolded_chrom <- unfolded_chrom[, -zero_var]
#' splsda_final <- mixOmics::splsda(unfolded_chrom, Y,
#'                                  ncomp = 2, keepX = list.keepX)
#' # Scores
#' scores <- as.data.frame(splsda_final$variates$X)
#' scores$Names <- rownames(scores)
#' scores <- merge(scores, metadata)
#' xyplot(comp2 ~ comp1, data = scores, groups = Type, pch = c(8, 1), cex = 2)
#' }
setMethod(f = "import_mixOmics", signature = "GCxGC",
          definition = function(chromatogram, model, mod_time,
                                time_range, sampling_rt) {
            
            if (!missing(chromatogram)) {
              if ( !(is(model, "mixo_splsda") | is(model, "mixo_plsda")) )
                stop("Please, provide a model builded by plsda or 
                splda function by mixOmics package")
              raw_loadings <- as.data.frame(model$loadings$X)
              chrom <- chromatogram@chromatogram
              n_row <- nrow(chrom)
              n_col <- ncol(chrom)
              loadings_2d <- lapply(raw_loadings,
                                    FUN = function(lds, n_row, n_col){
                matrix(lds, ncol = n_col, nrow = n_row, byrow = TRUE)
              }, n_row = n_row, n_col = n_col  )
              
            } else{
              raw_loadings <- as.data.frame(model$loadings$X)
              len_1d <- floor(sampling_rt * mod_time)
              len_2d <- floor(ncol(raw_loadings) / len_1d)
              loadings_2d <- lapply(raw_loadings,
                                    FUN = function(lds, n_row, n_col){
                matrix(lds,nrow = len_1d, ncol = len_2d, byrow = TRUE)
              }  )
            }
            mod_time <- chromatogram@mod_time
            time <- chromatogram@time
            foreign_model <- new("PLSDA", loadings = loadings_2d,
                                 mod_time = mod_time, time = time)
            foreign_model
            
          })


