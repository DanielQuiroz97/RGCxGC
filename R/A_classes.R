####  Parent class GCxGG ####
#' Class GCxGC
#' 
#' Class \emph{GCxGC} defines the superclass of bidimensional comprehensive
#' gas chromatography 
#' 
#' The validity function evaluates if the provied file can be readed in a
#' netCDF format. The validation employs the function 
#' \code{\link[RNetCDF]{open.nc}} to check if the provided file inherits to
#' NetCDF.
#'
#' @slot name the name of a netCDF file to which the data will be retrieved
#' @slot mod_time A integer with the modulation time for the second dimension.
#'   Note the integer should be provide with an \emph{L} at the end of 
#'   the number.
#' @exportClass GCxGC
#' @importFrom RNetCDF open.nc
setClass(Class = "GCxGC", 
         slots = c(name = "character", mod_time = "numeric"),
         validity =  function(object) {
           val_chrom <- RNetCDF::open.nc(object@name)
           if (inherits(val_chrom, "NetCDF")) T
           else paste(chrom_name, "is not a valid NetCDF file")
         })

#### inherited class  raw_GCxGC ####
#' Subclass raw_GCxGC
#' 
#' Subclass \emph{raw_GCxGC} are contained in \emph{GCxGC} super class. It
#' contains a dedicated slot to storage the folded bidimensional chromatogram.
#' 
#' In the first creation of a \emph{raw_GCxGC} object, the slot for the
#' chromatogram. To read and fold the chromatogram use the function
#' \code{\link{read_chrom}}.
#' 
#' @slot chromatogram a numeric matrix.
#' @slot time a vector of to elements with the range of the first dimenstion
#'   retention time
#' @exportClass raw_GCxGC
setClass(Class = "raw_GCxGC", slots = c(chromatogram = "matrix",
                                        time = "vector"),
         contains = "GCxGC")

#### inherited class  preproc_GCxGC ####
#' Subclass preproc_GCxGC
#' 
#' Subclass \emph{preproc_GCxGC} are contained in \emph{raw_GCxGC} super class.
#' It contains a dedicated slot to storage the preprocessed bidimensional
#' chromatogram.
#' 
#' After reading a bidimensional chromatogram, you can perform serveral
#'  preprocessing technicas as smothing, or baseline correction. It will
#' create an object of subclass preproc_GCxGC.
#' 
#' @exportClass preproc_GCxGC
setClass(Class = "preproc_GCxGC", contains = "raw_GCxGC")

#### inherited class  aligned_GCxGC ####
#' Subclass aligned_GCxGC
#' 
#' Subclass \emph{aligned_GCxGC} are contained in \emph{raw_GCxGC} super class.
#' It is not contained in the \emph{prepec_GCxGC} due to raw chromatograms can
#' be aligned without a previous preproccesing technique Although it can
#' improve the performance of the aligment, it is not mandatory.
#' 
#' You can perform the aligment after some preprocessing technique as:
#' baseline correction, or signal smoothing to imporve the performance of the
#' alignment function.
#' 
#' @exportClass aligned_GCxGC
setClass("aligned_GCxGC", contains = 'raw_GCxGC')

#### inherited class  batch_2DCOW ####
#' Subclass batch_2DCOW
#' 
#' Subclass \emph{batch_2DCOW} are contained in \emph{raw_GCxGC} super
#' class. \emph{batch_2DCOW} contains multiple aligned chromatograms, which the
#' first one is consider as the reference chromatogram.
#' 
#' You can perform the alignment after some preprocessing technic as:
#' baseline correction, or signal smothing to imporve the performance of the
#' aligment function, or with the raw chromatogram.
#' 
#' @exportClass batch_2DCOW
setClass("batch_2DCOW", slots = c(Batch_2DCOW = "list"),
         contains = 'raw_GCxGC')

####  Parent joined_chrom ####
#' Class joined_chrom
#' 
#' Class \emph{joined_chrom} defines the superclass to gather single
#' chromatogram as well batch cromatograms into a singe list previos
#' multiway principal compoment analysis
#'  
#' @slot chromatograms A named list with all chromatograms.
#' @slot time The time range of chromatographic run
#' @slot groups A data.frame containing the experiment metadata with
#'  a column named as \emph{Name}
#' @slot mod_time modulation time of the second dimension
#' @exportClass joined_chrom
setClass("joined_chrom", slots = c(chromatograms = "list",
                                   groups = "data.frame",
                                   time = "vector",
                                   mod_time = "numeric"),
         validity = function(object){
           if (is(objetc, "raw_GCxGC")) TRUE
           else print("A chromatogram of class raw_GCxGC or 
                      preproc_GCxGC is needed")
         })
####  Parent MPCA ####
#' Class MPCA
#' 
#' Class \emph{MPCA} defines the superclass of Multiway Principal Component
#' Analysis
#' 
#' @slot scores A matrix with the eigenvalues of projected chromatograms
#'  into principal components space.
#' @slot loadings The eigenvectors of each principal component.
#' @slot summary The summary of the multiway principal component analysis.
#' @slot groups A data.frame with the experiment metadata. It must have a column
#' \emph{Names} to join with chromatograms.
#' @slot time The time range of chromatographic run
#' @slot mod_time modulation time of the second dimension
#' @exportClass MPCA
setClass("MPCA", slots = c(scores = "matrix", loadings = "list",
                           time = "vector", mod_time = "numeric",
                           summary = "list", groups = "data.frame"))

#### Parent foring loadings ####
#' Class foreing_load
#' 
#' Class \emph{forein_load} defines the superclass for the resulting loadings
#'  from foreing classification models, especially for PLS-DA model from
#'  mixOmics package
#'
#' @slot loadings A matrix with the eigenvectors of each latent variable
#' @slot time A two numeric vector of leght two with the range of the
#'  chromatographic run time
#' @slot mod_time The modulation period of the chromatographic experiment
#' @slot sampling_rate The sampling rate of the mass analyzer
#' @exportClass foreign_load
setClass("foreign_load", slots = c(loadings = "list", time = "vector",
                                   mod_time = "numeric"))