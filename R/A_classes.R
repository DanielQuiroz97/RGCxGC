####  Parent class GCxGG ####
#' Class GCxGC
#' 
#' Class \emph{GCxGC} defines the superclass of two-dimensional comprehensive
#' gas chromatography 
#' 
#' The validity function evaluates if the provided file can be readed as a 
#' NetCDF file. The validation function employs the function 
#' \code{\link[RNetCDF]{open.nc}} to check if the provided file inherits to
#' NetCDF class.
#'
#' @slot name the name of a NetCDF file to where the data will be retrieved.
#' @slot mod_time a integer with the modulation time for the second dimension.
#' @exportClass GCxGC
#' @importFrom RNetCDF open.nc
setClass(Class = "GCxGC", 
         slots = c(name = "character", mod_time = "numeric"),
         validity =  function(object) {
           val_chrom <- RNetCDF::open.nc(object@name)
           if (inherits(val_chrom, "NetCDF")) {
             #close.nc(val_chrom)
             return(T)
             }
           else paste(chrom_name, "is not a valid NetCDF file")
         })

#### inherited class  raw_GCxGC ####
#' Subclass raw_GCxGC
#' 
#' Subclass \emph{raw_GCxGC} are contained in \emph{GCxGC} super class. It
#' contains a dedicated slot to storage the folded two-dimensional chromatogram.
#' 
#' In the first creation of a \emph{raw_GCxGC} object, the slot for the
#' chromatogram is not created yet. To read and fold the chromatogram 
#' use the \code{\link{read_chrom}}  function.
#' 
#' @slot chromatogram a numeric matrix.
#' @slot time a vector of two elements with the range of the first dimenstion
#'   run time
#' @exportClass raw_GCxGC
setClass(Class = "raw_GCxGC", slots = c(chromatogram = "matrix",
                                        time = "vector"),
         contains = "GCxGC")

#### inherited class  preproc_GCxGC ####
#' Subclass preproc_GCxGC
#' 
#' Subclass \emph{preproc_GCxGC} are contained in \emph{raw_GCxGC} super class.
#' It contains a dedicated slot to storage the preprocessed two-dimensional
#' chromatogram.
#' 
#' After reading a two-dimensional chromatogram, you can perform serveral
#' preprocessing techniques such as smoothing or baseline correction.
#' This action will create an object of a preproc_GCxGC subclass.
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
#' first one is the reference chromatogram.
#' 
#' You can perform the alignment after some preprocessing technic as:
#' baseline correction, or signal smothing to improve the performance of the
#' aligment function, or with the raw chromatogram.
#' 
#' @exportClass batch_2DCOW
setClass("batch_2DCOW", slots = c(Batch_2DCOW = "list"),
         contains = 'raw_GCxGC')

####  Parent joined_chrom ####
#' Class joined_chrom
#' 
#' Class \emph{joined_chrom} defines the superclass to gather single
#' chromatogram, as well as batch chromatograms into a single list, prior to
#' multiway principal compoment analysis or unfolding them.
#'  
#' @slot chromatograms a named list with all chromatograms.
#' @slot time the time range of the chromatographic run
#' @slot groups a data.frame containing the experiment metadata with
#'  a column named as \emph{Names}.
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
####  Parent class for projection methods such as PCA and PLS-DA ####
#' Class projected
#' 
#' The \emph{projected} class defines the superclass for projection methods,
#' specially for multiway principal component analysis and discriminant
#' analysis based on partial least squares. The class represents the 
#' convergence of in-package results (m_prcomp) and the foreing 
#' building model (PLS-DA) procedure.
#' 
#' @slot loadings The eigenvectors of each principal component.
#' @slot time The time range of chromatographic run
#' @slot mod_time modulation time of the second dimension
setClass("projected", slots = c(loadings = "list", time = "vector",
                                mod_time = "numeric"))


#' Subclass MPCA
#' \emph{MPCA} subclass for Multiway Principal Component Analysis object
#' @slot scores A matrix with the eigenvalues of projected chromatograms
#'  into principal components space.
setClass("MPCA", slots = c(scores = "matrix", summary = "list",
                           groups = "data.frame"), contains = "projected")

#' Subclass PLSDA
#' 
#' Class \emph{PLSDA} defines the class to import foreign results of 
#' partial least squares-discriminant analysis performed with
#' mixOmics package.
setClass("PLSDA", contains = "projected")