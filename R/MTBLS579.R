#' @title  Chromatograms from Dioagnostic Metabolite Biomarkers 
#'         of Chronic Typhoid Carraige study
#'         
#' 
#' @description   The dataset was retrieved from MetaboLights with the 
#' identifier number MTBLS79 \url{https://www.ebi.ac.uk/metabolights/MTBLS579}.
#' Two groups from the entire study was downloaded: control and \emph{S. typhi}
#' carriage. The name files of control group are: 08GB, 09GB, and 14GB,
#' which has the following native name 08_GB.cdf, 14_GB.cdf, and 09_GB.cdf in the
#' MetaboLights database. For the \emph{S. typhi} group the names are:
#' 34GB, 24GB, 29GB, which has the native name of 34_GB.cdf, 24_GB.cdf and
#' 29_GB.cdf in MetaboLights database.
#' 
#' Due to the large size of chromatograms, these data is a subset of the whole
#' chromatograms from 7 min to 18 min of chromatografic run. If you would
#' like to acces the whole formated chromatograms, please go to 
#' \url{https://github.com/DanielQuiroz97/MTBLS579}.
#' 
#' The original study was developed by
#'  \insertCite{Nasstrom2018;textual}{RGCxGC}.
#' 
#' @docType data
#' 
#' @usage data(MTBLS579)
#' 
#' @keywords datasets MetaboLights
#' 
#' @format A joined_chrom object containing four slots:
#' \describe{
#'   \item{chromatograms}{A named list with the two-dimensional chromatograms}
#'   \item{groups}{The metadata containing two varaibles and six observations}
#'   \item{time}{The retantion time range of the chromatographic run}
#'   \item{mod_time}{The modulation time}
#' }
#' 
#' @source \url{https://www.ebi.ac.uk/metabolights/MTBLS579}
#' @references 
#'     \insertAllCited{}
"MTBLS579"