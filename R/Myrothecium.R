#' @title  Microbial metabolism kinetics of \emph{Myrothecium sp.}
#' 
#' @description This object contains six chromatograms of the microbial
#' metabolic kinetics. \emph{Myrothecium sp.} were cultured in CMA.
#' Inoculation was made from Petri dishes with the fully-grown fungal cells and
#' sterile distilled water (to wash the surface of the plate). The plate was 
#' scraped with a sterile glass handle to obtain the spore suspension.
#' The suspension was liquated and the concentration of 2.4 x 105 spores/mL was
#' determined using a Neubauer chamber and an optical microscope. Then, 50 uL
#' of the suspension was inoculated in a flow chamber into the tubes containing
#' the culture medium. Tubes were kept at 25 celsius degree in a growth
#' chamber with 12h of photoperiod. 
#' 
#' A solid-phase microextraction (SPME) assay containing a
#' DVB / CAR / PDMS (Divinylbenzene / Carboxene / Polymethylsiloxane 50/30 mm)
#' fiber was placed into the tube headspace.
#' 
#' A set of columns consisting of HP-5MS 30m x 0.25mm x 0.25 um
#' connected to a Supelcowax 1m x 0.10mm × 0.10um 
#' with a 1m x 0.25mm deactivated silica capillary being allocated between
#' them. In these tests, a modulation period of 5s was used.
#' 
#' For GCxGC-QMS data acquisition, GCMSsolution version 5.3 software was used.
#' The temperature program were 60-165 celcius at 3 celcius/min;
#' 165 celcius - 260celcius at 20 celcius/min; 260 celcius (5 min);
#' flow rate was 0.6 mL/min (Helium 5.0 carrier gas); splitless injection mode,
#' ion source temperature 200 celcius, interface temperature
#' 260 celcius; voltage 0,9 kV; mass range 50-380 m/z; acquisition rate 25Hz and
#' electron ionization (70eV).
#' 
#' The original study was developed by
#' \insertCite{Quiroz-Moreno2020;textual}{RGCxGC}.
#' 
#' @docType data
#' 
#' @usage data(Myrothecium)
#' 
#' @keywords datasets antagonism microbial
#' @importFrom Rdpack reprompt
#' @format A joined_chrom object containing four slots:
#' \describe{
#'   \item{chromatograms}{A named list with the two-dimensional chromatograms}
#'   \item{groups}{The metadata containing two varaibles and six observations}
#'   \item{time}{The retantion time range of the chromatographic run}
#'   \item{mod_time}{The modulation time}
#' }
#' @references 
#'     \insertAllCited{}
"Myrothecium"