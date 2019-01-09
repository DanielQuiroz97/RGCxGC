#' @export
#' @docType methods
#' @rdname scores-methods
setGeneric(name = "scores",
           def = function(Object)
             setGeneric("scores"))

#' @title  Method plot_scores
#' @rdname scores-methods
#' @aliases scores,MPCA-method
#' @description `scores` exports the scores matrix of the previously MPCA
#'  performed.
#' 
#' @details  This function takes the scores of MPCA and retrieves the score
#'   matrix.
#' 
#' @param Object a MPCA object
#' @exportMethod scores
#' @examples 
#' 
#' data(MTBLS579)
#' # MPCA with mean-centered and scaled data
#' MTBLS579_mpca <- m_prcomp(MTBLS579, center = TRUE, scale = TRUE)
#' # Export scores matrix
#' scores(MTBLS579_mpca)
setMethod(f = "scores",
          signature = "MPCA",
          definition = function(Object){
            scores <- as.data.frame(Object@scores)
            scores$Names <- rownames(scores)
            if (nrow(scores) == nrow(Object@groups)){
              scores <- merge(scores, Object@groups)
              scores
            } else scores
          })