
#' Create stageR object
#'
#' Constructor function for \code{\link{stageRClass}}. A stageR class is a class used for stage-wise analysis in high throughput settings.
#' In its most basic form, it consists of a vector of p-values for the screening hypothesis and a matrix of p-values for the confirmation hypotheses.
#'
#' @param pScreen A vector of screening hypothesis p-values.
#' @param pConfirmation A matrix of confirmation hypothesis p-values. The number of rows should be equal to the length of \code{pScreen}.
#' @param pScreenAdjusted logical, indicating whether the supplied p-values for the screening hypothesis have already been adjusted for multiplicity according to the FDR.
#' @return An instance of an object of the \code{\link{stageRClass}}
#' @references Van den Berge K., Soneson C., Robinson M.D., Clement L.
#' @examples
#' # create a \code{\link{stageRClass}} object
#' stageRObj <- buildStageR(pScreen=runif(10), pConfirmation=matrix(runif(30),nrow=10,ncol=3))
#' @name buildStageR
#' @rdname buildStageR
#' @export
buildStageR <- function(pScreen, pConfirmation, pScreenAdjusted=FALSE)
{
  stageR <- new("stageR")
  stageR@pScreen <- pScreen
  stageR@pConfirmation <- pConfirmation
  stageR@pScreenAdjusted <- pScreenAdjusted
  return(stageR)
}
