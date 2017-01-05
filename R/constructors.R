
#' Create stageR object
#'
#' Constructor function for \code{\link{stageRClass}}. A stageR class is a class used for stage-wise analysis in high throughput settings.
#' In its most basic form, it consists of a vector of p-values for the screening hypothesis and a matrix of p-values for the confirmation hypotheses.
#'
#' @param pScreen A vector of screening hypothesis p-values.
#' @param pConfirmation A matrix of confirmation hypothesis p-values. The number of rows should be equal to the length of \code{pScreen}.
#' @param pScreenAdjusted logical, indicating whether the supplied p-values for the screening hypothesis have already been adjusted for multiplicity according to the FDR.
#' @param tx2gene Only applicable for transcript-level analysis. A \code{\link[base]{data.frame}} with transcript IDs in the first columns and gene IDs in the second column. The rownames from \code{pConfirmation} must be contained in the transcript IDs from \code{tx2gene}, and the names from \code{pScreen} must be contained in the gene IDs.
#' @return An instance of an object of the \code{\link{stageRClass}}
#' @references Van den Berge K., Soneson C., Robinson M.D., Clement L.
#' @examples
#' # create a \code{\link{stageRClass}} object
#' stageRObj <- stageR(pScreen=runif(10), pConfirmation=matrix(runif(30),nrow=10,ncol=3))
#' stageRObj <- stageRTx(pScreen=runif(10), pConfirmation=matrix(runif(30),nrow=10,ncol=3), tx2gene=data.frame(transcripts=paste0("transcript",1:10),genes=paste0("gene",rep(1:2,each=5))))
#' @name stageR
#' @rdname stageR
#' @export
stageR <- function(pScreen, pConfirmation, pScreenAdjusted=FALSE)
{
  if(length(pScreen)!=nrow(pConfirmation)) stop("The number of screening hypothesis p-values must be equal to the number of rows in pConfirmation.")
  if(!identical(as.character(names(pScreen)),as.character(rownames(pConfirmation)))) warning("The features (names) in pScreen are not identical to the features (rownames) in pConfirmation.")
  stageR <- new("stageR")
  stageR@pScreen <- pScreen
  stageR@pConfirmation <- pConfirmation
  stageR@pScreenAdjusted <- pScreenAdjusted
  stageR@adjusted <- FALSE
  return(stageR)
}

#' @name stageRTx
#' @rdname stageR
#' @export
stageRTx <- function(pScreen, pConfirmation, pScreenAdjusted=FALSE, tx2gene){
  if(any(is.na(match(rownames(pConfirmation),tx2gene[,1])))) stop("not all transcript names in pConfirmation match with a transcript ID from the tx2gene object.")
  if(any(is.na(match(names(pScreen),tx2gene[,2])))) stop("not all gene names in pScreen match with a gene ID from the tx2gene object.")
  stageR <- new("stageRTx")
  stageR@pScreen <- pScreen
  stageR@pConfirmation <- pConfirmation
  stageR@pScreenAdjusted <- pScreenAdjusted
  stageR@adjusted <- FALSE
  stageR@tx2gene <- tx2gene
  return(stageR)
}

