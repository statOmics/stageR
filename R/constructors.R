#' @include stageRClasses.R allGenerics.R

#' @title Create stageR object
#' @description
#' Constructor function for \code{\link{stageRClass}}. A stageR class is a class used for stage-wise analysis in high throughput settings.
#' In its most basic form, it consists of a vector of p-values for the screening hypothesis and a matrix of p-values for the confirmation hypotheses.
#' @param pScreen A vector of screening hypothesis p-values.
#' @param pConfirmation A matrix of confirmation hypothesis p-values. When constructing a \code{\link{stageRClass}} object, the number of rows should be equal to the length of \code{pScreen}. For a \code{\link{stageRTxClass}} object, the dimensions can be different.
#' @param pScreenAdjusted logical, indicating whether the supplied p-values for the screening hypothesis have already been adjusted for multiplicity according to the FDR.
#' @param ... Additional arguments.
#' @return An instance of an object of the \code{\link{stageRClass}}
#' @references
#' Van den Berge K., Soneson C., Robinson M.D., Clement L. (2017). stageR: a general stage-wise method for controlling the gene-level false discovery rate in differential expression and differential transcript usage. Genome Biology 18:151. https://doi.org/10.1186/s13059-017-1277-0
#' @examples
#' # create a \code{\link{stageRClass}} object
#' pScreen <- runif(10)
#' names(pScreen) <- paste0("gene",1:10)
#' pConfirmation <- matrix(runif(30),nrow=10,ncol=3)
#' rownames(pConfirmation) <-  paste0("gene",1:10)
#' stageRObj <- stageR(pScreen=pScreen, pConfirmation=pConfirmation)
#' pConfirmationTx <- matrix(runif(10),ncol=1)
#' names(pScreen) <- paste0("gene",rep(1:2,each=5))
#' stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmationTx, tx2gene=data.frame(transcripts=paste0("transcript",1:10),genes=paste0("gene",rep(1:2,each=5))))
#' @name stageR
#' @rdname stageR
#' @export
#setMethod("stageR", signature=signature(pScreen="numeric", pConfirmation="matrix"),
#         definition=function(pScreen, pConfirmation, pScreenAdjusted=FALSE){
stageR <- function(pScreen, pConfirmation, pScreenAdjusted=FALSE){
  if(length(pScreen)!=nrow(pConfirmation))
    stop("The number of screening hypothesis p-values must be equal to the number of rows in pConfirmation.")
  if(!identical(as.character(names(pScreen)),as.character(rownames(pConfirmation))))
    warning("The features (names) in pScreen are not identical to the features (rownames) in pConfirmation.")
  stageR <- new("stageR")
  stageR@pScreen <- pScreen
  stageR@pConfirmation <- pConfirmation
  stageR@pScreenAdjusted <- pScreenAdjusted
  stageR@adjusted <- FALSE
  return(stageR)
}
#)

#' @title Create stageRTx object.
#' @description
#' Constructor function for \code{\link{stageRTxClass}}. A stageR class is a class used for stage-wise analysis in high throughput settings.
#' In its most basic form, it consists of a vector of p-values for the screening hypothesis, a matrix of p-values for the confirmation hypotheses and a tx2gene object for linking genes to transcripts.
#' @param pScreen A vector of screening hypothesis p-values.
#' @param pConfirmation A matrix of confirmation hypothesis p-values. The number of rows should be equal to the length of \code{pScreen}.
#' @param pScreenAdjusted logical, indicating whether the supplied p-values for the screening hypothesis have already been adjusted for multiplicity according to the FDR.
#' @param tx2gene Only applicable for transcript-level analysis. A \code{\link[base]{data.frame}} with transcript IDs in the first columns and gene IDs in the second column. The rownames from \code{pConfirmation} must be contained in the transcript IDs from \code{tx2gene}, and the names from \code{pScreen} must be contained in the gene IDs.
#' @param ... Additional arguments.
#' @return An instance of an object of the \code{\link{stageRTxClass}}
#' @references
#' Van den Berge K., Soneson C., Robinson M.D., Clement L. (2017). stageR: a general stage-wise method for controlling the gene-level false discovery rate in differential expression and differential transcript usage. Genome Biology 18:151. https://doi.org/10.1186/s13059-017-1277-0
#' @examples
#' # create a \code{\link{stageRClass}} object
#' pScreen <- runif(10)
#' names(pScreen) <- paste0("gene",1:10)
#' pConfirmation <- matrix(runif(30),nrow=10,ncol=3)
#' rownames(pConfirmation) <-  paste0("gene",1:10)
#' stageRObj <- stageR(pScreen=pScreen, pConfirmation=pConfirmation)
#' pConfirmationTx <- matrix(runif(10),ncol=1)
#' names(pScreen) <- paste0("gene",rep(1:2,each=5))
#' stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmationTx, tx2gene=data.frame(transcripts=paste0("transcript",1:10),genes=paste0("gene",rep(1:2,each=5))))
#' @name stageRTx
#' @rdname stageRTx
#' @export
stageRTx <- function(pScreen, pConfirmation, pScreenAdjusted=FALSE, tx2gene){
  if(is.null(names(pScreen))) stop("pScreen does not have names, please set the names of the corresponding genes to the pScreen vector.")
  if(any(is.na(match(rownames(pConfirmation),tx2gene[,1]))))
    stop("not all transcript names in pConfirmation match with a transcript ID from the tx2gene object.")
  if(any(is.na(match(names(pScreen),tx2gene[,2]))))
    stop("not all gene names in pScreen match with a gene ID from the tx2gene object.")
  stageR <- new("stageRTx")
  stageR@pScreen <- pScreen
  stageR@pConfirmation <- pConfirmation
  stageR@pScreenAdjusted <- pScreenAdjusted
  stageR@adjusted <- FALSE
  stageR@tx2gene <- tx2gene
  return(stageR)
}

setValidity("stageR",function(object){
  if(length(pScreen)!=nrow(pConfirmation))
    message("The number of screening hypothesis p-values must be equal to the number of rows in pConfirmation.")

  if(!identical(as.character(names(pScreen)),as.character(rownames(pConfirmation))))
    message("The features (names) in pScreen are not identical to the features (rownames) in pConfirmation.")

  if(any(is.na(getPConfirmation(object))))
    message("NA confirmation stage p-values are not allowed.")

})

setValidity("stageRTx",function(object){
  if(any(is.na(match(rownames(pConfirmation),tx2gene[,1]))))
    message("not all transcript names in pConfirmation match with a transcript ID from the tx2gene object.")

  if(any(is.na(match(names(pScreen),tx2gene[,2]))))
    message("not all gene names in pScreen match with a gene ID from the tx2gene object.")

  if(any(is.na(getPConfirmation(object))))
    message("NA confirmation stage p-values are not allowed.")

})

