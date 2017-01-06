#' @include stageRClasses.R


print.stageR <- function(object, ...){
    cat("stageR object, containing: \n")
    if(!object@adjusted){
	cat(paste0(length(object@pScreen)," screening hypothesis p-values \n"))
	cat(paste0(ncol(object@pConfirmation)," confirmation hypotheses for ",nrow(object@pConfirmation)," genes \n"))
    }
    if(object@adjusted){
	cat(paste0(length(object@pScreen)," screening hypothesis p-values \n"))
	cat(paste0(ncol(object@pConfirmation)," confirmation hypotheses for",nrow(object@pConfirmation)," genes \n"))
	cat(paste0("adjusted p-values on a ",object@alpha*100,"% OFDR level with the following FWER correction method: ",object@method," \n"))

    }
}

setMethod("show","stageR",function(object) print.stageR(object))


print.stageRTx <- function(object, ...){
  cat("stageRTx object, containing: \n")
  if(!object@adjusted){
    cat(paste0(length(object@pScreen)," screening hypothesis p-values \n"))
    cat(paste0(ncol(object@pConfirmation)," confirmation hypotheses for ",nrow(object@pConfirmation)," genes \n"))
  }
  if(object@adjusted){
    cat(paste0(length(object@pScreen)," screening hypothesis p-values \n"))
    cat(paste0(ncol(object@pConfirmation)," confirmation hypotheses for",nrow(object@pConfirmation)," genes \n"))
    cat(paste0("adjusted p-values on a ",object@alpha*100,"% OFDR level with the following FWER correction method: ",object@method," \n"))

  }
}

setMethod("show","stageRTx",function(object) print.stageRTx(object))





