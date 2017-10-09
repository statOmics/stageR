#' @include stageRClasses.R

setMethod("show", "stageR", function(object){
  cat("stageR object, containing: \n")
  if (!isAdjusted(object)) {
    cat(paste0("- ",length(getPScreen(object)), " screening hypothesis p-values \n"))
    cat(paste0("- ",
      ncol(getPConfirmation(object)),
      " confirmation hypotheses for ",
      nrow(getPConfirmation(object)),
      " genes \n"
    ))
  }
  if (isAdjusted(object)) {
    cat(paste0("- ",length(getPScreen(object)), " screening hypothesis p-values \n"))
    cat(paste0("- ",
      ncol(getPConfirmation(object)),
      " confirmation hypotheses for ",
      nrow(getPConfirmation(object)),
      " genes \n"
    ))
    cat(
      paste0(
        "- adjusted p-values on a ",
        getAlpha(object) * 100,
        "% OFDR level with the following FWER correction method: ",
        getMethod(object),
        " \n"
      )
    )

  }
})

setMethod("show", "stageRTx", function(object){
  cat("stageRTx object, containing: \n")
  if (!isAdjusted(object)) {
    cat(paste0("- ", length(getPScreen(object)), " screening hypothesis p-values \n"))
    cat(paste0("- ",
      ncol(getPConfirmation(object)),
      " confirmation hypothesis for ",
      nrow(getPConfirmation(object)),
      " transcripts \n"
    ))
  }
  if (isAdjusted(object)) {
    cat("- ",paste0(length(getPScreen(object)), " screening hypothesis p-values \n"))
    cat(paste0("- ",
      ncol(getPConfirmation(object)),
      " confirmation hypothesis for ",
      nrow(getPConfirmation(object)),
      " transcripts \n"
    ))
    cat(
      paste0(
        "- adjusted p-values on a ",
        getAlpha(object) * 100,
        "% OFDR level with the following FWER correction method: ",
        getMethod(object),
        " \n"
      )
    )

  }
})
