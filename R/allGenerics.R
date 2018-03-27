#' @include stageRClasses.R

#' @rdname stageR
#' @export
setGeneric("stageR", function(pScreen, pConfirmation, ...) standardGeneric("stageR"))

#' @rdname getPScreen
#' @export
setGeneric("getPScreen", function(object, ...) standardGeneric("getPScreen"))

#' @rdname getPConfirmation
#' @export
setGeneric("getPConfirmation", function(object, ...) standardGeneric("getPConfirmation"))

#' @rdname stageWiseAdjustment
#' @export
setGeneric("stageWiseAdjustment", function(object, method, alpha, ...) standardGeneric("stageWiseAdjustment"))

#' @rdname getAdjustedPValues
#' @export
setGeneric("getAdjustedPValues", function(object, onlySignificantGenes, order, ...) standardGeneric("getAdjustedPValues"))

#' @rdname adjustedAlphaLevel
#' @export
setGeneric("adjustedAlphaLevel", function(object, ...) standardGeneric("adjustedAlphaLevel"))

#' @rdname getResults
#' @export
setGeneric("getResults", function(object, ...) standardGeneric("getResults"))

#' @rdname getSignificantGenes
#' @export
setGeneric("getSignificantGenes", function(object, ...) standardGeneric("getSignificantGenes"))

#' @rdname getSignificantTx
#' @export
setGeneric("getSignificantTx", function(object, ...) standardGeneric("getSignificantTx"))

#' @rdname getAlpha
#' @export
setGeneric("getAlpha", function(object, ...) standardGeneric("getAlpha"))

#' @rdname getTx2gene
#' @export
setGeneric("getTx2gene", function(object, ...) standardGeneric("getTx2gene"))

#' @rdname isPScreenAdjusted
#' @export
setGeneric("isPScreenAdjusted", function(object, ...) standardGeneric("isPScreenAdjusted"))

#' @rdname isAdjusted
#' @export
setGeneric("isAdjusted", function(object, ...) standardGeneric("isAdjusted"))

#' @rdname getMethod
#' @export
setGeneric("getMethod", function(object, ...) standardGeneric("getMethod"))
