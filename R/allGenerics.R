#' @include stageRClasses.R

setGeneric("stageR", function(pScreen, pConfirmation, ...) standardGeneric("stageR"))

setGeneric("getPScreen", function(object, ...) standardGeneric("getPScreen"))

setGeneric("getPConfirmation", function(object, ...) standardGeneric("getPConfirmation"))

setGeneric("stageWiseAdjustment", function(object, method, alpha, ...) standardGeneric("stageWiseAdjustment"))

setGeneric("getAdjustedPValues", function(object, onlySignificantGenes, order, ...) standardGeneric("getAdjustedPValues"))

setGeneric("adjustedAlphaLevel", function(object, ...) standardGeneric("adjustedAlphaLevel"))

setGeneric("getResults", function(object, ...) standardGeneric("getResults"))

setGeneric("getSignificantGenes", function(object, ...) standardGeneric("getSignificantGenes"))

setGeneric("getSignificantTx", function(object, ...) standardGeneric("getSignificantTx"))

setGeneric("getAlpha", function(object, ...) standardGeneric("getAlpha"))

setGeneric("getTx2gene", function(object, ...) standardGeneric("getTx2gene"))

setGeneric("isPScreenAdjusted", function(object, ...) standardGeneric("isPScreenAdjusted"))

setGeneric("isAdjusted", function(object, ...) standardGeneric("isAdjusted"))

setGeneric("getMethod", function(object, ...) standardGeneric("getMethod"))
