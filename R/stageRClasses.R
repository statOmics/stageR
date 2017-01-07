#' The stageR class
#'
#' This class is used for adjusting p-values with stage-wise testing for high-throughput studies.
#'
#' @slot pScreen A vector of p-values for the screening hypothesis.
#' @slot pConfirmation A matrix of p-values for the confirmation hypotheses.
#' @slot adjustedP A matrix of adjusted p-values. This slot should be accessed through \code{\link{getAdjustedP}}. Alternatively, significance results can be accessed through \code{\link{getResults}}.
#' @slot method Character string indicating the method used for FWER correction in the confirmation stage of the stage-wise analysis. Can be any of \code{"none"}, \code{"holm"}, \code{"dte"}, \code{"dtu"}, \code{"user"}. \code{"none"} will not adjust the p-values in the confirmation stage. \code{"holm"} is an adapted Holm procedure for a stage-wise analysis, where the method takes into account the fact that genes in the confirmation stage have already passed the screening stage, hence the procedure will be more powerful for the most significant p-value as compared to the standard Holm procedure. \code{"dte"} is the adjusted Holm-Shaffer procedure for differential transcript expression analysis. \code{"dtu"} is the adjusted Holm-Shaffer procedure for differential transcript usage. \code{"user"} indicates a user-defined adjustment that should be specified with the \code{adjustment} argument.
#' @slot alpha the OFDR level on which the stage-wise analysis should be controlled.
#' @slot alphaAdjusted the adjusted significance level to compare against FWER-adjusted p-values of the confirmation stage to decide on significance of the hypothesis test.
#' @slot pScreenAdjusted logical, indicating whether the supplied p-values for the screening hypothesis have already been adjusted for multiplicity according to the FDR.
#' @slot tx2gene matrix with transcript IDs in the first column and gene IDs in the second column to be used for DTE and DTU analysis. All rownames from \code{pConfirmation} should match with a transcript ID and all names from \code{pScreen} should match with a gene ID.
#' @references
#' Van den Berge K., Soneson C., Robinson M.D., Clement L. "A general and powerful stage-wise testing procedure for differential expression and differential transcript usage." Submitted, 2017.
#' R. Heller, E. Manduchi, G. R. Grant, and W. J. Ewens, “A flexible two-stage procedure for identifying gene sets that are differentially expressed.” Bioinformatics (Oxford, England), vol. 25, pp. 1019–25, 2009.
#' S. Holm, “A Simple Sequentially Rejective Multiple Test Procedure,” Scandinavian Journal of Statistics, vol. 6, no. 2, pp. 65–70, 1979.
#' J. P. Shaffer, “Modified Sequentially Rejective Multiple Test Procedures,” Journal of the American Statistical Asso- ciation, vol. 81, p. 826, 1986.
#' @name stageRClass
#' @rdname stageRClass
#' @exportClass stageR
setClass("stageR",slots=c(pScreen="numeric", pConfirmation="matrix", adjustedP="matrix", method="character", alpha="numeric", alphaAdjusted="numeric", adjusted="logical", pScreenAdjusted="logical"))

#' @name stageRClass
#' @rdname stageRClass
#' @exportClass stageRTx
setClass("stageRTx",slots=c(pScreen="numeric", pConfirmation="matrix", adjustedP="matrix", method="character", alpha="numeric", alphaAdjusted="numeric", adjusted="logical", pScreenAdjusted="logical", tx2gene="data.frame"))


