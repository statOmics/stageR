#' @include stageRClasses.R allGenerics.R constructors.R
#' @import tidyverse dplyr

.createGeneTibble <- function(pScreen=NULL, pConfirmation, tx2gene=NULL, weights=NULL){
  require(tidyverse) ; require(reshape2)

  if(is.null(tx2gene)){
    df <- data.frame(geneID=rownames(pConfirmation), pConfirmation)
    # convert to long format
    df <- reshape2::melt(df, id.vars="geneID", variable.name="hypothesis", value.name="pvalue")
    nestedDf <- df %>% group_by(geneID) %>% nest(.key="data")
  } else { #transcript data
    df <- data.frame(txID=tx2gene[,1], geneID=tx2gene[,2], pvalue=pConfirmation[as.character(tx2gene[,1]),])
    # note that weights are only applicable on tx-level data.
    if(!is.null(weights)) df$weights <- weights
    nestedDf <- df %>% group_by(geneID) %>% nest(.key="data")
  }

  if(!is.null(pScreen)){
    nestedDf$genePval <- pScreen[as.character(nestedDf$geneID)]
    nestedDf <- nestedDf[,c("geneID","genePval","data")]
  }
  return(nestedDf)
}


## no confirmation stage adjustment function
.noAdjustment <- function(df){
  df$padj <- df$pvalue
  return(df)
}

## holm adjustment in confirmation stage for complex DGE experiments
.holmAdjustment <- function(df){
  pval <- df$pvalue
  o <- order(pval)
  if(all(!is.na(pval))){ #if no NA's, standard Holm with screening correct.
    n <- length(pval)
  } else { #if NA's present, only correct for non NA p-values
    n <- length(pval[!is.na(pval)])
  }
  # Holm adjustment: passing screening stage implies 1 false hypothesis
  adjustment <- c(n-1,(n-1):1)
  if(length(adjustment)!=length(pval)){ #if NA values are present
    adjustment <- c(adjustment, rep(1,length(pval)-length(adjustment)))
  }
  rowAdjusted <- pval[o]*adjustment
  rowAdjusted <- pmin(rowAdjusted,1)
  rowAdjusted <- cummax(rowAdjusted)
  rowBack <- vector(length=length(pval))
  rowBack[o] <- rowAdjusted
  df$padj <- rowBack
  return(df)
}

## user specified adjustment in confirmation stage for complex DGE experiments
.userAdjustment <- function(df, adjustment){
  pval <- df$pvalue
  o <- order(pval)
  rowAdjusted <- pval[o]*adjustment
  rowAdjusted <- pmin(rowAdjusted,1)
  # check monotone increase of adjusted p-values
  rowAdjusted <- cummax(rowAdjusted)
  rowBack <- vector(length=length(pval))
  rowBack[o] <- rowAdjusted
  df$padj <- rowBack
  return(df)
}

## DTE adjustment in confirmation stage
.dteAdjustment <- function(df){
  pval <- df$pvalue
  o <- order(pval)
  n <- length(pval)
  # DTE adjustment: passing screening stage implies 1 false hypothesis
  if(n==1) adjustment=0 else adjustment=c(n-1,(n-1):1)
  rowAdjusted <- pval[o]*adjustment
  rowAdjusted <- pmin(rowAdjusted,1)
  rowAdjusted <- cummax(rowAdjusted)
  rowBack <- vector(length=length(pval))
  rowBack[o] <- rowAdjusted
  df$padj <- rowBack
  return(df)
}

## DTU adjustment in confirmation stage
.dtuAdjustment <- function(df){
  pval <- df$pvalue
  o <- order(pval)
  n <- length(pval)
  # DTU adjustment: passing screening stage implies 2 false hypotheses
  if(n==1) stop("Some genes have only one transcript; this is incompatible with DTU correction. Remove these transcripts.")
  if(n==2) adjustment=c(0,0) else adjustment=c(n-2,n-2,(n-2):1)
  rowAdjusted <- pval[o]*adjustment
  rowAdjusted <- pmin(rowAdjusted,1)
  rowAdjusted <- cummax(rowAdjusted)
  rowBack <- vector(length=length(pval))
  rowBack[o] <- rowAdjusted
  df$padj <- rowBack
  return(df)
}

# aggregating transcript p-values
.aggregatePvalues <- function(geneTibble, aggMethod){
  require(aggregation)
  if(aggMethod %in% c("sidak","fisher")){
    aggMethod <- get(aggMethod)
    pScreen <- unlist(map(geneTibble$data, function(x) aggMethod(x$pvalue)))
  } else if(aggMethod=="lancaster"){
    aggMethod <- get(aggMethod)
    pScreen <- unlist(map(geneTibble$data, function(x) aggMethod(x$pvalue, weights=x$weights)))
  }
  geneTibble$genePval <- pScreen
  return(geneTibble)
}


.stageWiseTest <- function(pScreen, pConfirmation, alpha,
                           method=c("none","holm","dte","dtu","user"),
                           adjustment=NULL, tx2gene=NULL, pScreenAdjusted,
                           allowNA=FALSE, aggMethod=NULL, weights=NULL){

  if(is.null(aggMethod) & is.null(pScreen)) stop("Neither pScreen nor aggMethod were specified. Please specify either screening stage p-values or an aggregation method for the confirmation stage p-values.")

  ## remove NA screening hypothesis p-values.
  if(allowNA){
    if(any(is.na(pScreen))){
      naFeatures <- which(is.na(pScreen))
      message(paste0("Removing ",length(naFeatures),
                     " features with NA screening hypothesis p-values. \n"))
      pScreen <- pScreen[-naFeatures]
      pConfirmation <- pConfirmation[-naFeatures,]
    }
  }

  ## check for NA values.
  if(!allowNA){
    if(any(is.na(pScreen)) | any(is.na(pConfirmation))){
      stop("NA p-values found in either the screening or confirmation tests.
           If you want to allow for NA p-values, set allowNA=TRUE.")
    }
  }
  method <- match.arg(method,c("none","holm","dte","dtu","user"))

  # create nested data frame
  geneTibble <- .createGeneTibble(pScreen=pScreen, pConfirmation=pConfirmation, tx2gene=tx2gene, weights=weights)

  # aggregate p-values to obtain screening stage p-value.
  aggMethod <- match.arg(aggMethod,c("sidak","fisher","lancaster"))
  if(is.null(pScreen)){
    geneTibble <- .aggregatePvalues(geneTibble=geneTibble, aggMethod=aggMethod)
  }

  # screening stage
  if(!pScreenAdjusted){
    pScreen <- geneTibble$genePval
    names(pScreen) <- geneTibble$geneID
    padjScreen <- p.adjust(pScreen,"BH")
  } else {
    padjScreen <- pScreen
  }
  geneTibble$genePadj <- padjScreen
  genesStageI <- padjScreen<=alpha

  # only select genes passing screening stage.
  geneTibbleStageI <- geneTibble[geneTibble$geneID%in%names(which(genesStageI)),]

  # confirmation stage adjustment
  if(method=="none"){

    geneTibbleStageII <- geneTibbleStageI
    geneTibbleStageII$data <- map(geneTibbleStageI$data, .noAdjustment)

  } else if(method=="holm"){

    geneTibbleStageII <- geneTibbleStageI
    geneTibbleStageII$data <- map(geneTibbleStageI$data, .holmAdjustment)

  } else if(method=="user"){
    if(length(adjustment)!=ncol(pConfirmation)){
      stop("the length of the adjustment vector is not equal to the number of confirmation hypotheses as defined by the number of columns in pConfirmation.")
    }

    geneTibbleStageII <- geneTibbleStageI
    geneTibbleStageII$data <- map(geneTibbleStageI$data, .userAdjustment,  adjustment=adjustment)

  } else if(method=="dte"){
    if(any(is.na(match(rownames(pConfirmation),tx2gene[,1])))){
      stop("not all transcript names in pConfirmation match with
           a transcript ID from the tx2gene object.")
    }
    if(any(is.na(match(names(pScreen),tx2gene[,2])))){
      stop("not all gene names in pScreen match with
           a gene ID from the tx2gene object.")
    }

    geneTibbleStageII <- geneTibbleStageI
    geneTibbleStageII$data <- map(geneTibbleStageI$data, .dteAdjustment)

  } else if(method=="dtu"){
    if(any(is.na(match(rownames(pConfirmation),tx2gene[,1])))){
      stop("not all transcript names in pConfirmation match with
           a transcript ID from the tx2gene object.")
    }
    if(any(is.na(match(names(pScreen),tx2gene[,2])))){
      stop("not all gene names in pScreen match with
           a gene ID from the tx2gene object.")
    }

    geneTibbleStageII <- geneTibbleStageI
    geneTibbleStageII$data <- map(geneTibbleStageI$data, .dtuAdjustment)

  } else stop("specify a valid method for the confirmation stage.")

  # calculate BH-adjusted s.l.
  G <- length(padjScreen) #nr of genes
  R <- sum(padjScreen<=alpha) #nr of rejections
  alphaAdjusted <- R/G*alpha

  # correct FWER-adjusted p-values acc. to BH-adjusted s.l.
  geneTibbleStageII$data <- map(geneTibbleStageII$data, function(x){
    naId <- is.na(x$padj)
    x <- x %>% mutate(padj_SW=padj)
    x$padj_SW[!naId] <- x$padj[!naId]*G/R
    x$padj_SW[!naId] <- pmin(x$padj_SW[!naId],1)
    return(x)
  })

  return(list(geneTibble=geneTibble, geneTibbleStageII=geneTibbleStageII, alphaAdjusted=alphaAdjusted))
}



#' adjust p-values in a two-stage analysis
#'
#' This function will adjust p-values according to a hierarchical two-stage testing paradigm.
#'
#' @param object an object of the \code{\link{stageRClass}} class.
#' @param method Character string indicating the method used for FWER correction in the confirmation stage of the stage-wise analysis. Can be any of \code{"none"}, \code{"holm"}, \code{"dte"}, \code{"dtu"}, \code{"user"}. \code{"none"} will not adjust the p-values in the confirmation stage. \code{"holm"} is an adapted Holm procedure for a stage-wise analysis, where the method takes into account the fact that genes in the confirmation stage have already passed the screening stage, hence the procedure will be more powerful for the most significant p-value as compared to the standard Holm procedure. \code{"dte"} is the adjusted Holm-Shaffer procedure for differential transcript expression analysis. \code{"dtu"} is the adjusted Holm-Shaffer procedure for differential transcript usage. \code{"user"} indicates a user-defined adjustment that should be specified with the \code{adjustment} argument.
#' @param alpha the OFDR on which to control the two-stage analysis.
#' @param tx2gene Only applicable when  \code{method} is \code{"dte"} or \code{"dtu"}.  A \code{\link[base]{data.frame}} with transcript IDs in the first column and gene IDs in the second column. The rownames from \code{pConfirmation} must be contained in the transcript IDs from \code{tx2gene}, and the names from \code{pScreen} must be contained in the gene IDs.
#' @param aggMethod The method to use to aggregate p-values across hypotheses. Only applicable when \code{pScreen} is not provided. Options are \code{"sidak"}, \code{"fisher"}, and \code{"lancaster"}, as provided by the \code{aggregation} R package from Yi et al. (2018).
#' @param weights The weights to use in p-value aggregation when \code{aggMethod} is \code{"lancaster"}. Note that this is only applicable for transcript-level data. Typically, these are set to the base mean expression of the transcript, as suggested in Yi et al. (2018).
#' @param adjustment a user-defined adjustment of the confirmation stage p-values. Only applicable when \code{method} is \code{"user"} and ignored otherwise.
#' @param ... Additional arguments passed to \code{.stageWiseTest}
#' @return
#' A stageR/stageRTx object with stage-wise adjusted p-values.
#' @references
#'
#' Van den Berge K., Soneson C., Robinson M.D., and Clement L. (2017). stageR: a general stage-wise method for controlling the gene-level false discovery rate in differential expression and differential transcript usage. Genome Biology 18:151. https://doi.org/10.1186/s13059-017-1277-0
#'
#' R. Heller, E. Manduchi, G. R. Grant, and W. J. Ewens, "A flexible two-stage procedure for identifying gene sets that are differentially expressed." Bioinformatics (Oxford, England), vol. 25, pp. 1019-25, 2009.
#'
#' S. Holm, "A Simple Sequentially Rejective Multiple Test Procedure," Scandinavian Journal of Statistics, vol. 6, no. 2, pp. 65-70, 1979.
#'
#' J. P. Shaffer, "Modified Sequentially Rejective Multiple Test Procedures," Journal of the American Statistical Association, vol. 81, p. 826, 1986.
#'
#' L. Yi, H. Pimentel, N.L. Bray, and L. Pachter, "Gene-level differential analysis at transcript-level resolution." Genome Biology 19:53, 2018.
#'
#' @examples
#'    # DGE
#'    pScreen=c(seq(1e-10,1e-2,length.out=100),seq(1e-2,.2,length.out=100),seq(.2,1,length.out=100))
#'    names(pScreen)=paste0("gene",1:300)
#'    set.seed(7)
#'    pConfirmation=matrix(runif(900),nrow=300,ncol=3)
#'    dimnames(pConfirmation)=list(paste0("gene",1:300),c("H1","H2","H3"))
#'    stageRObj <- stageR(pScreen=pScreen, pConfirmation=pConfirmation)
#'    stageRObj <- stageWiseAdjustment(stageRObj, method="holm", alpha=0.05)
#'
#'
#'    # tx with pScreen
#'    set.seed(1)
#'    genes=paste0("gene",sample(1:200,1000,replace=TRUE))
#'    nGenes=length(table(genes))
#'    transcripts=paste0("tx",1:1000)
#'    tx2gene=data.frame(transcripts,genes)
#'    #gene-wise p-values
#'    pScreen=c(seq(1e-10,1e-2,length.out=nGenes-100),seq(1e-2,.2,length.out=50),seq(.2,1,length.out=50))
#'    names(pScreen)=names(table(genes)) #discards genes that are not simulated
#'    pConfirmation=matrix(runif(1000),nrow=1000,ncol=1)
#'    rownames(pConfirmation)=transcripts
#'    stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation, tx2gene=tx2gene)
#'    stageRObj <- stageWiseAdjustment(stageRObj, method="dte", alpha=0.05)
#'
#'    # tx with aggMethod lancaster
#'    set.seed(1)
#'    genes=paste0("gene",sample(1:200,1000,replace=TRUE))
#'    nGenes=length(table(genes))
#'    transcripts=paste0("tx",1:1000)
#'    tx2gene=data.frame(transcripts,genes)
#'    pScreen=NULL
#'    pConfirmation=matrix(c(runif(100,1e-10,1e-4),runif(900)),nrow=1000,ncol=1)
#'    rownames(pConfirmation)=transcripts
#'    aggMethod="lancaster"
#'    weights=seq(0,1,length=nrow(pConfirmation))
#'    stageRObj <- stageRTx(pConfirmation=pConfirmation ,pScreenAdjusted=pScreenAdjusted, tx2gene=tx2gene)
#'    stageRObj <- stageWiseAdjustment(stageRObj, method="dte", alpha=0.05, aggMethod=aggMethod, weights=weights)
#'
#' @name stageWiseAdjustment
#' @rdname stageWiseAdjustment
#' @aliases stageWiseAdjustment stageWiseAdjustment,stageR stageWiseAdjustment,stageRTx
#' @export
setMethod("stageWiseAdjustment",
          signature=signature(object="stageR",
                              method="character",
                              alpha="numeric"),
          definition=function(object, method, alpha, adjustment=NULL, weights=NULL, ...){
            if(!is.null(weights)) stop("Weights can only provided for transcript-level data, hence with an object of class stageRTx.")
            pScreen <- getPScreen(object)
            if(length(pScreen)==0) pScreen <- NULL
            pConfirmation <- getPConfirmation(object)
            pScreenAdjusted <- isPScreenAdjusted(object)
            stageAdjPValues <- .stageWiseTest(pScreen=pScreen,
                                              pConfirmation=pConfirmation,
                                              alpha=alpha,
                                              method=method,
                                              pScreenAdjusted=pScreenAdjusted,
                                              adjustment=adjustment, ...)
            object@geneTibble <- stageAdjPValues[["geneTibble"]]
            object@adjustedP <- stageAdjPValues[["geneTibbleStageII"]]
            object@alphaAdjusted <- stageAdjPValues[["alphaAdjusted"]]
            if(length(getPScreen(object))==0) object@pScreen <- stageAdjPValues[["geneTibble"]]$genePval
            object@method <- method
            object@alpha <- alpha
            object@adjusted <- TRUE
            return(object)
          })
#' @describeIn stageWiseAdjustment Adjust p-values in a two-stage analysis
setMethod(
  "stageWiseAdjustment",
  signature = signature(
    object = "stageRTx",
    method = "character",
    alpha = "numeric"
  ),
  definition = function(object, method, alpha, tx2gene, ...) {
    pScreen <- getPScreen(object)
    if(length(pScreen)==0) pScreen <- NULL
    pConfirmation <- getPConfirmation(object)
    pScreenAdjusted <- isPScreenAdjusted(object)
    tx2gene <- getTx2gene(object)
    stageAdjPValues <-.stageWiseTest(
        pScreen = pScreen,
        pConfirmation = pConfirmation,
        alpha = alpha,
        method = method,
        pScreenAdjusted = pScreenAdjusted,
        tx2gene = tx2gene,
        ...
      )
    object@geneTibble <- stageAdjPValues[["geneTibble"]]
    object@adjustedP <- stageAdjPValues[["geneTibbleStageII"]]
    object@alphaAdjusted <- stageAdjPValues[["alphaAdjusted"]]
    if(length(getPScreen(object))==0) object@pScreen <- stageAdjPValues[["geneTibble"]]$genePval
    object@method <- method
    object@alpha <- alpha
    object@adjusted <- TRUE
    return(object)
  }
)

#' Return screening hypothesis p-values from a \code{\link{stageRClass}} object.
#'
#' @param object an object of the \code{\link{stageRClass}} class.
#' @param ... Additional arguments
#' @return
#' A vector of screening stage (aggregated) p-values.
#' @examples
#' pScreen=c(seq(1e-10,1e-2,length.out=100),seq(1e-2,.2,length.out=100),seq(.2,1,length.out=100))
#' names(pScreen)=paste0("gene",1:300)
#' pConfirmation=matrix(runif(900),nrow=300,ncol=3)
#' dimnames(pConfirmation)=list(paste0("gene",1:300),c("H1","H2","H3"))
#' stageRObj <- stageR(pScreen=pScreen, pConfirmation=pConfirmation)
#' getPScreen(stageRObj)
#' @name getPScreen
#' @rdname getPScreen
#' @aliases getPScreen getPScreen,stageR getPScreen,stageRTx
#' @export
setMethod("getPScreen",signature=signature(object="stageR"),
          definition=function(object){return(object@pScreen)})
#' @describeIn getPScreen Return screening hypothesis p-values from a \code{\link{stageRClass}} object.
setMethod("getPScreen",signature=signature(object="stageRTx"),
          definition=function(object){return(object@pScreen)})

#' Return unadjusted confirmation hypothesis p-values from a \code{\link{stageRClass}} object.
#'
#' @param object an object of the \code{\link{stageRClass}} class.
#' @return
#' A matrix of the unadjusted p-values to be used in the confirmation stage.
#' @examples
#' pScreen=c(seq(1e-10,1e-2,length.out=100),seq(1e-2,.2,length.out=100),seq(.2,1,length.out=100))
#' names(pScreen)=paste0("gene",1:300)
#' pConfirmation=matrix(runif(900),nrow=300,ncol=3)
#' dimnames(pConfirmation)=list(paste0("gene",1:300),c("H1","H2","H3"))
#' stageRObj <- stageR(pScreen=pScreen, pConfirmation=pConfirmation)
#' getPConfirmation(stageRObj)
#' @aliases getPConfirmation getPConfirmation,stageR getPConfirmation,stageRTx
#' @name getPConfirmation
#' @rdname getPConfirmation
#' @export
setMethod("getPConfirmation",signature=signature(object="stageR"),
          definition=function(object){return(object@pConfirmation)})
#' @describeIn getPConfirmation Return unadjusted confirmation hypothesis p-values from a \code{\link{stageRClass}} object.
setMethod("getPConfirmation",signature=signature(object="stageRTx"),
          definition=function(object){return(object@pConfirmation)})





.getAdjustedP <- function(object, onlySignificantGenes=FALSE, order=TRUE){
  ## this function is used in getAdjustedPValues
  ## to return the adjusted p-values for a stageR class.
  message(paste0("The returned adjusted p-values are based on a ",
                 "stage-wise testing approach and are only valid for ",
                 "the provided target OFDR level of ",
                 getAlpha(object)*100,
                 "%. If a different target OFDR level is of interest,",
                 "the entire adjustment should be re-run. \n"))
  if(onlySignificantGenes){ #significant genes
      if(order){
        padj <- object@adjustedP
        o <- order(padj$genePval)
        return(padj[o,])
      } else {
        return(object@adjustedP)
      }
  } else { #all genes
    tibble <- object@geneTibble
    padj <- object@adjustedP
    tibble$data[match(padj$geneID,tibble$geneID)] <- padj$data
    if(order){
      o <- order(tibble$genePval)
      return(tibble[o,])
    } else {
      return(tibble)
    }
  }
}

#' Retrieve the stage-wise adjusted p-values.
#'
#' This functions returns the stage-wise adjusted p-values for an object from the  \code{\link{stageRClass}} class. Note, that the p-values should have been adjusted with the \code{\link{stageWiseAdjustment,stageR,character,numeric-method}} function prior to calling this function.
#'
#' @param object an object of the \code{\link{stageRClass}} class.
#' @param onlySignificantGenes logical. If FALSE (default), all genes are returned. If TRUE, only the genes significant for the screening hypothesis are returned.
#' @param order logical. If TRUE (default), the returned matrix of adjusted p-values are ordered based on the screening hypothesis adjusted p-value.
#' @param ... Other arguments passed to .getAdjustedP or .getAdjustedPTx
#' @return
#' For complex DGE experiments (stageR object), a matrix of adjusted p-values where every row corresponds to a gene, and every column corresponds to a contrast. The first column will be the BH FDR adjusted p-value from the screening step.
#' For transcript-level experiments (stageRTx object), a matrix of adjusted p-values where every row corresponds to a transcript.
#' @details
#' The function returns FDR adjusted p-values for the screening hypothesis and stage-wise adjusted p-values for the confirmation hypothesis p-values. For features that were not significant in the screening hypothesis, the confirmation stage adjusted p-values are set to \code{NA}.
#' The adjusted p-values in the output of \code{getAdjustedPValues} can directly be compared to alpha, the OFDR level specified in \code{stageWiseAdjustment}, to flag significant features.
#' @examples
#' pScreen=c(seq(1e-10,1e-2,length.out=100),seq(1e-2,.2,length.out=100),seq(.2,1,length.out=100))
#' names(pScreen)=paste0("gene",1:300)
#' pConfirmation=matrix(runif(900),nrow=300,ncol=3)
#' dimnames(pConfirmation)=list(paste0("gene",1:300),c("H1","H2","H3"))
#' stageRObj <- stageR(pScreen=pScreen, pConfirmation=pConfirmation)
#' stageRObj <- stageWiseAdjustment(stageRObj, method="holm", alpha=0.05)
#' head(getAdjustedPValues(stageRObj, onlySignificantGenes=TRUE, order=TRUE))
#' @references
#' Van den Berge K., Soneson C., Robinson M.D., Clement L. (2017). stageR: a general stage-wise method for controlling the gene-level false discovery rate in differential expression and differential transcript usage. Genome Biology 18:151. https://doi.org/10.1186/s13059-017-1277-0
#'
#' @name getAdjustedPValues
#' @rdname getAdjustedPValues
#' @aliases getAdjustedPValues getAdjustedPValues,stageR getAdjustedPValues,stageRTx
#' @export
setMethod("getAdjustedPValues",
          signature=signature(object="stageR",
                              onlySignificantGenes="logical",
                              order="logical"),
          definition=function(object, onlySignificantGenes, order, ...){
            if(!isAdjusted(object)){
              stop("adjust p-values first using stageWiseAdjustment")}
            return(.getAdjustedP(object=object, onlySignificantGenes, order, ...))
          })
#' @describeIn getAdjustedPValues Retrieve the stage-wise adjusted p-values.
setMethod("getAdjustedPValues",
          signature=signature(object="stageRTx",
                              onlySignificantGenes="logical",
                              order="logical"),
          definition=function(object, onlySignificantGenes, order, ...){
            if(!isAdjusted(object)){
              stop("adjust p-values first using stageWiseAdjustment")}
            return(.getAdjustedP(object=object, onlySignificantGenes, order, ...))
          })

#' Get adjusted significance level from the screening stage.
#'
#' This functions returns the adjusted significance level from the screening stage that should be used to compare confirmation stage FWER adjusted p-values against.
#'
#' @param object an object of the \code{\link{stageRClass}} class.
#' @details
#' The adjusted significance level is calculated as the fraction of significant features in the screening stage multiplied the alpha level.
#' @return
#' Scalar, the adjusted significance level from the screening stage.
#' @examples
#' pScreen=c(seq(1e-10,1e-2,length.out=100),seq(1e-2,.2,length.out=100),seq(.2,1,length.out=100))
#' names(pScreen)=paste0("gene",1:300)
#' pConfirmation=matrix(runif(900),nrow=300,ncol=3)
#' dimnames(pConfirmation)=list(paste0("gene",1:300),c("H1","H2","H3"))
#' stageRObj <- stageR(pScreen=pScreen, pConfirmation=pConfirmation, pScreenAdjusted=FALSE)
#' stageRObj <- stageWiseAdjustment(stageRObj, method="holm", alpha=0.05)
#' adjustedAlphaLevel(stageRObj)
#'# @method stageR-method
#' @references
#' Van den Berge K., Soneson C., Robinson M.D., Clement L. (2017). stageR: a general stage-wise method for controlling the gene-level false discovery rate in differential expression and differential transcript usage. Genome Biology 18:151. https://doi.org/10.1186/s13059-017-1277-0
#'
#' R. Heller, E. Manduchi, G. R. Grant, and W. J. Ewens, "A flexible two-stage procedure for identifying gene sets that are differentially expressed." Bioinformatics (Oxford, England), vol. 25, pp. 1019-25, 2009.
#'
#' @seealso \code{\link{stageR}}, \code{\link{stageRClass}}
#' @name adjustedAlphaLevel
#' @rdname adjustedAlphaLevel
#' @aliases adjustedAlphaLevel adjustedAlphaLevel,stageR adjustedAlphaLevel,stageRTx
#' @export
setMethod("adjustedAlphaLevel",signature=signature(object="stageR"),
          definition=function(object){
            if(!isAdjusted(object)){
              stop("adjust p-values first using stageWiseAdjustment")}
            return(object@alphaAdjusted)
          })
#' @describeIn adjustedAlphaLevel Get adjusted significance level from the screening stage.
setMethod("adjustedAlphaLevel",signature=signature(object="stageRTx"),
          definition=function(object){
            if(!isAdjusted(object)){
              stop("adjust p-values first using stageWiseAdjustment")}
            return(object@alphaAdjusted)
          })



.getResults <- function(object){
  alpha <- getAlpha(object)
  res <- do.call(rbind, map(object@adjustedP$data,function(x) x$padj_SW<=alpha))+0
  res <- cbind(object@adjustedP$genePadj<=alpha,res)
  dimnames(res) <- list(object@adjustedP$geneID, c("screen",colnames(getPConfirmation(object))))
  return(res)
}

#' Get significance results according to a stage-wise analysis.
#'
#' This functions returns a matrix that indicates whether a specific feature is significant for a specific hypothesis of interest according to a stage-wise analysis. The function is not applicable to transcript-level analysis.
#'
#' @param object an object of the \code{\link{stageRClass}} class.
#' @details
#' The FDR adjusted screening hypothesis p-values are compared to the alpha level specified. The FWER adjusted confirmation stage p-values are compared to the adjusted significance level from the screening stage.
#' @return
#' A logical matrix with rows corresponding to genes and columns corresponding to contrasts, where the first column represents the screening stage on the aggregated p-values. A 0 represents a non-significant test, a 1 represents a significant test according to the stage-wise analysis.
#' @examples
#' pScreen=c(seq(1e-10,1e-2,length.out=100),seq(1e-2,.2,length.out=100),seq(.2,1,length.out=100))
#' names(pScreen)=paste0("gene",1:300)
#' pConfirmation=matrix(runif(900),nrow=300,ncol=3)
#' dimnames(pConfirmation)=list(paste0("gene",1:300),c("H1","H2","H3"))
#' stageRObj <- stageR(pScreen=pScreen, pConfirmation=pConfirmation)
#' stageRObj <- stageWiseAdjustment(stageRObj, method="holm", alpha=0.05)
#' head(getResults(stageRObj))
#' @references
#' Van den Berge K., Soneson C., Robinson M.D., Clement L. (2017). stageR: a general stage-wise method for controlling the gene-level false discovery rate in differential expression and differential transcript usage. Genome Biology 18:151. https://doi.org/10.1186/s13059-017-1277-0
#'
#' @name getResults
#' @rdname getResults
#' @aliases getResults getResults,stageR getResults,stageRTx
#' @export
setMethod("getResults",signature=signature(object="stageR"),
          definition=function(object){
            if(!isAdjusted(object)){
              stop("adjust p-values first using stageWiseAdjustment")}
            return(.getResults(object))
          })

#' Return significant genes when performing transcript-level analysis.
#'
#' This functions returns a matrix with significant genes by aggregated testing of its respective transcripts.
#'
#' @param object an object of the \code{\link{stageRClass}} class.
#' @return
#' A matrix with significant genes and their corresponding FDR-adjusted screening stage (aggregated) p-value.
#' @examples
#' #make identifiers linking transcripts to genes
#' set.seed(1)
#' genes=paste0("gene",sample(1:200,1000,replace=TRUE))
#' nGenes=length(table(genes))
#' transcripts=paste0("tx",1:1000)
#' tx2gene=data.frame(transcripts,genes)
#' #gene-wise q-values
#' pScreen=c(seq(1e-10,1e-2,length.out=nGenes-100),seq(1e-2,.2,length.out=50),seq(.2,1,length.out=50))
#' names(pScreen)=names(table(genes)) #discards genes that are not simulated
#' pConfirmation=matrix(runif(1000),nrow=1000,ncol=1)
#' rownames(pConfirmation)=transcripts
#' stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation ,pScreenAdjusted=TRUE, tx2gene=tx2gene)
#' stageRObj <- stageWiseAdjustment(stageRObj, method="dte", alpha=0.05)
#' head(getSignificantGenes(stageRObj))
#' @references
#' Van den Berge K., Soneson C., Robinson M.D., Clement L. (2017). stageR: a general stage-wise method for controlling the gene-level false discovery rate in differential expression and differential transcript usage. Genome Biology 18:151. https://doi.org/10.1186/s13059-017-1277-0
#'
#' @name getSignificantGenes
#' @rdname getSignificantGenes
#' @aliases getSignificantGenes getSignificantGenes,stageRTx
#' @export
setMethod("getSignificantGenes",signature=signature(object="stageRTx"),
          definition=function(object){
            if(!isAdjusted(object)){
              stop("adjust p-values first using stageWiseAdjustment")}
            if(class(object)!="stageRTx"){
              stop("this function only works on an object of class stageRTx")}
            padj <- getAdjustedPValues(object, onlySignificantGenes=TRUE, order=FALSE)
            padjGene <- padj[padj$genePadj<=getAlpha(object),]
            return(padjGene)
          })


#' Return significant transcripts when performing transcript-level analysis.
#'
#' This functions returns a matrix with significant transctripts according to a stage-wise analysis.
#'
#' @param object an object of the \code{\link{stageRClass}} class.
#' @return
#' A matrix of significant transcripts with their corresponding stage-wise adjusted p-value (i.e. FDR and FWER correction).
#' @references
#' Van den Berge K., Soneson C., Robinson M.D., Clement L. (2017). stageR: a general stage-wise method for controlling the gene-level false discovery rate in differential expression and differential transcript usage. Genome Biology 18:151. https://doi.org/10.1186/s13059-017-1277-0
#' @examples
#' #make identifiers linking transcripts to genes
#' set.seed(1)
#' genes=paste0("gene",sample(1:200,1000,replace=TRUE))
#' nGenes=length(table(genes))
#' transcripts=paste0("tx",1:1000)
#' tx2gene=data.frame(transcripts,genes)
#' #gene-wise q-values
#' pScreen=c(seq(1e-10,1e-2,length.out=nGenes-100),seq(1e-2,.2,length.out=50),seq(.2,1,length.out=50))
#' names(pScreen)=names(table(genes)) #discards genes that are not simulated
#' pConfirmation=matrix(runif(1000),nrow=1000,ncol=1)
#' rownames(pConfirmation)=transcripts
#' stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation ,pScreenAdjusted=TRUE, tx2gene=tx2gene)
#' stageRObj <- stageWiseAdjustment(stageRObj, method="dte", alpha=0.05)
#' head(getSignificantTx(stageRObj))
#' @name getSignificantTx
#' @rdname getSignificantTx
#' @aliases getSignificantTx getSignificantTx,stageR getSignificantTx,stageRTx
#' @export
setMethod("getSignificantTx",signature=signature(object="stageRTx"),
          definition=function(object){
            if(!isAdjusted(object)){
              stop("adjust p-values first using stageWiseAdjustment")}
            if(class(object)!="stageRTx"){
              stop("this function only works on an object of class stageRTx")}
            padj <- getAdjustedPValues(object, onlySignificantGenes=TRUE, order=FALSE)
            padjLong <- unnest(padj)
            return(padjLong[padjLong$padj_SW<=getAlpha(object),])
          })

#' Retrieve the significance level for the stage-wise adjustment.
#'
#' This functions returns the significance level on which the stage-wise adjustment is based.
#'
#' @param object an object of the \code{\link{stageRClass}} or \code{stageRTxClass} class.
#' @param ... Additional arguments
#' @return
#' Returns a calar vector with the OFDR alpha level that was specified by the user.
#' @examples
#' pScreen=c(seq(1e-10,1e-2,length.out=100),seq(1e-2,.2,length.out=100),seq(.2,1,length.out=100))
#' names(pScreen)=paste0("gene",1:300)
#' pConfirmation=matrix(runif(900),nrow=300,ncol=3)
#' dimnames(pConfirmation)=list(paste0("gene",1:300),c("H1","H2","H3"))
#' stageRObj <- stageR(pScreen=pScreen, pConfirmation=pConfirmation)
#' stageRObj <- stageWiseAdjustment(stageRObj, method="holm", alpha=0.05)
#' getAlpha(stageRObj)
#' @references
#' Van den Berge K., Soneson C., Robinson M.D., Clement L. (2017). stageR: a general stage-wise method for controlling the gene-level false discovery rate in differential expression and differential transcript usage. Genome Biology 18:151. https://doi.org/10.1186/s13059-017-1277-0
#'
#' @name getAlpha
#' @rdname getAlpha
#' @aliases getAlpha getAlpha,stageR getAlpha,stageRTx
#' @export
setMethod("getAlpha",signature=signature(object="stageR"),
          definition=function(object, ...){
            if(is.null(object@alpha)){
              stop("No significance level was specified yet.
                   Maybe you need to run stageWiseAdjustment first.")}
            return(object@alpha)
          })
#' @describeIn getAlpha Retrieve the significance level for the stage-wise adjustment.
setMethod("getAlpha",signature=signature(object="stageRTx"),
          definition=function(object, ...){
            if(is.null(object@alpha)){
              stop("No significance level was specified yet.
                   Maybe you need to run stageWiseAdjustment first.")}
            return(object@alpha)
          })

#' Retrieve the data frame linking genes to transcripts.
#'
#' This functions returns a data frame that links the genes with the transcripts being analysed.
#'
#' @param object an object of the \code{stageRTxClass} class.
#' @param ... Additional arguments
#' @return
#' A matrix linking gene to transcript identifiers.
#' @examples
#' #make identifiers linking transcripts to genes
#' set.seed(1)
#' genes=paste0("gene",sample(1:200,1000,replace=TRUE))
#' nGenes=length(table(genes))
#' transcripts=paste0("tx",1:1000)
#' tx2gene=data.frame(transcripts,genes)
#' #gene-wise q-values
#' pScreen=c(seq(1e-10,1e-2,length.out=nGenes-100),seq(1e-2,.2,length.out=50),seq(.2,1,length.out=50))
#' names(pScreen)=names(table(genes)) #discards genes that are not simulated
#' pConfirmation=matrix(runif(1000),nrow=1000,ncol=1)
#' rownames(pConfirmation)=transcripts
#' stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation ,pScreenAdjusted=TRUE, tx2gene=tx2gene)
#' getTx2gene(stageRObj)
#' @references
#' Van den Berge K., Soneson C., Robinson M.D., Clement L. (2017). stageR: a general stage-wise method for controlling the gene-level false discovery rate in differential expression and differential transcript usage. Genome Biology 18:151. https://doi.org/10.1186/s13059-017-1277-0
#'
#' @name getTx2gene
#' @rdname getTx2gene
#' @aliases getTx2gene getTx2gene,stageRTx
#' @export
setMethod("getTx2gene",signature=signature(object="stageRTx"),
          definition=function(object, ...){
            return(object@tx2gene)
          })

#' Are the screening p-values adjusted for multiplicity?
#'
#' This functions returns a logical stating whether the screening hypothesis p-values are already adjusted for multiple testing according to the BH FDR criterion.
#'
#' @param object an object of the \code{\link{stageRClass}} or \code{stageRTxClass} class.
#' @param ... Additional arguments
#' @return
#' A logical stating whether the screening hypothesis p-values are already adjusted for multiple testing according to the BH FDR criterion.
#' @examples
#' pScreen=c(seq(1e-10,1e-2,length.out=100),seq(1e-2,.2,length.out=100),seq(.2,1,length.out=100))
#' names(pScreen)=paste0("gene",1:300)
#' pConfirmation=matrix(runif(900),nrow=300,ncol=3)
#' dimnames(pConfirmation)=list(paste0("gene",1:300),c("H1","H2","H3"))
#' stageRObj <- stageR(pScreen=pScreen, pConfirmation=pConfirmation)
#' isPScreenAdjusted(stageRObj)
#' @references
#' Van den Berge K., Soneson C., Robinson M.D., Clement L. (2017). stageR: a general stage-wise method for controlling the gene-level false discovery rate in differential expression and differential transcript usage. Genome Biology 18:151. https://doi.org/10.1186/s13059-017-1277-0
#'
#' @name isPScreenAdjusted
#' @rdname isPScreenAdjusted
#' @aliases isPScreenAdjusted isPScreenAdjusted,stageR isPSCreenAdjusted,stageRTx
#' @export
setMethod("isPScreenAdjusted",signature=signature(object="stageR"),
          definition=function(object, ...){
            return(object@pScreenAdjusted)
          })
#' @describeIn isPScreenAdjusted Are the screening p-values adjusted for multiplicity?
setMethod("isPScreenAdjusted",signature=signature(object="stageRTx"),
          definition=function(object, ...){
            return(object@pScreenAdjusted)
          })

#' Has stage-wise adjustment already been performed on the object?
#'
#' This functions returns a logical stating whether the p-values have already been adjusted according to the stage-wise method.
#'
#' @param object an object of the \code{\link{stageRClass}} or \code{stageRTxClass} class.
#' @param ... Additional arguments
#' @return
#' A logical stating whether the p-values have already been adjusted according to the stage-wise method
#' @examples
#' pScreen=c(seq(1e-10,1e-2,length.out=100),seq(1e-2,.2,length.out=100),seq(.2,1,length.out=100))
#' names(pScreen)=paste0("gene",1:300)
#' pConfirmation=matrix(runif(900),nrow=300,ncol=3)
#' dimnames(pConfirmation)=list(paste0("gene",1:300),c("H1","H2","H3"))
#' stageRObj <- stageR(pScreen=pScreen, pConfirmation=pConfirmation)
#' isAdjusted(stageRObj)
#' stageRObj <- stageWiseAdjustment(stageRObj, method="holm", alpha=0.05)
#' isAdjusted(stageRObj)
#' @references
#' Van den Berge K., Soneson C., Robinson M.D., Clement L. (2017). stageR: a general stage-wise method for controlling the gene-level false discovery rate in differential expression and differential transcript usage. Genome Biology 18:151. https://doi.org/10.1186/s13059-017-1277-0
#'
#' @name isAdjusted
#' @rdname isAdjusted
#' @aliases isAdjusted isAdjusted,stageR isAdjusted,stageRTx
#' @export
setMethod("isAdjusted",signature=signature(object="stageR"),
          definition=function(object, ...){
            return(object@adjusted)
          })
#' @describeIn isAdjusted Has stage-wise adjustment already been performed on the object?
setMethod("isAdjusted",signature=signature(object="stageRTx"),
          definition=function(object, ...){
            return(object@adjusted)
          })

#' Retrieve FWER correction method.
#'
#' This functions retrieves the method used for FWER multiple testing correction in the confirmation stage of a stage-wise analysis.
#'
#' @param object an object of the \code{\link{stageRClass}} or \code{stageRTxClass} class.
#' @param ... Additional arguments
#' @return
#' Returns a character vector of length 1 specifying the FWER correction method that is used in the confirmation stage of the stage-wise analysis.
#' @examples
#' pScreen=c(seq(1e-10,1e-2,length.out=100),seq(1e-2,.2,length.out=100),seq(.2,1,length.out=100))
#' names(pScreen)=paste0("gene",1:300)
#' pConfirmation=matrix(runif(900),nrow=300,ncol=3)
#' dimnames(pConfirmation)=list(paste0("gene",1:300),c("H1","H2","H3"))
#' stageRObj <- stageR(pScreen=pScreen, pConfirmation=pConfirmation)
#' stageRObj <- stageWiseAdjustment(stageRObj, method="holm", alpha=0.05)
#' getMethod(stageRObj)
#' @references
#' Van den Berge K., Soneson C., Robinson M.D., Clement L. (2017). stageR: a general stage-wise method for controlling the gene-level false discovery rate in differential expression and differential transcript usage. Genome Biology 18:151. https://doi.org/10.1186/s13059-017-1277-0
#'
#' @name getMethod
#' @rdname getMethod
#' @aliases getMethod getMethod,stageR getMethod,stageRTx
#' @export
setMethod("getMethod",signature=signature(object="stageR"),
          definition=function(object, ...){
            return(object@method)
          })
#' @describeIn getMethod Retrieve FWER correction method.
setMethod("getMethod",signature=signature(object="stageRTx"),
          definition=function(object, ...){
            return(object@method)
          })
