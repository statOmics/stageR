#' @include stageRClasses.R allGenerics.R constructors.R

.stageWiseTest <- function(pScreen, pConfirmation, alpha, method=c("none","holm","dte","dtu","user"), adjustment=NULL, tx2gene=NULL, pScreenAdjusted, allowNA=FALSE){

  if(allowNA){
    if(any(is.na(pScreen))){
      naFeatures <- which(is.na(pScreen))
      message(paste0("Removing ",length(naFeatures)," features with NA screening hypothesis p-values. \n"))
      pScreen <- pScreen[-naFeatures]
      pConfirmation <- pConfirmation[-naFeatures,]
    }
  }

  ## check for NA values
  if(!allowNA){
    if(any(is.na(pScreen)) | any(is.na(pConfirmation)))
      stop("NA p-values found in either the screening or confirmation tests. If you want to allow for NA p-values, set allowNA=TRUE.")
  }
  method <- match.arg(method,c("none","holm","dte","dtu","user"))

  #screening stage
  if(!pScreenAdjusted)
    padjScreen <- p.adjust(pScreen,"BH") else
      padjScreen <- pScreen
  significanceOrdering <- order(padjScreen)
  genesStageI <- padjScreen<=alpha

  if(method=="none"){

    pAdjConfirmation <- matrix(nrow=nrow(pConfirmation),ncol=ncol(pConfirmation),
                               dimnames=list(c(rownames(pConfirmation)),colnames(pConfirmation)))
    pAdjConfirmation[genesStageI,] <- pConfirmation[genesStageI,]
    padjScreenReturn <- padjScreen

  } else if(method=="holm"){

    padjScreenReturn <- padjScreen
    ## only do correction for genes that passed the screening stage
    pAdjConfirmation <- matrix(nrow=nrow(pConfirmation),ncol=ncol(pConfirmation),
                               dimnames=list(c(rownames(pConfirmation)),colnames(pConfirmation)))

    for(k in seq_len(which(genesStageI))){
        row <- pConfirmation[which(genesStageI)[k],]
        # Holm correction conditional on passing the screening stage.
        o <- order(row)
        if(all(!is.na(row))){ #if no NA's, standard Holm with screening stage correction
          n <- length(row)
        } else { #if NA's present, only correct for non NA p-values
          n <- length(row[!is.na(row)])
        }
        # Holm adjustment: passing screening stage implies 1 false hypothesis
        adjustment <- c(n-1,(n-1):1)
        if(length(adjustment)!=length(row)) adjustment <- c(adjustment,
                                                           rep(1,length(row)-length(adjustment)))
        rowAdjusted <- row[o]*adjustment
        rowAdjusted <- pmin(rowAdjusted,1)
        rowAdjusted <- cummax(rowAdjusted)
        rowBack <- vector(length=length(row))
        rowBack[o] <- rowAdjusted
        pAdjConfirmation[genesStageI[k],] <- rowBack
    }

  } else if(method=="user"){
    if(length(adjustment)!=ncol(pConfirmation))
      stop("the length of the adjustment vector is not equal to the number of confirmation hypotheses as defined by the number of columns in pConfirmation.")
    padjScreenReturn=padjScreen
    pAdjConfirmation <- matrix(nrow=nrow(pConfirmation),ncol=ncol(pConfirmation),
                               dimnames=list(c(rownames(pConfirmation)),colnames(pConfirmation)))
    for(k in seq_len(which(genesStageI))){
        row <- pConfirmation[which(genesStageI)[k],]
        o <- order(row)
        rowAdjusted <- row[o]*adjustment
        rowAdjusted <- pmin(rowAdjusted,1)
        # check monotone increase of adjusted p-values
        rowAdjusted <- cummax(rowAdjusted)
        rowBack <- vector(length=length(row))
        rowBack[o] <- rowAdjusted
        rowBack
      pAdjConfirmation[genesStageI[k],] <- rowBack
    }

  } else if(method=="dte"){

    if(any(is.na(match(rownames(pConfirmation),tx2gene[,1]))))
      stop("not all transcript names in pConfirmation match with a transcript ID from the tx2gene object.")
    if(any(is.na(match(names(pScreen),tx2gene[,2]))))
      stop("not all gene names in pScreen match with a gene ID from the tx2gene object.")
    significantGenes <- names(padjScreen)[genesStageI]
    geneForEachTx <- tx2gene[match(rownames(pConfirmation),tx2gene[,1]),2]
    txLevelAdjustments <- sapply(significantGenes,function(gene){
      id <- which(geneForEachTx %in% gene)
      row <- pConfirmation[id,]
      #make sure names are passed along if only one tx
      if(length(id)==1) names(row)=rownames(pConfirmation)[id]
      o <- order(row)
      n <- length(row)
      # DTE adjustment: passing screening stage implies 1 false hypothesis
      if(n==1) adjustment=0 else adjustment=c(n-1,(n-1):1)
      rowAdjusted <- row[o]*adjustment
      rowAdjusted <- pmin(rowAdjusted,1)
      rowAdjusted <- cummax(rowAdjusted)
      rowBack <- vector(length=length(row))
      rowBack[o] <- rowAdjusted
      names(rowBack) <- names(row)
      rowBack
    })
    pAdjConfirmation <- matrix(nrow=nrow(pConfirmation),ncol=1)
    rownames(pAdjConfirmation) <- paste0(geneForEachTx,".",rownames(pConfirmation))
    # adjusted p-values for screening hypothesis
    padjScreenReturn <- padjScreen[geneForEachTx]
    # adjusted p-values for confirmation hypothesis
    pAdjConfirmation[names(unlist(txLevelAdjustments)),1] = unlist(txLevelAdjustments)

  } else if(method=="dtu"){

    if(any(is.na(match(rownames(pConfirmation),tx2gene[,1]))))
      stop("not all transcript names in pConfirmation match with a transcript ID from the tx2gene object.")
    if(any(is.na(match(names(pScreen),tx2gene[,2]))))
      stop("not all gene names in pScreen match with a gene ID from the tx2gene object.")
    # adjust screening
    significantGenes <- names(padjScreen)[genesStageI]
    geneForEachTx <- as.character(tx2gene[match(rownames(pConfirmation),tx2gene[,1]),2])
    txLevelAdjustments <- sapply(significantGenes,function(gene){
      id <- which(geneForEachTx %in% gene)
      row <- pConfirmation[id,]
      o <- order(row)
      n <- length(row)
      # DTU adjustment: passing screening stage implies 2 false hypotheses
      if(n==2) adjustment=c(0,0) else adjustment=c(n-2,n-2,(n-2):1)
      rowAdjusted <- row[o]*adjustment
      rowAdjusted <- pmin(rowAdjusted,1)
      rowAdjusted <- cummax(rowAdjusted)
      rowBack <- vector(length=length(row))
      rowBack[o] <- rowAdjusted
      names(rowBack) <- names(row)
      rowBack
    })
    pAdjConfirmation <- matrix(nrow=nrow(pConfirmation),ncol=1)
    rownames(pAdjConfirmation) <- paste0(geneForEachTx,".",rownames(pConfirmation))
    # adjusted p-values for screening hypothesis
    padjScreenReturn <- padjScreen[as.character(geneForEachTx)]
    # adjusted p-values for confirmation hypothesis
    pAdjConfirmation[names(unlist(txLevelAdjustments)),1] <- unlist(txLevelAdjustments)

  } else stop("method must be either one of 'holm' or ... ")

  #BH-adjusted s.l.
  alphaAdjusted <- sum(padjScreen<=alpha)/length(padjScreen)*alpha
  #Correct FWER-adjusted p-values acc. to BH-adjusted s.l.
  pAdjConfirmation[!is.na(pAdjConfirmation)] <- pmin(pAdjConfirmation[!is.na(pAdjConfirmation)]*length(padjScreen)/sum(padjScreen<=alpha),1)
  if(!(method %in% c("dte","dtu"))){
    pAdjStage <- cbind(padjScreenReturn,pAdjConfirmation)
    colnames(pAdjStage)[1] <- "padjScreen"
  }
  if(method %in% c("dte","dtu")){
    pAdjStage <- cbind(pAdjConfirmation,padjScreenReturn)[,2:1]
    colnames(pAdjStage) <- c("gene","transcript")
  }
  return(list(pAdjStage=pAdjStage, alphaAdjusted=alphaAdjusted))
}

.getAdjustedP <- function(object, onlySignificantGenes=FALSE, order=TRUE){
  ## this function is used in getAdjustedPValues to return the adjusted p-values for a stageR class.
  warning(paste0("The returned adjusted p-values are based on a stage-wise testing approach and are only valid for the provided target OFDR level of ",getAlpha(object)*100,"%. If a different target OFDR level is of interest, the entire adjustment should be re-run. \n"), call.=FALSE)
  if(onlySignificantGenes){ #significant genes
    genesStageI <- object@adjustedP[,"padjScreen"]<=getAlpha(object)
    if(sum(genesStageI)==0){
      message(paste0("No genes were found to be significant on a ",alpha*100,"% OFDR level."))
    } else {
      if(order){
        sigGenes <- object@adjustedP[genesStageI,]
        o <- order(sigGenes[,"padjScreen"])
        return(sigGenes[o,])
      } else {
        sigGenes <- object@adjustedP[genesStageI,]
        return(sigGenes)
      }
    }
  } else { #all genes
    if(order){
      o <- order(object@adjustedP[,"padjScreen"])
      return(object@adjustedP[o,])
    } else {
      return(object@adjustedP)
    }
  }
}

.getAdjustedPTx <- function(object, onlySignificantGenes=FALSE, order=TRUE){
  ## this function is used in getAdjustedPValues to return the adjusted p-values for a stageRTx class.
  warning(paste0("The returned adjusted p-values are based on a stage-wise testing approach and are only valid for the provided target OFDR level of ",getAlpha(object)*100,"%. If a different target OFDR level is of interest, the entire adjustment should be re-run. \n"), call.=FALSE)
  tx2gene <- getTx2gene(object)
  pConfirmation <- getPConfirmation(object)
  geneForEachTx <- tx2gene[match(rownames(pConfirmation),tx2gene[,1]),2]

  if(onlySignificantGenes){ #significant genes
    genesStageI <- which(object@adjustedP[,"gene"]<=getAlpha(object))
    if(sum(genesStageI)==0){
      message(paste0("No genes were found to be significant on a ",alpha*100,"% OFDR level."))
    } else {
      if(order){ #sort
        ordGenes <- order(object@adjustedP[genesStageI,1])
        sigGeneIDs <- unlist(lapply(strsplit(names(genesStageI),split=".",fixed=TRUE), function(x) x[1] ))
        #order acc to gene significance
        idList <- sapply(unique(sigGeneIDs[ordGenes]), function(gene) which(geneForEachTx%in%gene))
        #order tx within gene
        idListOrdTx <- lapply(idList, function(x) x[order(pConfirmation[x,])])
        outData <- object@adjustedP[unlist(idListOrdTx),]
        outData <- data.frame("geneID"=sigGeneIDs[ordGenes],"txID"=unlist(lapply(strsplit(rownames(outData),split=".",fixed=TRUE), function(x) x[2] )),outData, row.names=NULL)
        return(outData)
      } else { #dont sort
        outData <- object@adjustedP[genesStageI,]
        outData <- data.frame("geneID"=unlist(lapply(strsplit(rownames(outData),split=".",fixed=TRUE), function(x) x[1] )),"txID"=unlist(lapply(strsplit(rownames(outData),split=".",fixed=TRUE), function(x) x[2] )),outData, row.names=NULL)
        return(outData)
      }
    }
  } else { #all genes
    if(order){ #sort
      ordGenes <- order(object@adjustedP[,"gene"])
      sigGeneIDs <- unlist(lapply(strsplit(rownames(object@adjustedP),split=".",fixed=TRUE), function(x) x[1] ))
      #order acc to gene significance
      idList <- sapply(unique(sigGeneIDs[ordGenes]), function(gene) which(geneForEachTx%in%gene))
      #order tx within gene
      idListOrdTx <- lapply(idList, function(x) x[order(pConfirmation[x,])])
      outData <- object@adjustedP[unlist(idListOrdTx),]
      outData <- data.frame("geneID"=sigGeneIDs[ordGenes],"txID"=unlist(lapply(strsplit(rownames(outData),split=".",fixed=TRUE), function(x) x[2] )),outData, row.names=NULL)
      return(outData)
    } else { #dont sort
      outData <- object@adjustedP
      outData <- data.frame("geneID"=unlist(lapply(strsplit(rownames(outData),split=".",fixed=TRUE), function(x) x[1] )),"txID"=unlist(lapply(strsplit(rownames(outData),split=".",fixed=TRUE), function(x) x[2] )),outData, row.names=NULL)
      return(outData)
    }
  }
}

.getResults <- function(object){
  adjustedPValues <- getAdjustedPValues(object, onlySignificantGenes=FALSE, order=FALSE)
  results <- matrix(0,nrow=nrow(adjustedPValues),ncol=ncol(adjustedPValues), dimnames=dimnames(adjustedPValues))
  results[adjustedPValues<=getAlpha(object)] = 1
  return(results)
}


#' adjust p-values in a two-stage analysis
#'
#' This function will adjust p-values according to a hierarchical two-stage testing paradigm.
#'
#' @param object an object of the \code{\link{stageRClass}} class.
#' @param method Character string indicating the method used for FWER correction in the confirmation stage of the stage-wise analysis. Can be any of \code{"none"}, \code{"holm"}, \code{"dte"}, \code{"dtu"}, \code{"user"}. \code{"none"} will not adjust the p-values in the confirmation stage. \code{"holm"} is an adapted Holm procedure for a stage-wise analysis, where the method takes into account the fact that genes in the confirmation stage have already passed the screening stage, hence the procedure will be more powerful for the most significant p-value as compared to the standard Holm procedure. \code{"dte"} is the adjusted Holm-Shaffer procedure for differential transcript expression analysis. \code{"dtu"} is the adjusted Holm-Shaffer procedure for differential transcript usage. \code{"user"} indicates a user-defined adjustment that should be specified with the \code{adjustment} argument.
#' @param alpha the OFDR on which to control the two-stage analysis.
#' @param tx2gene Only applicable when  \code{method} is \code{"dte"} or \code{"dtu"}.  A \code{\link[base]{data.frame}} with transcript IDs in the first columns and gene IDs in the second column. The rownames from \code{pConfirmation} must be contained in the transcript IDs from \code{tx2gene}, and the names from \code{pScreen} must be contained in the gene IDs.
#' @param adjustment a user-defined adjustment of the confirmation stage p-values. Only applicable when \code{method} is \code{"user"} and ignored otherwise.
#' @param ... Additional arguments passed to \code{stageWiseTest}
#' @return
#' A stageR/stageRTx object with stage-wise adjusted p-values.
#' @references
#' Van den Berge K., Soneson C., Robinson M.D., Clement L. 2017. A general and powerful stage-wise testing procedure for differential expression and differential transcript usage. \url{http://biorxiv.org/content/early/2017/02/16/109082}
#' R. Heller, E. Manduchi, G. R. Grant, and W. J. Ewens, "A flexible two-stage procedure for identifying gene sets that are differentially expressed." Bioinformatics (Oxford, England), vol. 25, pp. 1019-25, 2009.
#'
#' S. Holm, "A Simple Sequentially Rejective Multiple Test Procedure," Scandinavian Journal of Statistics, vol. 6, no. 2, pp. 65-70, 1979.
#' J. P. Shaffer, "Modified Sequentially Rejective Multiple Test Procedures," Journal of the American Statistical Association, vol. 81, p. 826, 1986.
#' @examples
#' pScreen=c(seq(1e-10,1e-2,length.out=100),seq(1e-2,.2,length.out=100),seq(.2,1,length.out=100))
#' names(pScreen)=paste0("gene",1:300)
#' pConfirmation=matrix(runif(900),nrow=300,ncol=3)
#' dimnames(pConfirmation)=list(paste0("gene",1:300),c("H1","H2","H3"))
#' stageRObj <- stageR(pScreen=pScreen, pConfirmation=pConfirmation)
#' stageRObj <- stageWiseAdjustment(stageRObj, method="holm", alpha=0.05)
#' getAdjustedPValues(stageRObj, onlySignificantGenes=TRUE, order=TRUE)
#'# @name stageWiseAdjustment
#'# @rdname stageWiseAdjustment
#' @describeIn stageWiseAdjustment Adjust p-values in a two-stage analysis
#' @export
setMethod("stageWiseAdjustment",signature=signature(object="stageR", method="character", alpha="numeric"),
          definition=function(object, method, alpha, adjustment=NULL, ...){
            pScreen <- getPScreen(object)
            pConfirmation <- getPConfirmation(object)
            pScreenAdjusted <- isPScreenAdjusted(object)
            stageAdjPValues <- .stageWiseTest(pScreen=pScreen, pConfirmation=pConfirmation, alpha=alpha, method=method,  pScreenAdjusted=pScreenAdjusted, adjustment=adjustment, ...)
            object@adjustedP <- stageAdjPValues[["pAdjStage"]]
            object@alphaAdjusted <- stageAdjPValues[["alphaAdjusted"]]
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
    pConfirmation <- getPConfirmation(object)
    pScreenAdjusted <- isPScreenAdjusted(object)
    tx2gene <- getTx2gene(object)
    stageAdjPValues <-
      .stageWiseTest(
        pScreen = pScreen,
        pConfirmation = pConfirmation,
        alpha = alpha,
        method = method,
        pScreenAdjusted = pScreenAdjusted,
        tx2gene = tx2gene,
        ...
      )
    object@adjustedP <- stageAdjPValues[["pAdjStage"]]
    object@alphaAdjusted <-
      stageAdjPValues[["alphaAdjusted"]]
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
#'# @name getPScreen
#'# @rdname getPScreen
#' @describeIn getPScreen Return screening hypothesis p-values from a \code{\link{stageRClass}} object.
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
#'# @name getPConfirmation
#'# @rdname getPConfirmation
#' @describeIn getPConfirmation Return unadjusted confirmation hypothesis p-values from a \code{\link{stageRClass}} object.
#' @export
setMethod("getPConfirmation",signature=signature(object="stageR"),
          definition=function(object){return(object@pConfirmation)})
#' @describeIn getPConfirmation Return unadjusted confirmation hypothesis p-values from a \code{\link{stageRClass}} object.
setMethod("getPConfirmation",signature=signature(object="stageRTx"),
          definition=function(object){return(object@pConfirmation)})


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
#' @examples
#' pScreen=c(seq(1e-10,1e-2,length.out=100),seq(1e-2,.2,length.out=100),seq(.2,1,length.out=100))
#' names(pScreen)=paste0("gene",1:300)
#' pConfirmation=matrix(runif(900),nrow=300,ncol=3)
#' dimnames(pConfirmation)=list(paste0("gene",1:300),c("H1","H2","H3"))
#' stageRObj <- stageR(pScreen=pScreen, pConfirmation=pConfirmation)
#' stageRObj <- stageWiseAdjustment(stageRObj, method="holm", alpha=0.05)
#' head(getAdjustedPValues(stageRObj, onlySignificantGenes=TRUE, order=TRUE))
#' @references
#' Van den Berge K., Soneson C., Robinson M.D., Clement L. 2017. A general and powerful stage-wise testing procedure for differential expression and differential transcript usage. http://biorxiv.org/content/early/2017/02/16/109082
#'
#'# @name getAdjustedPValues
#'# @rdname getAdjustedPValues
#' @describeIn getAdjustedPValues Retrieve the stage-wise adjusted p-values.
#' @export
setMethod("getAdjustedPValues",signature=signature(object="stageR", onlySignificantGenes="logical", order="logical"),
          definition=function(object, onlySignificantGenes, order, ...){
            if(!isAdjusted(object)) stop("adjust p-values first using stageWiseAdjustment")
            return(.getAdjustedP(object=object, ...))
          })
#' @describeIn getAdjustedPValues Retrieve the stage-wise adjusted p-values.
setMethod("getAdjustedPValues",signature=signature(object="stageRTx", onlySignificantGenes="logical", order="logical"),
          definition=function(object, onlySignificantGenes, order, ...){
            if(!isAdjusted(object)) stop("adjust p-values first using stageWiseAdjustment")
            return(.getAdjustedPTx(object=object, ...))
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
#' Van den Berge K., Soneson C., Robinson M.D., Clement L. 2017. A general and powerful stage-wise testing procedure for differential expression and differential transcript usage. \url{http://biorxiv.org/content/early/2017/02/16/109082}
#'
#' R. Heller, E. Manduchi, G. R. Grant, and W. J. Ewens, "A flexible two-stage procedure for identifying gene sets that are differentially expressed." Bioinformatics (Oxford, England), vol. 25, pp. 1019-25, 2009.
#'
#' @seealso \code{\link{stageR}}, \code{\link{stageRClass}}
#'# @name adjustedAlphaLevel
#'# @rdname adjustedAlphaLevel
#' @export
#' @describeIn adjustedAlphaLevel Get adjusted significance level from the screening stage.
setMethod("adjustedAlphaLevel",signature=signature(object="stageR"),
          definition=function(object){
            if(!isAdjusted(object)) stop("adjust p-values first using stageWiseAdjustment")
            return(object@alphaAdjusted)
          })
#' @describeIn adjustedAlphaLevel Get adjusted significance level from the screening stage.
setMethod("adjustedAlphaLevel",signature=signature(object="stageRTx"),
          definition=function(object){
            if(!isAdjusted(object)) stop("adjust p-values first using stageWiseAdjustment")
            return(object@alphaAdjusted)
          })

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
#' Van den Berge K., Soneson C., Robinson M.D., Clement L. 2017. A general and powerful stage-wise testing procedure for differential expression and differential transcript usage. http://biorxiv.org/content/early/2017/02/16/109082
#'
#'# @name getResults
#'# @rdname getResults
#' @describeIn getResults Get significance results according to a stage-wise analysis.
#' @export
setMethod("getResults",signature=signature(object="stageR"),
          definition=function(object){
            if(!isAdjusted(object)) stop("adjust p-values first using stageWiseAdjustment")
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
#' pScreen=c(seq(1e-10,1e-2,length.out=nGenes-100),seq(1e-2,.2,length.out=50),seq(50))
#' names(pScreen)=names(table(genes)) #discards genes that are not simulated
#' pConfirmation=matrix(runif(1000),nrow=1000,ncol=1)
#' rownames(pConfirmation)=transcripts
#' stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation ,pScreenAdjusted=TRUE, tx2gene=tx2gene)
#' stageRObj <- stageWiseAdjustment(stageRObj, method="dte", alpha=0.05)
#' head(getSignificantGenes(stageRObj))
#' @references
#' Van den Berge K., Soneson C., Robinson M.D., Clement L. 2017. A general and powerful stage-wise testing procedure for differential expression and differential transcript usage. http://biorxiv.org/content/early/2017/02/16/109082
#'
#'# @name getSignificantGenes
#'# @rdname getSignificantGenes
#' @describeIn getSignificantGenes Return significant genes when performing transcript level analysis.
#' @export
setMethod("getSignificantGenes",signature=signature(object="stageRTx"),
          definition=function(object){
            if(!isAdjusted(object)) stop("adjust p-values first using stageWiseAdjustment")
            if(class(object)!="stageRTx") stop("this function only works on an object of class stageRTx")
            adjustedPValues <- getAdjustedPValues(object, onlySignificantGenes=FALSE, order=FALSE)
            geneIDs <- adjustedPValues$geneID
            pScreenAdjusted <- adjustedPValues[,"gene"]
            significantGeneIDs <- which(pScreenAdjusted<=getAlpha(object))
            significantGeneNames <- geneIDs[significantGeneIDs]
            geneAdjustedPValues <- adjustedPValues[significantGeneIDs,"gene"]
            dups <- duplicated(significantGeneNames)
            significantGenes <- matrix(geneAdjustedPValues[!dups],ncol=1,dimnames=list(significantGeneNames[!dups],"FDR adjusted p-value"))
            return(significantGenes)
          })


#' Return significant transcripts when performing transcript-level analysis.
#'
#' This functions returns a matrix with significant transctripts according to a stage-wise analysis.
#'
#' @param object an object of the \code{\link{stageRClass}} class.
#' @return
#' A matrix of significant transcripts with their corresponding stage-wise adjusted p-value (i.e. FDR and FWER correction).
#' @references
#' Van den Berge K., Soneson C., Robinson M.D., Clement L. 2017. A general and powerful stage-wise testing procedure for differential expression and differential transcript usage. Submitted.
#' @examples
#' #make identifiers linking transcripts to genes
#' set.seed(1)
#' genes=paste0("gene",sample(1:200,1000,replace=TRUE))
#' nGenes=length(table(genes))
#' transcripts=paste0("tx",1:1000)
#' tx2gene=data.frame(transcripts,genes)
#' #gene-wise q-values
#' pScreen=c(seq(1e-10,1e-2,length.out=nGenes-100),seq(1e-2,.2,length.out=50),seq(50))
#' names(pScreen)=names(table(genes)) #discards genes that are not simulated
#' pConfirmation=matrix(runif(1000),nrow=1000,ncol=1)
#' rownames(pConfirmation)=transcripts
#' stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation ,pScreenAdjusted=TRUE, tx2gene=tx2gene)
#' stageRObj <- stageWiseAdjustment(stageRObj, method="dte", alpha=0.05)
#' head(getSignificantTx(stageRObj))
#'# @name getSignificantTx
#'# @rdname getSignificantTx
#' @describeIn getSignificantTx Return significant transcripts when performing transcript-level analysis.
#' @export
setMethod("getSignificantTx",signature=signature(object="stageRTx"),
          definition=function(object){
            if(!isAdjusted(object)) stop("adjust p-values first using stageWiseAdjustment")
            if(class(object)!="stageRTx") stop("this function only works on an object of class stageRTx")
            adjustedPValues <- getAdjustedPValues(object, onlySignificantGenes=FALSE, order=FALSE)
            txIDs <- adjustedPValues$txID
            significantTxIDs <- which(adjustedPValues[,"transcript"]<=getAlpha(object))
            significantTxNames <- txIDs[significantTxIDs]
            significantTranscripts <- matrix(adjustedPValues[significantTxIDs,"transcript"],ncol=1,dimnames=list(significantTxNames,"stage-wise adjusted p-value"))
            return(significantTranscripts)
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
#' Van den Berge K., Soneson C., Robinson M.D., Clement L. 2017. A general and powerful stage-wise testing procedure for differential expression and differential transcript usage. http://biorxiv.org/content/early/2017/02/16/109082
#'
#'# @name getAlpha
#'# @rdname getAlpha
#' @describeIn getAlpha Retrieve the significance level for the stage-wise adjustment.
#' @export
setMethod("getAlpha",signature=signature(object="stageR"),
          definition=function(object, ...){
            if(is.null(object@alpha)) stop("No significance level was specified yet. Maybe you need to run stageWiseAdjustment first.")
            return(object@alpha)
          })
#' @describeIn getAlpha Retrieve the significance level for the stage-wise adjustment.
setMethod("getAlpha",signature=signature(object="stageRTx"),
          definition=function(object, ...){
            if(is.null(object@alpha)) stop("No significance level was specified yet. Maybe you need to run stageWiseAdjustment first.")
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
#' pScreen=c(seq(1e-10,1e-2,length.out=nGenes-100),seq(1e-2,.2,length.out=50),seq(50))
#' names(pScreen)=names(table(genes)) #discards genes that are not simulated
#' pConfirmation=matrix(runif(1000),nrow=1000,ncol=1)
#' rownames(pConfirmation)=transcripts
#' stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation ,pScreenAdjusted=TRUE, tx2gene=tx2gene)
#' getTx2gene(stageRObj)
#' @references
#' Van den Berge K., Soneson C., Robinson M.D., Clement L. 2017. A general and powerful stage-wise testing procedure for differential expression and differential transcript usage. http://biorxiv.org/content/early/2017/02/16/109082
#'
#'# @name getTx2gene
#'# @rdname getTx2gene
#' @describeIn getTx2gene Retrieve the data frame linking genes to transcripts.
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
#' Van den Berge K., Soneson C., Robinson M.D., Clement L. 2017. A general and powerful stage-wise testing procedure for differential expression and differential transcript usage. http://biorxiv.org/content/early/2017/02/16/109082
#'
#'# @name isPScreenAdjusted
#'# @rdname isPScreenAdjusted
#' @describeIn isPScreenAdjusted Are the screening p-values adjusted for multiplicity?
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
#' Van den Berge K., Soneson C., Robinson M.D., Clement L. 2017. A general and powerful stage-wise testing procedure for differential expression and differential transcript usage. http://biorxiv.org/content/early/2017/02/16/109082
#'
#'# @name isAdjusted
#'# @rdname isAdjusted
#' @describeIn isAdjusted Has stage-wise adjustment already been performed on the object?
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
#' Van den Berge K., Soneson C., Robinson M.D., Clement L. 2017. A general and powerful stage-wise testing procedure for differential expression and differential transcript usage. http://biorxiv.org/content/early/2017/02/16/109082
#'
#'# @name getMethod
#'# @rdname getMethod
#' @describeIn getMethod Retrieve FWER correction method.
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
