#' @include stageRClasses.R allGenerics.R

.stageWiseTest <- function(pScreen, pConfirmation, alpha, method=c("none","holm","dte","dtu","user"), adjustment=NULL, tx2gene=NULL, pScreenAdjusted){
    method <- match.arg(method,c("none","holm","dte","dtu","user"))

        if(method=="none"){

          if(!pScreenAdjusted) padjScreen <- p.adjust(pScreen,"BH") else padjScreen <- pScreen
          significanceOrdering <- order(padjScreen)
          genesStageI <- padjScreen<alpha
	        pAdjConfirmation <- matrix(nrow=nrow(pConfirmation),ncol=ncol(pConfirmation), dimnames=list(c(rownames(pConfirmation)),colnames(pConfirmation)))
	        pAdjConfirmation[genesStageI,] <- pConfirmation[genesStageI,]
	        padjScreenReturn=padjScreen

    } else if(method=="holm"){

      if(!pScreenAdjusted) padjScreen <- p.adjust(pScreen,"BH") else padjScreen <- pScreen
      significanceOrdering <- order(padjScreen)
      genesStageI <- padjScreen<alpha
      padjScreenReturn=padjScreen
	## only do correction for genes that passed the screening stage
      pAdjConfirmation <- matrix(nrow=nrow(pConfirmation),ncol=ncol(pConfirmation), dimnames=list(c(rownames(pConfirmation)),colnames(pConfirmation)))

      pAdjConfirmation[genesStageI,] <- t(sapply(1:length(which(genesStageI)), function(i){
      row <- pConfirmation[which(genesStageI)[i],]
        # Holm correction conditional on passing the screening stage.
        o <- order(row)
        n <- length(row)
        # Holm adjustment: passing screening stage implies 1 false hypothesis
        adjustment <- c(n-1,(n-1):1)
        rowAdjusted <- row[o]*adjustment
        rowAdjusted <- pmin(rowAdjusted,1)
        rowAdjusted <- cummax(rowAdjusted)
        rowBack <- vector(length=length(row))
        rowBack[o] <- rowAdjusted
        rowBack
      }))

    } else if(method=="user"){
      if(length(adjustment)!=ncol(pConfirmation)) stop("the length of the adjustment vector is not equal to the number of confirmation hypotheses as defined by the number of columns in pConfirmation.")
      if(!pScreenAdjusted) padjScreen <- p.adjust(pScreen,"BH") else padjScreen <- pScreen
      significanceOrdering <- order(padjScreen)
      genesStageI <- padjScreen<alpha
      padjScreenReturn=padjScreen
      pAdjConfirmation <- matrix(nrow=nrow(pConfirmation),ncol=ncol(pConfirmation), dimnames=list(c(rownames(pConfirmation)),colnames(pConfirmation)))
      pAdjConfirmation[genesStageI,] <- t(sapply(1:length(which(genesStageI)), function(i){
        row <- pConfirmation[which(genesStageI)[i],]
        o <- order(row)
		    rowAdjusted <- row[o]*adjustment
		    rowAdjusted <- pmin(rowAdjusted,1)
		    # check monotone increase of adjusted p-values
		    rowAdjusted <- cummax(rowAdjusted)
		    rowBack <- vector(length=length(row))
		    rowBack[o] <- rowAdjusted
		    rowBack
	  }))

    } else if(method=="dte"){

      if(any(is.na(match(rownames(pConfirmation),tx2gene[,1])))) stop("not all transcript names in pConfirmation match with a transcript ID from the tx2gene object.")
      if(any(is.na(match(names(pScreen),tx2gene[,2])))) stop("not all gene names in pScreen match with a gene ID from the tx2gene object.")
      # adjust screening
      if(!pScreenAdjusted) padjScreen <- p.adjust(pScreen,"BH") else padjScreen <- pScreen
      significanceOrdering <- order(padjScreen)
      genesStageI <- padjScreen<=alpha
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

      if(any(is.na(match(rownames(pConfirmation),tx2gene[,1])))) stop("not all transcript names in pConfirmation match with a transcript ID from the tx2gene object.")
      if(any(is.na(match(names(pScreen),tx2gene[,2])))) stop("not all gene names in pScreen match with a gene ID from the tx2gene object.")
      # adjust screening
      if(!pScreenAdjusted) padjScreen <- p.adjust(pScreen,"BH") else padjScreen <- pScreen
      significanceOrdering <- order(padjScreen)
      genesStageI <- padjScreen<=alpha
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
      pAdjConfirmation[names(unlist(txLevelAdjustments)),1] = unlist(txLevelAdjustments)

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
      pAdjStage=cbind(pAdjConfirmation,padjScreenReturn)[,2:1]
      colnames(pAdjStage) = c("gene","transcript")
    }
    return(list(pAdjStage=pAdjStage, alphaAdjusted=alphaAdjusted))
}

.getAdjustedP <- function(object, onlySignificantGenes=FALSE, order=TRUE){
  warning(paste0("The returned adjusted p-values are based on a stage-wise testing approach and are only valid for the provided target OFDR level of ",object@alpha*100,"%. If a different target OFDR level is of interest, the entire adjustment should be re-run."), call.=FALSE)
	if(onlySignificantGenes){ #significant genes
	    genesStageI <- object@adjustedP[,1]<=object@alpha
	    if(sum(genesStageI)==0){
	      message(paste0("No genes were found to be significant on a ",alpha*100,"% OFDR level."))
	      } else {
	      if(order){
	        sigGenes=object@adjustedP[genesStageI,]
	        o=order(sigGenes[,1])
	        return(sigGenes[o,])
	      } else {
	        sigGenes=object@adjustedP[genesStageI,]
	        return(sigGenes)
	      }
	      }
	} else { #all genes
	  if(order){
	    o=order(object@adjustedP[,1])
	    return(object@adjustedP[o,])
	  } else {
	    return(object@adjustedP)
	  }
	}
}

.getAdjustedPTx <- function(object, onlySignificantGenes=FALSE, order=TRUE){
  warning(paste0("The returned adjusted p-values are based on a stage-wise testing approach and are only valid for the provided target OFDR level of ",object@alpha*100,"%. If a different target OFDR level is of interest, the entire adjustment should be re-run."), call.=FALSE)
  tx2gene=object@tx2gene
  pConfirmation=getPConfirmation(object)
  geneForEachTx <- tx2gene[match(rownames(pConfirmation),tx2gene[,1]),2]

  if(onlySignificantGenes){ #significant genes
    genesStageI <- which(object@adjustedP[,1]<=object@alpha)
    if(sum(genesStageI)==0){
      message(paste0("No genes were found to be significant on a ",alpha*100,"% OFDR level."))
    } else {
      if(order){ #sort
        ordGenes = order(object@adjustedP[genesStageI,1])
        sigGeneIDs=unlist(lapply(strsplit(names(genesStageI),split=".",fixed=TRUE), function(x) x[1] ))
        #order acc to gene significance
        idList = sapply(unique(sigGeneIDs[ordGenes]), function(gene) which(geneForEachTx%in%gene))
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
      ordGenes = order(object@adjustedP[,1])
      sigGeneIDs=unlist(lapply(strsplit(rownames(object@adjustedP),split=".",fixed=TRUE), function(x) x[1] ))
      #order acc to gene significance
      idList = sapply(unique(sigGeneIDs[ordGenes]), function(gene) which(geneForEachTx%in%gene))
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
  if(class(object)!="stageR") stop("object should be from the stageR class.")
    results=matrix(0,nrow=nrow(object@adjustedP),ncol=ncol(object@adjustedP), dimnames=dimnames(object@adjustedP))
    results[object@adjustedP<=object@alpha] = 1
    return(results)
}


#' adjust p-values in a two-stage analysis
#'
#' This function will adjust p-values according to a hierarchical two-stage testing paradigm.
#'
#' @param object an object of the \code{\link{stageRClass}} class.
#' @param method Character string indicating the method used for FWER correction in the confirmation stage of the stage-wise analysis. Can be any of \code{"none"}, \code{"holm"}, \code{"dte"}, \code{"dtu"}, \code{"user"}. \code{"none"} will not adjust the p-values in the confirmation stage. \code{"holm"} is an adapted Holm procedure for a stage-wise analysis, where the method takes into account the fact that genes in the confirmation stage have already passed the screening stage, hence the procedure will be more powerful for the most significant p-value as compared to the standard Holm procedure. \code{"dte"} is the adjusted Holm-Shaffer procedure for differential transcript expression analysis. \code{"dtu"} is the adjusted Holm-Shaffer procedure for differential transcript usage. \code{"user"} indicates a user-defined adjustment that should be specified with the \code{adjustment} argument.
#' @param alpha the OFDR on which to control the two-stage analysis.
#' @param adjustment a user-defined adjustment of the confirmation stage p-values. Only applicable when \code{method} is \code{"none"} and ignored otherwise.
#' @param tx2gene Only applicable when  \code{method} is \code{"dte"} or \code{"dtu"}.  A \code{\link[base]{data.frame}} with transcript IDs in the first columns and gene IDs in the second column. The rownames from \code{pConfirmation} must be contained in the transcript IDs from \code{tx2gene}, and the names from \code{pScreen} must be contained in the gene IDs.
#' @references
#' Van den Berge K., Soneson C., Robinson M.D., Clement L. 2017. A general and powerful stage-wise testing procedure for differential expression and differential transcript usage. http://biorxiv.org/content/early/2017/02/16/109082
#'
#' R. Heller, E. Manduchi, G. R. Grant, and W. J. Ewens, “A flexible two-stage procedure for identifying gene sets that are differentially expressed.” Bioinformatics (Oxford, England), vol. 25, pp. 1019–25, 2009.
#'
#' S. Holm, “A Simple Sequentially Rejective Multiple Test Procedure,” Scandinavian Journal of Statistics, vol. 6, no. 2, pp. 65–70, 1979.
#'
#' J. P. Shaffer, “Modified Sequentially Rejective Multiple Test Procedures,” Journal of the American Statistical Asso- ciation, vol. 81, p. 826, 1986.
#'
#' @examples
#' pScreen=c(seq(1e-10,1e-2,length.out=100),seq(1e-2,.2,length.out=100),seq(.2,1,length.out=100))
#' names(pScreen)=paste0("gene",1:300)
#' pConfirmation=matrix(runif(900),nrow=300,ncol=3)
#' dimnames(pConfirmation)=list(paste0("gene",1:300),c("H1","H2","H3"))
#' stageRObj <- stageR(pScreen=pScreen, pConfirmation=pConfirmation)
#' stageRObj <- stageWiseAdjustment(stageRObj, method="holm", alpha=0.05)
#' getAdjustedPValues(stageRObj, onlySignificantGenes=TRUE)
#' @name stageWiseAdjustment
#' @rdname stageWiseAdjustment
#' @export
setMethod("stageWiseAdjustment",signature=signature(object="stageR", method="character", alpha="numeric"),
	  definition=function(object, method, alpha, ...){
	      pScreen=getPScreen(object)
	      pConfirmation=getPConfirmation(object)
	      pScreenAdjusted=object@pScreenAdjusted
	      stageAdjPValues <- .stageWiseTest(pScreen=pScreen, pConfirmation=pConfirmation, alpha=alpha, method=method,  pScreenAdjusted=pScreenAdjusted, ...)
	      object@adjustedP <- stageAdjPValues[["pAdjStage"]]
	      object@alphaAdjusted <- stageAdjPValues[["alphaAdjusted"]]
	      object@method <- method
	      object@alpha <- alpha
	      object@adjusted <- TRUE
	      return(object)
	  })
setMethod("stageWiseAdjustment",signature=signature(object="stageRTx", method="character", alpha="numeric"),
          definition=function(object, method, alpha, ...){
            pScreen=getPScreen(object)
            pConfirmation=getPConfirmation(object)
            pScreenAdjusted=object@pScreenAdjusted
            tx2gene=object@tx2gene
            stageAdjPValues <- .stageWiseTest(pScreen=pScreen, pConfirmation=pConfirmation, alpha=alpha, method=method,  pScreenAdjusted=pScreenAdjusted, tx2gene=tx2gene, ...)
          object@adjustedP <- stageAdjPValues[["pAdjStage"]]
          object@alphaAdjusted <- stageAdjPValues[["alphaAdjusted"]]
          object@method <- method
          object@alpha <- alpha
          object@adjusted <- TRUE
          return(object)
          })

#' Return screening hypothesis p-values from a \code{\link{stageRClass}} object.
#'
#' @param object an object of the \code{\link{stageRClass}} class.
#' @examples
#' pScreen=c(seq(1e-10,1e-2,length.out=100),seq(1e-2,.2,length.out=100),seq(.2,1,length.out=100))
#' names(pScreen)=paste0("gene",1:300)
#' pConfirmation=matrix(runif(900),nrow=300,ncol=3)
#' dimnames(pConfirmation)=list(paste0("gene",1:300),c("H1","H2","H3"))
#' stageRObj <- stageR(pScreen=pScreen, pConfirmation=pConfirmation)
#' getPScreen(stageRObj)
#' @name getPScreen
#' @rdname getPScreen
#' @export
setMethod("getPScreen",signature=signature(object="stageR"),
	  definition=function(object){return(object@pScreen)})
setMethod("getPScreen",signature=signature(object="stageRTx"),
          definition=function(object){return(object@pScreen)})

#' Return unadjusted confirmation hypothesis p-values from a \code{\link{stageRClass}} object.
#'
#' @param object an object of the \code{\link{stageRClass}} class.
#' @examples
#' pScreen=c(seq(1e-10,1e-2,length.out=100),seq(1e-2,.2,length.out=100),seq(.2,1,length.out=100))
#' names(pScreen)=paste0("gene",1:300)
#' pConfirmation=matrix(runif(900),nrow=300,ncol=3)
#' dimnames(pConfirmation)=list(paste0("gene",1:300),c("H1","H2","H3"))
#' stageRObj <- stageR(pScreen=pScreen, pConfirmation=pConfirmation)
#' getPConfirmation(stageRObj)
#' @name getPConfirmation
#' @rdname getPConfirmation
#' @export
setMethod("getPConfirmation",signature=signature(object="stageR"),
	  definition=function(object){return(object@pConfirmation)})
setMethod("getPConfirmation",signature=signature(object="stageRTx"),
          definition=function(object){return(object@pConfirmation)})


#' Retrieve the stage-wise adjusted p-values.
#'
#' This functions returns the stage-wise adjusted p-values for an object from the  \code{\link{stageRClass}} class. Note, that the p-values should have been adjusted with the \code{\link{stageWiseAdjustment}} function prior to calling this function.
#'
#' @param object an object of the \code{\link{stageRClass}} class.
#' @param onlySignificantGenes logical. If FALSE (default), all genes are returned. If TRUE, only the genes significant for the screening hypothesis are returned.
#' @param order logical. If TRUE (default), the returned matrix of adjusted p-values are ordered based on the screening hypothesis adjusted p-value.
#' @details
#' The function returns FDR adjusted p-values for the screening hypothesis and stage-wise adjusted p-values for the confirmation hypothesis p-values. For features that were not significant in the screening hypothesis, the confirmation stage adjusted p-values are set to 1.
#' @examples
#' pScreen=c(seq(1e-10,1e-2,length.out=100),seq(1e-2,.2,length.out=100),seq(.2,1,length.out=100))
#' names(pScreen)=paste0("gene",1:300)
#' pConfirmation=matrix(runif(900),nrow=300,ncol=3)
#' dimnames(pConfirmation)=list(paste0("gene",1:300),c("H1","H2","H3"))
#' stageRObj <- stageR(pScreen=pScreen, pConfirmation=pConfirmation)
#' stageRObj <- stageWiseAdjustment(stageRObj, method="holm", alpha=0.05)
#' head(getAdjustedPValues(stageRObj))
#' @references
#' Van den Berge K., Soneson C., Robinson M.D., Clement L. 2017. A general and powerful stage-wise testing procedure for differential expression and differential transcript usage. http://biorxiv.org/content/early/2017/02/16/109082
#'
#' @name getAdjustedPValues
#' @rdname getAdjustedPValues
#' @export
setMethod("getAdjustedPValues",signature=signature(object="stageR"),
	  definition=function(object, ...){
	    if(!object@adjusted) stop("adjust p-values first using stageWiseAdjustment")
	      return(.getAdjustedP(object=object, ...))
	  })
setMethod("getAdjustedPValues",signature=signature(object="stageRTx"),
          definition=function(object, ...){
            if(!object@adjusted) stop("adjust p-values first using stageWiseAdjustment")
            return(.getAdjustedPTx(object=object, ...))
          })

#' Get adjusted significance level from the screening stage.
#'
#' This functions returns the adjusted significance level from the screening stage that should be used to compare confirmation stage FWER adjusted p-values to.
#'
#' @param object an object of the \code{\link{stageRClass}} class.
#' @details
#' The adjusted significance level is calculated as the fraction of significant features in the screening stage times the alpha level.
#' @examples
#' pScreen=c(seq(1e-10,1e-2,length.out=100),seq(1e-2,.2,length.out=100),seq(.2,1,length.out=100))
#' names(pScreen)=paste0("gene",1:300)
#' pConfirmation=matrix(runif(900),nrow=300,ncol=3)
#' dimnames(pConfirmation)=list(paste0("gene",1:300),c("H1","H2","H3"))
#' stageRObj <- stageR(pScreen=pScreen, pConfirmation=pConfirmation)
#' stageRObj <- stageWiseAdjustment(stageRObj, method="holm", alpha=0.05)
#' adjustedAlphaLevel(stageRObj)
#' @references
#' Van den Berge K., Soneson C., Robinson M.D., Clement L. 2017. A general and powerful stage-wise testing procedure for differential expression and differential transcript usage. http://biorxiv.org/content/early/2017/02/16/109082
#'
#' R. Heller, E. Manduchi, G. R. Grant, and W. J. Ewens, “A flexible two-stage procedure for identifying gene sets that are differentially expressed.” Bioinformatics (Oxford, England), vol. 25, pp. 1019–25, 2009.
#'
#' @name adjustedAlphaLevel
#' @rdname adjustedAlphaLevel
#' @export
setMethod("adjustedAlphaLevel",signature=signature(object="stageR"),
	  definition=function(object){
	    if(!object@adjusted) stop("adjust p-values first using stageWiseAdjustment")
	    return(object@alphaAdjusted)
	    })
setMethod("adjustedAlphaLevel",signature=signature(object="stageRTx"),
          definition=function(object){
            if(!object@adjusted) stop("adjust p-values first using stageWiseAdjustment")
            return(object@alphaAdjusted)
            })

#' Get significance results according to a stage-wise analysis.
#'
#' This functions returns a matrix that indicates whether a specific feature is significant for a specific hypothesis of interest according to a stage-wise analysis.
#'
#' @param object an object of the \code{\link{stageRClass}} class.
#' @details
#' The FDR adjusted screening hypothesis p-values are compared to the alpha level specified. The FWER adjusted confirmation stage p-values are compared to the adjusted significance level from the screening stage.
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
#' @name getResults
#' @rdname getResults
#' @export
setMethod("getResults",signature=signature(object="stageR"),
	  definition=function(object){
	    if(!object@adjusted) stop("adjust p-values first using stageWiseAdjustment")
	    return(.getResults(object))
	    })

#' Return significant genes when performing transcript level analysis.
#'
#' This functions returns a matrix with significant genes by aggregated testing of its respective transcripts.
#'
#' @param object an object of the \code{\link{stageRClass}} class.
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
#' @name getSignificantGenes
#' @rdname getSignificantGenes
#' @export
setMethod("getSignificantGenes",signature=signature(object="stageRTx"),
          definition=function(object){
            if(!object@adjusted) stop("adjust p-values first using stageWiseAdjustment")
            if(class(object)!="stageRTx") stop("this function only works on an object of class stageRTx")
            adjustedPValues=getAdjustedPValues(object, onlySignificantGenes=FALSE, order=FALSE)
            geneIDs=adjustedPValues$geneID
            pScreenAdjusted=adjustedPValues[,"gene"]
            significantGeneIDs=which(pScreenAdjusted<=object@alpha)
            significantGeneNames=geneIDs[significantGeneIDs]
            geneAdjustedPValues=adjustedPValues[significantGeneIDs,"gene"]
            dups=duplicated(significantGeneNames)
            significantGenes=matrix(geneAdjustedPValues[!dups],ncol=1,dimnames=list(significantGeneNames[!dups],"FDR adjusted p-value"))
            return(significantGenes)
          })


#' Return significant transcripts when performing transcript level analysis.
#'
#' This functions returns a matrix with significant transctripts according to a stage-wise analysis.
#'
#' @param object an object of the \code{\link{stageRClass}} class.
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
#' @name getSignificantTx
#' @rdname getSignificantTx
#' @export
setMethod("getSignificantTx",signature=signature(object="stageRTx"),
          definition=function(object){
            if(!object@adjusted) stop("adjust p-values first using stageWiseAdjustment")
            if(class(object)!="stageRTx") stop("this function only works on an object of class stageRTx")
            adjustedPValues=getAdjustedPValues(object, onlySignificantGenes=FALSE, order=FALSE)
            txIDs=adjustedPValues$txID
            significantTxIDs=which(adjustedPValues[,"transcript"]<=object@alpha)
            significantTxNames=txIDs[significantTxIDs]
            significantTranscripts=matrix(adjustedPValues[significantTxIDs,"transcript"],ncol=1,dimnames=list(significantTxNames,"stage-wise adjusted p-value"))
            return(significantTranscripts)
          })

