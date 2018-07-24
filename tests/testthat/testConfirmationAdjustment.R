library(testthat)
library(stageR)

context("Test the confirmation p-value correction methods.")

## note that there is no screening selection (all genes pass).
set.seed(723)
pTest = matrix(runif(150,1e-10,1e-1),nrow=15,ncol=10, dimnames=list(paste0("gene",1:15), paste0("H",1:10)))
pScreen = rep(1e-5,15) ; names(pScreen)=rownames(pTest)
stageRObj = stageR(pScreen=pScreen, pConfirmation=pTest)

## Holm
pTestHolm = t(apply(pTest,1,function(row){
  o=order(row)
  row=cummax(pmin(row[o]*c(9,9:1),1))
  rowBack=vector(length=length(row))
  rowBack[o]=row
  rowBack
}))
colnames(pTestHolm)=paste0("H",1:10)

hlp <- getAdjustedPValues(stageWiseAdjustment(stageRObj, method="holm",alpha=0.05), FALSE, FALSE)
pHolmStageR <- do.call(rbind,lapply(map(hlp$data,t), function(x) as.numeric(x["padj_SW",])))
dimnames(pHolmStageR) <- dimnames(pTestHolm)

test_that("Test that Holm correction is correct",{
  expect_equal(pHolmStageR, pTestHolm, tolerance=1e-6)
})

## user
adjustment=1:10
pTestUser = t(apply(pTest,1,function(row){
  o=order(row)
  row=cummax(pmin(row[o]*adjustment,1))
  rowBack=vector(length=length(row))
  rowBack[o]=row
  rowBack
}))
colnames(pTestUser)=paste0("H",1:10)

hlp <- getAdjustedPValues(stageWiseAdjustment(stageRObj, method="user",alpha=0.05,adjustment=1:10), FALSE, FALSE)
pUserStageR <- do.call(rbind,lapply(map(hlp$data,t), function(x) as.numeric(x["padj_SW",])))
dimnames(pUserStageR) <- dimnames(pTestUser)


test_that("Test that user correction is correct",{
  expect_equal(pUserStageR, pTestUser, tolerance=1e-6)
})

## none
adjustment=rep(1,10)
pTestNone = t(apply(pTest,1,function(row){
  o=order(row)
  row=cummax(pmin(row[o]*adjustment,1))
  rowBack=vector(length=length(row))
  rowBack[o]=row
  rowBack
}))
colnames(pTestNone)=paste0("H",1:10)

hlp <- getAdjustedPValues(stageWiseAdjustment(stageRObj, method="none",alpha=0.05), FALSE, FALSE)
pNoneStageR <- do.call(rbind,lapply(map(hlp$data,t), function(x) as.numeric(x["padj_SW",])))
dimnames(pNoneStageR) <- dimnames(pTestNone)


test_that("Test that none correction is correct",{
  expect_equal(pNoneStageR, pTestNone, tolerance=1e-6)
})

## DTE
pScreen=rep(1e-5,9)
names(pScreen)=paste0("gene",1:9)
pTx=pTest[,1,drop=FALSE]
rownames(pTx)=paste0("transcript",1:15)
tx2gene = data.frame(transcript=paste0("transcript",1:15), gene=paste0("gene",c(rep(c(1,2),each=4),3:9)))
stageRTxObj = stageRTx(pScreen=pScreen, pConfirmation=pTx, tx2gene=tx2gene)
# adjust manually
#gene1
gene1P = pTx[1:4,]
o=order(gene1P)
gene1PAdj <- vector(length=length(gene1P))
gene1PAdj[o] = cummax(pmin(gene1P[o]*c(3,3,2,1),1))
#gene2
gene2P = pTx[5:8,]
o=order(gene2P)
gene2PAdj <- vector(length=length(gene2P))
gene2PAdj[o] = cummax(pmin(gene2P[o]*c(3,3,2,1),1))
#others
geneOthers=rep(0,7)
allAdjP=unname(c(gene1PAdj, gene2PAdj, geneOthers))

#stageR
stageRTxObjDTE <- stageWiseAdjustment(stageRTxObj, method="dte", alpha=0.05)
hlp <- getAdjustedPValues(stageRTxObjDTE, FALSE, FALSE)
pDTEStageR <- unlist(lapply(map(hlp$data,t), function(x) as.numeric(x["padj_SW",])))

test_that("Test that DTE correction is correct",{
  expect_equal(pDTEStageR, allAdjP, tolerance=1e-6)
})

## DTU
pScreen=rep(1e-5,2)
names(pScreen)=paste0("gene",1:2)
pTx=pTest[1:8,1,drop=FALSE]/10
rownames(pTx)=paste0("transcript",1:8)
tx2gene = data.frame(transcript=paste0("transcript",1:8), gene=paste0("gene",rep(c(1,2),each=4)))
stageRTxObj = stageRTx(pScreen=pScreen, pConfirmation=pTx, tx2gene=tx2gene)
# adjust manually
#gene1
gene1P = pTx[1:4,]
o=order(gene1P)
gene1PAdj = cummax(pmin(gene1P[o]*c(2,2,2,1),1))
gene1Back=vector(length=length(gene1PAdj))
gene1Back[o] = gene1PAdj
#gene2
gene2P = pTx[5:8,]
o=order(gene2P)
gene2PAdj = cummax(pmin(gene2P[o]*c(2,2,2,1),1))
gene2Back=vector(length=length(gene2PAdj))
gene2Back[o] = gene2PAdj

allAdjP=unname(c(gene1Back, gene2Back))

#stageR
stageRTxObjDTU <- stageWiseAdjustment(stageRTxObj, method="dtu", alpha=0.05)
hlp <- getAdjustedPValues(stageRTxObjDTU, FALSE, FALSE)
pDTUStageR <- unlist(lapply(map(hlp$data,t), function(x) as.numeric(x["padj_SW",])))

test_that("Test that DTU correction is correct",{
  expect_equal(pDTUStageR, allAdjP, tolerance=1e-6)
})

rm(pTest, pScreen, stageRObj, pTestHolm, pTestUser, pTestNone, pTx, tx2gene, stageRTxObj, gene1P, o, gene1PAdj, gene2P, gene2PAdj, geneOthers, allAdjP, gene1Back, gene2Back)