#' create a report on a filtered SummarizedExperiment
#' @param se SummarizedExperiment
#' @param sampvar character(1) name of a variable in colData to be regarded as a factor 
#' @param bspar instance of BiocSingularParam, defaults to `IrlbaParam()`
#' @param npc numeric(1) number of PCs to be computed (passed to runPCA as `rank`)
#' @export
rept_filtse = function(se, sampvar = "Sample_Group", bspar=IrlbaParam(), npc=6) {
  options(BiocSingularParam.default=bspar)
  xxt = runPCA(t(assay(se)), rank=npc, BSPARAM=bsparam())
  newdf = data.frame(y=factor(colData(se)[[sampvar]]), as.data.frame(t(assay(se))))
  rf1 = randomForest(y~., data=newdf, importance=TRUE)
  ans = list(pcaout=xxt, rfout=rf1, xlabs=colData(se)[[sampvar]])
  class(ans) = "reportObj"
  ans
}

#' summary for report
#' @param object instance of reportObj
#' @param \dots not used
#' @export
summary.reportObj = function(object, ...) {
  xxt = object$pcaout
  par(mfrow=c(1,2), mar=c(4,3,2,2))
  stats:::biplot.prcomp(xxt, xlabs=object$xlabs, cex=c(.8,.4))
  stats:::biplot.prcomp(xxt, xlabs=object$xlabs, choices=2:3, cex=c(.8,.4))
  print(object$rfout)
  varImpPlot(object$rfout, n.var=15)
}
