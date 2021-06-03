#' make a data.frame tailored to COPDgene EPIC arrays, for fitting nlme::gls
#' @import SummarizedExperiment
#' @import nlme
#' @param se SummarizedExperiment with one CpG
#' @param idvar character(1) ID in colData(se)
#' @param gendvar character(1)
#' @param racevar character(1)
#' @param phasevar character(1)
#' @param smkvar character(1)
#' @param agevar character(1)
#' @param sampname character(1)
#' @note Designed to work with the HDF5SummarizedExperiments made from COPDGene EPIC arrays
#' @export
makedf_copdgene = function(se, idvar="SID", gendvar="Gender", racevar="race", phasevar="Phase_study",
   smkvar="smoking_status", agevar="age_visit", sampname="Sample.name" ) {
 stopifnot(nrow(se)==1)
 cd = colData(se)
 data.frame(beta=as.numeric(assay(se)),id=cd[[idvar]],
    ph=cd[[phasevar]], gender=cd[[gendvar]], race=cd[[racevar]], smkstat=cd[[smkvar]],
   age=cd[[agevar]], sname=cd[[sampname]])
}

#df1 = makedf(zz[1,])
#cs1 = corCompSymm(.5, ~1|id)
#csmod1 = gls(beta~race+gender+smkstat, correlation=cs1, data=df1, na.action=na.omit, subset=gender!="U")
#get_icc(csmod1)
#

#' for a vector of indices into a SummarizedExperiment, extract ICCs for a GLS model
#' @param se SummarizedExperiment
#' @param inds numeric() row indices
#' @param makedf.in function that generates a data.frame for use by nlme:gls, consumes se and dots,
#' must return a data.frame with beta and id as columns, fixedfmla will use beta and any
#' columns other than id to build model
#' @param fixedfmla formula with beta as dependent variable
#' @param \dots passed to makedf.in, within sapply
#' @examples
#' requireNamespace("HDF5Array")
#' requireNamespace("SummarizedExperiment")
#' se6 = HDF5Array::loadHDF5SummarizedExperiment(system.file("chr6_hdf5", package="meffilTools"))
#' makedf_simple = function(se, ...) {
#'   stopifnot(nrow(se)==1)
#'   data.frame(beta=as.numeric(SummarizedExperiment::assay(se)),
#'        id=as.numeric(factor(se$Sample_Group)))
#' }
#' gg = get_iccs(se6, inds=1:5, makedf.in=makedf_simple, fixedfmla=beta~1)
#' gg
#' @export
get_iccs = function(se, inds=1:10, makedf.in=makedf_copdgene, fixedfmla=beta~race+gender+smkstat, ...) {
 get_icc = function(fit) corMatrix(fit$modelStruct$corStruct)[[2]][1,2]
 ans = sapply(inds, function(x, ...) {
   df1 = makedf.in(se[x,], ...)
   cs1 = corCompSymm(.5, ~1|id)
   mod = gls(fixedfmla, correlation=cs1, data=df1, na.action=na.omit)
   get_icc(mod)
   })
 names(ans) = rownames(se)[inds]
 ans
}
