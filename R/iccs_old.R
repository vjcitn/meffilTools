#' make a data.frame tailored COPDgene EPIC arrays, for fitting nlme::gls
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
makedf_old= function(se, idvar="SID", gendvar="Gender", racevar="race", phasevar="Phase_study",
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
#' @param fixedfmla formula with beta as dependent variable
get_iccs_old = function(se, inds=1:10, fixedfmla=beta~race+gender+smkstat) {
 get_icc = function(fit) corMatrix(fit$modelStruct$corStruct)[[2]][1,2]
 ans = sapply(inds, function(x) {
   df1 = makedf(se[x,])
   cs1 = corCompSymm(.5, ~1|id)
   mod = gls(fixedfmla, correlation=cs1, data=df1, na.action=na.omit, subset=gender!="U")
   get_icc(mod)
   })
 names(ans) = rownames(se)[inds]
 ans
}
