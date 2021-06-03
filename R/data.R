
#' epic annotation of CpGs
#' @docType data
#' @format data.frame
#' @note Derived from IlluminaHumanMethylationEPICanno.ilm10b4.hg19.
"epic_cpgs"

#' epic annotation of CpGs
#' @docType data
#' @format GRanges
#' @note Derived from IlluminaHumanMethylationEPICanno.ilm10b4.hg19.  Includes
#' in_prom in mcols, which gives the symbol for gene g when the associated CpG
#' is located in (TSS(g)-2000, TSS(g)+200)
"epic_granges"

#' subset of Horvath clock CpGs available on EPIC
#' @docType data
#' @format GRanges
#' @note derived from 13059_2013_3156_MOESM23_ESM.csv, supplemental to `https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115#Sec40`
"egage"

#' 353 CpGs identified in Horvath S Genome Biology 2013
#' @docType data
#' @format character vector
#' @note not all are present on EPIC array
"horvathCpGs"
