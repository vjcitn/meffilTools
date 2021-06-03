
#' use epic_granges annotation to filter an EPIC RangedSummarizedExperiment to promoter regions
#' and add associated gene symbols
#' @param se RangedSummarizedExperiment with EPIC data
#' @export
subset_to_promoters = function(se) {
 stopifnot(is(se, "RangedSummarizedExperiment"))
 peg = epic_granges[-which(is.na(epic_granges$in_prom))]
 okpeg = peg[names(peg) %in% rownames(se)]
 SE_prom = se[ names(okpeg),]
 pegsym = peg$in_prom
 names(pegsym) = names(peg)
 rowRanges(SE_prom)$symbol = pegsym[names(rowRanges(SE_prom))]
 SE_prom
}
