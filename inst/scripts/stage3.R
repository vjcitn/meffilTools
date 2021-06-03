library(meffilTools)
data(epic_granges)
library(SummarizedExperiment)
library(GenomicRanges)
library(gdsfmt)
library(HDF5Array)

m = getmat()
ii = index.gdsn(m$obj, "matrix")
mat = read.gdsn(ii)
cpgs = read.gdsn(index.gdsn(m$obj, "row.names"))
samps = read.gdsn(index.gdsn(m$obj, "col.names"))
dimnames(mat) = list(cpgs, samps)

mat2se = function(mat) {
  ss = read.csv(meffilTools::sample_sheet_path())
  sss = ss[-which(duplicated(ss$Sample.name)),]
  rownames(sss) = sss$Sample.name
  skp = sss[colnames(mat),]
  ini = SummarizedExperiment::SummarizedExperiment(mat)
  colData(ini) = S4Vectors::DataFrame(skp)
  assayNames(ini) = "betas"
  ini
}
 


makeses = function(inmat, tag="chr20") {
 cur = epic_granges[which(seqnames(epic_granges)==tag)]
 curmat = inmat[names(cur),]
 se = mat2se(curmat)
 rowRanges(se) = epic_granges[rownames(se)]
 tag = gsub("chr", "", tag)
 save(se, file=paste0("se", tag, ".rda"))
 saveHDF5SummarizedExperiment(se, paste0("se", tag, "_hdf5"))
 se
 }


#makeses(mat, "chr11")
#makeses(mat, "chr10")
#makeses(mat, "chr9")
#makeses(mat, "chr8")
#makeses(mat, "chr7")
#gc()
#makeses(mat, "chr6")
#gc()
#makeses(mat, "chr5")
#gc()
#makeses(mat, "chr4")
#gc()
#makeses(mat, "chr3")
#gc()
#makeses(mat, "chr2")
#gc()
#makeses(mat, "chr1")
#gc()
#makeses(mat, "chrX")
makeses(mat, "chrY")


