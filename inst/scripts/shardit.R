
library(meffilTools)
data(epic_granges)
library(SummarizedExperiment)
library(GenomicRanges)
library(gdsfmt)
library(HDF5Array)

get_full_mat = function(gdspath) {
  ref = openfn.gds(gdspath)
  on.exit(closefn.gds(ref))
  ans = read.gdsn(index.gdsn(ref, "matrix"))
  rn = read.gdsn(index.gdsn(ref, "row.names"))
  cn = read.gdsn(index.gdsn(ref, "col.names"))
  dimnames(ans) = list(rn, cn)
  ans
}

shard_to_hdf5 = function(fullmat, coldata, chrname="chr22",
    cpginfo = epic_granges, prefix="", overwrite=FALSE, assayname="betas") {
  stopifnot(is(cpginfo, "GenomicRanges"))
  stopifnot(is(coldata, "DataFrame"))
  stopifnot(all(rownames(coldata)==colnames(fullmat)))
  fn = paste0(prefix, chrname, "_hdf5")
  if (dir.exists(fn) && overwrite==FALSE) stop(paste("won't overwrite", fn))
  cur = cpginfo[which(seqnames(cpginfo)==chrname)]
  di = setdiff( names(cur), rownames(fullmat))
  if (length(di)>0) message(sprintf("%d unmeasured EPIC CpG", length(di)))
  okn = intersect(rownames(fullmat), names(cur))
  curmat = fullmat[okn,]
  se = SummarizedExperiment(curmat)
  assayNames(se) = assayname
  colData(se) = coldata
  rowRanges(se) = cur[okn]
  se = se[order(start(rowRanges(se))),]
  saveHDF5SummarizedExperiment(se, fn, replace=overwrite)
}


fm = get_full_mat("../betas.gds")
ssheet = meffil.read.samplesheet(base="../../data-epic-demo", pattern="Demo_SampleSheet.csv")
cd = DataFrame(ssheet)
#shard_to_hdf5(fm, cd, chrname="chr22", overwrite=TRUE)
mychr = paste("chr", 1:22, sep="")
lapply(mychr, function(x) {shard_to_hdf5(fm, cd, chrname=x, overwrite=TRUE); cat(x)})

