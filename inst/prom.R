
library(meffilTools)
requireNamespace("HDF5Array")
library(SummarizedExperiment)
se6 = HDF5Array::loadHDF5SummarizedExperiment(system.file("chr6_hdf5", package="meffilTools"))
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
pp = promoters(TxDb.Hsapiens.UCSC.hg19.knownGene)
g19 = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
syms = mapIds(org.Hs.eg.db, keys=g19$gene_id, keytype="ENTREZID", column="SYMBOL") 
g19$symbol = syms
p19 = promoters(g19)
pom = findOverlaps(p19, rowRanges(se6))
