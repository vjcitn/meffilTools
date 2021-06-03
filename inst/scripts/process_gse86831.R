path="."

samplesheet <- meffil.read.samplesheet(base=path, pattern="newss.csv")


qc.file <- "./ncbi.qc-report.html"
author <- "Illumina, et al."
study <- "EPIC demo dataset"
number.pcs <- 2
norm.file <- "./ncbi.normalization-report.html"
#Generate QC objects.

qc.objects <- meffil.qc(samplesheet, cell.type.reference="blood gse35069 complete", verbose=T)
#QC report.

qc.summary <- meffil.qc.summary(qc.objects, verbose=T)

norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs=number.pcs, verbose=T)

meffil.normalize.samples(norm.objects, gds.filename="bigdemo.gds", verbose=TRUE)
