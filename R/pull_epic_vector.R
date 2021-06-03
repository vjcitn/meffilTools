#' grab a vector of betas assuming epic dimensions
#' @import gdsfmt
#' @param colid character(1) or NULL column id to select, assumes 'col.names' present
#' @param cpg character(1) or NULL CpG to select, assumes 'row.names' present
#' @param base character(1) folder in which .gds file is stored
#' @param gds character(1) local name of .gds file
#' @param ncpg integer(1) number of cpgs, the 'row' extent of 'matrix', assumed present in gds
#' @param nsamps integer(1) number of columns, the 'col' extent of 'matrix'
#' @examples
#' \dontrun{
#' x = pull_epic_vector()
#' stopifnot(nrow(x) == 865859)
#' head(x)
#' y = pull_epic_vector(cpg="cg18478105", colid=NULL)
#' head(y)
#' }
#' @note This function uses openfn.gds and on.exit(closefn.gds(...)), extracts a
#' matrix and assigns dimnames.  We might want to make dimnames assignment optional.
#" At most one of col
#' @export
pull_epic_vector = function(colid = "TOE440084_204234730121_R06C01_379130", cpg=NULL,
    base="/udd/resiq/meffil/QC/COPD.blood",
   gds = "beta.copd.blood.gds", ncpg=865859L, nsamps=12147L) {
  chk = c(colid, cpg)
  stopifnot(length(chk)==1)
  bb = openfn.gds(paste(base, gds, sep="/"))
  on.exit(closefn.gds(bb))
  cpgnames <- index.gdsn(bb, "row.names")
  allcpgs = read.gdsn(cpgnames, 1, ncpg)
  rm(cpgnames)
  if (!is.null(cpg)) stopifnot(cpg %in% allcpgs)
  sampleIDS <- index.gdsn(bb, "col.names")
  allids = read.gdsn(sampleIDS, 1, nsamps)
  stopifnot(colid %in% allids)
  rm(sampleIDS)  # avoid confusion
  BETAS <- index.gdsn(bb, "matrix")
  if (!is.null(colid)) {
   sampind = match(colid, allids)
   ans = read.gdsn(BETAS, start=c(1, sampind), count=c(length(allcpgs), 1), simplify="none")
   dimnames(ans) = list(allcpgs, colid)
   }
  else {
   cpgind = match(cpg, allcpgs)
   ans = read.gdsn(BETAS, start=c(cpgind, 1), count=c(1, length(allids)), simplify="none")
   dimnames(ans) = list(cpg, allids)
   }
  ans
}

#' build a forkable meffil GDS reference
#' @note hard codes paths for CDNM project
#' @export
getmat = function() {
  ncpg=865859L 
  nsamps=12147L
  base="/udd/resiq/meffil/QC/COPD.blood"
  gds = "beta.copd.blood.gds"
  bb = openfn.gds(paste(base, gds, sep="/"), allow.fork=TRUE)
  cpgnames <- index.gdsn(bb, "row.names")
  allcpgs = read.gdsn(cpgnames, 1, ncpg)
  rm(cpgnames)
  sampleIDS <- index.gdsn(bb, "col.names")
  allids = read.gdsn(sampleIDS, 1, nsamps)
  BETAS <- index.gdsn(bb, "matrix")
  ans = list(obj=bb, BETAS=BETAS, cpgs=allcpgs, sampids=allids)
  class(ans) = "forkable.meffilref"
  ans
}

#' manage a meffilref
#' @import gdsfmt
#' @param con forkable.meffilref instance
#' @param \dots not used
#' @export
close.forkable.meffilref = function(con, ...)
  closefn.gds(con$obj)

#' print meth
#' @param x instance of forkable.meffilref
#' @param \dots not used
#' @export
print.forkable.meffilref = function(x, ...) {
 cat("forkable.meffilref instance, wraps\n")
 print(x$obj)
}


#' get a named vector of betas from meffil GDS instance
#' @param x instance of S3 class forkable.meffilref
#' @param cpg character(1) feature ID
#' @param ncpu integer(1) number of CPUs to use to pick through the GDS row data
#' @note Row access is slow, so we use parallel ingestion of chunks.
#' @examples
#' \dontrun{
#' tst1 = getmat()
#' tst1
#' head(get_cpg(tst1))
#' close(tst1)
#' }
#' @export
get_cpg = function(x, cpg="cg18478105", ncpu=1) {
 stopifnot(inherits(x, "forkable.meffilref"))
 desc = objdesp.gdsn(x$BETAS)
 dims = desc$dim
 nsamp = dims[2]
 ch = BBmisc::chunk(seq_len(nsamp), n.chunks=ncpu)
 lens = vapply(ch, length, integer(1))
 starts = c(1L, cumsum(lens))[1:ncpu]
 options(mc.cores=ncpu)
 cind = match(cpg, x$cpgs)
 ans = as.numeric(unlist(mclapply(seq_len(length(ch)), 
    function(z) read.gdsn(x$BETAS, start=c(cind, starts[z]), count=c(1, lens[z])))))
 names(ans) = x$sampids
 ans
}

