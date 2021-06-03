find_anticorr = function(se, inuse) {
 stopifnot(length(inuse)==1)
 if (nrow(se)==1) return(inuse)
 if (nrow(se)==2) return(setdiff(rownames(se), inuse))
 ct = cor(t(as.matrix(assay(se))))
 att = try(which.min(ct[inuse,]))
 if (!inherits(att, "try-error")) return(names(att))
 att
}
library(parallel)
## 0/0 packages newly attached/loaded, see sessionInfo() for details.
options(mc.cores=6)
addcpg = mclapply(seq_len(length(gg)),
   function(x) try(find_anticorr(gg[[x]], topp[x])))
augmented = unique(c(topp, unlist(addcpg)))
newse = pidsSE_prom[augmented,]
remm2 = rowMads(assay(newse))
orem = order(remm2, decreasing=TRUE)
top2 = newse[ orem[seq_len(90)], ]
set.seed(1234)
rownames(top2) = rowRanges(top2)$symbol
o2p = rept_filtse(top2)
summary(o2p)
