## ---- warning=FALSE, message=FALSE--------------------------------------------
library(mirtargetbenchmark)

## ---- eval=FALSE--------------------------------------------------------------
#  head(gene_expr)

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(gene_expr[1:5,1:8])

## ---- eval=FALSE--------------------------------------------------------------
#  head(mir_expr)

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(mir_expr[1:5,1:5])

