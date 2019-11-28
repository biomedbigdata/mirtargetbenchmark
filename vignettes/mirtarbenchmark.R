## ---- warning=FALSE, message=FALSE--------------------------------------------
library(mirtargetbenchmark)

## ---- eval=FALSE--------------------------------------------------------------
#  head(gene_expr)

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(gene_expr[1:5,1:5])

## ---- eval=FALSE--------------------------------------------------------------
#  head(mir_expr)

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(mir_expr[1:5,1:5])

## ----echo = FALSE-------------------------------------------------------------
sprintf("The number of genes left after the filtering are: %s",ncol(gene_expr))


## ---- echo=FALSE--------------------------------------------------------------
sprintf("The number of miRNAs left after the filtering are: %s",ncol(mir_expr))

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(results[1:5,1:5])

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(targetscan[1:5,1:5])

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(converted_regression_matrix[2:6,1:5])

## ----warning=FALSE,cache=TRUE, message=FALSE----------------------------------

tool_list <- list(sample_targetscan, sample_mirwalk, sample_pita)
common_list <- common_data(gene_miRNA_id_conversion(results,cgid = 'ensembl_gene_id_version', mir_bool = FALSE)
                           , tool_list, tool_names = c('targetscan','mirwalk','pita'))

