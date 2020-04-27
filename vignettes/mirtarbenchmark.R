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

## ----results='hide', message=FALSE--------------------------------------------
gene_expr <- expression_preprocessing(gene_expr,missing_value = -9.9658,missing_value_filter = 90 ,
                                      variance_filter = 0.2, mean_filter = 0.2)

## ----echo = FALSE-------------------------------------------------------------
sprintf("The number of genes left after the filtering are: %s",ncol(gene_expr))


## ---- results='hide',message=FALSE--------------------------------------------
mir_expr <- expression_preprocessing(mir_expr,missing_value = 0,missing_value_filter = 90 ,
                                      variance_filter = 0.2, mean_filter = 0.2)

## ---- echo=FALSE--------------------------------------------------------------
sprintf("The number of miRNAs left after the filtering are: %s",ncol(mir_expr))

## ----results = 'hide', warning=FALSE, message=FALSE---------------------------
results <- regression_results(dependent_expr = gene_expr, independent_expr = mir_expr, 
                              model_choice = 's', rsq_filter = 0.2)

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(results[1:5,1:5])

## ---- warning=FALSE, message=FALSE--------------------------------------------
targetscan <- convert_tool_data_into_matrix(pred_data = targetscan_pred[1:1000,],miRNA_colname = "miRNA",
                                            gene_colname = "Gene.ID",score_colname = "context...score")

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(targetscan[1:5,1:5])

## ----results = 'hide', warning=FALSE, message=FALSE---------------------------

converted_regression_matrix <- gene_miRNA_id_conversion(mat_data = results,cgid = "ENSEMBL_VERSION",
                                                        mir_bool = FALSE) 


## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(converted_regression_matrix[2:6,1:5])

## ----warning=FALSE, message=FALSE---------------------------------------------

tool_list <- list(sample_targetscan, sample_mirwalk, sample_pita)

#extracting common genes and miRNAs
common_list <- common_data(sample_simple_regression_results
                           , tool_list, tool_names = c('targetscan','mirwalk','pita'))

#using miRNAs with more than 500 validated target genes in miRTarBase
common_validated_targets <- validated_miRNA_filter(common_list, matrix_names = c("targetscan","mirwalk","pita",
                                  "regression"),mirtarbase7, threshold_valmir = 100)

#list of tool matrices
common_list_scores <- list(common_validated_targets[["targetscan"]], common_validated_targets[["mirwalk"]],                                  common_validated_targets[["pita"]])

#name the elements (to be used in ensembl methods)
names(common_list_scores) <- c("targetscan", "mirwalk", "pita")

#results and analysis
output <- tool_vs_regression(common_validated_targets[["regression"]],common_list_scores, tool_names = c('targetscan','mirwalk','pita'))


## ---- echo=FALSE, results='asis'----------------------------------------------
print("The green, red and pink lines represent TargetScan , PITA and miRWalk respectively")

## ---- echo=FALSE, results='asis'----------------------------------------------

output$targetscan$threshold <- (output$targetscan$threshold - min(output$targetscan$threshold))/(max(output$targetscan$threshold) - min(output$targetscan$threshold))
output$mirwalk$threshold <- (output$mirwalk$threshold - min(output$mirwalk$threshold))/(max(output$mirwalk$threshold) - min(output$mirwalk$threshold))
output$pita$threshold <- (output$pita$threshold - min(output$pita$threshold))/(max(output$pita$threshold) - min(output$pita$threshold))

plot(output$targetscan$threshold, output$targetscan$PPV, type = 'l',ylim = c(0,1), xlim = c(0,1), col = "green", xlab = "thresholds", ylab = "precision")
lines(output$pita$threshold, output$pita$PPV, type = 'l', col = "red")
lines(output$mirwalk$threshold, output$mirwalk$PPV, type = 'l', col = "pink")

plot(output$targetscan$threshold, output$targetscan$SENS, type = 'l',ylim = c(0,1), xlim = c(0,1), col = "green", xlab = "thresholds", ylab = "recall")
lines(output$pita$threshold, output$pita$SENS, type = 'l', col = "red")
lines(output$mirwalk$threshold, output$mirwalk$SENS, type = 'l', col = "pink")


## ----warning=FALSE, message=FALSE---------------------------------------------

#covert binding sites into binary
mircode[mircode > 1] <- 1
mirtarbase[mirtarbase > 1] <- 1
mirtarbase7[mirtarbase7 > 1] <- 1

tool_list <- list(mircode, mirtarbase, mirtarbase7)

#extracting common genes and miRNAs
common_list <- common_data(sample_simple_regression_results
                           , tool_list, tool_names = c('mircode','mirtarbase','mirtarbase7'))

#list of tool matrices
common_list_1 <- list(common_list[["mircode"]], common_list[["mirtarbase"]], common_list[["mirtarbase7"]])

#results and analysis
output_1 <- regression_vs_tool(common_list[["regression"]],common_list_1, tool_names = c('mircode','mirtarbase','mirtarbase7'))


## ---- echo=FALSE, results='asis'----------------------------------------------

   print("The metrics shown below were obtained using a threshold of zero for coefficients")
   precision <- cbind(output_1$mircode$PPV[11], output_1$mirtarbase$PPV[11], output_1$mirtarbase7$PPV[11])
   recall <- cbind(output_1$mircode$SENS[11], output_1$mirtarbase$SENS[11], output_1$mirtarbase7$SENS[11])
   metrics <- rbind(precision, recall)
   rownames(metrics) <- c("precision", "recall")
   colnames(metrics) <- c("mircode", "mirtarbase", "mirtarbase7")
  knitr::kable(metrics)
  

## ----warning=FALSE, message=FALSE---------------------------------------------

ppv_performance <- c(output$targetscan$PPV[1], output$mirwalk$PPV[1], output$pita$PPV[1])

names(ppv_performance) <- c("targetscan", "mirwalk", "pita")

ppv_performance <- sort(ppv_performance, decreasing = TRUE)

new_common_list <- list()

library(rlist)

for (x in 1:length(ppv_performance)) {
  new_common_list <- list.append(new_common_list, common_list_scores[[names(ppv_performance)[x]]])
}

ensemble_predictions <- ensemble_scores(new_common_list, tool_names = names(ppv_performance), ppv_performance)


## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(head(ensemble_predictions[,1:5]))

## ----warning=FALSE, message=FALSE---------------------------------------------

#compare the performance with other tools
common_list_scores <- list.append(common_list_scores, ensemble_predictions)

output_with_ensembl <- tool_vs_regression(common_validated_targets[["regression"]],
                                          common_list_scores, tool_names = c(c('targetscan','mirwalk','pita',
                                                                               'ensembl')))




## ---- echo=FALSE, results='asis'----------------------------------------------
print("The green, red, pink and black lines represent TargetScan , PITA, miRWalk and ensemble respectively")

## ---- echo=FALSE, results='asis'----------------------------------------------

#scale all the scores to have range of 0-1
output_with_ensembl$targetscan$threshold <- (output_with_ensembl$targetscan$threshold - min(output_with_ensembl$targetscan$threshold))/(max(output_with_ensembl$targetscan$threshold) - min(output_with_ensembl$targetscan$threshold))
output_with_ensembl$mirwalk$threshold <- (output_with_ensembl$mirwalk$threshold - min(output_with_ensembl$mirwalk$threshold))/(max(output_with_ensembl$mirwalk$threshold) - min(output_with_ensembl$mirwalk$threshold))
output_with_ensembl$pita$threshold <- (output_with_ensembl$pita$threshold - min(output_with_ensembl$pita$threshold))/(max(output_with_ensembl$pita$threshold) - min(output_with_ensembl$pita$threshold))


plot(output_with_ensembl$targetscan$threshold, output_with_ensembl$targetscan$PPV, type = 'l',ylim = c(0,1), xlim = c(0,1), col = "green", xlab = "thresholds", ylab = "precision")
lines(output_with_ensembl$pita$threshold, output_with_ensembl$pita$PPV, type = 'l', col = "red")
lines(output_with_ensembl$mirwalk$threshold, output_with_ensembl$mirwalk$PPV, type = 'l', col = "pink")
lines(output_with_ensembl$ensembl$threshold, output_with_ensembl$ensembl$PPV, type = 'l', col = "black")


plot(output_with_ensembl$targetscan$threshold, output_with_ensembl$targetscan$SENS, type = 'l',ylim = c(0,1), xlim = c(0,1), col = "green", xlab = "thresholds", ylab = "recall")
lines(output_with_ensembl$pita$threshold, output_with_ensembl$pita$SENS, type = 'l', col = "red")
lines(output_with_ensembl$mirwalk$threshold, output_with_ensembl$mirwalk$SENS, type = 'l', col = "pink")
lines(output_with_ensembl$ensembl$threshold, output_with_ensembl$ensembl$SENS, type = 'l', col = "black")

