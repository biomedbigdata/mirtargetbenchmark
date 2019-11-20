#' Convert tool prediction dataframe into a matrix with rows as miRNAs, columns as genes and scores as values
#' @export
#' @import dplyr
#' @import foreach
#' @import parallel
#' @import doSNOW
#' @param pred_data The predictions in the form of a data frame
#' @param miRNA_index The column index of miRNA in the dataframe
#' @param gene_index The column index of gene in the dataframe
#' @param score_index The column index of scores in your dataframe
convert_tool_data_into_matrix <- function(pred_data, miRNA_index, gene_index, score_index){
                                 print(colnames(pred_data))
                                 pred_data <- pred_data[,c(miRNA_index,gene_index,score_index)]
                                 print("okay")
                                 genes <- unique(pred_data[,gene_index])
                                 miRNAs <- unique(pred_data[,miRNA_index])
                                 print("okay2")
                                 tool_mat <- data.frame(matrix(0, ncol = length(genes), nrow = length(miRNAs)))
                                 colnames(tool_mat) <- genes
                                 rownames(tool_mat) <- miRNAs
                                 print("okay3")
                                 tool_pred <- foreach(x = seq(1,length(genes),1), .combine = cbind, .inorder = TRUE,
                                                         .packages = c('doSNOW', 'dplyr'))%dopar% {

                                                           col <- genes[x]


                                                           vector <- tryCatch({
                                                             a <- pred_data%>%filter(UQ(as.symbol(colnames(pred_data)[2])) == col)
                                                             row_indices <- match(rownames(tool_mat), a[,1])
                                                             vec <- numeric(length(miRNAs))
                                                             vec <- a[,3][row_indices]
                                                             print("success")
                                                             vec
                                                           },
                                                           error = function(e)
                                                           {
                                                             print(e)
                                                             return(NA)
                                                           }

                                                           )


                                                         }
                                 colnames(tool_pred) <- genes
                                 rownames(tool_pred) <- miRNAs
                                 return(tool_pred)
}
