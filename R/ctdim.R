#' @title dataframe to matrix
#' @description
#' Convert tool prediction dataframe into a matrix with rows as miRNAs, columns as genes and scores as values
#' @importFrom dplyr filter
#' @importFrom graphics hist
#' @import foreach
#' @import parallel
#' @import doSNOW
#' @importFrom rlang UQ
#' @param pred_data The predictions in the form of a data frame
#' @param miRNA_colname The column name of miRNA in the dataframe
#' @param gene_colname The column name of gene in the dataframe
#' @param score_colname The column name of scores in your dataframe
#' @export
convert_tool_data_into_matrix <- function(pred_data, miRNA_colname, gene_colname, score_colname){
                                 pred_data <- pred_data[,c(miRNA_colname,gene_colname,score_colname)]

                                 #remove the NAs from all the columns
                                 list_indices <- 1:nrow(pred_data)
                                 if (length(which(is.na(pred_data[miRNA_colname])))>0){
                                      list_indices <- list_indices[-which(is.na(pred_data[miRNA_colname]))]
                                 }
                                 if (length(which(is.na(pred_data[miRNA_colname])))>0){
                                   list_indices <- list_indices[-which(is.na(pred_data[gene_colname]))]
                                 }
                                 if (length(which(is.na(pred_data[miRNA_colname])))>0){
                                   list_indices <- list_indices[-which(is.na(pred_data[score_colname]))]
                                 }
                                 pred_data <- pred_data[list_indices,]

                                 # create vectors for gene names and miRNA names
                                 genes <- unique(pred_data[,gene_colname])
                                 miRNAs <- unique(pred_data[,miRNA_colname])

                                 tool_mat <- data.frame(matrix(0, ncol = length(genes), nrow = length(miRNAs)))
                                 colnames(tool_mat) <- genes
                                 rownames(tool_mat) <- miRNAs

                                 if(length(genes) < 1){
                                   stop("The dataframe is empty or all values are NA")
                                 }
                                 x <- NULL
                                 tool_pred <- foreach(x = seq(1, length(genes),1), .combine = cbind, .inorder = TRUE,
                                                         .packages = c('doSNOW', 'dplyr'))%dopar% {

                                                           col <- genes[x]


                                                           vector <- tryCatch({
                                                             a <- pred_data%>%filter(UQ(as.symbol(colnames(pred_data)[2])) == col)
                                                             row_indices <- match(rownames(tool_mat), a[,1])
                                                             vec <- numeric(length(miRNAs))
                                                             vec <- a[,3][row_indices]

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
                                 # convert NA's to 0's
                                 tool_pred[is.na(tool_pred)] <- 0
                                 return(tool_pred)
}
