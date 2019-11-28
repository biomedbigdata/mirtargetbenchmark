#' Filter the tool and regression matrixes using experimentally validated data
#' @export
#' @param matrix_list A list containing all the tool data and regression matrix
#' @param matrix_names A list of names of all the data sets in the list
#' @param val_mat A matrix containg experimentally validated gene-miRNA interactions.
#' @param threshold_valmir A threshold for filtering miRNAs. Only the miRNAs which have more than threshold_valmir
#' gene targets will be kept. It is 100 by default.

validated_miRNA_filter <- function(matrix_list, matrix_names, val_mat, threshold_valmir = 100){
                          keep <- c()

                          for (x in seq(1,nrow(val_mat),1)) {

                            if(length(which(val_mat[x,] > 0 )) >= threshold_valmir){
                              keep <- c(keep,rownames(val_mat)[x])
                            }
                          }

                          print(length(keep))

                          val_mat <- val_mat[keep,]

                          keep1 <- intersect(rownames(matrix_list[[matrix_names[1]]]),rownames(val_mat))



                          for (x in seq(1, length(matrix_list),1)) {
                            matrix_list[[matrix_names[x]]] <- matrix_list[[matrix_names[[x]]]][keep1,]
                          }
                          print(rownames(matrix_list[[matrix_names[1]]]))
                          names(matrix_list) <- matrix_names
                          return(matrix_list)
}
