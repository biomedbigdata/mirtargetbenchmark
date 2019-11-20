#' extract common genes and miRNAs from all the tool prediction matrices and regression matrix
#' @export
#' @import rlist
#' @param reg_mat A matrix containg the regression coefficients
#' @param tools A list containing prediction matrices/dataframes from different tools.Each element in the list is a matrix/dataframe
#' @param filter Genes/miRNAs with more than 95\% values equal to filter will be removed for all the tools. NA by default.

common_data <- function(reg_mat,tools,filter = NA){


              common <- list.append(tools,reg_mat)
              row_ids <- Reduce(intersect,lapply(common, rownames))
              col_ids <- Reduce(intersect,lapply(common, colnames))

              for (x in seq(1,length(common)-1,1)) {
                mat <- common[x][[1]]
                mat <- mat[row_ids,col_ids]
                if(is.na(filter)){
                  mat <- mat[-which(rowMeans(is.na(mat)) > 0.95),-which(colMeans(is.na(mat)) > 0.95)]
                }
                else{
                  mat <- mat[-which(rowMeans(mat == filter) > 0.95),-which(colMeans(mat == filter) > 0.95)]
                }
                common[x][[1]] <- mat
              }
                row_ids <- Reduce(intersect,lapply(common, rownames))
                col_ids <- Reduce(intersect,lapply(common, colnames))

                for (x in seq(1,length(common),1)) {
                  mat <- common[x][[1]]
                  mat <- mat[row_ids,col_ids]
                  common[x][[1]] <- mat
                }
                return(common)
}
