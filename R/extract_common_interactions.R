#' extract common genes and miRNAs from all the tool prediction matrices and regression matrix
#' @export
#' @import rlist
#' @param reg_mat A matrix containg the regression coefficients
#' @param tools A list containing prediction matrices/dataframes from different tools.Each element in the list is a matrix/dataframe
#' @param tool_names A vector of the names of the tools in the same sequence as tools.
#' @param filter_bool
#' @param filter_threshold
#' @param filter Genes/miRNAs with more than filter_threshold values equal to filter will be removed for all the tools. NA by default.

common_data <- function(reg_mat,tools,tool_names,filter_bool = FALSE,filter_value = NA, filter_threshold = 0.95 ){

              #append the regression matrix to the list of tools
              common <- list.append(tools,reg_mat)

              names(common) <- c(tool_names,"regression")

              #extract common rows and columns from all the matrices
              row_ids <- Reduce(intersect,lapply(common, rownames))
              col_ids <- Reduce(intersect,lapply(common, colnames))

              if(filter_bool){
                #remove the rows and columns with very few values
                for (x in seq(1,length(common)-1,1)) {
                  mat <- common[[tool_names[x]]]
                  mat <- mat[row_ids,col_ids]

                  if(is.na(filter_value)){
                    mat <- mat[-which(rowMeans(is.na(mat)) > filter_threshold),]
                    mat <- mat[,-which(colMeans(is.na(mat)) > filter_threshold)]
                  }
                  else{
                    mat <- mat[-which(rowMeans(mat == filter_value) > filter_threshold),]
                    mat <- mat[,-which(colMeans(mat == filter_value) > filter_threshold)]
                  }
                  common[[tool_names[x]]] <- mat
                }
              }
              tool_names <- c(tool_names,"regression")

                #extract the common interactions again
                row_ids <- Reduce(intersect,lapply(common, rownames))
                col_ids <- Reduce(intersect,lapply(common, colnames))

                #keep only the common columns and rows for all the matrices
                for (x in seq(1,length(common),1)) {

                  mat <- common[[tool_names[x]]]
                  mat <- mat[row_ids,col_ids]
                  common[[tool_names[x]]] <- mat
                }
                names(common) <- tool_names
                return(common)
}
