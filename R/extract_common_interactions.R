#' extract common genes and miRNAs from all the tool prediction matrices and regression matrix
#' @export
#' @import rlist
#' @param reg_mat A matrix containg the regression coefficients
#' @param tools A list containing prediction matrices/dataframes from different tools.Each element in the list is a matrix/dataframe
#' @param tool_names A vector of the names of the tools in the same sequence as tools.

common_data <- function(reg_mat,tools,tool_names){

                #append the regression matrix to the list of tools
                common <- list.append(tools,reg_mat)

                names(common) <- c(tool_names,"regression")

                #extract common rows and columns from all the matrices
                row_ids <- Reduce(intersect,lapply(common, rownames))
                col_ids <- Reduce(intersect,lapply(common, colnames))

                if(is.null(row_ids)){
                  stop("There are no common miRNAs between the matrices")
                }
                else if(is.null(col_ids)){
                  stop("There are no common genes between the matrices")
                }
                else{

                    tool_names <- c(tool_names,"regression")

                    #keep only the common columns and rows for all the matrices
                    for (x in seq(1,length(common),1)) {

                      mat <- common[[tool_names[x]]]
                      mat <- mat[row_ids,col_ids]
                      common[[tool_names[x]]] <- mat
                    }
                    names(common) <- tool_names
              }
                return(common)
}
