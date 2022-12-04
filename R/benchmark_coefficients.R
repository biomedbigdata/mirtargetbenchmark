#' @title Compare with Regresssion
#' @description
#' Benchmark tools by comparing regression against tools
#' @import rlist
#' @param reg_mat A matrix containing the regression coefficients in the form of a matrix/dataframe.
#' @param tool_mat_list A list containing dataframes from different tools.
#' @param tool_names A list of the names of tools in the same order as tool_mat_list
#' @export
regression_vs_tool <- function(reg_mat,tool_mat_list,tool_names){

                      names(tool_mat_list) <- tool_names

                      #convert the regression matrix into a vector
                      reg_coeff_vector <- as.vector(unlist(reg_mat))

                      vector_tools <- list()

                      #convert the dataframes in the list into vectors and store in vector_tools list
                      for (x in seq(1,length(tool_mat_list),1)) {
                        vector_tools <- list.append(vector_tools,as.vector(unlist(tool_mat_list[[tool_names[x]]])))
                      }

                      names(vector_tools) <- tool_names


                      analysis <- list()

                      min_coeff <- min(reg_coeff_vector)

                      for (tool in seq(1,length(vector_tools),1)) {


                        vector_tools1 <- vector_tools[[tool_names[tool]]]
                        vector_results <- reg_coeff_vector


                        something <- foreach (x = seq(min_coeff,0,-min_coeff/10), .combine =rbind,.inorder = TRUE) %do% {

                          vector_results1 <- vector_results



                          #threshold the regression coefficients based on the iterator x
                          vector_results1[which(vector_results1 > x) ] <- 0

                          vector_results1[which(vector_results1 != 0)] <- 1

                          #calculate various metrics to perform analysis
                          class_1 <- length(which(vector_results1 == 1))
                          class_0 <- length(which(vector_results1 == 0))

                          tp <- length(intersect(which(vector_results1 == 1), which(vector_tools1 == 1)))
                          tn <- length(intersect(which(vector_results1 == 0), which(vector_tools1 == 0)))
                          fp <- length(intersect(which(vector_results1 == 0), which(vector_tools1 == 1)))
                          fn <- length(intersect(which(vector_results1 == 1), which(vector_tools1 == 0)))
                          sens <- tp/(tp+fn)
                          spec <- tn/(tn + fp)
                          ppv <- tp/(tp+fp)
                          data.frame(ones = class_1,zeros = class_0,threshold = x, TP = tp, TN = tn, FP = fp, FN = fn,SENS = sens, FPR = 1 - spec,PPV = ppv)
                        }


                        analysis <- list.append(analysis,something)

                      }

                      names(analysis) <- tool_names
                      return(analysis)
}
