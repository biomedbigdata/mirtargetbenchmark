#' Benchmark tools by comparing regression against tools
#' @export
#' @param reg_mat A matrix containing the regression coefficients in the form of a matrix/dataframe.
#' @param vector_tools A list containing dataframes from different tools.
#' @param val_mat A validation dataset containing the number of binding sites to validate the miRNAs.
#' @param threshold_valmir A threshold value to filter the miRNAs.miRNAs with less than threshold gene targets will be removed.
#' The default is 100.
regression_vs_tool <- function(reg_mat,tool_mat_list, val_mat, threshold_valmir = 100){
                      print(threshold_valmir)
                      keep <- c()
                      for (x in seq(1,nrow(val_mat),1)) {
                        if(length(which(val_mat[x,] > 0 )) >= threshold_valmir){
                          keep <- c(keep,rownames(val_mat[x,]))
                        }
                      }

                      keep <- intersect(rownames(tool_mat_list[1][[1]]),rownames(val_mat))

                      for (x in seq(1, length(tool_mat_list),1)) {
                        tool_mat_list[x][[1]] <- tool_mat_list[x][[1]][keep,]
                      }

                      reg_mat <- reg_mat[keep,]

                      reg_coeff_vector <- as.vector(unlist(reg_mat))
                      vector_tools <- list()
                      for (x in seq(1,length(tool_mat_list),1)) {
                        vector_tools <- list.append(vector_tools,as.vector(unlist(tool_mat_list[x][[1]])))
                      }

                      precision <- list()

                      for (tool in seq(1,length(vector_tools),1)) {

                        something <- foreach (x = seq(-1,0,0.01), .combine =rbind) %do% {
                          vector_results1 <- reg_coeff_vector
                          vector_tools1 <- vector_tools[tool][[1]]

                          vector_results1[vector_results1 > x ] <- 0
                          vector_results1[vector_results1 < x] <- 1

                          class_1 <- length(which(vector_results1 == 1))
                          class_0 <- length(which(vector_results1 == 0))

                          tp <- length(intersect(which(vector_results1 == 1), which(vector_tools1 == 1)))
                          tn <- length(intersect(which(vector_results1 == 0),which(vector_tools1 == 0)))
                          fp <- length(intersect(which(vector_results1 == 0), which(vector_tools1 == 1)))
                          fn <- length(intersect(which(vector_results1 == 1), which(vector_tools1 == 0)))
                          sens <- tp/(tp+fn)
                          spec <- tn/(tn + fp)
                          ppv <- tp/(tp+fp)
                          data.frame(ones = class_1,zeros = class_0,threshold = x, TP = tp, TN = tn, FP = fp, FN = fn,SENS = sens, FPR = 1 - spec,PPV = ppv)
                        }



                        precision <- list.append(precision,something$PPV)

                      }
                      return(precision)
}
