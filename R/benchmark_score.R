#' @title Regression vs Scores
#' @description
#' Benchmark various tools against regression using the scores.
#' @import rlist
#' @import foreach
#' @param reg_mat A matrix containing the regression coefficients in the form of a matrix/dataframe.
#' @param tool_mat_list A list containing dataframes from different tools.
#' @param tool_names A list of the names of tools in the same order as tool_mat_list
#' @param threshold A threshold coefficient value for coefficients.All the coefficients less than the threshold
#' will be changed to 1 and others will be 0. threshold = 0 by default.
#' @export
tool_vs_regression <- function(reg_mat, tool_mat_list, tool_names, threshold =  0){

                    names(tool_mat_list) <- tool_names
                    vector_tools <- list()

                    #convert the dataframes in the list into vectors and store in vector_tools list
                    for (x in seq(1,length(tool_mat_list),1)) {
                      vector_tools <- list.append(vector_tools, as.vector(unlist(tool_mat_list[[tool_names[x]]])))
                  }

                    names(vector_tools) <- tool_names

                    #convert the regression matrix into a vector
                    reg_coeff_vector <- as.vector(unlist(reg_mat))

                    #convert the vector into binary vector using the threshold value
                    reg_coeff_vector[which(reg_coeff_vector > threshold) ] <- 0
                    reg_coeff_vector[which(reg_coeff_vector != 0) ] <- 1


                    all_tool_analysis = list()


                    for (elem in seq(1, length(vector_tools),1)) {

                      vector_tools1 <- vector_tools[[tool_names[elem]]]
                      reg_coeff_vector1 <- reg_coeff_vector


                      #keep non zero scores for tool vectors and keep the corresponding values in regression vector
                      keep_indices <- which(vector_tools1 != 0)
                      #keep_indices2 <- which(reg_coeff_vector1 != 0)
                      #keep_indices <- c(keep_indices1,keep_indices2)

                      vector_tools1 <- vector_tools1[keep_indices]
                      reg_coeff_vector1 <- reg_coeff_vector1[keep_indices]


                        analysis <- foreach(x = seq(0,1,0.01),.combine = rbind,.inorder = TRUE)%do% {

                        vector_tools2 <- vector_tools1

                        #calculate the x quantile
                        quantiles <- quantile(vector_tools2, x)

                        #threshold based on the quantiles
                        vector_tools2[which(vector_tools2 >  quantiles)] <- 1
                        vector_tools2[which(vector_tools2 != 1) ] <- 0

                        #calculate various metrics to perform the analysis
                        tp <- length(intersect(which(vector_tools2 == 1),which(reg_coeff_vector1 == 1)))
                        tn <- length(intersect(which(vector_tools2 == 0),which(reg_coeff_vector1 == 0)))
                        fp <- length(intersect(which(vector_tools2 == 1),which(reg_coeff_vector1 == 0)))
                        fn <- length(intersect(which(vector_tools2 == 0),which(reg_coeff_vector1 == 1)))
                        sens <- tp/(tp+fn)
                        spec <- tn/(tn + fp)
                        ppv <- tp/(tp+fp)
                        data.frame(threshold = quantiles,TP = tp, TN = tn, FP = fp, FN = fn,SENS = sens, FPR = 1 - spec,PPV = ppv)

                      }

                     all_tool_analysis <- list.append(all_tool_analysis,analysis)
                   }

                  names(all_tool_analysis) <- tool_names
                  return(all_tool_analysis)
}
