#' Benchmark various tools against regression using the scores.
#' @export
#' @import ROCR
#' @import pracma
#' @import foreach
#' @param reg_mat A matrix containing the regression coefficients in the form of a matrix/dataframe.
#' @param vector_tools A list containing dataframes from different tools.
#' @param threshold A threshold coefficient value for coefficients.All the coefficients less than the threshold
#' will be changed to 1 and others will be 0. threshold = 0 by default.
tool_vs_regression <- function(reg_mat, tool_mat_list, threshold =  0){
                    vector_tools <- list()
                    for (x in seq(1,length(tool_mat_list),1)) {
                      vector_tools <- list.append(vector_tools, as.vector(unlist(tool_mat_list[x][[1]])))
                    }
                    reg_coeff_vector <- as.vector(unlist(reg_mat))
                    reg_coeff_vector[reg_coeff_vector > threshold ] <- 0
                    reg_coeff_vector[reg_coeff_vector < threshold ] <- 1

                    predictions <- c()
                    performances <- c()
                    all_tool_analysis <- list()
                    roc <- c()
                    pr <- c()

                    for (elem in seq(1, length(vector_tools),1)) {
                      predictions <- c(predictions,prediction(vector_tools[elem][[1]], reg_coeff_vector))
                    }

                    for (elem in seq(1, length(predictions),1)) {
                      performances <- c(performances,performance(predictions[elem][[1]],"prec","rec"))
                      roc <- cbind(roc,performance(predictions[elem][[1]], measure = "auc")@y.values[[1]])
                    }
                    for (elem in seq(1, length(performances),1)) {
                      if(elem == 1){
                        plot(performances[elem][[1]],colorize = TRUE)
                      }
                      else{
                        plot(performances[elem][[1]], add = TRUE, colorize = TRUE)
                      }
                      x <- performances[elem][[1]]@x.values[[1]]
                      y <- performances[elem][[1]]@y.values[[1]]

                      pr <- cbind(pr, trapz(x[-1],y[-1]))
                    }

                    AUC <- rbind(roc,pr)
                    rownames(AUC) <- c("roc", "pr")



                    for (elem in seq(1, length(vector_tools),1)) {

                      quantiles <- tapply(vector_tools[elem][[1]], cut(vector_tools[elem][[1]], 100), max)

                      analysis <- foreach(x = seq(0,length(quantiles),1),.combine = rbind,.inorder = TRUE)%do% {

                        vector_tools1 <- vector_tools[elem][[1]]


                        vector_tools1[vector_tools1 >  quantiles[x]] <- 1

                        vector_tools1[vector_tools1 < quantiles[x]] <- 0

                        tp <- length(intersect(which(vector_tools1 == 1),which(reg_coeff_vector == 1)))
                        tn <- length(intersect(which(vector_tools1 == 0),which(reg_coeff_vector == 0)))
                        fp <- length(intersect(which(vector_tools1 == 1),which(reg_coeff_vector == 0)))
                        fn <- length(intersect(which(vector_tools1 == 0),which(reg_coeff_vector == 1)))
                        sens <- tp/(tp+fn)
                        spec <- tn/(tn + fp)
                        ppv <- tp/(tp+fp)
                        data.frame(threshold = x,TP = tp, TN = tn, FP = fp, FN = fn,SENS = sens, FPR = 1 - spec,PPV = ppv)
                      }
                      plot(analysis$SENS, analysis$PPV, type = "l", lty = 1, col = "blue"
                          ,xlim = c(0,0.5), ylim = c(0.1,0.6), main = "Precision-Recall Curve")
                      par(new=TRUE)
                      all_tool_analysis <- list.append(all_tool_analysis,analysis)
                   }


                  return(list(AUC,all_tool_analysis))
}
