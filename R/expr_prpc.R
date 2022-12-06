#' @title preprocessing expression data
#' @description
#' Filter the expression data
#' @importFrom stats var
#' @importFrom stats quantile
#' @importFrom graphics hist
#' @importFrom magrittr %>%
#' @param expr_data The input expression matrix containing expression values
#' @param missing_value_filter Genes/miRNAs with more than missing_value\% will be removed.
#' @param missing_value The expression value which represents the missing value in your expression matrix.
#' @param variance_filter Genes/miRNAs with variance less than variance_filter quantile will be removed.Default is 0.
#' @param mean_filter Genes/miRNAs with mean less than mean_filter quantile will be removed.Default is 0.
#' @export
expression_preprocessing <- function(expr_data,missing_value_filter = 0,missing_value,
                                     variance_filter = 0,mean_filter = 0){

                              hist(as.numeric(unlist(as.vector(expr_data), use.names = FALSE)), border = "blue", xlab = "initial expression values",
                                   ylab = "initial frequency", main = "initial histogram of expression data")

                              #missing value filter
                              if(missing_value_filter){
                              res <- colSums(expr_data == missing_value)/nrow(expr_data)*100
                              # if there are more than mentioned number of missing values
                              if(length(which(res>missing_value_filter)) > 0){
                                expr_data <- expr_data[,-which(res>missing_value_filter)]
                              }
                              }

                              #variance filter
                              if(variance_filter){

                                l <- c()
                                expr_data <- tryCatch({
                                  for(x in seq(1,ncol(expr_data), 1)){
                                    l <- c(l, var(expr_data[,x]))
                                  }
                                  u <- which(l < quantile(l,variance_filter))

                                  expr_data <- expr_data[,-u]
                                  expr_data
                                  }
                                  ,
                                  error = function(e){
                                      stop("The variance_filter value is very high, please try a lower value ")
                                  }

                                  )
                              }

                              if(mean_filter){
                                l <- c()
                                expr_data <- tryCatch({
                                  for(x in seq(1,ncol(expr_data), 1)){
                                    l <- c(l, mean(expr_data[,x]))
                                  }

                                  u <- which(l < quantile(l,mean_filter))

                                  expr_data <- expr_data[,-u]
                                  expr_data
                                },
                                  error = function(e){
                                    stop("The mean_filter value is very high, please try a lower value ")
                                  }
                                )
                              }

                              hist(as.numeric(unlist(as.vector(expr_data), use.names = FALSE)), border = "blue", xlab = "final expression values",
                                    ylab = "final frequency", main = "final histogram of expression data")

                              return(as.data.frame(expr_data))
}


