#' Filter the expression data
#' @export
#' @param expr_data The input expression matrix containing expression values
#' @param missing_value_filter Genes/miRNAs with more than missing_value\% will be removed.
#' @param missing_value The expression value which represents the missing value in your expression matrix.
#' @param variance_filter Genes/miRNAs with variance less than variance_filter quantile will be removed.Default is 0.
#' @param mean_filter Genes/miRNAs with mean less than mean_filter quantile will be removed.Default is 0.
expression_preprocessing <- function(expr_data,missing_value_filter = 0,missing_value,
                                     variance_filter = 0,mean_filter = 0){

                              plot(hist(expr_data), border = "blue", xlab = "initial expression values",
                                   ylab = "initial frequency")

                              #missing value filter
                              if(missing_value_filter){
                              res <- colSums(expr_data == missing_value)/nrow(expr_data)*100
                              expr_data <- expr_data[,-which(res>missing_value_filter)]
                              }
                              #variance filter
                              if(variance_filter){
                              l <- c()
                              for(x in seq(1,ncol(expr_data), 1)){
                                l <- c(l, var(expr_data[,x]))
                              }
                              u <- which(l < quantile(l,variance_filter))

                              expr_data <- expr_data[,-u]
                              }


                              if(mean_filter){
                                l <- c()
                                for(x in seq(1,ncol(expr_data), 1)){
                                  l <- c(l, mean(expr_data[,x]))
                                }
                                print(l)
                                u <- which(l < quantile(l,mean_filter))

                                expr_data <- expr_data[,-u]
                              }

                              plot(hist(expr_data), border = "blue", xlab = "final expression values",
                                    ylab = "final frequency")
                              return(as.data.frame(expr_data))
}


