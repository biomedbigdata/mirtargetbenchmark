#' Regression models on the expression data
#' @export
#' @import glmnet
#' @import foreach
#' @import doSNOW
#' @import parallel
#' @param gene_expr input  preprocessed gene expression matrix
#' @param mir_expr  input preprocessed miRNA expression matrix
#' @param num String input that specifies the regression model to be used.
#' 'e','l','r' and 's' for elastic net, lasso, ridge and simple linear regression respectievly.
#' @param rsq_filter A filter to select models/genes using the r squared value.All genes with rsq < rsq_filter will be removed.
regression_results <- function(gene_expr , mir_expr , num, rsq_filter){
                      switch(num,
                      e = {
                        results <- foreach(x = seq(1,ncol(gene_expr),1),.combine = cbind,.inorder = TRUE,
                                .packages = c('glmnet','doSNOW'))%dopar%{

                                  y <- as.matrix(gene_expr[, x])
                                  X <- as.matrix(mir_expr)
                                  print("started")

                                  models <- foreach(alpha_value = seq(0, 1, 0.2)) %do% {
                                    cv.glmnet(X, y, alpha = alpha_value)
                                  }


                                  models.cvm <- sapply(models, function(model) {
                                    min(model$cvm)
                                  })

                                  #return model with smallest residual sum of squares
                                  models <- models[[which.min(models.cvm)]]


                                  lambda_value = models$lambda.min
                                  coeff = as.matrix(coef(models, s = lambda_value))

                                  coefficients <- as.vector(coeff)
                                  r2 <- models$glmnet.fit$dev.ratio[which(models$glmnet.fit$lambda == models$lambda.min)]
                                  coefficients[1] <- r2

                                  coeff_data <- as.data.frame(coefficients)
                                  return(coeff_data)

                        }},
                      l = {
                        results <- foreach(x = seq(1,ncol(gene_expr),1),.combine = cbind,.inorder = TRUE,
                                             .packages = c('glmnet','doSNOW'))%dopar%{

                                    y <- as.matrix(gene_expr[, x ])
                                    X <- as.matrix(mir_expr)

                                    print("started")

                                    models <- cv.glmnet(X, y, alpha = 1)

                                    lambda_value = models$lambda.min

                                    coeff = as.matrix(coef(models, s = lambda_value))


                                    coefficients <- as.vector(coeff)
                                    r2 <- models$glmnet.fit$dev.ratio[which(models$glmnet.fit$lambda == models$lambda.min)]
                                    coefficients[1] <- r2


                                    coeff_data <- as.data.frame(coefficients)
                                    print("coef_enter")

                                    return(coeff_data)
                                             }
                      },
                      r = {
                        results <- foreach(x = seq(1,ncol(gene_expr),1),.combine = cbind,.inorder = TRUE,
                                           .packages = c('glmnet','doSNOW'))%dopar%{

                                             y <- as.matrix(gene_expr[, x ])
                                             X <- as.matrix(mir_expr)

                                             print("started")

                                             models <- cv.glmnet(X, y, alpha = 0)

                                             lambda_value = models$lambda.min

                                             coeff = as.matrix(coef(models, s = lambda_value))


                                             coefficients <- as.vector(coeff)
                                             r2 <- models$glmnet.fit$dev.ratio[which(models$glmnet.fit$lambda == models$lambda.min)]
                                             coefficients[1] <- r2


                                             coeff_data <- as.data.frame(coefficients)
                                             print("coef_enter")

                                             return(coeff_data)
                                           }

                      },
                      s = {
                        results <- foreach(x = seq(1,ncol(gene_expr),1), .combine = cbind, .inorder = TRUE)%dopar%{

                          working_data <- as.data.frame(cbind(gene_expr[,x], mir_expr))
                          colnames(working_data)[1] = "gene"
                          model <- lm(gene~. ,data  = working_data)

                          coefficients <- as.vector(model$coefficients)
                          coefficients[1] <- summary(model)$r.squared
                          coeff_data <- as.data.frame(coefficients)


                          return(coeff_data)
                        }
                        },
                      stop("Incorrect input!!")
                      )
                      rows <- c("rsq", colnames(mir_expr))
                      colnames(results) <- colnames(gene_expr)
                      rownames(results) <- rows

                      keep <- which(results["rsq",] > rsq_filter)

                      results <- results[,keep]

                      return(results)

}


