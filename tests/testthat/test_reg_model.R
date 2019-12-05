library(mirtargetbenchmark)

test_that("check if the result is a matrix or dataframe",{
  expect_equal(class(regression_results(gene_expr, mir_expr,'s',0.2)), "data.frame")

})


