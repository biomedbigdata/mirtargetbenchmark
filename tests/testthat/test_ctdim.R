library(mirtargetbenchmark)

test_that("check if the tool predictions are converted into a matrix with values",{
  expect_equal(class(convert_tool_data_into_matrix(targetscan_pred[1:1000,],miRNA_colname = "miRNA",
                                                   gene_colname = "Gene.ID",score_colname = "context...score")),"matrix")

})

test_that("check if the number of rows of the output is equal to the number of unique miRNAs in input predictions",{
  expect_equal(nrow(convert_tool_data_into_matrix(targetscan_pred[1:1000,],miRNA_colname = "miRNA",
  gene_colname = "Gene.ID",score_colname = "context...score")),length(unique(targetscan_pred[1:1000,]$miRNA)))
})


test_that("check if the number of columns of the output is equal to the number of unique genes in input predictions",{
  expect_equal(ncol(convert_tool_data_into_matrix(targetscan_pred[1:1000,], miRNA_colname = "miRNA",
  gene_colname = "Gene.ID",score_colname = "context...score")), length(unique(targetscan_pred[1:1000,]$Gene.ID)))
})

test_that("check if the class of scores in predictions is equal to the class of elements of output matrix",{
  expect_equal(class(convert_tool_data_into_matrix(targetscan_pred[1:1000,],miRNA_colname = "miRNA",
  gene_colname = "Gene.ID",score_colname = "context...score")[1,1]), class(targetscan_pred[1:1000,]$context...score[1]))

})
