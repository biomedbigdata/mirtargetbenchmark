library(mirtargetbenchmark)

test_that("check if there are any duplicate or NA gene ids",{
    expect_equal(anyNA(colnames(gene_miRNA_id_conversion(convert_tool_data_into_matrix(targetscan_pred[1:1000,],
    miRNA_colname = "miRNA",gene_colname = "Gene.ID",score_colname = "context...score"),
    cgid = 'ensembl_gene_id_version', cmirid = 'refseq_mrna'))), FALSE)
})

test_that("check if there are any duplicate or NA miRNA ids",{
    expect_equal(anyNA(rownames(gene_miRNA_id_conversion(convert_tool_data_into_matrix(targetscan_pred[1:1000,],
    miRNA_colname = "miRNA",gene_colname = "Gene.ID",score_colname = "context...score"),
    cgid = 'ensembl_gene_id_version', cmirid = 'refseq_mrna'))), FALSE)
})
