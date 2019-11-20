#' Convert gene names and miRNA names into ensembl ids and accession ids respectively
#' @export
#' @importFrom biomaRt getBM
#' @importFrom biomaRt useEnsembl
#' @import miRBaseConverter
#' @param mat_data expression matrix or tool prediction matrix with columns as genes and rows as miRNAs
#' @param cgid A string containing current gene id format.
#' @param dgid A string containing desired gene id format.'ensembl_gene_id' by default.
#' @param gene_bool A boolean value. If it is TRUE, gene ids will be converted.TRUE by default
#' @param mir_bool A boolean value. If it is TRUE, miRNA ids will be converted.TRUE by default
#' @param cmirid A string containing current miRNA id format.empty by default.version will be automatically detected
#' @param dmirid A string containing desired miRNA id format.
gene_miRNA_id_conversion <- function(mat_data,gene_bool = TRUE,cgid,
                                     dgid = 'ensembl_gene_id',mir_bool = TRUE,cmirid = ' ', dmirid = 'Accession'){

                              if(gene_bool){
                              if(cgid == 'ensembl_gene_id_version'){
                                colnames(mat_data) <- sub("\\.\\d+$", "", colnames(mat_data))
                              }
                              else{
                                ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",host="www.ensembl.org")
                                ensembl_id <- getBM(attributes=c(dgid,cgid),filters = cgid,
                                                    values = colnames(mat_data), mart = ensembl)

                                sequenc <- match(colnames(mat_data), ensembl_id[,cgid])
                                colnames(mat_data) <- ensembl_id[,dgid][sequenc]

                                mat_data <- mat_data[, !duplicated(colnames(mat_data))]
                              }
                            }

                            if(mir_bool){
                              if(cmirid == "Accession"){
                                res <- as.data.frame(miRNA_AccessionToName(rownames(mat_data),targetVersion = dmirid))
                              }
                              else if(dmirid == "Accession"){
                                cmirid <- ver <- checkMiRNAVersion(rownames(mat_data), verbose = TRUE)
                                res <- as.data.frame(miRNA_NameToAccession(rownames(mat_data),version = ver))
                              }
                              else{
                                cmirid <- ver <- checkMiRNAVersion(rownames(mat_data), verbose = TRUE)
                                miRNAVersionConvert(rownames(mat_data), targetVersion = dmirid, exact = TRUE,
                                                    verbose = TRUE)
                              }

                              print("done")
                              nas <- which(is.na(res[,dmirid]))
                              if(length(nas) == 0){
                                rownames(mat_data) <- res[,dmirid]
                              }
                              else{
                              mat_data <- mat_data[-which(is.na(res[,dmirid])),]
                              print("done2")
                              rownames(mat_data) <- na.omit(res[,dmirid])

                              }
                            }
                            return(mat_data)
}
