#' @title names to ids
#' @description
#' Convert gene names and miRNA names into ensembl ids and accession ids respectively
#' @import org.Hs.eg.db
#' @import AnnotationDbi
#' @import miRBaseConverter
#' @importFrom stats aggregate
#' @importFrom stats formula
#' @param mat_data expression matrix or tool prediction matrix with columns as genes and rows as miRNAs
#' @param cgid A string containing current gene id format.'ACCNUM' by default
#' @param dgid A string containing desired gene id format.'ENSEMBL' by default.
#' @param gene_bool A boolean value. If it is TRUE, gene ids will be converted.TRUE by default
#' @param mir_bool A boolean value. If it is TRUE, miRNA ids will be converted.TRUE by default
#' @param cmirid A string containing current miRNA id format.empty by default.version will be automatically detected
#' @param dmirid A string containing desired miRNA id format.
#' @export
gene_miRNA_id_conversion <- function(mat_data,gene_bool = TRUE,cgid='ACCNUM',
                                     dgid = 'ENSEMBL',mir_bool = TRUE,cmirid = ' ', dmirid = 'Accession'){

                              # first convert the tool data into a dataframe
                              mat_data <- as.data.frame(mat_data)
                              if(gene_bool){

                              if(cgid == 'ENSEMBL_VERSION'){
                                #convert the gene ids into ensemble gene ids
                                colnames(mat_data) <- sub("\\.\\d+$", "", colnames(mat_data))
                              }
                              else{
                                # remove the versions from the colnames to just keep the gene ids
                                colnames(mat_data) <- sub("\\.\\d+$", "", colnames(mat_data))
                                # remove the multiple ids for on gene to keep the first id (eg. pita)
                                colnames(mat_data) <- sub(";.*","", colnames(mat_data))
                                # convert the gene ids
                                colens <- select(org.Hs.eg.db, keys=colnames(mat_data),
                                          columns=c(dgid, cgid), keytype=cgid)
                                # merge the 1:many mappings in one row with comma separated ensembl mappings
                                colens_1 <- aggregate(formula(paste0(dgid,"~",cgid)), colens, paste, collapse=",")
                                # find the indices of current ids in colens inside the colnames of mat_data
                                updated_colids <- match(colens_1[[cgid]], colnames(mat_data))
                                # update the arrangement of column names
                                mat_data <- mat_data[,updated_colids]
                                # update the column names to ensembl id's
                                colnames(mat_data) <- colens_1[[dgid]]

                              }
                                # remove the NA's
                                mat_data <- mat_data[,!is.na(colnames(mat_data))]

                                # uniqueify the column names
                                colnames(mat_data) <- make.unique(colnames(mat_data))
                            }

                            if(mir_bool){
                              # remove the versions from the end of miRNA ids, if it results in duplicates, then make unique
                              rownames(mat_data) <- make.unique(sub("\\.\\d+$", "", rownames(mat_data)))

                              if(cmirid == "Accession"){
                                res <- as.data.frame(miRNA_AccessionToName(rownames(mat_data),targetVersion = dmirid))
                              }
                              else if(dmirid == "Accession"){
                                ver <- checkMiRNAVersion(rownames(mat_data), verbose = TRUE)
                                res <- as.data.frame(miRNA_NameToAccession(rownames(mat_data),version = ver))
                              }
                              else{
                                ver <- checkMiRNAVersion(rownames(mat_data), verbose = TRUE)
                                res <- as.data.frame(miRNAVersionConvert(rownames(mat_data), targetVersion = dmirid, exact = TRUE,
                                                    verbose = TRUE))
                              }

                              # remove the NA's if there are any
                              nas <- which(is.na(res[,dmirid]))
                              if(length(nas) == 0){
                                rownames(mat_data) <- res[,dmirid]
                              }
                              else{

                              mat_data <- mat_data[-which(is.na(res[,dmirid])),]

                              rownames(mat_data) <- res[,dmirid][-which(is.na(res[,dmirid]))]

                              }
                            # make unique if duplicates exist
                            rownames(mat_data) <- make.unique(rownames(mat_data))
                            }
                            return(mat_data)
}
