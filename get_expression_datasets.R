###fetal expression
fetal_expr <- read_csv('fetal_single_cell_expression/gene_fraction_celltype.txt') #this is from Jay Shendure's lab: https://descartes.brotmanbaty.org/bbi/human-gene-expression-during-development/
fetal_expr$RowID <- gsub("[.].*", "", fetal_expr$RowID)

fetal_expr_gene_ids <- fetal_expr$RowID
fetal_expr_all_cell_types <- as.matrix(fetal_expr[, -1])
rownames(fetal_expr_all_cell_types) <- fetal_expr_gene_ids

fetal_expr_hypometh_1 <- fetal_expr[proms_good_power$gene_id, ]
fetal_expr_hypometh_2 <- fetal_expr[tss_with_hypometh$gene_id, ]
fetal_expr_hypometh_3 <- fetal_expr[proms_all[df_k36me3_vs_expr$tx_id]$gene_id, ]


###germline(testis) expression
load(file = "median_expr_testis.rda") ###from GTEx (same dataset used in Boukas et al., 2020)

###esc expression
library(tximport)
library(EnsDb.Hsapiens.v75)
files <- paste0("h1_esc_rnaseq/quants/", 
                list.files("h1_esc_rnaseq/quants"), "/quant.sf") ###from ENCODE after downloading fastq files and mapping/quantifying transcripts with Salmon

names(files) <- paste0("healthy", 1:2)

edb <- EnsDb.Hsapiens.v75
#k <- keys(edb, keytype = "GENEID")
#df <- select(edb, keys = k, keytype = "GENEID", columns = "TXNAME")

df <- transcripts(edb, return.type = "DataFrame")
tx2gene <- df[, c(8, 7)]
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, 
                countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = TRUE)


expr_h1_esc <- log2(txi$counts + 1)