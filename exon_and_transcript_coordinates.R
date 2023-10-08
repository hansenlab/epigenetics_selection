###obtain loeuf estimates of canonical transcripts from gnomAD
constraint <- read_tsv('constraint.txt')
constraint <- as(constraint, "DataFrame")
rownames(constraint) <- constraint$transcript

#obtain coordinates of coding exons and full transcripts
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
all_coding_exons <- cdsBy(edb, filter = SeqNameFilter(c(1:22)), columns = "tx_id")
all_coding_exons <- unlist(all_coding_exons)
genome(seqinfo(all_coding_exons)) <- "hg19"
seqlevelsStyle(all_coding_exons) <- "ucsc"
exons_to_use <- all_coding_exons[which(all_coding_exons$tx_id %in% constraint$transcript[which(constraint$canonical == TRUE)])]


tx_all <- transcripts(edb, filter = SeqNameFilter(c(1:22)), columns = c("tx_id", "gene_id"))
seqlevelsStyle(tx_all) <- "ucsc"
genome(tx_all) <- "hg19"
tx_all <- tx_all[which(tx_all$tx_id %in% constraint$transcript[which(constraint$canonical == TRUE)])]

