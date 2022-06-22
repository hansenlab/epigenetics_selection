seqlevelsStyle(proms_good_power) <- "NCBI"
df <- data.frame(chr=seqnames(proms_good_power), start = start(proms_good_power)-1, end = end(proms_good_power))
write.table(df, file="/dcl01/hansen/data/db151_snps_hg19/all_ranges_in_proms_good_power.bed", 
            quote=F, sep="\t", row.names=F, col.names=F)

#bedtools intersect -a all_snps_db151_hg19.vcf -b all_ranges_in_proms_good_power.bed -header > all_snps_in_proms_good_power.vcf
#


#######parsing of the vcf file using the VariationAnnotation package
library(VariantAnnotation)
snps <- readVcf('/dcl01/hansen/data/db151_snps_hg19/all_snps_in_proms_good_power.vcf', "hg19")

snps_granges <- rowRanges(snps)
metadata_df <- info(snps)
topmed_af <- metadata_df$TOPMED
number_of_alleles <- lengths(topmed_af)


snps_granges_hg19 <- snps_granges[which(number_of_alleles == 2)]
topmed_af_hg19 <- topmed_af[which(number_of_alleles == 2)]
topmed_af_hg19 <- unlist(topmed_af_hg19)


snps_granges_hg19$REF_AF <- topmed_af_hg19[seq(1, length(topmed_af_hg19), by = 2)]
snps_granges_hg19$ALT_AF <- topmed_af_hg19[seq(2, length(topmed_af_hg19), by = 2)]

snps_granges_hg19_df <- data.frame(chr = seqnames(snps_granges_hg19), 
                                   start = start(snps_granges_hg19), 
                                   end = end(snps_granges_hg19), 
                                   strand = strand(snps_granges_hg19), 
                                   ref = snps_granges_hg19$REF, 
                                   alt = snps_granges_hg19$ALT, 
                                   ref_af = snps_granges_hg19$REF_AF, 
                                   alt_af = snps_granges_hg19$ALT_AF)

snps_granges_hg19 <- makeGRangesFromDataFrame(snps_granges_hg19_df, keep.extra.columns = TRUE)
seqlevelsStyle(snps_granges_hg19) <- "ucsc"

all_cpgs <- cpgs_with_meth_level
genome(all_cpgs) <- "hg19"
end(all_cpgs) <- start(all_cpgs)+1
overlaps <- findOverlaps(all_cpgs, proms_good_power)
cpgs_proms_good_power <- all_cpgs[unique(queryHits(overlaps))]

overlaps <- findOverlaps(snps_granges_hg19, cpgs_proms_good_power)
cg_snps <- snps_granges_hg19[unique(queryHits(overlaps))]
non_cg_snps <- snps_granges_hg19[-unique(queryHits(overlaps))]

cg_df <- data.frame(chr = seqnames(cg_snps), 
                    start = start(cg_snps), 
                    end = end(cg_snps), 
                    strand = strand(cg_snps), 
                    ref = cg_snps$ref, 
                    alt = cg_snps$alt.value, 
                    ref_af = cg_snps$ref_af, 
                    alt_af = cg_snps$alt_af)

non_cg_df <- data.frame(chr = seqnames(non_cg_snps), 
                        start = start(non_cg_snps), 
                        end = end(non_cg_snps), 
                        strand = strand(non_cg_snps), 
                        ref = non_cg_snps$ref, 
                        alt = non_cg_snps$alt.value, 
                        ref_af = non_cg_snps$ref_af, 
                        alt_af = non_cg_snps$alt_af)

save(cg_df, non_cg_df, file = "/dcl01/hansen/data/db151_snps_hg19/all_snps_in_proms_good_power_granges_df.rda")


