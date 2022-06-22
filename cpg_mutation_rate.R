###
chimp_good_cov_cpgs <- cpgs_with_meth_level_chimp[which(cpgs_with_meth_level_chimp$coverage >= 10)]
rhesus_good_cov_cpgs <- cpgs_with_meth_level_rhesus[which(cpgs_with_meth_level_rhesus$coverage >= 10)]

human_proms_chimp_cpgs_overlaps <- findOverlaps(chimp_good_cov_cpgs, proms_good_power - 1500)
chimp_good_cov_cpgs_proms <- chimp_good_cov_cpgs[unique(queryHits(human_proms_chimp_cpgs_overlaps))]
chimp_good_cov_cpgs_proms$tx_id <- proms_good_power$tx_id[subjectHits(human_proms_chimp_cpgs_overlaps)]
chimp_good_cov_cpgs_proms$oe_lof_upper <- proms_good_power$oe_lof_upper[subjectHits(human_proms_chimp_cpgs_overlaps)]

chimp_good_cov_cpgs_proms_rhesus_cpgs_overlaps <- findOverlaps(chimp_good_cov_cpgs_proms, rhesus_good_cov_cpgs)
chimp_good_cov_cpgs_proms <- chimp_good_cov_cpgs_proms[unique(queryHits(chimp_good_cov_cpgs_proms_rhesus_cpgs_overlaps))]
chimp_good_cov_cpgs_proms$meth_level_rhesus <- rhesus_good_cov_cpgs$meth_level[subjectHits(chimp_good_cov_cpgs_proms_rhesus_cpgs_overlaps)]

chimp_good_cov_cpgs_proms$cpg_in_human <- "yes"
human_overlaps <- findOverlaps(chimp_good_cov_cpgs_proms, cpgs_proms_good_power)
chimp_good_cov_cpgs_proms$cpg_in_human[-queryHits(human_overlaps)] <- "no"

chimp_good_cov_cpgs_proms$meth_level_chimp <- chimp_good_cov_cpgs_proms$meth_level
chimp_good_cov_cpgs_proms$meth_level[which(chimp_good_cov_cpgs_proms$cpg_in_human == "yes")] <- cpgs_proms_good_power$meth_level[
  subjectHits(human_overlaps)]
chimp_good_cov_cpgs_proms$meth_level[which(chimp_good_cov_cpgs_proms$cpg_in_human == "no")] <- NA 


##
meth_cpgs <- chimp_good_cov_cpgs_proms[which(chimp_good_cov_cpgs_proms$meth_level_chimp >= 0.8 & 
                                               chimp_good_cov_cpgs_proms$meth_level_rhesus >= 0.8)]
hypometh_cpgs <- chimp_good_cov_cpgs_proms[which(chimp_good_cov_cpgs_proms$meth_level_chimp <= 0.2 & 
                                                   chimp_good_cov_cpgs_proms$meth_level_rhesus <= 0.2)]

meth_cpgs$loeuf_quartile <- 1
meth_cpgs$loeuf_quartile[which(meth_cpgs$oe_lof_upper > quantile(proms_good_power$oe_lof_upper, 0.25) & 
                                 meth_cpgs$oe_lof_upper <= quantile(proms_good_power$oe_lof_upper, 0.5))] <- 2
meth_cpgs$loeuf_quartile[which(meth_cpgs$oe_lof_upper > quantile(proms_good_power$oe_lof_upper, 0.5) & 
                                 meth_cpgs$oe_lof_upper <= quantile(proms_good_power$oe_lof_upper, 0.75))] <- 3
meth_cpgs$loeuf_quartile[which(meth_cpgs$oe_lof_upper > quantile(proms_good_power$oe_lof_upper, 0.75))] <- 4


hypometh_cpgs$loeuf_quartile <- 1
hypometh_cpgs$loeuf_quartile[which(hypometh_cpgs$oe_lof_upper > quantile(proms_good_power$oe_lof_upper, 0.25) & 
                                     hypometh_cpgs$oe_lof_upper <= quantile(proms_good_power$oe_lof_upper, 0.5))] <- 2
hypometh_cpgs$loeuf_quartile[which(hypometh_cpgs$oe_lof_upper > quantile(proms_good_power$oe_lof_upper, 0.5) & 
                                     hypometh_cpgs$oe_lof_upper <= quantile(proms_good_power$oe_lof_upper, 0.75))] <- 3
hypometh_cpgs$loeuf_quartile[which(hypometh_cpgs$oe_lof_upper > quantile(proms_good_power$oe_lof_upper, 0.75))] <- 4

prop.table(table(as.factor(is.na(meth_cpgs$meth_level))))[2]
prop.table(table(as.factor(is.na(hypometh_cpgs$meth_level))))[2]

prop.table(table(as.factor(meth_cpgs$meth_level <= 0.2)))[2]
prop.table(table(as.factor(hypometh_cpgs$meth_level >= 0.8)))[2]

table(as.factor(is.na(meth_cpgs$meth_level[which(meth_cpgs$loeuf_quartile == 4)])))
table(as.factor(is.na(hypometh_cpgs$meth_level[which(hypometh_cpgs$loeuf_quartile == 4)])))
##
####mutation rate per generation estimation
118/320000
0.00036875/(118 + 1372)

meth_mut_rate <- 2.474832e-07
hypometh_mut_rate <- 7.971193e-08






