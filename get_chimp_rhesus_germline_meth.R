#####
library(rtracklayer)
cpgs_with_meth_level_rhesus <- import("~/Downloads/germline_methylation/rhesus_germline_meth.bw", 
                                      format = "BigWig")
cpgs_with_coverage_rhesus <- import("~/Downloads/germline_methylation/rhesus_germline_meth_coverage.bw", 
                                    format = "BigWig")

cpgs_with_meth_level_rhesus$meth_level <- cpgs_with_meth_level_rhesus$score
cpgs_with_meth_level_rhesus$coverage <- cpgs_with_coverage_rhesus$score

#####
cpgs_with_meth_level_chimp <- import("~/Downloads/germline_methylation/chimp_germline_meth.bw", 
                                     format = "BigWig")
cpgs_with_coverage_chimp <- import("~/Downloads/germline_methylation/chimp_germline_meth_coverage.bw", 
                                   format = "BigWig")

cpgs_with_meth_level_chimp$meth_level <- cpgs_with_meth_level_chimp$score
cpgs_with_meth_level_chimp$coverage <- cpgs_with_coverage_chimp$score
