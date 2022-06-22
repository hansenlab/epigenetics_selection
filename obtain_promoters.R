####
load(file = "proms_to_use_and_good_power.rda") #proms_to_use contains the high-confidence promoters from Boukas et al., 2020
proms_to_use$low_power_likely_intolerant <- FALSE
proms_to_use$low_power_likely_intolerant[which(proms_to_use$oe_lof < 0.5 & proms_to_use$exp_lof < 10)] <- TRUE

proms_good_power <- proms_to_use[which(proms_to_use$exp_lof >= 10)]
proms_good_power <- keepNonOverlappingRanges(proms_good_power)
proms_good_power <- keepNonOverlappingRanges(proms_good_power) #for some reason I need to run it twice

overlapping_low_power_intolerant <- proms_good_power$tx_id[
  unique(queryHits(findOverlaps(proms_good_power, proms_to_use[which(proms_to_use$low_power_likely_intolerant == TRUE)], 
                                ignore.strand = TRUE)))]


proms_good_power <- proms_good_power[-which(proms_good_power$oe_lof_upper > 0.5 & 
                                              proms_good_power$tx_id %in% overlapping_low_power_intolerant)]
save(proms_good_power, file = "proms_good_power_exp_lof_more_than_10.rda")