###three functions used to make the decile plots for DNAm (binary), 
###H3K4me3 chip-seq data (binary: presence or absence of mark), and H3K36me3 chip-seq data (presence/absence as well)
###

getPercentOfMethPromoters <- function(meth_vector, loeuf_vector, quantile_values, what = c("less", "greater", "in between"), cutoff){
  if (what == "less"){
    length(which(meth_vector[which(loeuf_vector <= quantile(loeuf_vector, quantile_values[1], na.rm = TRUE))] >= cutoff)) / 
      length(meth_vector[which(loeuf_vector <= quantile(loeuf_vector, quantile_values[1], na.rm = TRUE))])
  } else if (what == "in between"){
    length(which(meth_vector[which(loeuf_vector <= quantile(loeuf_vector, quantile_values[2], na.rm = TRUE) & 
                                     loeuf_vector > quantile(loeuf_vector, quantile_values[1], na.rm = TRUE))] >= cutoff)) / 
      length(meth_vector[which(loeuf_vector <= quantile(loeuf_vector, quantile_values[2], na.rm = TRUE) & 
                                 loeuf_vector > quantile(loeuf_vector, quantile_values[1], na.rm = TRUE))])
  } else {
    length(which(meth_vector[which(loeuf_vector > quantile(loeuf_vector, quantile_values[1], na.rm = TRUE))] >= cutoff)) / 
      length(meth_vector[which(loeuf_vector > quantile(loeuf_vector, quantile_values[1], na.rm = TRUE))])
  }
}

getPercentOfPromotersWithSignal <- function(promoter_granges, metric_1, metric_2, quantile_values, what = c("less", "greater", "in between")){
  if (what == "less"){
    length(which(promoter_granges$signal_strength[which(metric_1 <= 
                                                          quantile(metric_2, quantile_values[1], na.rm = TRUE))] > 0)) / 
      length(promoter_granges$signal_strength[which(metric_1 <= 
                                                      quantile(metric_2, quantile_values[1], na.rm = TRUE))])
  } else if (what == "in between"){
    length(which(promoter_granges$signal_strength[which(metric_1 <= 
                                                          quantile(metric_2, quantile_values[2], na.rm = TRUE) & 
                                                          metric_1 > 
                                                          quantile(metric_2, quantile_values[1], na.rm = TRUE))] > 0)) / 
      length(promoter_granges$signal_strength[which(metric_1 <= 
                                                      quantile(metric_2, quantile_values[2], na.rm = TRUE) & 
                                                      metric_1 > 
                                                      quantile(metric_2, quantile_values[1], na.rm = TRUE))])
  } else {
    length(which(promoter_granges$signal_strength[which(metric_1 > 
                                                          quantile(metric_2, quantile_values[1], na.rm = TRUE))] > 0)) / 
      length(promoter_granges$signal_strength[which(metric_1 > 
                                                      quantile(metric_2, quantile_values[1], na.rm = TRUE))])
  }
}

getPercentOfGenesWithCodingPeak <- function(mark_vector, loeuf_vector, quantile_values, what = c("less", "greater", "in between")){
  if (what == "less"){
    length(which(mark_vector[which(loeuf_vector <= quantile(loeuf_vector, quantile_values[1], na.rm = TRUE))] == "yes")) / 
      length(mark_vector[which(loeuf_vector <= quantile(loeuf_vector, quantile_values[1], na.rm = TRUE))])
  } else if (what == "in between"){
    length(which(mark_vector[which(loeuf_vector <= quantile(loeuf_vector, quantile_values[2], na.rm = TRUE) & 
                                     loeuf_vector > quantile(loeuf_vector, quantile_values[1], na.rm = TRUE))] == "yes")) / 
      length(mark_vector[which(loeuf_vector <= quantile(loeuf_vector, quantile_values[2], na.rm = TRUE) & 
                                 loeuf_vector > quantile(loeuf_vector, quantile_values[1], na.rm = TRUE))])
  } else {
    length(which(mark_vector[which(loeuf_vector > quantile(loeuf_vector, quantile_values[1], na.rm = TRUE))] == "yes")) / 
      length(mark_vector[which(loeuf_vector > quantile(loeuf_vector, quantile_values[1], na.rm = TRUE))])
  }
}


