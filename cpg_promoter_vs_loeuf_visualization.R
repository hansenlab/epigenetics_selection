meth_vs_loeuf_vec <- c(getPercentOfMethPromoters(meth_percentage, human_loeufs, 0.1, "less", 0.8), 
         getPercentOfMethPromoters(meth_percentage, human_loeufs, c(0.1, 0.2), "in between", 0.8),
         getPercentOfMethPromoters(meth_percentage, human_loeufs, c(0.2, 0.3), "in between", 0.8),
         getPercentOfMethPromoters(meth_percentage, human_loeufs, c(0.3, 0.4), "in between", 0.8),
         getPercentOfMethPromoters(meth_percentage, human_loeufs, c(0.4, 0.5), "in between", 0.8),
         getPercentOfMethPromoters(meth_percentage, human_loeufs, c(0.5, 0.6), "in between", 0.8),
         getPercentOfMethPromoters(meth_percentage, human_loeufs, c(0.6, 0.7), "in between", 0.8),
         getPercentOfMethPromoters(meth_percentage, human_loeufs, c(0.7, 0.8), "in between", 0.8),
         getPercentOfMethPromoters(meth_percentage, human_loeufs, c(0.8, 0.9), "in between", 0.8),
         getPercentOfMethPromoters(meth_percentage, human_loeufs, 0.9, "greater", 0.8))



##
quartz(file = "methylation_vs_loeuf_promoter.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
plot(meth_vs_loeuf_vec, pch = 19, cex = 1.2, col = alpha("red", 0.75), bty = 'l', ylim = c(0, 0.22),
     main = "proximal promoter\n(+/- 500bp from TSS)", font.main = 1, ylab = "% genes w/ methylation", xlab = "LOEUF decile", 
     xaxt = 'n', yaxt = 'n')
axis(1, at = 1:10, cex.axis = 0.7)
axis(2, at = c(0, 0.20))
dev.off()



