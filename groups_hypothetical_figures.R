######
quartz(file = "methylation_groups.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
plot(density(rnorm(500, 0.01, 0.1), from = 0, to = 1), bty = 'l', xlab = "s_het", xlim = c(0, 1), main = "", 
     lwd = 3, col = alpha("red", 0.75), ylim = c(0, 4.2), yaxt = "n", xaxt = "n")
lines(density(rnorm(500, 0.2, 0.1), from = 0, to = 1), bty = 'l', xlab = "s_het", xlim = c(0, 1), main = "", 
      lwd = 3, col = "cornflowerblue")
lines(density(rnorm(500, 0.5, 0.1), from = 0, to = 1), bty = 'l', xlab = "s_het", xlim = c(0, 1), main = "", 
      lwd = 3, col = "orange")
#legend <- legend("topright", legend =  c("group 1", "group 2", "group 3"), col = c(alpha("red", 0.75), 
#                                "cornflowerblue", "orange"), lwd = 2.5, lty = "solid", bty = 'n')
axis(2, at = c(0, 2, 4))
axis(1, at = c(0, 0.4, 0.8), labels = c(0, 0.005, 0.01))
dev.off()

quartz(file = "methylation_groups_2.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
plot(density(rnorm(500, 0.2, 0.1), from = 0, to = 1), bty = 'l', xlab = "s_het", xlim = c(0, 0.5), main = "", 
     lwd = 3, col = alpha("red", 0.75), ylim = c(0, 4.5), yaxt = "n", xaxt = "n")
lines(density(rnorm(500, 0.2, 0.1), from = 0, to = 1), bty = 'l', xlab = "s_het", xlim = c(0, 0.5), main = "", 
      lwd = 3, col = "cornflowerblue")
lines(density(rnorm(500, 0.2, 0.1), from = 0, to = 1), bty = 'l', xlab = "s_het", xlim = c(0, 0.5), main = "", 
      lwd = 3, col = "orange")
legend <- legend("bottom", legend =  c("group 1", "group 2", "group 3"), col = c(alpha("red", 0.75), 
                                                                                 "cornflowerblue", "orange"), lwd = 2.5, lty = "solid", bty = 'n', cex = 0.62)
axis(2, at = c(0, 2, 4))
axis(1, at = c(0, 0.25, 0.5), labels = c(0, 0.005, 0.01))
dev.off()

quartz(file = "methylation_groups_3.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
par(mar = c(3,5.5,2,2))
plot(1, 0.5, bty = 'l', main = "", ylab = "methylation frequency\nin population",
     pch = 19, cex = 1.5, col = rgb(0,0,0,0.7), ylim = c(0, 1), xlim = c(0.8, 3.2), yaxt = "n", xaxt = "n", cex.lab = 1.25)
points(2, 0.57, bty = 'l', pch = 19, cex = 1.5, col = rgb(0,0,0,0.7))
points(3, 0.42, bty = 'l', pch = 19, cex = 1.5, col = rgb(0,0,0,0.7))
#legend <- legend("topright", legend =  c("group 1", "group 2", "group 3"), col = c(alpha("red", 0.75), 
#                                "cornflowerblue", "orange"), lwd = 2.5, lty = "solid", bty = 'n')
dev.off()

quartz(file = "methylation_groups_4.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
par(mar = c(3,5.5,2,2))
plot(1, 0.25, bty = 'l', main = "", ylab = "methylation frequency\nin population",
     pch = 19, cex = 1.5, col = rgb(0,0,0,0.7), ylim = c(0, 1), xlim = c(0.8, 3.2), yaxt = "n", xaxt = "n", cex.lab = 1.25)
points(2, 0.5, bty = 'l', pch = 19, cex = 1.5, col = rgb(0,0,0,0.7))
points(3, 0.75, bty = 'l', pch = 19, cex = 1.5, col = rgb(0,0,0,0.7))
#legend <- legend("topright", legend =  c("group 1", "group 2", "group 3"), col = c(alpha("red", 0.75), 
#                                "cornflowerblue", "orange"), lwd = 2.5, lty = "solid", bty = 'n')
dev.off()


quartz(file = "methylation_groups_3.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1, 2))
vec <- runif(10, 0.3, 0.6)
plot(vec, bty = 'l', main = "", xlab = "LOEUF decile", ylab = "",
     pch = 19, cex = 1.5, col = rgb(0,0,0,0.7), ylim = c(0, 1), xlim = c(0.8, 10.2), yaxt = "n", xaxt = "n", cex.lab = 1.25)
axis(1, at = c(1:10), cex.axis = 0.62)
axis(2, at = c(0, 0.5, 1))

vec <- seq(0.15, 0.85, by = 0.1) + rnorm(8, 0, 0.1)
vec <- c(0.05, vec, 0.95)
plot(vec, bty = 'l', main = "", xlab = "LOEUF decile", ylab = "",
     pch = 19, cex = 1.5, col = rgb(0,0,0,0.7), ylim = c(0, 1), xlim = c(0.8, 10.2), yaxt = "n", xaxt = "n", cex.lab = 1.25)
axis(1, at = c(1:10), cex.axis = 0.62)
axis(2, at = c(0, 0.5, 1))
#legend <- legend("topright", legend =  c("group 1", "group 2", "group 3"), col = c(alpha("red", 0.75), 
#                                "cornflowerblue", "orange"), lwd = 2.5, lty = "solid", bty = 'n')
dev.off()

quartz(file = "methylation_groups_4.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
par(mar = c(3,5.5,2,2))
plot(1, 0.25, bty = 'l', main = "", ylab = "methylation frequency\nin population",
     pch = 19, cex = 1.5, col = rgb(0,0,0,0.7), ylim = c(0, 1), xlim = c(0.8, 3.2), yaxt = "n", xaxt = "n", cex.lab = 1.25)
points(2, 0.5, bty = 'l', pch = 19, cex = 1.5, col = rgb(0,0,0,0.7))
points(3, 0.75, bty = 'l', pch = 19, cex = 1.5, col = rgb(0,0,0,0.7))
#legend <- legend("topright", legend =  c("group 1", "group 2", "group 3"), col = c(alpha("red", 0.75), 
#                                "cornflowerblue", "orange"), lwd = 2.5, lty = "solid", bty = 'n')
dev.off()


quartz(file = "methylation_groups_5.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
plot(density(rnorm(500, 0.2, 0.2), from = 0, to = 1), bty = 'l', xlab = "LOEUF", xlim = c(0, 1), main = "", 
     lwd = 3, col = alpha("red", 0.75), ylim = c(0, 3), yaxt = "n", xaxt = "n")
lines(density(rnorm(500, 0.5, 0.2), from = 0, to = 1), bty = 'l', xlab = "LOEUF", xlim = c(0, 1), main = "", 
      lwd = 3, col = "orange")
legend <- legend("top", legend =  c("group 1", "group 2"), col = c(alpha("red", 0.75), 
                                                                   "orange"), lwd = 2.5, lty = "solid", bty = 'n')
axis(2, at = c(0, 1, 2))
axis(1, at = c(0, 0.4, 0.8),labels = c("0", "0.5", "1"))
dev.off()