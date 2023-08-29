###
###x --> y --> z
###y <- x + noise
###z <- y + noise
###y* <- y + m_noise
###the question is how m_noise affects inferences about the association between x and z when controlling for y

mediation_noise <- sapply(seq(0, 0.5, by = 0.05), function(xx) {
  p_vals <- replicate(10000, {
    x <- rnorm(10000)
    y <- x + rnorm(10000)
    y_star <- y + rnorm(10000, 0, xx)
    z <- y + rnorm(10000)
    
    res1 <- residuals(lm(x ~ y_star))
    res2 <- residuals(lm(z ~ y_star))
    perm_test <- permcor(res1, res2)
    perm_test_no_mediation <- permcor(x, z)
    c(perm_test, perm_test[1]/perm_test_no_mediation[1])
  })
  p_vals
})

quartz(file = "mediation_noise.pdf", width = 11, height = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(2, 5))
hist(mediation_noise[1, seq(2, 30000, by = 3)], font.main = 1, lty = 0, breaks = 40, 
     freq = FALSE, xaxt = 'n', yaxt = 'n', col = alpha("red", 0.57), xlab = "p value", main = "")
axis(1, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"))
axis(2, at = c(0, 1.2), labels = c("0", "1.2"))
hist(mediation_noise[2, seq(2, 30000, by = 3)], font.main = 1, lty = 0, col = alpha("red", 0.57), breaks = 40, 
     freq = FALSE, xaxt = 'n', yaxt = 'n', xlab = "", main = "")
axis(1, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"))
axis(2, at = c(0, 1), labels = c("0", "1"))
hist(mediation_noise[3, seq(2, 30000, by = 3)], font.main = 1, lty = 0, col = alpha("red", 0.57), breaks = 40, 
     freq = FALSE, xaxt = 'n', yaxt = 'n', xlab = "", main = "")
axis(1, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"))
axis(2, at = c(0, 2.5), labels = c("0", "2.5"))
hist(mediation_noise[4, seq(2, 30000, by = 3)], font.main = 1, lty = 0, col = alpha("red", 0.57), breaks = 40, 
     freq = FALSE, xaxt = 'n', xlim = c(0, 1), yaxt = 'n', xlab = "", main = "")
axis(1, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"))
axis(2, at = c(0, 10), labels = c("0", "10"))
hist(mediation_noise[5, seq(2, 30000, by = 3)], font.main = 1, lty = 0, col = alpha("red", 0.57), breaks = 40, 
     freq = FALSE, xaxt = 'n', xlim = c(0, 1), yaxt = 'n', xlab = "", main = "")
axis(1, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"))
axis(2, at = c(0, 30), labels = c("0", "30"))


hist(mediation_noise[1, seq(3, 30000, by = 3)], font.main = 1, lty = 0, breaks = 40, 
     freq = FALSE, xaxt = 'n', yaxt = 'n', col = alpha("red", 0.57), main = "", xlab = "corr. coef. ratio")
axis(1, at = c(-0.06, 0, 0.06), labels = c("-0.06", "0", "0.06"))
axis(2, at = c(0, 25), labels = c("0", "25"))
abline(v = 0, lty = "longdash", col = rgb(0,0,0,0.7))
hist(mediation_noise[2, seq(3, 30000, by = 3)], font.main = 1, lty = 0, breaks = 40, 
     freq = FALSE, xaxt = 'n', yaxt = 'n', col = alpha("red", 0.57), main = "", xlab = "")
axis(1, at = c(-0.06, 0, 0.06), labels = c("-0.06", "0", "0.06"))
axis(2, at = c(0, 22), labels = c("0", "22"))
abline(v = 0, lty = "longdash", col = rgb(0,0,0,0.7))
hist(mediation_noise[3, seq(3, 30000, by = 3)], font.main = 1, lty = 0, breaks = 40, 
     freq = FALSE, xaxt = 'n', yaxt = 'n', col = alpha("red", 0.57), main = "", xlab = "")
axis(1, at = c(-0.05, 0, 0.05), labels = c("-0.05", "0", "0.05"))
axis(2, at = c(0, 22), labels = c("0", "22"))
abline(v = 0, lty = "longdash", col = rgb(0,0,0,0.7))
hist(mediation_noise[4, seq(3, 30000, by = 3)], font.main = 1, lty = 0, breaks = 40, 
     freq = FALSE, xaxt = 'n', yaxt = 'n', col = alpha("red", 0.57), main = "", xlab = "")
axis(1, at = c(-0.04, 0.08), labels = c("-0.04", "0.08"))
axis(2, at = c(0, 25), labels = c("0", "25"))
abline(v = 0, lty = "longdash", col = rgb(0,0,0,0.7))
hist(mediation_noise[5, seq(3, 30000, by = 3)], font.main = 1, lty = 0, breaks = 40, 
     freq = FALSE, xaxt = 'n', yaxt = 'n', col = alpha("red", 0.57), main = "", xlab = "")
axis(1, at = c(-0.02, 0.10), labels = c("-0.02", "0.10"))
axis(2, at = c(0, 22), labels = c("0", "22"))
abline(v = 0, lty = "longdash", col = rgb(0,0,0,0.7))
dev.off()




