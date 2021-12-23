library(gdata)
dat <- read.xls("../data/Mesh Bias Trial.xlsx", sheet = 5)

head(dat)

dat$bin <- ifelse(dat$Retained.Escaped == "R", 1, 0)
fit <- glm(bin ~ Length, family = "binomial", data = dat)


pred_df <- data.frame(Length = seq(min(dat$Length), max(dat$Length), length = 1e3))
pred <- predict(fit, newdata = pred_df, se.fit = TRUE)
pred_df$mean <- plogis(pred$fit)
pred_df$lwr <- plogis(pred$fit - 1.96 * pred$se.fit)
pred_df$upr <- plogis(pred$fit + 1.96 * pred$se.fit)

l50 <- -coef(fit)[1] / coef(fit)[2]

pred_df$Length[which.min((pred_df$mean - 0.75)^2)] - 
pred_df$Length[which.min((pred_df$mean - 0.25)^2)]


pdf("../tex/sel_paper/figures/Selectivity_trial_14mm.pdf", height = 7, width = 8)
with(dat, plot(Length, bin, bty = "l", xlab = "Length", ylab = "Retained (1:yes, 0:no)", yaxt = "n"))
axis(side = 2, at = seq(0, 1, by = 0.1))
with(pred_df, lines(Length, mean))
with(pred_df, lines(Length, lwr, lty = 2))
with(pred_df, lines(Length, upr, lty = 2))
lines(c(l50, l50), c(0, 0.5), col = "blue")
lines(c(0, l50), c(0.5, 0.5), col = "blue")
text(37, 0.05, label = expression(L[50] *"= 35.6"))
dev.off()
