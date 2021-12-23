##----------------------------
## Prob(capture|contact)
## CM: 13/9/2021
## want the first passage time, which is a function of length
## use this to get the time of contact
##----------------------------
library(reshape)
library(ggplot2); theme_set(theme_bw())
source("selectivity_functions.R")

## time
T <- 100
t <- 1:T
l <- seq(0.001, 60, length = 300)
L <- length(l)

load("length_fit_bb.RData")
theta <- fit_bb$par
vcov <- solve(I_bb)
## visualise
## parameters
fixed <- c("L50w" = 25, "SRw" = 1, "pw" = 0.01, "pc" = 0.01)

## contact probability
pc <- as.numeric(fixed["pc"]) ## corner
pf <- as.numeric(plogis(theta["logitpf"])) ## fyke
pw <- as.numeric(fixed["pw"]) ## wall 
## contact selection
L50c <- as.numeric(exp(theta["lnL50c"]))
L50f <- as.numeric(exp(theta["lnL50f"]))
L50w <- as.numeric(fixed["L50w"]) 
SRc <- as.numeric(exp(theta["lnSRc"]))
SRf <- as.numeric(exp(theta["lnSRf"]))
SRw <- as.numeric(fixed["SRw"])

## simulate to see what P(capture|contact) should look like

nsim <- 100
res <- expand.grid(l = seq(23, 32, length = 20),
                   sim = 1:nsim,
                   f_contact = NA,
                   f_capture = NA,
                   f_time_contact = NA,
                   c_contact = NA,
                   c_capture = NA,
                   c_time_contact = NA)
state <- matrix(0, nrow = 4, ncol = length(t))

m <- nrow(res)

for(j in 1:nrow(res)){
    print(round(j/m * 100))
    state <- state * 0
    state[1,1] <- 1 ## free
    l <- res$l[j]
    for(i in 2:length(t)){
        Q <- matrix(0, nrow = 4, ncol = 4)
        ## free
        Q[1,2] <- pf * rl(l, L50f, SRf)
        Q[1,3] <- pc * rl(l, L50c, SRc)
        Q[1,4] <- pw * (1 - rl(l, L50w, SRw))
        Q[1,1] <- -pf * rl(l, L50f, SRf) - pc * rl(l, L50c, SRc) - pw * (1 - rl(l, L50w, SRw))
        ## fyke
        Q[2,1] <- 1 - rl(l, L50f, SRf)
        Q[2,2] <- rl(l, L50f, SRf) - 1
        ## corner
        Q[3,3] <- rl(l, L50c, SRc) - 1
        Q[3,4] <- 1 - rl(l, L50c, SRc)
        ## free all stay in free so all zero intensities
        ## use the difference in time here
        Phi <- state[,i-1] %*% expm::expm((t[i] - t[i-1]) * Q)
        state[,i] <- rmultinom(1, size = 1, prob = Phi)
    }
    ## individual state occupied vector
    state_ind <- apply(state, 2, function(z){which(z == 1)})
    res$f_contact[j] <- 2 %in% state_ind
    res$f_capture[j] <- state_ind[length(t)] == 2
    if(res$f_contact[j]){
        res$f_time_contact[j] <- t[which(state_ind == 2)][1] ## ignoring re-contacts for now
    }
    res$c_contact[j] <- 3 %in% state_ind
    res$c_capture[j] <- state_ind[length(t)] == 3
    if(res$c_contact[j]){
        res$c_time_contact[j] <- t[which(state_ind == 3)][1] ## ignoring re-contacts for now
    }
}

f_res <- aggregate(f_capture ~ l, function(z){sum(z)/length(z)}, data = subset(res, f_contact))
c_res <- aggregate(c_capture ~ l, function(z){sum(z)/length(z)}, data = subset(res, c_contact))

## Analytical
## lengths to test
lvec <- seq(10, 70, length = 1000)
Q <- array(0, dim = c(4, 4, length(lvec)))
## free
Q[1,2,] <- pf * rl(lvec, L50f, SRf)
Q[1,3,] <- pc * rl(lvec, L50c, SRc)
Q[1,4,] <- pw * (1 - rl(lvec, L50w, SRw))
Q[1,1,] <- -colSums(Q[1,2:4,])
## fyke
Q[2,1,] <- 1 - rl(lvec, L50f, SRf)
Q[2,2,] <- -Q[2,1,]
## corner
Q[3,4,] <- 1 - rl(lvec, L50c, SRc)
Q[3,3,] <- -Q[3,4,]

## probability of capture given contact
pcgc <- apply(Q, 3, function(z){
    ## full system
    P_full <- expm::expm(T * z)
    ## assume fyke and corner are absorbing
    ## fyke
    zf_absorb <- z
    zf_absorb[2, ] <- 0
    Pf_absorb <- expm::expm(T * zf_absorb)
    ## as all start in state 1 (free)
    pfgcf <- P_full[1, 2] / Pf_absorb[1,2]
    ## corner
    zc_absorb <- z
    zc_absorb[3, ] <- 0
    Pc_absorb <- expm::expm(T * zc_absorb)
    ## as all start in state 1 (free)
    pcgcc <- P_full[1, 3] / Pc_absorb[1,3]
    return(c(pfgcf = pfgcf, pcgcc = pcgcc))
})

pcgc <- as.data.frame(t(pcgc))
pcgc$l <- lvec

## contact selections
## fyke
idx <- which.min((pcgc$pfgcf - 0.5)^2)
L50f_48 <- pcgc$l[idx]
idx25 <- which.min((pcgc$pfgcf - 0.25)^2)
idx75 <- which.min((pcgc$pfgcf - 0.75)^2)
SRf_48 <- pcgc$l[idx75] - pcgc$l[idx25]

## corner
idx <- which.min((pcgc$pcgcc - 0.5)^2)
L50c_48 <- pcgc$l[idx]
idx25 <- which.min((pcgc$pcgcc - 0.25)^2)
idx75 <- which.min((pcgc$pcgcc - 0.75)^2)
SRc_48 <- pcgc$l[idx75] - pcgc$l[idx25]

par(mfrow = c(2, 1), mar = c(2, 2, 1, 1))
with(f_res, plot(l, f_capture, xlim = c(20, 35), pch = 19, col = "cadetblue"))
curve(rl(x, L50f, SRf), add = TRUE, col = "cadetblue", n = 1e3, lty = 2)
with(pcgc, lines(l, pfgcf, col = "cadetblue"))
##
with(c_res, plot(l, c_capture, xlim = c(20, 35), pch = 19, col = "indianred1"))
curve(rl(x, L50c, SRc), add = TRUE, col = "indianred1", n = 1e3, lty = 2)
with(pcgc, lines(l, pcgcc, col = "indianred1"))
 
## figure for Russell
## Bevacqua parameters
library(gdata)
dat <- read.xls(xls = "../data/Yellow Eel Survey Details.xlsx", sheet = "Total stage and corrections")

enc <- subset(dat, Method == "Enclosure" & Location %in% c("Feeagh", "Furnace"))
enc$lmm <- enc$Length * 10
## trunk area
## trunk area
## fish density: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4297919/
## 1.09 kg/L
## 0.00109 g/mm^3
rho <- 0.001
ab_fit <- lm(log(Weight) ~ log(lmm), data = enc)
a <- exp(coef(ab_fit)[1])
b <- coef(ab_fit)[2]

enc$A <- a * (enc$Length * 10)^(b - 1)/ rho

estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

## Fyke
Lmin <- min(enc$lmm[enc$Method.Type == "Fyke"])
Lmax <- estimate_mode(enc$lmm[enc$Method.Type == "Fyke"]) ## modal length
Amin <- a * Lmin^(b - 1) / rho
Amax <- a * Lmax^(b - 1) / rho
A50 <- (Amax + Amin) / 2
alpha <- 0.95
eta <- log((1 - 0.95)/(1 + 0.95)) / (Amax - A50)

A <- a * (lvec * 10)^(b - 1)/ rho
sel_f <- plogis(-eta * (A - A50))

## Corner
Lmin <- min(enc$lmm[enc$Method.Type == "Corner"])
Lmax <- estimate_mode(enc$lmm[enc$Method.Type == "Corner"]) ## modal length
Amin <- a * Lmin^(b - 1) / rho
Amax <- a * Lmax^(b - 1) / rho
A50 <- (Amax + Amin) / 2
alpha <- 0.95
eta <- log((1 - 0.95)/(1 + 0.95)) / (Amax - A50)
A <- a * (lvec * 10)^(b - 1)/ rho
sel_c <- plogis(-eta * (A - A50))

##---------------------------
## PROPORTION AT LENGTH PLOT
##---------------------------
## include proportion at length plot here
enc$count <- 1
brks <- seq(10.5, 80.5, by = 1)

midpoints <- function(x){
    lower <- as.numeric(gsub(",.*","",gsub("\\(|\\[|\\)|\\]","", x)))
    upper <- as.numeric(gsub(".*,","",gsub("\\(|\\[|\\)|\\]","", x)))
    return(lower+(upper-lower)/2)
}
enc$len2 <- cut(enc$Length, breaks = brks, right = FALSE)
enc$len2 <-  midpoints(enc$len2)

enc_lf <- aggregate(count ~ len2 + Method.Type, FUN = sum, data = enc)
enc_lf_cast <- cast(enc_lf, len2 ~ Method.Type, value = "count")
enc_lf_cast$Corner[is.na(enc_lf_cast$Corner)] <- 0
enc_lf_cast$Fyke[is.na(enc_lf_cast$Fyke)] <- 0
enc_lf_cast$Total <- with(enc_lf_cast, Corner + Fyke)
enc_lf_cast$pFyke <- with(enc_lf_cast, Fyke / Total)


library(mvtnorm)

ns <- 1e3
thetas <- rmvnorm(ns, theta, vcov)

pf_hat <- matrix(NA, nrow = ns, ncol = length(lvec))

for(i in 1:ns){
    L50c <- exp(thetas[i,"lnL50c"])
    L50f <- exp(thetas[i,"lnL50f"])
    L50w <- fixed["L50w"]
    SRc <- exp(thetas[i,"lnSRc"])
    SRf <- exp(thetas[i,"lnSRf"])
    SRw <- fixed["SRw"]
    pc <- fixed["pc"]
    pf <- plogis(thetas[i,"logitpf"])
    pw <- fixed["pw"]
    ##
    Q <- array(0, dim = c(4, 4, length(lvec)))
    ## free
    Q[1,2,] <- pf * rl(lvec, L50f, SRf)
    Q[1,3,] <- pc * rl(lvec, L50c, SRc)
    Q[1,4,] <- pw * (1 - rl(lvec, L50w, SRw))
    Q[1,1,] <- -colSums(Q[1,2:4,])
    ## fyke
    Q[2,1,] <- 1 - rl(lvec, L50f, SRf)
    Q[2,2,] <- -Q[2,1,]
    ## corner
    Q[3,4,] <- 1 - rl(lvec, L50c, SRc)
    Q[3,3,] <- -Q[3,4,]
    ## escaped all stay in escaped so all zero intensities
    rho <- apply(Q, 3, function(z){
        Phi <- matrix(c(1, 0, 0, 0), nr = 1) %*% expm::expm(48 * z)
        rho <- Phi[2] / sum(Phi[2:3])
        return(rho)
    })
    pf_hat[i,] <- rho
}

matplot(lvec, t(pf_hat), type = "l")

pred <- data.frame(l = lvec,
                   mean = apply(pf_hat, 2, mean),
                   lwr = apply(pf_hat, 2, quantile, p = 0.025),
                   upr = apply(pf_hat, 2, quantile, p = 0.975))


xlim <- c(15, 60)

png("../tex/sel_paper/figures/FigX_Contact_selection.png", height = 8, width = 6, units = "in", res = 300)
par(mfrow = c(3, 1), mar = c(2, 2, 1, 2), oma = c(2, 2, 1, 2))
## proportion retained
with(enc_lf_cast, plot(len2, pFyke, cex = log(Total + 1), bty = "l", xlim = xlim, xlab = "", ylab = ""))
mtext(side = 2, line = 2.5, text = "Proportion in Fyke")
with(pred, matlines(l, cbind(lwr, mean, upr), lty = c(2, 1, 2), col = 1))
abline(h = 0.5, lty = 3)
legend("topright", legend = c("Mean", "95% CI"), lty = c(1,2), bty = "n", cex = 1.2)
legend("topleft", legend = c("A"), bty = "n", cex = 1.2)
##
hist(enc$Length[enc$Method.Type == "Fyke"], xlim = xlim,
     right = FALSE, breaks = brks, border = "white", main = "")
box(bty = "u")
legend("topleft", legend = "B", bty = "n", cex = 1.2)
mtext(side = 2, text = "Count", line = 2.5)
par(new = TRUE)
##plot(lvec, sel_f, col = "darkslategray4", type = "l", lwd = 1, xlim = xlim, axes = FALSE)
curve(rl(x, L50f, SRf), from = 0, to = 80, n = 1e3, xlim = xlim, axes = FALSE, col = "darkblue", lty = 2)
with(pcgc, lines(l, pfgcf, col = "darkblue"))
## curve(rl(x, L50f_48, SRf_48), col = "red", lty = 2, add = TRUE)
legend("topright", legend = c("Instantaneous", "@48 hrs"),
       lty = c(2, 1), col = c("darkblue", "darkblue"), bty = "n", cex = 1.2)
axis(side = 4)
mtext(side = 4, text = "Selection probability", line = 0.5, outer = TRUE)
##
hist(enc$Length[enc$Method.Type == "Corner"], xlim = xlim,
     right = FALSE, main = "", breaks = brks, border = "white")
legend("topleft", legend = "C", bty = "n", cex = 1.2)
par(new = TRUE)
curve(rl(x, L50c, SRc), from = 0, to = 80, n = 1e3, xlim = xlim, axes = FALSE, col = "darkblue", lty = 2)
##plot(lvec, sel_f, col = "darkslategray4", type = "l", lwd = 1, xlim = xlim, axes = FALSE)
with(pcgc, lines(l, pcgcc, col = "darkblue"))
## curve(rl(x, L50c_48, SRc_48), col = "red", lty = 2, add = TRUE, n = 2e3)
axis(side = 4)
mtext(side = 2, text = "Count", line = 2.5)
mtext(side = 1, text = "Length (cm)", line = 2.5)
legend("topright", legend = c("Instantaneous", "@48 hrs"),
       lty = c(2, 1), col = c("darkblue", "darkblue"), bty = "n", cex = 1.2)
mtext(side = 4, text = "Selection probability", line = 2.5, outer = FALSE)
dev.off()


##---------
## SANDBOX
##---------
## comparison with simulation
par(mfrow = c(2, 1), mar = c(0, 0, 0, 0), oma = c(4, 4, 2, 2))
with(f_res, plot(l, f_capture, xlim = c(20, 35), ylim = c(0, 1), type = "n", bty = "l", xaxt = "n", yaxs = "i", xaxs = "i"))
curve(rl(x, L50f, SRf), add = TRUE, col = "cadetblue", n = 1e3, lty = 2)
with(pcgc, lines(l, pfgcf, col = "cadetblue"))
lines(lvec, sel_f, col = "grey")
legend("bottomright", legend = c("Instantaneous", "Time: 48 hours"),
       lty = c(2, 1), col = "cadetblue", bty = "n")
legend("topleft", legend = "Internal fyke net", bty = "n")
##
with(c_res, plot(l, c_capture, xlim = c(20, 35), ylim = c(0, 1), type = "n", bty = "l", yaxt = "n", yaxs = "i", xaxs = "i"))
axis(2, at = seq(0, 0.8, by = 0.2))
curve(rl(x, L50c, SRc), add = TRUE, col = "indianred1", n = 1e3, lty = 2)
with(pcgc, lines(l, pcgcc, col = "indianred1"))
lines(lbev/10, sel_c, col = "grey")
legend("bottomright", legend = c("Instantaneous", "Time: 48 hours"),
       lty = c(2, 1), col = "indianred1", bty = "n")
legend("topleft", legend = "Corner net", bty = "n")
lines(lbev/10, sel, col = "red")
mtext(side = 1, text = "Length (cm)", line = 2.5)
mtext(side = 2, text = "Selection probability", line = 2.5, outer = TRUE)


## time of contact
par(mfrow = c(2, 1))
hist(subset(res, f_contact)$f_time_contact, breaks = t, main = "Fyke first contact time", probability = TRUE)
##
hist(subset(res, c_contact)$c_time_contact, breaks = t, main = "Corner first contact time", probability = TRUE)

## other data
non_enc <- subset(dat, Method != "Enclosure" & Location %in% c("Feeagh", "Furnace") & Year %in% 2015:2017)

hist(non_enc$Length,
     main = "Fyke surveys", breaks = seq(10, 100), border = "white", xlim = c(10, 100))
par(new = TRUE)
plot(lvec, sel_f, col = "darkslategray4", type = "l", lwd = 1, xlim = c(10, 100), axes = FALSE)
with(pcgc, lines(l, pfgcf, col = "darkblue"))

par(mfrow = c(2, 1))
hist(enc$Length[enc$Method.Type == "Fyke"],
     main = "Enclosure surveys: 2015-2017", breaks = seq(10, 80), border = "white", xlim = c(10, 80), xlab = "Length (cm)")
par(new = TRUE)
plot(lvec, sel_f, col = "darkslategray4", type = "l", lwd = 1, xlim = c(10, 80), axes = FALSE, xlab = "")
with(pcgc, lines(l, pfgcf, col = "darkblue"))

##
hist(non_enc$Length,
     main = "All fyke surveys: 2015-2017", breaks = seq(10, 80), border = "white", xlim = c(10, 80), xlab = "")
par(new = TRUE)
plot(lvec, sel_f, col = "darkslategray4", type = "l", lwd = 1, xlim = c(10, 80), axes = FALSE, xlab = "Length (cm)")
with(pcgc, lines(l, pfgcf, col = "darkblue"))

non_enc <- subset(dat, Method != "Enclosure" & Location %in% c("Feeagh", "Furnace"))

ggplot(non_enc, aes(x = Length)) +
    geom_histogram(, breaks = 0:100) +
    facet_wrap(~Year,ncol = 2, scales = "free_y") +
    geom_vline(xintercept = c(26.21, 32.1), linetype = 2)

enc$floor_len <- floor(enc$Length)

tmp <- hist(subset(enc, Method.Type == "Fyke")$Length, breaks = 0:100, right = FALSE, main = "right = FALSE", ylim = c(0, 20))
x11()
tmp <- hist(subset(enc, Method.Type == "Fyke")$Length, breaks = 0:100, right = TRUE, main = "right = TRUE", ylim = c(0, 20))

aggregate(Sex ~ floor_len + Method.Type, FUN = sum, data = enc)


