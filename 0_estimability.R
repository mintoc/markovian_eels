##------------------------------------------
## Explore how estimable the parameters are
## CM: using Gill and King 2004
##
##------------------------------------------

source("selectivity_functions.R")

library(expm)

## time
T <- 48
nl <- rlnorm(1000, meanlog = log(15), sdlog = 0.4) ## one eel per 10m2
## closest mm
nl <- round(nl, 1)
nl_tab <- table(nl)
l <- as.numeric(names(nl_tab))
n <- as.numeric(nl_tab)
L <- length(l)

## parameters
## contact rate
pc <- 0.01 ## corner
pf <- 0.01 ## fyke
pw <- 0.01 ## wall 
## contact selection
L50c <- L50w <- 15
SRc <- SRw <- 2
L50f <- 20
SRf <- 4

set.seed(102)
dat <- NULL

for(k in 1:L){
    Q <- get_Q(l = l[k], L50f = L50f, L50c = L50c, L50w = L50w, SRf = SRf, SRc = SRc, SRw = SRw, pf = pf, pc = pc, pw = pw)
    P <- expm(T * Q)
    pi_l <- matrix(c(1, 0, 0, 0), nr = 1) %*% P
    nl <- rmultinom(1, size = n[k], prob = pi_l)
    tmp <- data.frame(l = l[k], fyke = nl[2], corner = nl[3])
    dat <- rbind(dat, tmp)
}

dat <- subset(dat, (fyke + corner) > 0)

dat$total <- with(dat, fyke + corner)

dat$pfyke <- with(dat, fyke / total)

with(dat, plot(l, pfyke, cex = log(total + 1), bty = "l"))

##------------
## ESTIMATION 
##------------
## ALL PARAMETERS
##----------------
## starting parameters
start <- c(rep(log(10), 3), rep(log(2), 3), rep(qlogis(0.01), 3))
fit0 <- nlminb(start = start,
              objective = get_sel_nll_all_continuous_faster,
              data = dat,
              T = T)

H <- optimHess(par = fit0$par, fn = get_sel_nll_all_continuous_faster, data = dat, T = T)
## solve(H)
## zeros in columns and rows concerning wall parameters, therefore fix these

##-----------------------
## WALL PARAMETERS FIXED
##-----------------------
## wall parameters fixed at correct values here
start <- c(rep(log(10), 2), rep(log(2), 2), rep(qlogis(0.01), 2))
fixed <- c(L50w = L50w, SRw = SRw, pw = pw)
fit1 <- nlminb(start = start,
              objective = get_sel_nll_wall_continuous_faster,
              data = dat,
              fixed = fixed,
              T = T,
              control = list(eval.max = 1e3, iter.max = 1e3))

H <- optimHess(par = fit1$par, fn = get_sel_nll_wall_continuous_faster, data = dat, fixed = fixed, T = T)

solve(H)

##------------------------------
## WALL PARAMETERS AND PC FIXED
##------------------------------
## wall parameters and pc fixed at correct values here
start <- c(rep(log(10), 2), rep(log(2), 2), 0.01)
fixed <- c(L50w = L50w, SRw = SRw, pw = pw, pc = pc)
fit2 <- nlminb(start = start,
               objective = get_sel_nll_wall_pc_continuous_faster,
               data = dat,
               fixed = fixed,
               T = T,
               control = list(eval.max = 1e3, iter.max = 1e3))

H <- optimHess(par = fit2$par, fn = get_sel_nll_wall_pc_continuous_faster, data = dat, fixed = fixed, T = T)

C <- solve(H)

cov2cor(C)

## check out the jacobian wrt wall pars of the proportion in the fyke over the total in the fyke and corner

##Q <- get_Q(l = 30, L50f = 25, L50c = 20, L50w = 20, SRf = 2, SRc = 3, SRw = 2, pf = 0.1, pc = 0.1, pw = 0.1)

Q <- get_Q(l = 30, L50f = 25, L50c = 20, L50w = 20, SRf = 2, SRc = 2, SRw = 2, pf = 0.1, pc = 0.1, pw = 0.2)
Phi <- matrix(c(1, 0, 0, 0), nr = 1) %*% expm(20 * Q)
(rho <- Phi[2] / sum(Phi[2:3]))
     

##--------------------
## likelihood surface
##--------------------
## wall parameters and pc fixed at correct values here
fixed <- c(L50w = L50w, SRw = SRw, pw = pw, pc = pc)

library(doSNOW)

## grid out the parameters and look at the likelihood surface
m <- 20 ## number of points to evaulate at per parameter
start <- expand.grid(lnL50c = seq(log(10), log(30), length = m),
                     lnL50f = seq(log(10), log(30), length = m),
                     lnSRc = seq(log(1), log(6), length = m),
                     lnSRf = seq(log(1), log(6), length = m),
                     logitpf = seq(qlogis(0.005), qlogis(0.02), length = m))

ll_grid <- rep(NA, nrow(start))

## close if open
close(pb)
stopCluster(cl)

## register 10 sessions
cl <- makeSOCKcluster(10)
registerDoSNOW(cl)

pb <- txtProgressBar(max = nrow(start), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

start$ll <- NA

ll_grid <- foreach(i = 1:nrow(start), .combine = rbind, .options.snow = opts) %dopar% {
    library(expm)    
    ll <- -get_sel_nll_wall_pc_continuous_faster(unlist(start[i,]), data = dat, fixed = fixed, T = T)
    tmp <- data.frame(i = i, ll = ll)
    tmp
}
close(pb)
stopCluster(cl)

##save(list = ls(), file = "profiles.RData")

## check all ordered by i
all(diff(ll_grid$i) == 1)

start$ll <- ll_grid$ll

## inspect control parameters close to maximum
idx <- which.max(start$ll)
theta_hat <- start[idx, ]

pars <- c("lnL50c", "lnSRc", "lnL50f", "lnSRf", "logitpf")

true_df <- data.frame(par = pars,
                      true = c(log(L50c), log(SRc),
                               log(L50f), log(SRf),
                               qlogis(pf)))

## profile log likelihood of each parameter
pdf("../tex/sel_paper/figures/Fig_0_profile_likelihood.pdf", height = 8, width = 7)
par(mfrow = c(3, 2), mar = c(2, 2, 1, 1), oma = c(2, 2, 1, 1))
pars <- c("lnL50c", "lnSRc", "lnL50f", "lnSRf", "logitpf")
for(j in 1:length(pars)){
    tmp <- merge(theta_hat[, pars[pars != pars[j]]], start)
    tmp <- tmp[order(tmp[, pars[j]]), ]
    plot(tmp[, pars[j]], tmp$ll, type = "l", bty = "l", xlab = "", ylab = "", col = "slategrey")
    ##points(tmp[, j], tmp$ll)
    legend("topleft", legend = paste0(letters[j], ". ", pars[j]), bty = "n", cex = 1.3)
    abline(v = subset(true_df, par == pars[j])$true, lty = 2)
}
mtext(side = 2, line = 0.5, text = "Profile log-likelihood", outer = TRUE)
mtext(side = 1, line = 0.5, text = "Parameter value", outer = TRUE)
dev.off()
