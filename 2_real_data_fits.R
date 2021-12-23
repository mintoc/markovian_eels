##----------------------------------------------
## Fit continuous Markov model to the real data
## CM, RP: 2/7/2021
##
##----------------------------------------------

library(gdata)
library(reshape)
library(expm)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
library(doSNOW)
library(metR)

dat <- read.xls(xls = "../data/Yellow Eel Survey Details.xlsx", sheet = "Total stage and corrections")

enc <- subset(dat, Method == "Enclosure" & Location %in% c("Feeagh", "Furnace"))

## NB! Remove Furnace 14/07/2017 due to jellyfish

enc$loc2 <- with(enc, paste0(Location, "_", Date))

## each row is a fish
enc$count <- 1

source("selectivity_functions.R")
##---------
## LENGTH
##---------
enc_lf <- aggregate(count ~ Length + loc2 + Method.Type, FUN = sum, data = enc)
enc_cast <- cast(enc_lf, Length + loc2 ~ Method.Type, value = "count")
enc_cast$Corner[is.na(enc_cast$Corner)] <- 0
enc_cast$Fyke[is.na(enc_cast$Fyke)] <- 0

enc_cast$Total <- with(enc_cast, Corner + Fyke)
enc_cast$pFyke <- with(enc_cast, Fyke / Total)

ggplot(enc_cast, aes(x = Length, y = pFyke)) +
    geom_point(aes(colour = loc2, size = Total))

names(enc_cast) <- tolower(names(enc_cast))
names(enc_cast)[names(enc_cast) == "length"] <- "l"

## grid of starting values
m <- 20 ## may need to adjust after running but try this out
## if it doesn't work - start likelihood from a few different places and if converges to
## same place then profile from there
start <- expand.grid(lnL50c = seq(log(20), log(30), length = m),
                     lnL50f = seq(log(20), log(30), length = m),
                     lnSRc = seq(log(0.5), log(7), length = m),
                     lnSRf = seq(log(0.5), log(7), length = m),
                     logitpf = seq(qlogis(0.006), qlogis(0.01), length = m))

## fixed parameter values
fixed <- c("L50w" = 25, "SRw" = 1, "pw" = 0.01, "pc" = 0.01)

Time <- 48

run_grid <- FALSE

if(run_grid){
    ## close if open
    close(pb)
    stopCluster(cl)
    ## register 10 sessions
    cl <- makeSOCKcluster(10)
    registerDoSNOW(cl)
    pb <- txtProgressBar(max = nrow(start), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    ll_grid <- foreach(i = 1:nrow(start), .combine = rbind, .options.snow = opts) %dopar% {
        ll <- -get_sel_nll_wall_pc_continuous_faster(unlist(start[i,]), data = enc_cast, fixed = fixed, Time = Time)
        tmp <- cbind(start[i,], ll = ll)
        tmp
    }
    close(pb)
    stopCluster(cl)
    ## save(list = "ll_grid", file = "real_length_profile_15_07_2021.RData")
}
load("real_length_profile_15_07_2021.RData")

idx <- which.max(ll_grid$ll)

ll_grid[idx, ]

## plot bivariate likelihood surfaces
get_surface_vals <- function(pars){
    ## pars <- c("lnL50c", "lnSRc")
    sub <- ll_grid[idx, !names(ll_grid) %in% c(pars, "ll")]
    prof <- merge(ll_grid, sub)
    return(prof)
}

surf_c <- get_surface_vals(c("lnL50c", "lnSRc"))
cutoff <- -500
surf_c[surf_c < cutoff] <- NA
surf_f <- get_surface_vals(c("lnL50f", "lnSRf"))
surf_f[surf_f < cutoff] <- NA

surf_pf <- get_surface_vals(c("logitpf"))
surf_pf <- surf_pf[order(surf_pf$logitpf),]
with(surf_pf, plot(logitpf, ll, type = "l", ylim = c(-258, -255)))


ll_lc <-
    ggplot(surf_c, aes(x = lnL50c, y = lnSRc)) +
    geom_raster(aes(fill = ll)) +
    geom_contour(aes(z = ll), colour = "white", breaks = c(seq(-270, -257, by = 1), -256.6)) +
    geom_text_contour(aes(z = ll), breaks = seq(-270, -257, by = 1), col = "black", nudge_y = 0.01) +
    scale_fill_distiller( palette="Spectral", guide = "colorbar", direction = -1) +
    xlab(expression(ln(L50[c]))) +
    ylab(expression(ln(SR[c])))    

ll_lf <-
    ggplot(surf_f, aes(x = lnL50f, y = lnSRf)) +
    geom_raster(aes(fill = ll)) +
    geom_contour(aes(z = ll), colour = "white", breaks = c(seq(-270, -257, by = 1), -256.6)) +
    geom_text_contour(aes(z = ll), breaks = seq(-270, -257, by = 1), col = "black", nudge_y = 0.01) +
    ##scale_fill_distiller( palette="PuBuGn", guide = "colorbar", direction = 1)
    scale_fill_distiller( palette="Spectral", guide = "colorbar", direction = -1) +
    xlab(expression(ln(L50[f]))) +
    ylab(expression(ln(SR[f])))

pdf("../tex/sel_paper/figures/Fig_2_real_profile.pdf", height = 10, width = 6)
grid.arrange(ll_lc, ll_lf, ncol = 1)
dev.off()

start_best <- ll_grid[which.max(ll_grid$ll),]
start_best <- start_best[, names(start_best) != "ll"]

fit_l <- nlminb(start = unlist(start_best),
              objective = get_sel_nll_wall_pc_continuous_faster,
              data = enc_cast,
              ##lower = c(rep(log(5), 2), rep(log(0.5), 2), log(0.0001)),
              ##upper = c(rep(log(40), 2), rep(log(5), 2), log(0.1)),
              fixed = fixed,
              Time = Time,
              control = list(eval.max = 1e3, iter.max = 1e3)
              )

I_l <- optimHess(fit_l$par, fn = get_sel_nll_wall_pc_continuous_faster, data = enc_cast, fixed = fixed, Time = Time)

solve(I_l)

save(list = c("fit_l", "I_l"), file = "length_fit.RData")

load("length_fit.RData")

## compare fits
## l50s constrained
start <- unlist(start_best)
start["lnL50"] <- as.numeric((start["lnL50c"] + start["lnL50f"])/2)
start <- start[!names(start) %in% c("lnL50c", "lnL50f")]

fit_l50cons <- nlminb(start = start,
                      objective = get_sel_nll_wall_pc_continuous_faster,
                      data = enc_cast,
                      fixed = fixed,
                      Time = Time,
                      constrain = c(l50 = 1, sr = 0),
                      control = list(eval.max = 1e3, iter.max = 1e3)
                      )

## SRs constrained
start <- unlist(start_best)
start["lnSR"] <- as.numeric((start["lnSRc"] + start["lnSRf"])/2)
start <- start[!names(start) %in% c("lnSRc", "lnSRf")]

fit_SRcons <- nlminb(start = start,
                      objective = get_sel_nll_wall_pc_continuous_faster,
                      data = enc_cast,
                      fixed = fixed,
                      Time = Time,
                      constrain = c(l50 = 0, sr = 1),
                      control = list(eval.max = 1e3, iter.max = 1e3)
                      )

## both constrained
start <- unlist(start_best)
start["lnSR"] <- as.numeric((start["lnSRc"] + start["lnSRf"])/2)
start["lnL50"] <- as.numeric((start["lnL50c"] + start["lnL50f"])/2)
start <- start[!names(start) %in% c("lnSRc", "lnSRf", "lnL50c", "lnL50f")]

fit_bothcons <- nlminb(start = start,
                      objective = get_sel_nll_wall_pc_continuous_faster,
                      data = enc_cast,
                      fixed = fixed,
                      Time = Time,
                      constrain = c(l50 = 1, sr = 1),
                      control = list(eval.max = 1e3, iter.max = 1e3)
                      )

## all free versus all constrained
D <- 2 * (-fit_l$objective + fit_bothcons$objective)
pchisq(D, 2, lower.tail = FALSE)
## highly significant difference overall

## SR effect
D <- 2 * (-fit_l$objective + fit_SRcons$objective)
pchisq(D, 1, lower.tail = FALSE)

## L50 effect
D <- 2 * (-fit_l$objective + fit_l50cons$objective)
pchisq(D, 1, lower.tail = FALSE)

## wants to believe the SR differs but the L50s are the same
## but parameters become highly unlikely

##---------------
## BETA-BINOMIAL
##---------------
fit_bb <- nlminb(start = c(unlist(start_best), logbeta = log(2)),
              objective = bb_nll,
              data = enc_cast,
              ##lower = c(rep(log(5), 2), rep(log(0.5), 2), log(0.0001)),
              ##upper = c(rep(log(40), 2), rep(log(5), 2), log(0.1)),
              fixed = fixed,
              Time = Time,
              control = list(eval.max = 1e3, iter.max = 1e3)
              )

I_bb <- optimHess(fit_bb$par, fn = bb_nll, data = enc_cast, fixed = fixed, Time = Time)

sqrt(diag(solve(I_l)))
sqrt(diag(solve(I_bb)))

save(list = c("fit_bb", "I_bb"), file = "length_fit_bb.RData")

## test overdispersion
D <- 2 * (-fit_bb$objective + fit_l$objective)

library(emdbook)
pchibarsq(D, df = 1, mix = 0.5, lower.tail=FALSE)

## BETA-BINOMIAL CONSTRAINTS

## compare fits
## l50s constrained
start <- unlist(start_best)
start["lnL50"] <- as.numeric((start["lnL50c"] + start["lnL50f"])/2)
start <- start[!names(start) %in% c("lnL50c", "lnL50f")]

fit_l50cons <- nlminb(start = c(unlist(start), logbeta = log(2)),
                      objective = bb_nll,
                      data = enc_cast,
                      fixed = fixed,
                      Time = Time,
                      constrain = c(l50 = 1, sr = 0),
                      control = list(eval.max = 1e3, iter.max = 1e3)
                      )

## SRs constrained
start <- unlist(start_best)
start["lnSR"] <- as.numeric((start["lnSRc"] + start["lnSRf"])/2)
start <- start[!names(start) %in% c("lnSRc", "lnSRf")]

fit_SRcons <- nlminb(start = c(unlist(start), logbeta = log(2)),
                      objective = bb_nll,
                      data = enc_cast,
                      fixed = fixed,
                      Time = Time,
                      constrain = c(l50 = 0, sr = 1),
                      control = list(eval.max = 1e3, iter.max = 1e3)
                      )

## both constrained
start <- unlist(start_best)
start["lnSR"] <- as.numeric((start["lnSRc"] + start["lnSRf"])/2)
start["lnL50"] <- as.numeric((start["lnL50c"] + start["lnL50f"])/2)
start <- start[!names(start) %in% c("lnSRc", "lnSRf", "lnL50c", "lnL50f")]

fit_bothcons <- nlminb(start = c(unlist(start), logbeta = log(2)),
                      objective = bb_nll,
                      data = enc_cast,
                      fixed = fixed,
                      Time = Time,
                      constrain = c(l50 = 1, sr = 1),
                      control = list(eval.max = 1e3, iter.max = 1e3)
                      )

aic <- function(mod){
    2 * mod$objective + 2 * length(mod$par)
}

aic(fit_bb)
aic(fit_l50cons)
aic(fit_SRcons)
aic(fit_bothcons)


## all free versus all constrained
D <- 2 * (-fit_bb$objective + fit_bothcons$objective)
pchisq(D, 2, lower.tail = FALSE)
## highly significant difference overall

## SR effect
D <- 2 * (-fit_bb$objective + fit_SRcons$objective)
pchisq(D, 1, lower.tail = FALSE)

## L50 effect
D <- 2 * (-fit_bb$objective + fit_l50cons$objective)
pchisq(D, 1, lower.tail = FALSE)

