##-----------------------
## Simulation experiment
## CM: 6/5/21
## Notes: 
## 
##-----------------------

library(expm)
library(doSNOW)
library(ggplot2); theme_set(theme_bw())
library(viridis)
library(reshape)
library(xtable)

source("selectivity_functions.R")

## simulation experimental design
sim_df <- expand.grid(L50c = seq(20, 30, by = 5),
                      L50f = seq(20, 30, by = 5),
                      L50w = 20,
                      SRc = c(2, 4),
                      SRf = c(2, 4),
                      SRw = 2,
                      pc = c(0.01, 0.02),
                      pf = c(0.01, 0.02),
                      pw = c(0.01, 0.02),
                      iter = 1:50)

############################
## WALL PARS FIXED PARALLEL
############################

## close if open
if("pb" %in% ls()){
    close(pb)
}
if("cl" %in% ls()){
    stopCluster(cl)
}

## register 10 sessions
cl <- makeSOCKcluster(10)
registerDoSNOW(cl)

pb <- txtProgressBar(max = nrow(sim_df), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## fixed parameter values
## change these
fixed <- c(L50w = 20, SRw = 2, pw = 0.01, pc = 0.02)

start <- expand.grid(lnL50c = seq(log(15), log(30), length = 3),
                     lnL50f = seq(log(15), log(30), length = 3),
                     lnSRc = seq(log(0.5), log(4), length = 3),
                     lnSRf = seq(log(0.5), log(4), length = 3),
                     logitpf = seq(log(0.005), log(0.05), length = 3))

## start <- expand.grid(lnL50c = log(25),
##                      lnL50f = log(25),
##                      lnSRc = log(4),
##                      lnSRf = log(4),
##                      logitpf = qlogis(0.02))

Time <- 48 ## hours

## parallel loop
all_sim_res <- foreach(i = 1:nrow(sim_df), .combine = rbind, .options.snow = opts) %dopar% {
##for(i in 1:nrow(sim_df)){
    ##print(i)
    set.seed(i)
    library(expm)
    ##library(DEoptim)
    ## lengths in population
    ##nl <- rlnorm(1000, meanlog = log(20) - 0.4^2/2, sdlog = 0.4) ## one eel per 10m2
    nl <- rlnorm(2000, meanlog = log(20) - 0.4^2/2, sdlog = 0.4) ## one eel per 10m2
    ## closest mm
    nl <- round(nl, 1)
    nl_tab <- table(nl)
    l <- as.numeric(names(nl_tab))
    n <- as.numeric(nl_tab)
    dat <- sim_datc(pars = unlist(sim_df[i,]), n = n, lvec = l, Time = Time)
    ## observed catches only
    dat <- subset(dat, (corner + fyke) > 0)
    ## with(dat, plot(l, fyke/(fyke + corner), cex = fyke + corner))
    ## starting values
    nll_grid <- apply(start, 1, function(x){
        get_sel_nll_wall_pc_continuous_faster(unlist(x), data = dat, fixed = fixed, Time = Time)
    })
    start_best <- start[which.min(nll_grid)[1],]
    fit <- nlminb(start = unlist(start_best),
                  objective = get_sel_nll_wall_pc_continuous_faster,
                  data = dat,
                  lower = c(rep(log(5), 2), rep(log(0.5), 2), log(0.0001)),
                  upper = c(rep(log(40), 2), rep(log(5), 2), log(0.1)),
                  fixed = fixed,
                  Time = Time,
                  control = list(eval.max = 1e3, iter.max = 1e3)
                  )
    ##
    if(fit$convergence == 1){
        theta <- rep(NA, 6)
    }else{
        theta <- fit$par
        names(theta) <- c("lnL50c", "lnL50f", "lnSRc", "lnSRf", "logitpf")
    }
    par_hat <- data.frame(
        i = i,
        L50c_hat = exp(theta["lnL50c"]),
        L50f_hat = exp(theta["lnL50f"]),
        SRc_hat = exp(theta["lnSRc"]),
        SRf_hat = exp(theta["lnSRf"]),
        pf_hat = plogis(theta["logitpf"]))
    res <- cbind(sim_df[i,], par_hat)
    ##
    res
}

close(pb)
stopCluster(cl)

## make sure they line up
all(diff(all_sim_res$i) == 1)

all <- cbind(sim_df, all_sim_res)

## proportion of non-converged
prop.table(table(is.na(rowSums(all_sim_res))))

##save(all, file = "sim_results/sim_fixed_21_07_2021.RData")
load("sim_results/sim_fixed_21_07_2021.RData")
## make a dataframe with id variables and estimated and true values
## difference between Fyke and Corner
all$dL50 <- with(all, L50f - L50c)
all$dSR <- with(all, SRf - SRc)

##idvars <- c("L50c", "L50w", "SRc", "SRw", "pc", "pf", "pw", "iter", "L50f", "SRf")
idvars <- c("L50c", "L50f", "SRc", "SRf", "dL50", "dSR", "pw", "pc", "pf", "iter")

estvars <- c("L50c_hat", "L50f_hat", "SRc_hat", "SRf_hat", "pf_hat")

all_cast <- NULL

for(v in estvars){
    tmp <- all[, idvars]
    tmp$compartment <- ifelse(length(grep("f", v)) > 0, "Fyke", "Corner")
    tmp$par <- gsub("(c|f|_hat)", "", v)
    tmp$true <- all[, gsub("_hat", "", v)]
    tmp$est <- all[, v]
    tmp$error <- with(tmp, est - true)
    tmp$perror <- with(tmp, (est - true)/true)    
    all_cast <- rbind(all_cast, tmp)
    rm(tmp)
}

## both
all_cast$label <- with(all_cast, paste(par, compartment, sep = "\n"))

all_cast$fpc <- factor(all_cast$pc)
all_cast$fpw <- factor(all_cast$pw)
all_cast$fpf <- factor(all_cast$pf)

pdf("../tex/sel_paper/figures/Fig_S0_sim_pc_pw_effects.pdf", height = 10, width = 8)
ggplot(all_cast, aes(x = label, y = perror, fill = interaction(fpc, fpw, sep = ":"))) +
    geom_boxplot() +
    facet_grid(dL50 ~ dSR, labeller = "label_both", scale = "free_y") +
    geom_hline(yintercept = 0, linetype = 2) +
    theme(legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    scale_fill_manual("pc:pw", values = c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c")) +
    ylab("Proportional error (estimate - true)/true") +
    xlab("")
dev.off()

sub_cast <- subset(all_cast, par != "p")

pdf("../tex/sel_paper/figures/Fig_1_simulations.pdf", height = 8, width = 7)
ggplot(sub_cast, aes(x = label, y = perror, fill = label)) +
    geom_boxplot() +
    facet_grid(dL50 ~ dSR, labeller = "label_both", scale = "free_y") +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_fill_manual(values = c("lightgrey", "slategrey", "lightgrey", "slategrey")) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Proportional error (estimate - true)/true") +
    xlab("")
dev.off()



ggplot(all_cast, aes(x = label, y = perror, fill = interaction(fpc, fpw, fpf, sep = ":"))) +
    geom_boxplot() +
    facet_grid(dL50 ~ dSR, labeller = "label_both", scale = "free_y") +
    geom_hline(yintercept = 0, linetype = 2) +
    theme(legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    scale_fill_manual("pc:pw", values = c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c")) +
    ylab("Proportional error (estimate - true)/true") +
    xlab("")

head(all_cast)

rmse <- function(error){
    sqrt(mean(error^2))
}

## RMSE table
tmp0 <- aggregate(error ~ par + compartment + L50c + L50f + SRc + SRf, FUN = rmse, data = subset(all_cast, par != "p"))
names(tmp0)[names(tmp0) == "error"] <- "rmse"
tmp1 <- cast(tmp0, compartment +  par + L50c + SRc ~ L50f + SRf, value = "rmse")

tmp1[,-(1:4)] <- round(tmp1[,-(1:4)], 2)

str(tmp1)
tmp1$L50c <- as.integer(tmp1$L50c)
tmp1$SRc <- as.integer(tmp1$SRc)

rmse_xtab <- xtable(tmp1)

print(rmse_xtab, include.rownames = FALSE, file = "../tex/sel_paper/sim_rmse0.tex")

write.table(tmp1, "../tex/sel_paper/sim_rmse.csv")
