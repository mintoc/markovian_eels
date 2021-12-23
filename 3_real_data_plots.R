##----------------------------------------------
## Plot the real data
## CM, RP: 2/7/2021
##
##----------------------------------------------

library(gdata)
library(reshape)
library(expm)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)

dat <- read.xls(xls = "../data/Yellow Eel Survey Details.xlsx", sheet = "Total stage and corrections")

enc <- subset(dat, Method == "Enclosure" & Location %in% c("Feeagh", "Furnace"))

## look at the length frequencies
enc$loc2 <- with(enc, paste0(Location, "_", Date))

levels <- c("Feeagh_16/07/2015", "Feeagh_13/08/2015", 
            "Feeagh_22/07/2016", "Feeagh_18/08/2016", 
            "Feeagh_05/07/2017", "Feeagh_17/08/2017",
            "Furnace_11/06/2015", "Furnace_24/07/2015",
            "Furnace_30/06/2016", "Furnace_14/07/2016",
            "Furnace_22/06/2017", "Furnace_14/07/2017"
            )

enc$loc2 <- factor(enc$loc2, levels = levels)
##
enc$count <- 1

pdf("../tex/sel_paper/figures/Fig_S1_location_length.pdf", height = 9, width = 8)
ggplot(enc, aes(x = Length)) +
    geom_histogram(aes(fill = Method.Type), alpha = 0.6, position = "identity") +
    ##geom_density(aes(colour = Method.Type), bw = 3) +
    ##geom_point(aes(colour = Method.Type)) +
    facet_wrap(~loc2, scales = "free_y", ncol = 3) +
    scale_fill_manual("Compartment", values = c("#F8766D", "#619CFF")) +
    scale_colour_manual("Compartment", values = c("#F8766D", "#619CFF")) +
    theme(legend.position = "bottom", panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ylab("Count") +
    xlab("Length (cm)")
dev.off()

## all data pooled
p1l <- ggplot(enc, aes(x = Length)) +
    geom_histogram(aes(fill = Method.Type), alpha = 0.6, position = "identity", binwidth = 1, breaks = seq(10, 75)) +
    scale_fill_manual("Compartment", values = c("#F8766D", "#619CFF")) +
    scale_colour_manual("Compartment", values = c("#F8766D", "#619CFF")) +
    theme(legend.position = c(0.8, 0.8), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ylab("Count") +
    xlab("")

## length proportions plot
enc$len2 <- round(enc$Length)
enc_lf <- aggregate(count ~ len2 + Method.Type, FUN = sum, data = enc)
enc_lf_cast <- cast(enc_lf, len2 ~ Method.Type, value = "count")
enc_lf_cast$Corner[is.na(enc_lf_cast$Corner)] <- 0
enc_lf_cast$Fyke[is.na(enc_lf_cast$Fyke)] <- 0
enc_lf_cast$Total <- with(enc_lf_cast, Corner + Fyke)
enc_lf_cast$pFyke <- with(enc_lf_cast, Fyke / Total)

p2l <- ggplot(enc_lf_cast, aes(x = len2, y = pFyke)) +
    geom_point(aes(size = Total), pch = 1, colour = "slategrey") +
    scale_size_continuous(range = c(2,10)) +
    xlab("") +
    geom_hline(yintercept = 0.5, linetype = 3) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none") +
    ylab("Proportion in Fyke nets")

## length proportions plot
enc$area2 <- round(enc$area_cm2, 1)
enc_tf <- aggregate(count ~ area2 + Method.Type, FUN = sum, data = enc)
enc_tf_cast <- cast(enc_tf, area2 ~ Method.Type, value = "count")
enc_tf_cast$Corner[is.na(enc_tf_cast$Corner)] <- 0
enc_tf_cast$Fyke[is.na(enc_tf_cast$Fyke)] <- 0
enc_tf_cast$Total <- with(enc_tf_cast, Corner + Fyke)
enc_tf_cast$pFyke <- with(enc_tf_cast, Fyke / Total)

##---------------
## include fits!
##---------------
source("selectivity_functions.R")

##---------
## LENGTH
##---------
## get the counts by sample and length by method
## pad out zeros
enc_lf <- aggregate(count ~ Length + loc2 + Method.Type, FUN = sum, data = enc)
enc_cast <- cast(enc_lf, Length + loc2 ~ Method.Type, value = "count")

enc_cast$Corner[is.na(enc_cast$Corner)] <- 0
enc_cast$Fyke[is.na(enc_cast$Fyke)] <- 0

enc_cast$Total <- with(enc_cast, Corner + Fyke)
enc_cast$pFyke <- with(enc_cast, Fyke / Total)

names(enc_cast) <- tolower(names(enc_cast))
names(enc_cast)[names(enc_cast) == "length"] <- "l"

load("length_fit_bb.RData")
theta <- fit_bb$par
vcov <- solve(I_bb)
se <- sqrt(diag(vcov))
## parameter estimate table
res <- data.frame(Par = c("L50c", "L50f", "SRc", "SRf", "pf", "beta"), Est0 = NA, Est1 = NA)
for(i in 1:nrow(res)){
    idx <- grep(res$Par[i], names(theta))
    par_name <- names(theta)[idx]
    thetai <- round(as.numeric(theta[idx]), 2)
    sei <- round(as.numeric(se[par_name]), 3)
    res$Est0[i] <- paste0(thetai, " (", sei, ")")
    ## original scale
    if(res$Par[i] != "pf"){
        lthetai <- as.numeric(theta[idx])
        thetai <- round(exp(lthetai), 2)
        lwr <- round(exp(lthetai - 1.96 * as.numeric(se[idx])), 2)
        upr <- round(exp(lthetai + 1.96 * as.numeric(se[idx])), 2)
        res$Est1[i] <- paste0(thetai, " (", lwr,"-", upr, ")")
    }else{
        lthetai <- as.numeric(theta[idx])
        thetai <- round(plogis(lthetai), 3)
        lwr <- round(plogis(lthetai - 1.96 * as.numeric(se[idx])), 3)
        upr <- round(plogis(lthetai + 1.96 * as.numeric(se[idx])), 3)
        res$Est1[i] <- paste0(thetai, " (", lwr,"-", upr, ")")
    }
    
}

library(xtable)
resx <- xtable(res, row.names = FALSE)

print(resx, include.rownames = FALSE, file = "../tex/sel_paper/par_table0.tex")

## draw parametric samples
library(mvtnorm)

ns <- 1e3
thetas <- rmvnorm(ns, theta, vcov)

lpred <- seq(min(enc_cast$l), max(enc_cast$l), length = 100)

pf_hat <- matrix(NA, nrow = ns, ncol = length(lpred))

Time <- 48

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
    Q <- array(0, dim = c(4, 4, length(lpred)))
    ## free
    Q[1,2,] <- pf * rl(lpred, L50f, SRf)
    Q[1,3,] <- pc * rl(lpred, L50c, SRc)
    Q[1,4,] <- pw * (1 - rl(lpred, L50w, SRw))
    Q[1,1,] <- -colSums(Q[1,2:4,])
    ## fyke
    Q[2,1,] <- 1 - rl(lpred, L50f, SRf)
    Q[2,2,] <- -Q[2,1,]
    ## corner
    Q[3,4,] <- 1 - rl(lpred, L50c, SRc)
    Q[3,3,] <- -Q[3,4,]
    ## escaped all stay in escaped so all zero intensities
    rho <- apply(Q, 3, function(z){
        Phi <- matrix(c(1, 0, 0, 0), nr = 1) %*% expm(Time * z)
        rho <- Phi[2] / sum(Phi[2:3])
        return(rho)
    })
    pf_hat[i,] <- rho
}

matplot(lpred, t(pf_hat), type = "l")

pred <- data.frame(l = lpred,
                   mean = apply(pf_hat, 2, mean),
                   lwr = apply(pf_hat, 2, quantile, p = 0.025),
                   upr = apply(pf_hat, 2, quantile, p = 0.975))


p2l_final <- p2l + geom_line(data = pred, aes(x = l, y = mean)) +
    geom_line(data = pred, aes(x = l, y = upr), linetype = 2) +
    geom_line(data = pred, aes(x = l, y = lwr), linetype = 2)

## contact selection curves
rfl_hat <- matrix(NA, nrow = ns, ncol = length(lpred))
rcl_hat <- matrix(NA, nrow = ns, ncol = length(lpred))

for(i in 1:ns){
    L50c <- exp(thetas[i,"lnL50c"])
    L50f <- exp(thetas[i,"lnL50f"])
    L50w <- fixed["L50w"]
    SRc <- exp(thetas[i,"lnSRc"])
    SRf <- exp(thetas[i,"lnSRf"])
    rfl_hat[i,] <- rl(lpred, L50f, SRf)
    rcl_hat[i,] <- rl(lpred, L50c, SRc)
}

rl_pred <- rbind(
    data.frame(
        Compartment = "Corner",
        l = lpred,
        mean = apply(rcl_hat, 2, mean),
        lwr = apply(rcl_hat, 2, quantile, p = 0.025),
        upr = apply(rcl_hat, 2, quantile, p = 0.975)),
    data.frame(
        Compartment = "Fyke",
        l = lpred,
        mean = apply(rfl_hat, 2, mean),
        lwr = apply(rfl_hat, 2, quantile, p = 0.025),
        upr = apply(rfl_hat, 2, quantile, p = 0.975))    
)
    
p3l <- ggplot(rl_pred, aes(x = l, y = mean, group = Compartment, colour = Compartment)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Compartment), alpha = 0.4, colour = NA) +
    geom_line() +
    ##xlim(15, 45) +
    scale_fill_manual("Compartment", values = c("#F8766D", "#619CFF")) +
    scale_colour_manual("Compartment", values = c("#F8766D", "#619CFF")) +    
    xlab("Length (cm)") +
    ylab("Contact selection probability") +
    theme(legend.position = "none", panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())

pdf("../tex/sel_paper/figures/Fig_3_real_length.pdf", height = 10, width = 6)
grid.arrange(p1l, p2l_final, p3l, ncol = 1)
dev.off()

