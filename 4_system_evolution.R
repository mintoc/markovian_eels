##----------------------------
## Plots of how the estimated system evolves over time
## CM:
##
##----------------------------
library(reshape)
library(ggplot2); theme_set(theme_bw())
source("selectivity_functions.R")

## time
T <- 48
t <- seq(0, T, by = 6)
l <- seq(0.001, 60, length = 300)
L <- length(l)

load("length_fit_bb.RData")
theta <- fit_bb$par

## visualise
## parameters
fixed <- c("L50w" = 25, "SRw" = 1, "pw" = 0.01, "pc" = 0.01)

## contact probability
## note corner and internal contact probs worked out below to give 39% caught after 48 hours
##pc <- as.numeric(fixed["pc"]) ## corner
##pf <- as.numeric(plogis(theta["logitpf"])) ## fyke
pf <- 0.003989465
pc <- 1.581258 * pf
pw <- as.numeric(fixed["pw"]) ## wall



## contact selection
L50c <- as.numeric(exp(theta["lnL50c"]))
L50f <- as.numeric(exp(theta["lnL50f"]))
L50w <- as.numeric(fixed["L50w"]) 
SRc <- as.numeric(exp(theta["lnSRc"]))
SRf <- as.numeric(exp(theta["lnSRf"]))
SRw <- as.numeric(fixed["SRw"])


phi_t <- NULL

for(j in t){
    Q <- array(0, dim = c(4, 4, length(l)))
    ## free
    Q[1,2,] <- pf * rl(l, L50f, SRf)
    Q[1,3,] <- pc * rl(l, L50c, SRc)
    Q[1,4,] <- pw * (1 - rl(l, L50w, SRw))
    Q[1,1,] <- -pf * rl(l, L50f, SRf) - pc * rl(l, L50c, SRc) - pw * (1 - rl(l, L50w, SRw))
    ## fyke
    Q[2,1,] <- 1 - rl(l, L50f, SRf)
    Q[2,2,] <- rl(l, L50f, SRf) - 1
    ## corner
    Q[3,3,] <- rl(l, L50c, SRc) - 1
    Q[3,4,] <- 1 - rl(l, L50c, SRc)
    ## free all stay in free so all zero intensities
    Phi <- apply(Q, 3, function(z){
        Phi <- matrix(c(1, 0, 0, 0), nr = 1) %*% expm::expm(j * z)
        return(Phi)
    })
    Phi <- as.data.frame(t(Phi))
    names(Phi) <- c("Free", "Fyke", "Corner", "Escape")
    Phi$l <- l
    Phi$Time <- j
    phi_t <- rbind(phi_t, Phi)
}

phi_long <- melt(phi_t, id.vars = c("l", "Time"))


pdf("../tex/sel_paper/figures/Fig_4_system_evolution.pdf", height = 10, width = 7)
ggplot(phi_long, aes(x = l, y = value, fill = variable)) + 
    geom_area(colour="slategrey", size = 0.15) +
    facet_wrap(~ Time, ncol = 2, labeller = "label_both") +
    scale_fill_manual("Compartment", values = c("cornsilk", "lightgreen", "cadetblue1", "cornflowerblue")) +
    coord_cartesian(xlim = c(0, 60), ylim = c(0, 1)) +
    theme(legend.position = c(0.75, 0.1), legend.direction = "vertical") +
    xlab("Length (cm)") +
    ylab("Proportion")
dev.off()


t48 <- subset(phi_t, Time == 48)


with(t48, plot(l, Fyke, bty = "l", col = "#619CFF", type = "l", lwd = 1.5, xlim = c(20, 40)))
with(t48, lines(l, Corner, bty = "l", col = "#F8766D"))


## get the contact rates so that 39% caught in total after 48 hours (Dorow 2020)
pc <- as.numeric(fixed["pc"]) ## corner
pf <- as.numeric(plogis(theta["logitpf"])) ## fyke

scalar <- pc/pf
getprop <- function(pftest){
    pctest <- scalar * pftest
    Q <- matrix(0, nrow = 4, ncol = 4)
    ## free
    l <- 40
    Q[1,2] <- pftest * rl(l, L50f, SRf)
    Q[1,3] <- pctest * rl(l, L50c, SRc)
    Q[1,4] <- pw * (1 - rl(l, L50w, SRw))
    Q[1,1] <- -pftest * rl(l, L50f, SRf) - pctest * rl(l, L50c, SRc) - pw * (1 - rl(l, L50w, SRw))
    ## fyke
    Q[2,1] <- 1 - rl(l, L50f, SRf)
    Q[2,2] <- rl(l, L50f, SRf) - 1
    ## corner
    Q[3,3] <- rl(l, L50c, SRc) - 1
    Q[3,4] <- 1 - rl(l, L50c, SRc)
    ## free all stay in free so all zero intensities
    Phi <- matrix(c(1, 0, 0, 0), nr = 1) %*% expm::expm(48 * Q)
    (sum(Phi[1,2:3]) - 0.39)^2
}

fit <- optim(par = 0.1, fn = getprop, method = "Brent", lower = 0, upper = 1)




