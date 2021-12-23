
## contact selection
rl <- function(l, L50, SR){plogis(2 * log(3)/SR * (l - L50))}

## Markov transition intensity matrix
get_Q <- function(l, L50f, L50c, L50w, SRf, SRc, SRw, pf, pc, pw){
    q <- matrix(0, nr = 4, nc = 4)
    q[1,] <- c(0,
               pf * rl(l, L50f, SRf),
               pc * rl(l, L50c, SRc),
               pw * (1 - rl(l, L50w, SRw)))
    q[2,] <- c(1 - rl(l, L50f, SRf), 0, 0, 0)
    q[3,] <- c(0, 0, 0, 1 - rl(l, L50c, SRc))
    q[4,] <- c(0, 0, 0, 0)
    q[1,1] <- -sum(q[1,])
    q[2,2] <- -sum(q[2,])
    q[3,3] <- -sum(q[3,])
    return(q)
}

## data simulation (c)ontinuous
sim_datc <- function(pars, n, lvec, Time){
    ## n is the number of individuals per length class
    ## lvec are the length classes
    L50c <- pars["L50c"]
    L50f <- pars["L50f"]
    L50w <- pars["L50w"]
    SRc <- pars["SRc"]
    SRf <- pars["SRf"]
    SRw <- pars["SRw"]
    pf <- pars["pf"]
    pc <- pars["pc"]
    pw <- pars["pw"]
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
    dat <- NULL
    Phi <- apply(Q, 3, function(z){
        Phi <- matrix(c(1, 0, 0, 0), nr = 1) %*% expm::expm(Time * z)
        return(Phi)
    })
    for(i in 1:length(lvec)){
        n_comp <- rmultinom(1, size = n[i], prob = Phi[,i])
        dat <- rbind(dat, as.data.frame(cbind(lvec[i], t(n_comp))))
    }
    names(dat) <- c("l", "free", "fyke", "corner", "escape")
    return(dat)
}

get_sel_nll_wall_continuous_faster <- function(theta, data, fixed, T){
    L50c <- exp(theta[1])
    L50f <- exp(theta[2])
    L50w <- fixed["L50w"]
    SRc <- exp(theta[3])
    SRf <- exp(theta[4])
    SRw <- fixed["SRw"]
    pc <- plogis(theta[5])
    pf <- plogis(theta[6])
    pw <- fixed["pw"]
    ll <- 0
    Q <- array(0, dim = c(4, 4, nrow(data)))
    ## free
    Q[1,2,] <- pf * rl(data$l, L50f, SRf)
    Q[1,3,] <- pc * rl(data$l, L50c, SRc)
    Q[1,4,] <- pw * (1 - rl(data$l, L50w, SRw))
    Q[1,1,] <- -colSums(Q[1,2:4,])
    ## fyke
    Q[2,1,] <- 1 - rl(data$l, L50f, SRf)
    Q[2,2,] <- -Q[2,1,]
    ## corner
    Q[3,4,] <- 1 - rl(data$l, L50c, SRc)
    Q[3,3,] <- -Q[3,4,]
    ## free all stay in free so all zero intensities
    rho <- apply(Q, 3, function(z){
        Phi <- matrix(c(1, 0, 0, 0), nr = 1) %*% expm::expm(T * z)
        rho <- Phi[2] / sum(Phi[2:3])
        rho <- min(rho, 1-1e-9) ## to avoid dbinom numerical issues
        rho <- max(rho, 1e-9)
        return(rho)
    })
    ll <- dbinom(x = data$fyke, size = data$fyke + data$corner, prob = rho, log = TRUE)
    return(-sum(ll))
}

get_sel_nll_wall_pc_continuous_faster <- function(theta, data, fixed, Time, constrain = c(l50 = 0, sr = 0)){
    ## L50
    if(constrain["l50"] == 0){
        L50c <- exp(theta["lnL50c"])
        L50f <- exp(theta["lnL50f"])
    }else{
        L50c <- L50f <- exp(theta["lnL50"])     
    }
    ## SR
    if(constrain["sr"] == 0){
        SRc <- exp(theta["lnSRc"])
        SRf <- exp(theta["lnSRf"])
    }else{
        SRc <- SRf <- exp(theta["lnSR"])     
    }
    pf <- plogis(theta["logitpf"])
    ## fixed pars
    L50w <- fixed["L50w"]
    SRw <- fixed["SRw"]
    pc <- fixed["pc"]
    pw <- fixed["pw"]
    ## transition intensity
    Q <- array(0, dim = c(4, 4, nrow(data)))
    ## free
    Q[1,2,] <- pf * rl(data$l, L50f, SRf)
    Q[1,3,] <- pc * rl(data$l, L50c, SRc)
    Q[1,4,] <- pw * (1 - rl(data$l, L50w, SRw))
    Q[1,1,] <- -pf * rl(data$l, L50f, SRf) - pc * rl(data$l, L50c, SRc) - pw * (1 - rl(data$l, L50w, SRw))
    ## fyke
    Q[2,1,] <- 1 - rl(data$l, L50f, SRf)
    Q[2,2,] <- rl(data$l, L50f, SRf) - 1
    ## corner
    Q[3,3,] <- rl(data$l, L50c, SRc) - 1
    Q[3,4,] <- 1 - rl(data$l, L50c, SRc)
    ## free all stay in free so all zero intensities
    rho <- apply(Q, 3, function(z){
        Phi <- matrix(c(1, 0, 0, 0), nr = 1) %*% expm::expm(Time * z)
        rho <- Phi[2] / sum(Phi[2:3])
        rho <- min(rho, 1-1e-9) ## to avoid dbinom numerical issues
        rho <- max(rho, 1e-9)
        return(rho)
    })
    ll <- dbinom(x = data$fyke, size = data$fyke + data$corner, prob = rho, log = TRUE)
    return(-sum(ll))
}


bb_nll <- function(theta, data, fixed, Time, constrain = c(l50 = 0, sr = 0)){
    ## beta-binomial
    ## L50
    if(constrain["l50"] == 0){
        L50c <- exp(theta["lnL50c"])
        L50f <- exp(theta["lnL50f"])
    }else{
        L50c <- L50f <- exp(theta["lnL50"])     
    }
    ## SR
    if(constrain["sr"] == 0){
        SRc <- exp(theta["lnSRc"])
        SRf <- exp(theta["lnSRf"])
    }else{
        SRc <- SRf <- exp(theta["lnSR"])     
    }
    pf <- plogis(theta["logitpf"])
    beta <- exp(theta["logbeta"])
    ## fixed pars
    L50w <- fixed["L50w"]
    SRw <- fixed["SRw"]
    pc <- fixed["pc"]
    pw <- fixed["pw"]
    ## transition intensity
    Q <- array(0, dim = c(4, 4, nrow(data)))
    ## free
    Q[1,2,] <- pf * rl(data$l, L50f, SRf)
    Q[1,3,] <- pc * rl(data$l, L50c, SRc)
    Q[1,4,] <- pw * (1 - rl(data$l, L50w, SRw))
    Q[1,1,] <- -pf * rl(data$l, L50f, SRf) - pc * rl(data$l, L50c, SRc) - pw * (1 - rl(data$l, L50w, SRw))
    ## fyke
    Q[2,1,] <- 1 - rl(data$l, L50f, SRf)
    Q[2,2,] <- rl(data$l, L50f, SRf) - 1
    ## corner
    Q[3,3,] <- rl(data$l, L50c, SRc) - 1
    Q[3,4,] <- 1 - rl(data$l, L50c, SRc)
    ## free all stay in free so all zero intensities
    rho <- apply(Q, 3, function(z){
        Phi <- matrix(c(1, 0, 0, 0), nr = 1) %*% expm::expm(Time * z)
        rho <- Phi[2] / sum(Phi[2:3])
        rho <- min(rho, 1-1e-9) ## to avoid dbinom numerical issues
        rho <- max(rho, 1e-9)
        return(rho)
    })
    alpha <- rho * beta / (1 - rho)
    ki <- data$fyke
    ni <- with(data, fyke + corner)
    ll <- lchoose(ni, ki) + lbeta(ki + alpha, ni - ki + beta) - lbeta(alpha, beta)
    ##ll <- dbinom(x = data$fyke, size = data$fyke + data$corner, prob = rho, log = TRUE)
    return(-sum(ll))
}



get_sel_nll_wall_pc_corner_continuous_faster <- function(theta, data, fixed, Time){
    L50c <- fixed["L50c"]
    L50f <- exp(theta["lnL50f"])
    L50w <- fixed["L50w"]
    SRc <- SRw <- fixed["SRc"]
    SRf <- exp(theta["lnSRf"])
    SRw <- fixed["SRw"]
    pc <- fixed["pc"]
    pf <- plogis(theta["logitpf"])
    pw <- fixed["pw"]
    ll <- 0
    Q <- array(0, dim = c(4, 4, nrow(data)))
    ## free
    Q[1,2,] <- pf * rl(data$l, L50f, SRf)
    Q[1,3,] <- pc * rl(data$l, L50c, SRc)
    Q[1,4,] <- pw * (1 - rl(data$l, L50w, SRw))
    Q[1,1,] <- -pf * rl(data$l, L50f, SRf) - pc * rl(data$l, L50c, SRc) - pw * (1 - rl(data$l, L50w, SRw))
    ## fyke
    Q[2,1,] <- 1 - rl(data$l, L50f, SRf)
    Q[2,2,] <- rl(data$l, L50f, SRf) - 1
    ## corner
    Q[3,3,] <- rl(data$l, L50c, SRc) - 1
    Q[3,4,] <- 1 - rl(data$l, L50c, SRc)
    ## free all stay in free so all zero intensities
    rho <- apply(Q, 3, function(z){
        Phi <- matrix(c(1, 0, 0, 0), nr = 1) %*% expm::expm(Time * z)
        rho <- Phi[2] / sum(Phi[2:3])
        rho <- min(rho, 1-1e-9) ## to avoid dbinom numerical issues
        rho <- max(rho, 1e-9)
        return(rho)
    })
    ll <- dbinom(x = data$fyke, size = data$fyke + data$corner, prob = rho, log = TRUE)
    return(-sum(ll))
}


get_sel_nll_wall_continuous_faster2 <- function(theta, data, fixed, T){
    ##print(as.character(theta))
    ## fixed parameters fixed at given values
    L50c <- exp(theta[1])
    ##L50f <- exp(theta[1]) + exp(theta[2])
    L50f <- exp(theta[2])
    L50w <- fixed["L50w"]
    SRc <- exp(theta[3])
    ##SRf <- exp(theta[3]) + exp(theta[4])
    SRf <- exp(theta[4])
    SRw <- fixed["SRw"]
    pc <- fixed["pc"]
    pf <- plogis(theta[5])
    pw <- fixed["pw"]
    ll <- 0
    Q <- array(0, dim = c(4, 4, nrow(data)))
    ## free
    Q[1,2,] <- pf * rl(data$l, L50f, SRf) ## consolidate rl and rl in one
    Q[1,3,] <- pc * rl(data$l, L50c, SRc)
    Q[1,4,] <- pc * (1 - rl(data$l, L50c, SRc)) + pw * (1 - rl(data$l, L50w, SRw))
    Q[1,1,] <- -colSums(Q[1,2:4,])
    ## fyke
    Q[2,1,] <- 1 - rl(data$l, L50f, SRf)
    Q[2,2,] <- -Q[2,1,]
    ## corner
    Q[3,4,] <- 1 - rl(data$l, L50c, SRc)
    Q[3,3,] <- -Q[3,4,]
    ## free all stay in free so all zero intensities
    rho <- apply(Q, 3, function(z){
        Phi <- matrix(c(1, 0, 0, 0), nr = 1) %*% expm::expm(T * z)
        rho <- Phi[2] / sum(Phi[2:3])
        rho <- min(rho, 1-1e-9) ## to avoid dbinom numerical issues
        rho <- max(rho, 1e-9)
        return(rho)
    })
    ll <- dbinom(x = data$fyke, size = data$fyke + data$corner, prob = rho, log = TRUE)
    return(-sum(ll))
}




get_sel_nll_fixed <- function(theta, data, fixed, T){
    ## fixed parameters fixed at given values
    L50c <- exp(theta[1])
    L50f <- exp(theta[1]) + exp(theta[2])
    L50w <- fixed["L50w"]
    SRc <- exp(theta[3])
    SRf <- exp(theta[3]) + exp(theta[4])
    SRw <- fixed["SRw"]
    pc <- fixed["pc"]
    pf <- plogis(theta[5])
    pw <- fixed["pw"]
    ll <- 0
    for(i in 1:nrow(data)){
        P <- get_markov(data$l[i],
                        L50f = L50f,
                        L50c = L50c,
                        L50w = L50w,
                        SRf = SRf,
                        SRc = SRc,
                        SRw = SRw,
                        pf = pf,
                        pc = pc,
                        pw = pw)
        Phi <- matrix(c(1, 0, 0, 0), nr = 1) %*% (P %^% T)
        rho <- Phi[2] / sum(Phi[2:3])
        lli <- dbinom(x = data$fyke[i], size = data$fyke[i] + data$corner[i], prob = rho, log = TRUE)
        ll <- ll + lli
    }
    return(-ll)
}

get_sel_nll_fixed_continuous <- function(theta, data, fixed, T){
    ## fixed parameters fixed at given values
    L50c <- exp(theta[1])
    L50f <- exp(theta[1]) + exp(theta[2])
    L50w <- fixed["L50w"]
    SRc <- exp(theta[3])
    SRf <- exp(theta[3]) + exp(theta[4])
    SRw <- fixed["SRw"]
    pc <- fixed["pc"]
    pf <- plogis(theta[5])
    pw <- fixed["pw"]
    ll <- 0
    for(i in 1:nrow(data)){
        Q <- get_Q(data$l[i],
                   L50f = L50f,
                   L50c = L50c,
                   L50w = L50w,
                   SRf = SRf,
                   SRc = SRc,
                   SRw = SRw,
                   pf = pf,
                   pc = pc,
                   pw = pw)
        Phi <- matrix(c(1, 0, 0, 0), nr = 1) %*% expm::expm(T * Q)
        rho <- Phi[2] / sum(Phi[2:3])
        lli <- dbinom(x = data$fyke[i], size = data$fyke[i] + data$corner[i], prob = rho, log = TRUE)
        ll <- ll + lli
    }
    return(-ll)
}


get_sel_nll_all_continuous_faster <- function(theta, data, T){
    L50c <- exp(theta[1])
    L50f <- exp(theta[2])
    L50w <- exp(theta[3])
    SRc <- exp(theta[4])
    SRf <- exp(theta[5])
    SRw <- exp(theta[6])
    pc <- plogis(theta[7])
    pf <- plogis(theta[8])
    pw <- plogis(theta[9])
    ll <- 0
    Q <- array(0, dim = c(4, 4, nrow(data)))
    ## free
    Q[1,2,] <- pf * rl(data$l, L50f, SRf)
    Q[1,3,] <- pc * rl(data$l, L50c, SRc)
    Q[1,4,] <- pw * (1 - rl(data$l, L50w, SRw))
    Q[1,1,] <- -colSums(Q[1,2:4,])
    ## fyke
    Q[2,1,] <- 1 - rl(data$l, L50f, SRf)
    Q[2,2,] <- -Q[2,1,]
    ## corner
    Q[3,4,] <- 1 - rl(data$l, L50c, SRc)
    Q[3,3,] <- -Q[3,4,]
    ## free all stay in free so all zero intensities
    rho <- apply(Q, 3, function(z){
        Phi <- matrix(c(1, 0, 0, 0), nr = 1) %*% expm::expm(T * z)
        rho <- Phi[2] / sum(Phi[2:3])
        rho <- min(rho, 1-1e-9) ## to avoid dbinom numerical issues
        rho <- max(rho, 1e-9)
        return(rho)
    })
    ll <- dbinom(x = data$fyke, size = data$fyke + data$corner, prob = rho, log = TRUE)
    return(-sum(ll))
}
