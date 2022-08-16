cat_sim = function(harv.prop.init, harv.prop.maint=0.001, iter=100, tt=10){
    #  Kathryn Venning, Corey Bradshaw, Frédérik Saltré
    # Global Ecology, Flinders University — globalecologyflinders.com
    # feral cat reduction on Kangaroo Island
    # requires library - Plotly
    ### update 07/02/2021
    ### modified by James Maino 15/8/2022
    ## update includes: first year fertility, final year survival, predator reduction feedback, removed quasi extinction, previous version 'OFFICIAL cat eradication models GitHub'


    ## function arguments
    # harv.prop.init <- 0.9 # initial harvest rate
    # harv.prop.maint <-  # maintenance harvest rate
    # iter <- 100 # replicate simulations
    # tt <-10 # time horizon


    ## functions
    # beta distribution shape parameter estimator function
    estBetaParams <- function(mu, var) {
    alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
    beta <- alpha * (1 / mu - 1)
    return(params = list(alpha = alpha, beta = beta))
    }

    ## source/matrix operators
    source("matrix_tools.R")

    # create Leslie matrix
    age.max = 7

    ## create vectors 
    #fertility 
    m.vec <- c((0.745/3), 0.745, 2.52, 2.52, 2.52, 2.52, 1.98) ## KI cat birth rates matrix, data for female offsping produced each year. Data from Budke, C & Slater, M (2009)

    # fertility errors based on Budke & Slater
    juv.m.sd <- mean(c(((0.745/3-0.352/3)/2),((1.58/3-0.745/3)/2))) #mean and standard deviations, juvenile fertility
    fy.m.sd <- mean(c(((0.745-0.352)/2),((1.58-0.745)/2))) #mean and standard deviations, juvenile fertility
    A.m.sd <- mean(c(((2.52-1.98)/2),((3.78-2.52)/2))) #mean and standard deviations, adult fertility
    m.sd.vec <- c(0.18*m.vec[1],0.18*m.vec[2],A.m.sd,A.m.sd,A.m.sd,A.m.sd,A.m.sd) #mean and standard deviations vector, juvenile and adult fertility 

    #survival
    s.vec <- c(0.46, 0.46, 0.7, 0.7, 0.7, 0.7) ##KI cat survival # probability of surviving from one year to the next. e.g surviving fourth year of life

    # survival errors based on Budke & Slater
    y1.2.S.sd <- mean(c(((0.46-0.27)/2),((0.73-0.46)/2))) #mean and standard deviations, juvenile survival
    A.S.sd <- mean(c(((0.7-0.55)/2),((0.78-0.7)/2))) #mean and standard deviations, adult survival
    s.sd.vec <- c(y1.2.S.sd,y1.2.S.sd,A.S.sd,A.S.sd,A.S.sd,A.S.sd) #mean and standard deviations vector, juvenile and adult survival

    # create matrix
    popmat <- matrix(data = 0, nrow=age.max, ncol=age.max)
    diag(popmat[2:age.max,]) <- s.vec
    popmat[age.max,age.max] <- 0
    popmat[1,] <- m.vec
    popmat.orig <- popmat ## save original matrix

    ## matrix properties
    max.lambda(popmat) ## 1-yr lambda
    max.r(popmat) # rate of population change, 1-yr
    stable.stage.dist(popmat) ## stable stage distribution
    R.val(popmat, age.max) # reproductive value
    gen.l <- G.val(popmat, age.max) # mean generation length

    ## initial population vector
    pop.found <- 1629 # +/- 661 founding population size Hohnen et al 2020 
    ssd <- stable.stage.dist(popmat)
    init.vec <- ssd * pop.found #initial population vector

    #################
    ## project
    ## set time limit for projection in 1-yr increments
    # yr.now <- 2020 # update if more data available post-2010
    #************************
    # yr.end <- 2030 # set projection end date
    #************************
    # t <- (yr.end - yr.now) #timeframe

    tot.F <- sum(popmat.orig[1,])
    popmat <- popmat.orig #resets matrix 

    ## set population storage matrices
    n.mat <- matrix(0, nrow=age.max,ncol=(tt+1)) #empty matrix
    n.mat[,1] <- init.vec #fill first matrix column with initial population vector

    ## set up projection loop
    for (i in 1:tt) {
    n.mat[,i+1] <- popmat %*% n.mat[,i]
    }

    n.pred <- colSums(n.mat) #number of predators - cats - through time period, no density reduction treatment, no carry capacity
    # yrs <- seq(yr.now, yr.end, 1)
    # plot(yrs, n.pred,type="b",lty=2,pch=19,xlab="year",ylab="N")

    # compensatory density feedback # K = carry capacity
    #population rate of increase relative to carry capacity. Larger distance between populationa and K = faster population growth
    K.max <- 2*pop.found
    K.min <- 1 #not used
    K.vec <- c(1,pop.found/2,pop.found,0.75*K.max,K.max) #1= K.min, .75 = red.thresh??
    red.thresh <- 0.75 #not used
    red.vec <- c(1,0.965,0.89,0.79,0.71)
    plot(K.vec,red.vec,pch=19,type="b")
    Kred.dat <- data.frame(K.vec,red.vec)

    # logistic power function a/(1+(x/b)^c) #fits logistic power function to population relative to carry capacity, K
    param.init <- c(1, 15000, 2.5)
    fit.lp <- nls(red.vec ~ a/(1+(K.vec/b)^c), 
                data = Kred.dat,
                algorithm = "port",
                start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                trace = TRUE,      
                nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
    fit.lp.summ <- summary(fit.lp)
    plot(K.vec,red.vec,pch=19,xlab="N",ylab="reduction factor")
    K.vec.cont <- seq(1,2*pop.found,1)
    pred.lp.fx <- coef(fit.lp)[1]/(1+(K.vec.cont/coef(fit.lp)[2])^coef(fit.lp)[3])
    lines(K.vec.cont,pred.lp.fx,lty=2,lwd=3,col="red")

    a.lp <- coef(fit.lp)[1]
    b.lp <- coef(fit.lp)[2]
    c.lp <- coef(fit.lp)[3]


    # ## compensatory density-feedback deterministic model
    # ## set population storage matrices
    # n.mat <- matrix(0, nrow=age.max,ncol=(t+1))
    # n.mat[,1] <- init.vec
    # popmat <- popmat.orig

    # ## set up projection loop
    # for (i in 1:t) {
    #   totN.i <- sum(n.mat[,i])
    #   pred.red <- a.lp/(1+(totN.i/b.lp)^c.lp)
    #   diag(popmat[2:age.max,]) <- s.vec*pred.red
    #   popmat[age.max,age.max] <- 0
    #   n.mat[,i+1] <- popmat %*% n.mat[,i]
    # }

    # n.pred <- colSums(n.mat)
    # plot(yrs, n.pred,type="b",lty=2,pch=19,xlab="year",ylab="N",ylim=c(0,1.05*K.max)) #untreated population increases, rate of increase relative to K, no stochastic sampling
    # abline(h=K.max, lty=2, col="red") #carry capacity

    #################################################### 
    ## iterations and quasi ext for each following model
    ####################################################
    # iter <- 10000 #final model run at 10 000
    itdiv <- iter/1000 #final model rate at iter/1000


    ##################################################################################################################################################
    ## high harvest for first 2 years, constant proportional harvest in remaining years
    #####################################################################################################################################################

    # storage
    minn.med.mat <- minn.lo.mat <- minn.up.mat <- pmin.med.mat <- pmin.lo.mat <- pmin.up.mat <- matrix(data=NA, ncol=length(harv.prop.maint), nrow=length(harv.prop.init)) #storage matrices

    m <- 1
    n <- 1

    # storage
    n.sums.mat <- p.sums.mat <- matrix(data=NA, nrow=iter, ncol=(tt+1))

    for (e in 1:iter) {
        
        popmat <- popmat.orig
        
        n.mat <- matrix(0, nrow=age.max,ncol=(tt+1))
        n.mat[,1] <- init.vec
        
        for (i in 1:tt) {
        # stochastic survival values
        s.alpha <- estBetaParams(s.vec, s.sd.vec^2)$alpha
        s.beta <- estBetaParams(s.vec, s.sd.vec^2)$beta
        s.stoch <- rbeta(length(s.alpha), s.alpha, s.beta)
        
        # stochastic fertilty sampler (gaussian)
        fert.stch <- rnorm(length(popmat[,1]), popmat[1,], m.sd.vec)
        fert.stoch <- ifelse(fert.stch < 0, 0, fert.stch)
        
        totN.i <- sum(n.mat[,i])
        pred.red <- a.lp/(1+(totN.i/b.lp)^c.lp)
        
        popmat[1,] <- fert.stoch
        diag(popmat[2:age.max,]) <- s.stoch*pred.red
        #popmat[age.max,age.max] <- 0
        
        n.mat[,i+1] <- popmat %*% n.mat[,i]
        
        ## harvest 
        if (i < 3) {
            n.mat[,i+1] <- n.mat[,i+1] - round(stable.stage.dist(popmat) * round(sum(n.mat[,i+1])*harv.prop.init[n], 0), 0)
        } else {
            n.mat[,i+1] <- n.mat[,i+1] - round(stable.stage.dist(popmat) * round(sum(n.mat[,i+1])*harv.prop.maint[m], 0), 0)
        }
        
        if (length(which(n.mat[,i+1] < 0)) > 0) {
            n.mat[which(n.mat[,i+1] < 0), i+1] <- 0
        }
        
        } # end i loop
        
        n.sums.mat[e,] <- as.vector(colSums(n.mat))
        p.sums.mat[e,] <- n.sums.mat[e,] / pop.found
        
        if (e %% itdiv==0) print(e) 
    } # end e loop (stochastic iterations)

    return(n.sums.mat)
}