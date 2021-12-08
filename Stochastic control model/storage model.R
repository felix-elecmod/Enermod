library(mlOSP)
library(tidyverse)

#'-----------------------------------------------------------------------------
#'@title Storage Optimisation
#'@param model list of model parameters
#'@param input.domain tba
#'@param method choose regression/approximation method
#'@description This function is an adaptation of the swing.fixed.design function
#'from the mlOSP package built by Mike Ludkovski
#'-----------------------------------------------------------------------------
storage.fixed.design <-function (model, input.domain = NULL, method = "km", inTheMoney.thresh = 0){
  # Determine time steps
  M <- model$T/model$dt
  # Determine storage capacity discretisation
  discr.size <- (model$capa.max - model$capa.min)/model$capa.steps
  # Record start time
  t.start <- Sys.time()
  # refraction time - not needed for storage - set to 0
  #refractN <- model$refract/model$dt
  # set up matrix to hold model fits
  # dimensions: time steps x discretised storage levels x regime
  #fits <- array(list(), dim=c(M,model$capa.steps,3))
  #fits <- rep(list(rep(list(),model$capa.steps)),M)
  fits.0 <- matrix(rep(list(),M*model$capa.steps),ncol=model$capa.steps,nrow=M)
  fits.plus <- matrix(rep(list(),M*model$capa.steps),ncol=model$capa.steps,nrow=M)
  fits.minus <- matrix(rep(list(),M*model$capa.steps),ncol=model$capa.steps,nrow=M)
  #fits <- array(list(),dim=3,dimnames = c("Time","Steps","Regime"))
  #fits <- list(matrix(rep(list(),M*model$capa.steps),nrow = M,ncol=model$capa.steps),
              # matrix(rep(list(),M*model$capa.steps),nrow = M,ncol=model$capa.steps),
              # matrix(rep(list(),M*model$capa.steps),nrow = M,ncol=model$capa.steps))
  # list to hold exogenous factor simulations
  grids <- list()
  cur.sim <- 0
  # create pilot simulations on the basis of which we create one step
  # one step simulations at each time steps
  if(model$pilot.nsims > 0) {
    grids[[1]] <- model$sim.func(matrix(rep(model$x0[1:model$dim], 
                                            model$pilot.nsims), nrow = model$pilot.nsims, byrow = T), 
                                 model, model$dt)
    for (i in 2:(M - 1)) grids[[i]] <- model$sim.func(grids[[i - 
                                                               1]], model, model$dt)
    grids[[1]] <- grids[[3]]
    cur.sim <- model$pilot.nsims
  }
  # for each time step we define a design size (number of elements from which
  # to simulate one step ahead)
  design.size <- rep(0, M)
  
  # Start loop backward through time
  for (i in (M - 1):1) {
    # for each time period, we loop through all storage levels
    for (kk in 1:model$capa.steps) {
      # Determine design size for time step i
      if (length(model$N) == 1) {design.size[i] <- model$N}else{design.size[i] <- model$N[i]}
      # determine the number of batches
      if (is.null(model$batch.nrep)) 
        n.reps <- 1
      else if (length(model$batch.nrep) == 1) 
        n.reps <- model$batch.nrep
      else n.reps <- model$batch.nrep[i]
      # determine the input domain
      if (is.null(input.domain)) {
        init.grid <- grids[[i]]
        init.grid <- init.grid[sample(1:min(design.size[i], 
                                            dim(init.grid)[1]), design.size[i], rep = F), 
                               , drop = F]
      }
      else if (length(input.domain) == 2 * model$dim | 
               length(input.domain) == 1) {
        if (length(input.domain) == 1) {
          my.domain <- matrix(rep(0, 2 * model$dim), 
                              ncol = 2)
          if (input.domain > 0) {
            for (jj in 1:model$dim) my.domain[jj, ] <- quantile(grids[[i]][, 
                                                                           jj], c(input.domain, 1 - input.domain))
          }
          else {
            for (jj in 1:model$dim) my.domain[jj, ] <- range(grids[[i]][, 
                                                                        jj])
          }
        }
        else my.domain <- input.domain
        # determine quasiMC method
        if (is.null(model$qmc.method)) {
          init.grid <- tgp::lhs(design.size[i], my.domain)
        }
        else {
          init.grid <- model$qmc.method(design.size[i], 
                                        dim = model$dim)
          for (jj in 1:model$dim) init.grid[, jj] <- my.domain[jj, 
                                                               1] + (my.domain[jj, 2] - my.domain[jj, 1]) * 
              init.grid[, jj]
        }
      }
      else {
        init.grid <- matrix(input.domain, nrow = length(input.domain)/model$dim)
        design.size[i] <- nrow(init.grid)
      }
      # not used in storage setting due to symmetric payoff structure
      #init.grid <- init.grid[model$swing.payoff(init.grid, 
      #                                          model) > inTheMoney.thresh, , drop = F]
      #design.size[i] <- dim(init.grid)[1]
      
      # set up matrix of independent variables
      all.X <- matrix(rep(0, (model$dim + 2) * design.size[i]), 
                      ncol = model$dim + 2)
      # set up big grid of S(t) to base simulations of S(t+1) to S(T) on
      big.grid <- matrix(rep(t(init.grid), n.reps), ncol = ncol(init.grid), 
                         byrow = TRUE)
      # simulate future payoffs depending on chosen regime in t
      
      # determine current storage level
      curr.level = model$capa.min +discr.size*kk
      
        # do nothing
        if(i==(M-1)){
          xtemp <- model$sim.func(big.grid, model, model$dt)
          qVal.inactive <- exp(-model$dt*model$r)*terminal.cond(xtemp,curr.level,(model$capa.max-model$capa.min)/2)
          }else{
        fsim <- storage.policy(big.grid,M-i,fits.0[i:(M-1),],fits.plus[i:(M-1),],fits.minus[i:(M-1),],model,offset = 0,
                               level = curr.level)
        
        qVal.inactive <- 0 + exp(-model$dt*model$r)*fsim$totPayoff
          }
        
        # charge
        new.level <- pmin(curr.level + storage.charging(curr.level,model),model$capa.max-curr.level)
        if(new.level>model$capa.max){qVal.charge = qVal.inactive}else{
          
        if(i==(M-1)){
          qVal.charge <- exp(-model$dt*model$r)*terminal.cond(xtemp,new.level,(model$capa.max-model$capa.min)/2)
        }else{
        fsim <- storage.policy(big.grid,M-i,fits.0[i:(M-1),],fits.plus[i:(M-1),],fits.minus[i:(M-1),],model,offset = 0,
                               level = new.level)
        qVal.charge <- storage.payoff(big.grid,model,regime=1,level = curr.level
                                        )+exp(-model$dt*model$r)*fsim$totPayoff
        }
        }
        #discharge
        new.level <- pmin(curr.level -storage.discharging(curr.level,model),curr.level-model$capa.min)
        if(new.level<model$capa.min){qVal.discharge = qVal.inactive}else{
          
        if(i==(M-1)){
          qVal.discharge <- exp(-model$dt*model$r)*terminal.cond(xtemp,new.level,(model$capa.max-model$capa.min)/2)
        }else{
        fsim <- storage.policy(big.grid,M-i,fits.0[i:(M-1),],fits.plus[i:(M-1),],fits.minus[i:(M-1),],model,offset = 0,
                               level = new.level)
        qVal.discharge <- storage.payoff(big.grid,model,regime=-1,level=curr.level)+exp(-model$dt*model$r)*fsim$totPayoff
        }
        }
        # collect qvalues of different regimes
        qValue <- list(qVal.inactive,qVal.charge,qVal.discharge)
      
        

        #cur.sim <- cur.sim + delayedPayoff$nsims
      #qValue = fsim$totPayoff - immPayoff
      #cur.sim <- cur.sim + fsim$nsims
    
      # estimate qValues for every regime 
      temp <- list()
      for (zz in 1:3){
      
      for (jj in 1:design.size[i]) {
        all.X[jj, model$dim + 1] <- mean(qValue[[zz]][jj + 
                                                  seq(from = 0, len = n.reps, by = design.size[i])])
        all.X[jj, model$dim + 2] <- var(qValue[[zz]][jj + 
                                                 seq(from = 0, len = n.reps, by = design.size[i])])
      }
      all.X[, 1:model$dim] <- init.grid
      
      
      
      
      if (method == "km") 
        temp[[zz]] <- DiceKriging::km(y ~ 0, design = data.frame(x = init.grid), 
                                         response = data.frame(y = all.X[, model$dim + 
                                                                           1]), noise.var = all.X[, model$dim + 2]/n.reps, 
                                         control = list(trace = F), coef.trend = 0, 
                                         coef.cov = model$km.cov, coef.var = model$km.var, 
                                         covtype = model$kernel.family)
      else if (method == "trainkm") 
        temp[[zz]] <- DiceKriging::km(y ~ 0, design = data.frame(x = init.grid), 
                                         response = data.frame(y = all.X[, model$dim + 
                                                                           1]), control = list(trace = F), lower = model$min.lengthscale, 
                                         upper = model$max.lengthscale, noise.var = all.X[, 
                                                                                          model$dim + 2]/n.reps, covtype = model$kernel.family)
      else if (n.reps < 10 & method == "mlegp") 
        temp[[zz]] <- laGP::newGP(X = init.grid, 
                                     Z = all.X[, model$dim + 1], d = list(mle = FALSE, 
                                                                          start = model$km.cov), g = list(start = 1, 
                                                                                                          mle = TRUE))
      else if (method == "hetgp") {
        hetData <- hetGP::find_reps(big.grid, qValue)
        temp[[zz]] <- hetGP::mleHetGP(X = list(X0 = hetData$X0, 
                                                  Z0 = hetData$Z0, mult = hetData$mult), Z = hetData$Z, 
                                         lower = model$min.lengthscale, upper = model$max.lengthscale, 
                                         noiseControl = list(g_bounds = c(1e-04, 100)), 
                                         covtype = model$kernel.family)
      }
      else if (method == "homgp") {
        hetData <- hetGP::find_reps(big.grid, qValue)
        temp[[zz]] <- hetGP::mleHomGP(X = list(X0 = hetData$X0, 
                                                  Z0 = hetData$Z0, mult = hetData$mult), Z = hetData$Z, 
                                         lower = model$min.lengthscale, upper = model$max.lengthscale, 
                                         covtype = model$kernel.family)
      }
      else if (model$dim == 1 & method == "spline") 
        temp[[zz]] <- smooth.spline(x = init.grid, 
                                       y = all.X[, 2], nknots = model$nk)
      else if (method == "rvm") {
        if (is.null(model$rvm.kernel)) 
          rvmk <- "rbfdot"
        else rvmk <- model$rvm.kernel
        temp[[zz]] <- kernlab::rvm(x = init.grid, 
                                      y = all.X[, model$dim + 1], kernel = rvmk)
      }
      else if (method == "npreg") {
        if (is.null(model$np.kertype)) 
          kertype = "gaussian"
        else kertype <- model$np.kertype
        if (is.null(model$np.regtype)) 
          regtype <- "lc"
        else regtype <- model$np.regtype
        if (is.null(model$np.kerorder)) 
          kerorder <- 2
        else kerorder <- model$np.kerorder
        temp[[zz]] <- np::npreg(txdat = init.grid, 
                                   tydat = all.X[, model$dim + 1], regtype = regtype, 
                                   ckertype = kertype, ckerorder = kerorder)
      }
      }
      fits.0[[i,kk]] <- temp[1]
      fits.plus[[i,kk]] <- temp[2]
      fits.minus[[i,kk]] <- temp[3]
    }
  }
  return(list(fit = list(fits.0,fits.plus,fits.minus), timeElapsed = Sys.time() - t.start)) 
              #nsims = cur.sim))
  }



#'-----------------------------------------------------------------------------
#'@title storage forward simulation
#'@description Adaption of swing.policy from MLOSP
#'-----------------------------------------------------------------------------
storage.policy <- function (x, MM, fit.0,fit.plus,fit.minus, model, offset = 1, use.qv = FALSE, level, 
                            verbose = FALSE, undiscounted=FALSE
                            ) 
{
  
  fit <- list(fit.0,fit.plus,fit.minus)
  nsim <- 0
  if (is(x, "matrix") | is(x, "numeric")) {
    curX <- model$sim.func(x, model, model$dt)
    nsim <- nrow(x)
  }
  if (is(x, "list")) 
    curX <- x[[1]]
  
  
  # set up results matrices
  payoff <- matrix(0, nrow=nrow(curX), ncol=MM)
  undiscounted.mat <-  matrix(0, nrow=nrow(curX), ncol=MM)
  action.mat <- matrix(0,nrow=nrow(curX),ncol=MM)
  cur.level <- matrix(0,ncol = MM+1,nrow = nrow(curX))
  cur.level[,1] <- rep(level,nrow(curX))

  i <- 1 # start in period 1
  #zz <- regime # determine running regime in period 0
  kk <- level #determine storage level at the end of period 0

  while (i < (MM + (use.qv == TRUE))) {
  rule <- matrix(0,ncol=3,nrow=nrow(curX))
    for (zz in 1:3){
      myFit <- fit[[zz]][i + 1 - offset,kk]
      myx <- curX[curX, , drop = F]
      if (is(myFit, "earth")) {
        rule[[zz]] <- predict(myFit, myx)
      }
      if (is(myFit, "deepnet")) {
        rule[[zz]] <- nn.predict(myFit, myx)
      }
      if (is(myFit, "nnet")) {
        rule[[zz]] <- predict(myFit, myx, type = "raw")
      }
      if (is(myFit, "smooth.spline")) {
        rule[[zz]] <- predict(myFit, myx)$y
      }
      if (is(myFit, "randomForest")) {
        obj <- predict(myFit, myx, predict.all = T)$individual
        rule[[zz]] <- apply(obj, 1, median)
      }
      if (is(myFit, "dynaTree")) {
        obj <- predict(myFit, myx, quants = F)
        rule[[zz]] <- obj$mean
      }
      if (is(myFit, "km")) {
        rule[[zz]] <- predict(myFit, data.frame(x = myx), 
                        type = "UK")$mean
      }
      if (is(myFit, "gpi")) {
        rule[[zz]] <- predGP(myFit, XX = myx, lite = TRUE)$mean
      }
      if (is(myFit, "lm")) {
        lenn <- length(myFit$coefficients)
        rule[[zz]] <- myFit$coefficients[1] + model$bases(myx) %*% 
          myFit$coefficients[2:lenn]
      }
      if (class(myFit) == "homGP" | class(myFit) == "hetGP" | 
          class(myFit) == "homTP") {
        rule[[zz]] <- predict(x = myx, object = myFit)$mean
      }
      if (class(myFit) == "rvm") {
        rule[[zz]] <- predict(myFit, new = myx)
      }
      if (class(myFit) == "npregression") {
        rule[[zz]] <- predict(myFit, new = myx)
      }
      
    }
      
    # choose optimal action 
      action <- 
        matrix(apply(rule,1,which.max),nrow=nrow(curX)) %>% as_tibble() %>% 
        mutate(V2=case_when(V1==1 ~ 0,
                            V1==2 ~ 1,
                            V1==3 ~ -1)) %>%select(-V1) %>%  as.matrix()
      
    # determine new storage level
      cur.level[,i+1] <- cur.level[,i]+
        (action==1)*storage.charging(cur.level[,i],model)+
        (action==-1)*storage.discharging(cur.level[,i],model)
      
    # impose physical boundary constraints
      cur.level[which(cur.level[,i+1]>model$capa.max),i+1] <- model$capa.max
      cur.level[which(cur.level[,i+1]<model$capa.min),i+1] <- model$capa.min
      
    # adjust actions and store to output matrix
      action[which(cur.level[,i]==cur.level[,i+1])] <- 0
      action.mat[,i] <- action
      
    # compute discounted payoffs  
      payoff[,i] <- exp(-(i) * model$dt * model$r) * storage.payoff(myx,model,
                                                             regime=action,
                                                             level=cur.level[,i])
    # compute undiscounted cash flows
      if(undiscounted){undiscounted.mat[,i] <- storage.payoff(myx,model,
                                                          regime=action,
                                                          level=cur.level[,i])}
    
    i <- i + 1
    #if (verbose == TRUE) 
    #  browser()
    if (is(x, "matrix") | is(x, "numeric")) {
      curX <- model$sim.func(curX, model, model$dt)
    }
    if (is(x, "list")) 
      curX <- x[[i]]
  }
  
  # Add last period
  action.mat[,MM] <- rep(0,nrow(curX)) # no actions in last period
  payoff[,MM] <- exp(-(i) * model$dt * model$r) *terminal.cond(curX,cur.level[,MM-1],(model$capa.max-model$capa.min)/2)
  if(undiscounted){undiscounted.mat[,MM] <- terminal.cond(curX,cur.level[,MM-1],(model$capa.max-model$capa.min)/2)}
  cur.level[,MM+1] <- cur.level[,MM]
  
  
  
  totPayoff = apply(payoff, 1, sum)

  return(list(payoff = payoff,
              action = action.mat,
              totPayoff = totPayoff,
              level = cur.level[,-1],
              if(undiscounted){cashflows = undiscounted.mat}))
}


#'-----------------------------------------------------------------------------
#'@title Storage payoff
#'@param x state vector/matrix
#'@param model storage model parameters
#'@param regime current regime choice
#'@description We define the payoff according to Boogert and De Jong, 2008:
#' "Gas Storage Valuation Using a Monte Carlo Method", Journal of Derivatives
#' We adjust the technical parameters such that they mirror those of an electic
#' battery.
#' In this we assume a bang-bang property
#'-----------------------------------------------------------------------------
storage.payoff <- function (x, model,regime,level){
  if(regime==1){
    # add bid-ask spread as a function of current price
    ask <- (x+model$spread(x)/2)/model$charging.efficiency
    payoff <- -ask*pmin(storage.charging(level,model),model$capa.max-level)
  }else if(regime==-1){
    bid <- (x-model$spread(x)/2)*model$discharging.efficiency
    payoff <- bid*pmin(storage.discharging(level,model),level-model$capa.min)
  }else if(regime==0){
    payoff = 0
  }
  return(payoff)
}


storage.charging <- function(level,model,degr.threshold = 0.1,reduced.rate=0.5){
  if(pmin((model$capa.max-level),(level-model$capa.min)<degr.threshold*(model$capa.max-model$capa.min))){
    temp <- reduced.rate*model$charging.rate}else{temp <- model$charging.rate}
  return(temp)
}

storage.discharging <- function(level,model,degr.threshold = 0.1,reduced.rate=0.5){
  if(pmin((model$capa.max-level),(level-model$capa.min)<degr.threshold*(model$capa.max-model$capa.min))){
    temp <- reduced.rate*model$discharging.rate
  }else{temp <- model$discharging.rate}
  return(temp)
}


terminal.cond <- function(X,ter.level,target.level){
  payoff <- -2*X*abs(ter.level - target.level)
}

sim.gbm2 <- function (x0, model, dt = model$dt) 
{
  len <- nrow(x0)
  if(is.null(len)){len=1}
  newX <- x0
  for (j in 1:model$dim) newX[, j] <- x0[, j, drop = F] * 
    exp(rnorm(len) * model$sigma[j] * sqrt(dt) + (model$r - 
                                                    model$div - model$sigma[j]^2/2) * dt)
  return(newX)
}

