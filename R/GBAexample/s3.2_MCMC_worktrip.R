#### Setting up functions ####
source("./s1.0_packages.R", echo=F)
sourceCpp("cpp_gibbs.cpp")


RGibbs <- function(N, burnin, thn, draw, targets, vars.int, condex, conds.vec) {
  print("Running Gibbs sampler")
  # Setup receiving matrix
  mcmc_mat <- matrix(integer(), nrow = N+burnin, ncol = length(draw), byrow=T, dimnames = list(NULL, names(vars.int)))
  # Full Gibbs Sampler Loop
  for(i in 1:(N+burnin)) {
    if((i/(N+burnin)) %% 0.05 == 0) cat(paste(round(100*i/(N+burnin)),"%|",sep=""))
    for(t in thn) {
      #Going through the conditionals
      it=0
      stop=F
      while(it < 100){
        for(c in 1:length(conds)){
          #Going through the target variables within the conditional
          for(j in 1:length(targets[[c]])) {
            #Making a "spread", varying the target variable to extract probabilities with 'matkey()'
            spread = draw[condex[[c]]]
            dims = var.dims[condex[[c]]]
            sprdtarg = which(targets[[c]][j] == condex[[c]])
            keys = sapply(vars.int[[ targets[[c]][j] ]], function(x) {
              spread[sprdtarg] = x
              return(matkey(spread,dims)+1)
            })
            #Check if stuck in dead end, if so, go back to last good draw and break from within conditional
            if(sum(conds.vec[[c]][keys]) == 0 ) {
              stop = T
              draw = mcmc_mat[i-1,]
              break
            } else {
              #Choose a new value for target variable
              draw[ targets[[c]][j] ] = sample(vars.int[[ targets[[c]][j] ]], 1, prob = conds.vec[[c]][keys])
            }
          }
          if(stop == T) break
        }
        if(stop == F & c == length(conds)) break
        it = it + 1
      }
    }
    if(stop==F) mcmc_mat[i,] <- draw
  }
  mcmc_mat <- mcmc_mat[-(1:burnin),]
  return(mcmc_mat)
}
RMH <- function(N, burnin, thn, draw, targets, vars.int, condex, conds.vec) {
  print("Running Metropolis Hastings")
  # Setup receiving matrix
  mcmc_mat <- matrix(integer(), nrow = N+burnin, ncol = length(draw), byrow=T, dimnames = list(NULL, names(vars.int)))
  # Full Gibbs Sampler Loop
  for(i in 1:(N+burnin)) {
    if((i/(N+burnin)) %% 0.05 == 0) cat(paste(round(100*i/(N+burnin)),"%|",sep=""))
    stop = F
    for(t in thn) {
      #Going through the conditionals
      for(c in 1:length(conds)){
        #Going through the target variables within the conditional
        for(j in 1:length(targets[[c]])) {
          #Making a "spread", varying the target variable to extract probabilities with 'matkey()'
          spread = draw[condex[[c]]]
          dims = var.dims[condex[[c]]]
          sprdtarg = which(targets[[c]][j] == condex[[c]])
          keys = sapply(vars.int[[ targets[[c]][j] ]], function(x) {
            spread[sprdtarg] = x
            return(matkey(spread,dims)+1)
          })
          #Check if stuck in dead end, if so, go back to last good draw and break from within conditional
          if(sum(conds.vec[[c]][keys]) == 0 ) {
            i = i - 2
            stop = T
            break
          } else {
            for(k in 1:1000) {
              #Randomly choose a new value for target variable
              draw[ targets[[c]][j] ] = sample(vars.int[[ targets[[c]][j] ]], 1, prob = conds.vec[[c]][keys])
              if(i<burnin) break
              #freq <- data.table(as.character(vars.int[[ targets[[c]][j] ]]))
              #freq <- merge(freq, data.table(table(mcmc_mat[ , targets[[c]][j] ])), by = 'V1')
              lr_old <- MLR(mcmc_mat[ , targets[[c]][j] ], conds.vec[[c]][keys], vars.int[[ targets[[c]][j] ]])
              lr_new <- MLR(c(mcmc_mat[ , targets[[c]][j] ], draw[ targets[[c]][j] ]), conds.vec[[c]][keys], vars.int[[ targets[[c]][j] ]])
              a = lr_new / lr_old
              if( a >= 1) {
                break
              } else {
                if( sample(c(0,1), 1, TRUE, c(1-a, a)) == 1 ) break 
              }
            }
          }
        }
        #Break outer conditional loop too
        if(stop == T) break
      }
      #Break from thinning loop too
      if(stop == T) break
    }
    mcmc_mat[i,] <- draw
  }
  mcmc_mat <- mcmc_mat[-(1:burnin),]
  return(mcmc_mat)
}
fun_MH <- function(N=10000, burnin=1000, thn=1, starts, vars, targets, condex, conds, lang='C') {
  #Setting up input parameters
  vars.fac = lapply(vars, function(x) factor(x, levels = x))
  vars.int = lapply(vars.fac, function(x) as.integer(x))
  var.dims = sapply(vars.int, function(x) length(x))
  conds.vec = lapply(conds, as.vector)
  condex.0 = lapply(condex, function(x) x-1)
  targets.0 = lapply(targets, function(x) x-1)
  
  #Convert characters into integers
  for(n in names(vars.fac)) starts[ , (n) := factor(starts[[n]], levels = vars.fac[[n]]),]
  starts.int = starts[ , lapply(.SD, as.integer)]
  starts = as.matrix(starts.int)[1, ]
  
  #Check consistency
  for(i in 1:length(conds)) {
    for( j in 1:length(condex[[i]]))
      if( !(dim(conds[[i]])[j] == var.dims[condex[[i]]][j]) ) {
        stop(paste("Conditionals on", names(var.dims[condex[[i]]][j]), "are not consistent, or incorrect variables supplied"))
      }
  }
  
  #Which language?
  if( lang == 'R' ) {
      mcmc_mat <- RMH(N, burnin, thn, starts, targets, vars.int, condex, conds.vec)
  } else if( lang == 'C') {
      mcmc_mat <- RcppMH(N, burnin, thn, starts, targets.0, vars.int, condex.0, conds.vec)
  } else {
    print("Invalid coding specified, use 'R' or 'C'")
  }
  colnames(mcmc_mat) <- names(draw)
  mcmc_mat <- data.table(mcmc_mat)
  #Convert integers back into characters
  for(n in names(vars.fac)) {
    mcmc_mat[ , (n) := factor(mcmc_mat[[n]], levels = as.integer(vars.fac[[n]]), labels = vars.fac[[n]]),]
  }
  
  
  return(combomat)
}
fun_gibbs <- function(N=10000, burnin=1000, thn=1, starts, vars, targets, condex, conds, lang='R') {
  #Setting up input parameters
  vars.fac = lapply(vars, function(x) factor(x, levels = x))
  vars.int = lapply(vars.fac, function(x) as.integer(x))
  var.dims = sapply(vars.int, function(x) length(x))
  conds.vec = lapply(conds, as.vector)
  condex.0 = lapply(condex, function(x) x-1)
  targets.0 = lapply(targets, function(x) x-1)
  
  #Convert characters into integers
  start = as.data.table(starts)
  for(n in names(vars.fac)) start[ , (n) := factor(starts[[n]], levels = vars.fac[[n]]),]
  starts.int = start[ , lapply(.SD, as.integer)]
  start = as.matrix(starts.int)[1, ]
  
  #Check consistency
  for(i in 1:length(conds)) {
    for( j in 1:length(condex[[i]]))
      if( !(dim(conds[[i]])[j] == var.dims[condex[[i]]][j]) ) {
        stop(paste("Conditionals on", names(var.dims[condex[[i]]][j]), "are not consistent, or incorrect variables supplied"))
      }
  }
  
  print("Running Gibbs sampler")
  #Which language?
  if( lang == 'R' ) {
    mcmc_mat <- RGibbs(N, burnin, thn, start, targets, vars.int, condex, conds.vec)
  } else if( lang == 'C') {
    mcmc_mat <- RcppGibbs(N, burnin, thn, start, targets.0, vars.int, condex.0, conds.vec)
  } else {
    print("Invalid coding specified, use 'R' or 'C'")
  }
  colnames(mcmc_mat) <- names(start)
  mcmc_mat <- data.table(mcmc_mat)
  #Convert integers back into characters
  for(n in names(vars.fac)) {
    mcmc_mat[ , (n) := factor(mcmc_mat[[n]], levels = as.integer(vars.fac[[n]]), labels = vars.fac[[n]]),]
  }
 
  return(mcmc_mat)
}
fun_gibbsmulti <- function(N=10000, burnin=1000, thn=1, starts, vars, targets, condex, conds, lang='R') {
  #Setting up input parameters
  vars.fac = lapply(vars, function(x) factor(x, levels = x))
  vars.int = lapply(vars.fac, function(x) as.integer(x))
  var.dims = sapply(vars.int, function(x) length(x))
  conds.vec = lapply(conds, as.vector)
  condex.0 = lapply(condex, function(x) x-1)
  targets.0 = lapply(targets, function(x) x-1)
  
  #Convert characters into integers
  starts = as.data.table(starts)
  for(n in names(vars.fac)) starts[ , (n) := factor(starts[[n]], levels = vars.fac[[n]]),]
  starts.int = starts[ , lapply(.SD, as.integer)]
  
  #Check consistency
  for(i in 1:length(conds)) {
    for( j in 1:length(condex[[i]]))
      if( !(dim(conds[[i]])[j] == var.dims[condex[[i]]][j]) ) {
        stop(paste("Conditionals on", names(var.dims[condex[[i]]][j]), "are not consistent, or incorrect variables supplied"))
      }
  }
  
  #Which language?
  if( lang == 'R' ) {
    # Parallel start points
    combomat = rbindlist(lapply(1:nrow(starts), function(d) {
      draw = as.matrix(starts.int)[d, ]
      # Setup receiving matrix
      mcmc_mat <- matrix(integer(), nrow = N+burnin, ncol = length(draw), byrow=T, dimnames = list(NULL, names(vars)))
      print(paste("Running Gibbs sampler on start #",d,"of",nrow(starts)))
      # Full Gibbs Sampler Loop
      for(i in 1:(N+burnin)) {
        if((i/(N+burnin)) %% 0.05 == 0) cat(paste(round(100*i/(N+burnin)),"%|",sep=""))
        stop = F
        for(t in thn) {
          #Going through the conditionals
          for(c in length(conds):1){
            #Going through the target variables within the conditional
            for(j in 1:length(targets[[c]])) {
              #Making a "spread", varying the target variable to extract probabilities with 'matkey()'
              spread = draw[condex[[c]]]
              dims = var.dims[condex[[c]]]
              sprdtarg = which(targets[[c]][j] == condex[[c]])
              keys = sapply(vars.int[[ targets[[c]][j] ]], function(x) {
                spread[sprdtarg] = x
                return(matkey(spread,dims)+1)
              })
              #Check if stuck in dead end, if so, go back to last good draw and break from within conditional
              if(sum(conds.vec[[c]][keys]) == 0 ) {
                i = i - 2
                stop = T
                break
              } else {
                #Randomly choose a new value for target variable
                draw[ targets[[c]][j] ] = sample(vars.int[[ targets[[c]][j] ]], 1, prob = conds.vec[[c]][keys])
              }
            }
            #Break outer conditional loop too
            if(stop == T) break
          }
          #Break from thinning loop too
          if(stop == T) break
        }
        mcmc_mat[i,] <- draw
      }
      print("")
      mcmc_mat <- mcmc_mat[-(1:burnin),]
      mcmc_mat <- data.table(mcmc_mat,"draw" = d)
      #Convert integers back into characters
      for(n in names(vars.fac)) {
        mcmc_mat[ , (n) := factor(mcmc_mat[[n]], levels = as.integer(vars.fac[[n]]), labels = vars.fac[[n]]),]
      }
      return(mcmc_mat)
    }))
  } else if( lang == 'C') {
    # Parallel start points in C
    combomat = rbindlist(lapply(1:nrow(starts.int), function(d) {
      draw = as.matrix(starts.int)[d, ]
      print(paste("Running Gibbs sampler on start #",d,"of",nrow(starts.int)))
      mcmc_mat <- RcppGibbs(N, burnin, thn, draw, targets.0, vars.int, condex.0, conds.vec)
      colnames(mcmc_mat) <- names(draw)
      mcmc_mat <- data.table(mcmc_mat,"draw" = d)
      #Convert integers back into characters
      for(n in names(vars.fac)) {
        mcmc_mat[ , (n) := factor(mcmc_mat[[n]], levels = as.integer(vars.fac[[n]]), labels = vars.fac[[n]]),]
      }
      return(mcmc_mat)
    }))
  } else {
    print("Invalid coding specified, use 'R' or 'C'")
  }
  
  
  return(combomat)
}
fun_gibbsleapfrog <- function(N=100000, warmup=10000, chnsize=1000, starts, vars, targets, condex, conds, marg, micro) {
 
  #Setting up input parameters
  vars.fac = lapply(vars, function(x) factor(x, levels = x))
  vars.int = lapply(vars.fac, function(x) as.integer(x))
  var.dims = sapply(vars.int, function(x) length(x))
  conds.vec = lapply(conds, as.vector)
  condex.0 = lapply(condex, function(x) x-1)
  targets.0 = lapply(targets, function(x) x-1)
  burn = ceiling((chnsize^2)/N)
  srmse = 999
  
  #Check consistency
  for(i in 1:length(conds)) {
    for( j in 1:length(condex[[i]]))
      if( !(dim(conds[[i]])[j] == var.dims[condex[[i]]][j]) ) {
        stop(paste("Conditionals on", names(var.dims[condex[[i]]][j]), "are not consistent, or incorrect variables supplied"))
      }
  }
  
  mcmc_mat <- starts[0,]
  
  # #Convert characters into integers
  # starts.int = rbind(sapply(names(vars.fac), function(n) as.integer(factor(starts[[n]][1], levels = vars.fac[[n]]))),
  #                    sapply(names(vars.fac), function(n) as.integer(factor(starts[[n]][2], levels = vars.fac[[n]]))))
  # 
  # # Run gibbs sampler
  # mcmc_mat <- rbind(RcppGibbs(warmup, burn, thn=1, starts.int[1,], targets.0, vars.int, condex.0, conds.vec),
  #                   RcppGibbs(warmup, burn, thn=1, starts.int[2,], targets.0, vars.int, condex.0, conds.vec))
  # 
  # #Formatting
  # mcmc_mat = as.data.table(mcmc_mat)
  # colnames(mcmc_mat) <- names(starts)
  # 
  # # Convert integers back into characters
  # for(n in names(vars.fac)) {
  #   mcmc_mat[ , (n) := factor(mcmc_mat[[n]], levels = as.integer(vars.fac[[n]]), labels = vars.fac[[n]]),]
  # }
  # 
  # # Check results
  # vn <- names(mcmc_mat[,!c("otract","dtract")])
  # #Tally up frequencies
  # res = rbindlist(lapply(vn, function(n) setNames(mcmc_mat[,.N, by = c('otract',n) ], c("otract","Variable","MCMC"))))
  # res[ , MCMC := MCMC/sum(MCMC)]
  # #Merging with marginals to compare
  # res = merge(marg, res, by = c('otract','Variable'), all.x = T)
  # res[is.na(MCMC), MCMC := 0]
  # 
  # #Check improvement
  # srmse = res[ , sqrt(sum((Marginals - MCMC)^2)/.N)]/mean(res$Marginals)
  # print(paste("First SRMSE = ",srmse))
  # 
  # burn = ceiling((chnsize^2)/N)
  # #Which var?
  # var = res[which.max(abs(Marginals - MCMC)),.(Var,Variable,otract)]$Var
  # #Get the ID of best one
  # idx = which(micro[[var]] == res[which.max(abs(Marginals - MCMC)),]$Variable)
  # starts = micro[idx, ][which.max(N),!"N"]
  # #Add OD
  # od = setNames(as.data.table(melt(conds$ODI[-1,-1,])),c("otract","dtract","INDUS","value"))
  # od = od[value>0 & INDUS == starts$INDUS,][which.max(value),.(otract,dtract)]
  # starts = cbind(starts, od)
  
  #Main loop
  for(i in 1:(N/chnsize)) {
    #Convert characters into integers
    starts.int = sapply(names(vars.fac), function(n) as.integer(factor(starts[[n]], levels = vars.fac[[n]])))
    
    # Run gibbs sampler
    mcmc_draw = RcppGibbs(chnsize, burn, thn=1, starts.int, targets.0, vars.int, condex.0, conds.vec)
    
    #Formatting
    mcmc_draw = as.data.table(mcmc_draw)
    colnames(mcmc_draw) <- names(starts)
    
    # Convert integers back into characters
    for(n in names(vars.fac)) {
      mcmc_draw[ , (n) := factor(mcmc_draw[[n]], levels = as.integer(vars.fac[[n]]), labels = vars.fac[[n]]),]
    }
    
    # Adding to end of simulated data
    mcmc_mat = rbind(mcmc_mat, mcmc_draw)
    
    # Check results
    vn <- names(mcmc_mat[,!c("otract","dtract")])
    #Tally up frequencies
    res <- rbindlist(lapply(vn, function(n) setNames(mcmc_mat[,.N, by = c('otract',n) ], c("otract","Variable","MCMC"))))
    res[ , MCMC := MCMC/sum(MCMC)]
    #Merging with marginals to compare
    res <- merge(marg, res, by = c('otract','Variable'), all.x = T)
    res[is.na(MCMC), MCMC := 0]
    
    #Check improvement
    srmse_new = res[ , sqrt(sum((Marginals - MCMC)^2)/.N)]/mean(res$Marginals)
    
    if((srmse - srmse_new) > 1e-2) {
      print(paste("Continuing Gibbs chain", i, "of", (N/chnsize),"at SRMSE =",srmse))
      starts = as.matrix(mcmc_mat[nrow(mcmc_mat),])[1,]
      burn = 0
    } else {
      print(paste("Hopping Gibbs chain", i, "of", (N/chnsize),"at SRMSE =",srmse))
      burn = ceiling((chnsize^2)/N)
      #Which var?
      var = res[which.max(abs(Marginals - MCMC)),.(Var,Variable,otract)]$Var
      #Get the ID of best one
      idx = which(micro[[var]] == res[which.max(abs(Marginals - MCMC)),]$Variable)
      starts = micro[idx, ][which.max(N),!"N"]
      
      #OD
      odi <- as.data.table(prop.table(table(mcmc_mat[,.(otract,dtract,INDUS)])))
      lodes <- setNames(as.data.table(conds$ODI)[V2!="00000000000" & value>0,],c("otract","dtract","INDUS","LODES"))
      odi <- merge(odi,lodes, by = c("otract","dtract","INDUS"))
      
      #Add OD
      od <- odi[INDUS==starts$INDUS,][which.max(abs(N-LODES)),.(otract,dtract)]
      starts <- cbind(starts, od)
      
      # #Marginal plot w/ tracts
      # ggplot(data=res) +
      #   geom_point(aes(x=Marginals, y=MCMC, color=Variable), alpha = 0.1) +
      #   theme_bw() + coord_fixed()
      # 
      # #OD plot
      # ggplot(data=odi[sample(nrow(odi),10000,T),]) + 
      #   geom_point(aes(x=N, y=LODES, color=INDUS), alpha = 0.1) +
      #   theme_bw() + coord_fixed()
      
    }
    srmse = srmse_new
    #Looping back over...
  }
  
  return(mcmc_mat)
}
#Raking function
fun_raketract.singledim <- function(tr, freq.dat, marg.dat, bnd=c(0,Inf), trm=NULL, forcestrict=F, skip=F) {
  #Pulling the marginal
  nm <- names(marg.dat)
  levs.og <- lapply(nm, function(n) levels(marg.dat[[n]][[n]]))
  names(levs.og) <- nm
  
  #Pulling the marginal & data for current tract
  marg.dat <- lapply(marg.dat, function(x) x[otract == tr, !"otract"])
  if('otract' %in% colnames(freq.dat)) freq.dat <- freq.dat[otract == tr, !"otract"]
  
  #Removing margin zeros from data
  zeros <- lapply(marg.dat, function(x) x[Freq==0, ])
  for(n in names(zeros)) freq.dat[['N']][ freq.dat[[n]] %in% zeros[[n]][[n]] ] <- 0
  freq.dat <- freq.dat[N>0,]
  
  #Removing the zeros from the margins
  marg.dat <- lapply(marg.dat, function(x) x[Freq!=0, ])
  
  #Relevel the margin factors
  for(n in names(marg.dat)) {
    levs = levels(marg.dat[[n]][[n]])[levels(marg.dat[[n]][[n]]) %in% marg.dat[[n]][[n]]]
    marg.dat[[n]][[n]] <- factor(marg.dat[[n]][[n]], levels = levs)
    }
  #Matching factors in data to margin factors
  for(n in names(marg.dat)) freq.dat[[n]] <- factor(freq.dat[[n]], levels = levels(marg.dat[[n]][[n]]))
  
  #Scaling the weights
  total <- sum(marg.dat[[1]]$Freq)
  colnames(freq.dat)[which("N" == colnames(freq.dat))] <- "weights.sim"
  freq.dat[ , weights.sim := total*weights.sim/sum(weights.sim)]
  
  #Check if only 1 factor level remaining, aka only one household type
  singles <- sapply(marg.dat, function(x) ifelse(nrow(x[Freq>0,])<=1, x[[1]],NA))
  
  if(sum(is.na(singles)) > 0 ) {
    #Removing the singles temporarily
    freq.dat.stripped <- freq.dat[ , c(names(singles[is.na(singles)]),"weights.sim"), with=F]
    #Setting up the raking data
    design.dat <- svydesign(ids = ~1, weights=~weights.sim, data=freq.dat.stripped)
    #Formula
    form <- as.formula(paste("~",paste(names(singles[is.na(singles)]),collapse = "+"),sep=""))
    #Population marginals
    pop <- c(`(Intercept)` = total, unlist(lapply(names(singles[is.na(singles)]), function(n) {
      setNames(marg.dat[[n]]$Freq[-1], paste(n,marg.dat[[n]][[n]][-1],sep=""))
    })))
    raked.dat = 'character'
    #Raking the data, strictest standards
    if(skip==F) {
      raked.dat <- suppressWarnings(
        tryCatch(
          calibrate(design.dat, formula = form, population = pop, calfun = "raking", bounds = bnd, trim = trm, force=forcestrict),
                 error = function(e) "error"))
    }
    #Try again with relaxed constraints
    if(all(class(raked.dat) == "character")) {
      print("Relaxed boundaries")
      raked.dat <- suppressWarnings(
        tryCatch(
          calibrate(design.dat, formula = form, population = pop, calfun = "linear", bounds = c(0,1000), trim=trm, eps=1, force=T),
                 error = function(e) "error"))
    }
    
    #Putting adding the weights back in
    if(all(class(raked.dat) == "character")) {
      freq.dat.raked <- cbind("otract" = tr, freq.dat, "weights.raked" = freq.dat$weights.sim)
    } else {
      freq.dat.raked <- cbind("otract" = tr, freq.dat, "weights.raked" = weights(raked.dat))
    }
  } else {
    #No change possible, just return the original, with eps = 0
    freq.dat.raked <- cbind("otract" = tr, freq.dat, "weights.raked" = freq.dat$weights.sim)
  }
  
  #Convert back to original factors
  for(n in nm) freq.dat.raked[[n]] <- as.character(freq.dat.raked[[n]])
  for(n in names(marg.dat)) freq.dat.raked[[n]] <- factor(freq.dat.raked[[n]], levels = levs.og[[n]])
  #Out
  return(freq.dat.raked)
}

#Load base data
set.seed(12345)
load("./db/d2.3_microdata.RData")
load("./db/d2.5_margins.RData")


#### Gibbs sample the origin-destination LODES to create conditional for individuals ####
if(!file.exists("./db/d3.2_MCMClodes.RData")) {
  
  load("./db/d2.4_lodes.RData")
  #### Setup conditionals for LODES ####
  mtxOI <- lodes$mtxOI
  mtxDI <- lodes$mtxDI
  mtxOD <- lodes$mtxOD
  
  #Adjust the Origins to match marginals
  #cbind(margins$indindus[tracts,'NATRES'],mtxOI[tracts,'NATRES'])
  mtxOI[tracts,] <- (margins$indindus[tracts,-1] + mtxOI[tracts,])/2
  
  ### Need to run Gibbs sampler without NOIND first to build conditional for population
  conds = list("OI" = prop.table(mtxOI),
               "DI" = prop.table(mtxDI),
               "OD" = prop.table(mtxOD))
  
  condex = list(c(1,3),c(2,3),c(1,2))
  targets = list(3,2,1)
  vars = list("otract" = dimnames(mtxOI)$otract,
              "dtract" = dimnames(mtxOD)$dtract,
              "INDUS" = dimnames(mtxDI)$INDUS)
  
  indus = melt(as.data.table(mtxDI[-1,],keep.rownames = T), id.vars='rn')[value>0, ][which.min(value),variable]
  origin = melt(as.data.table(mtxOI[-1,],keep.rownames = T), id.vars='rn')[variable == indus & value>0, ][which.min(value),rn]
  dest = melt(as.data.table(mtxOD[-1,-1],keep.rownames = T), id.vars='rn')[rn == origin & value>0, ][which.min(value),rn]

  starts = data.table('otract' = origin, 'dtract' = dest, 'INDUS' = indus)[,lapply(.SD, as.character)]
  
  
  #Estimating joint probabilities using Gibbs sampler in C++ with Rcpp
  mcmc_odi <- fun_gibbs(N=1e8, burnin=10000, thn=1, starts, vars, targets, condex, conds, lang='C')
  
  #### Setup LODES Gibbs output ####
  #Making a matrix
  mtxODI <- prop.table(table(mcmc_odi[,!"draw"]))
  #Making a NOIND matrix and setting 0s and 1s accordingly
  noind.ratio <- as.vector(colSums(margins$indindus)['NOIND']/sum(colSums(margins$indindus)[-1]))
  noind.mtx <- diag(c("00000000000" = sum(mtxOD[1,])*noind.ratio, margins$indindus[,'NOIND']))
  noind.mtx <- prop.table(noind.mtx)
  mcmc.odi <- prop.table(abind("NOIND" = noind.mtx, mtxODI, along=3))
  names(dimnames(mcmc.odi)) <- c("otract","dtract","INDUS")
  
  #Getting proportions
  res.od = as.data.table(prop.table(table(mcmc_odi[,.(otract,dtract)])))
  res.oi = as.data.table(prop.table(table(mcmc_odi[,.(otract,INDUS)])))
  res.di = as.data.table(prop.table(table(mcmc_odi[,.(dtract,INDUS)])))
  res.oi[ , INDUS := as.character(INDUS)]
  res.di[ , INDUS := as.character(INDUS)]
  
  dat.od = setNames(melt(as.data.table(prop.table(mtxOD),keep.rownames = T), id.vars='rn'),c("otract","dtract","LODES"))[,  dtract := as.character(dtract)]
  dat.oi = setNames(melt(as.data.table(prop.table(mtxOI),keep.rownames = T), id.vars='rn'),c("otract","INDUS","LODES"))
  dat.di = setNames(melt(as.data.table(prop.table(mtxDI),keep.rownames = T), id.vars='rn'),c("dtract","INDUS","LODES"))
  dat.oi[ , INDUS := as.character(INDUS)]
  dat.di[ , INDUS := as.character(INDUS)]
  
  #merging proportions to gibbs results
  res.od = merge(res.od, dat.od, by=c("otract","dtract"), all=T)
  res.oi = merge(res.oi, dat.oi, by=c("otract","INDUS"), all=T)
  res.di = merge(res.di, dat.di, by=c("dtract","INDUS"), all=T)
  #Checking for zeros
  # res.od = res.od[ , (c('N','LODES')) := lapply(.SD, function(x) ifelse(is.na(x),0,x)), .SDcols = c('N','LODES')]
  # res.oi = res.oi[ , (c('N','LODES')) := lapply(.SD, function(x) ifelse(is.na(x),0,x)), .SDcols = c('N','LODES')]
  # res.di = res.di[ , (c('N','LODES')) := lapply(.SD, function(x) ifelse(is.na(x),0,x)), .SDcols = c('N','LODES')]
    
  #### Output ####
  mcmc.odi.res <- list("OD" = res.od, "OI" = res.oi, "DI" = res.di)
  
  #Adding tiny value for NATRES b/c it's bullshit
  
  # #Finding origins with no possible destination
  # missing <- lapply(dimnames(mcmc.odi)$INDUS, function(x) which(rowSums(mcmc.odi[,,x])==0))
  # names(missing) <- dimnames(mcmc.odi)$INDUS
  # 
  # #The missing industry
  # for(indus in names(missing)) {
  #   if(length(missing[[indus]])>0) {
  #     #The origin tracts with no destinations
  #     for(otract in missing[[indus]]) {
  #       #Possible destinations
  #       dtract = mtxOD[otract,][mtxOD[otract,]>0]/sum(mtxOD)
  #       #Pushing values in
  #       mcmc.odi[otract,names(dtract),indus] <- dtract
  #     }
  #   }
  # }
  # #Readjusting to 1.
  # mcmc.odi = prop.table(mcmc.odi)
  
  #Check results
  ggplot() +
    geom_point(data = mcmc.odi.res$OD[sample(mcmc.odi.res$OD[,.I],10000,replace = F), ], aes(x=N, y=LODES, color = "OD"), alpha=0.5, size=0.1) +
    #geom_point(data = mcmc.odi.res$OD, aes(x=N, y=LODES, color = "OD")) +
    geom_point(data = mcmc.odi.res$OI, aes(x=N, y=LODES, color = "OI"), alpha=0.5, size=0.1) +
    geom_point(data = mcmc.odi.res$DI, aes(x=N, y=LODES, color = "DI"), alpha=0.5, size=0.1) +
    #xlim(0,1e-3) + ylim(0,1e-3) +
    theme_bw() + coord_fixed()
  
  mcmc.odi.res$OD[LODES>0 & N==0,]
  mcmc.odi.res$DI[LODES>0 & N==0,]
  mcmc.odi.res$OI[LODES>0 & N==0,]
  mcmc.odi.res$OD[LODES==0 & N>0,]
  mcmc.odi.res$DI[LODES==0 & N>0,]
  mcmc.odi.res$OI[LODES==0 & N>0,]
  
  #
  save(mcmc_odi, file = "./db/d3.2_MCMClodes_raw.RData")
  save(mcmc.odi, mcmc.odi.res, file = "./db/d3.2_MCMClodes.RData")
  #rm(list=setdiff(ls(),c("fun_gibbs","matkey","RcppGibbs","tracts","tables","microdata","margins","LODES")))
  #
  rm(list = grepl("mcmc.",ls()))
}

#### Gibbs sample Individuals with OD ####
if(!file.exists("./db/d3.2_MCMCindod.RData")) {
  
  load("./db/d2.4_lodes.RData")
  load("./db/d3.2_MCMClodes.RData")
  
  #### Setup LODES conditionals ####
  mtxOI <- lodes$mtxOI[c("00000000000",tracts),]
  mtxDI <- lodes$mtxDI[c("00000000000",tracts),]
  mtxOD <- lodes$mtxOD[c("00000000000",tracts),c("00000000000",tracts)]
  
  #OI
  #Average out the lodes OI with the marginal OI
  mtxOI[tracts,] <- (margins$indindus[tracts,-1] + mtxOI[tracts,])/2
  #Determine proportion of unemployed outside region, approximately
  noind.ratio <- colSums(margins$indindus)[1] / sum(colSums(margins$indindus)[-1])
  mtxOI <- prop.table(cbind("NOIND" = c("00000000000" = as.vector(sum(mtxOI[1,]) * noind.ratio), margins$indindus[tracts,'NOIND']), mtxOI))
  names(dimnames(mtxOI)) <- c("otract","INDUS")
  
  #DI
  mtxDI <- prop.table(cbind( "NOIND" = c("00000000000" = as.vector(sum(mtxDI[1,]) * noind.ratio), margins$indindus[tracts,'NOIND']), mtxDI))
  names(dimnames(mtxDI)) <- c("dtract","INDUS")
  
  #OD
  mtxOD <- prop.table(diag(c("00000000000"=0, prop.table(margins$indindus[tracts,'NOIND']))) + prop.table(mtxOD))
  mtxOD <- mtxOD[tracts,c("00000000000",tracts)]
  
  #### MCMC for individuals with OD ####
  # Population conditional data
  data.ind = microdata$pums[,.(AGE,GEND,RELATE,OCC,INDUS)]
  pjoint.int = prop.table(table(data.ind[ , lapply(.SD,as.integer)]))
  
  #Setup variables
  vars = c(dimnames(prop.table(table(data.ind))),
           list("otract" = c("00000000000",tracts),
                "dtract" = c("00000000000",tracts)))
  
  #Starts. Initial values, we know the two isolated optima!!
  starts <- rbind(as.data.table(table(data.ind))[N>0 & INDUS != "NOIND",][which.max(N),!"N"],
                  as.data.table(table(data.ind))[N>0 & INDUS == "NOIND",][which.max(N),!"N"])
  #Add OD
  od.occ <- setNames(as.data.table(melt(mcmc.odi[-1,-1,])),c("otract","dtract","INDUS","value"))
  od.occ <- od.occ[value>0 & INDUS == starts[1,INDUS],][which.max(value),.(otract,dtract)]
  
  od.noocc <- setNames(as.data.table(melt(mcmc.odi[-1,-1,])),c("otract","dtract","INDUS","value"))
  od.noocc <- od.noocc[value>0 & INDUS == starts[2,INDUS],][which.max(value),.(otract,dtract)]
  
  #Combining
  starts <- cbind(starts, rbind(od.occ,od.noocc))
  
  #Conditionals
  conds <- list("IND" = pjoint.int,
                "OI" = mtxOI[c("00000000000",tracts),],
                "ODI" = mcmc.odi[c("00000000000",tracts),c("00000000000",tracts),])
  
  #age.gend.relate.occ.Indus -> origin.dest.industry
  condex = list(c(1,2,3,4,5), c(6,5), c(6,7,5))
  targets = list(c(1,2,3,4,5), 6, 7)
  
  #Marginals for checking
  marg <- list("AGE" = apply(margins$indagesex, c(1,3), sum),
               "GEND" = apply(margins$indagesex, c(1,2), sum),
               "RELATE" = as.matrix(margins$indrelate),
               "OCC"  = margins$indocc,
               "INDUS"  = margins$indindus)
  
  marg <- rbindlist(lapply(names(marg), function(n) {
    out <- setNames(melt(marg[[n]]), c("otract","Variable","Marginals"))
    out$Var <- n
    return(out)
  }))
  
  marg[ , otract := factor(otract, levels = tracts)]
  marg[, Marginals := Marginals/sum(Marginals)]
  
  #Microdata for checking
  micro <- as.data.table(prop.table(table(data.ind)))
  
  #Generate individuals with OD
  ratio = marg[Var == 'OCC' & Variable != 'NOOCC',sum(Marginals)] / marg[Var == 'OCC',sum(Marginals)]
  
  # mcmc.indod <- fun_gibbsleapfrog(N=1e5, warmup=1e4, chnsize=1e4, starts[1,], vars, targets, condex, conds, marg, micro)
  # mcmc.occindod <- fun_MH(N=10000*ratio, burnin=1000, thn=1, starts[1,], vars, targets, condex, conds, lang='R')
  # mcmc.nooccindod <- fun_MH(N=10000*(1-ratio), burnin=1000, thn=1, starts[2,], vars, targets, condex, conds, lang='R')
  n = 1e9
  
  #Run and Combine
  mcmc.indod <- rbind(fun_gibbs(N=n*ratio, burnin=1000, thn=1, starts[1,], vars, targets, condex, conds, lang='C'),
                      fun_gibbs(N=n*(1-ratio), burnin=1000, thn=1, starts[2,], vars, targets, condex, conds, lang='C'))
  
  # Check results ####
  #results list
  mcmc.indod.res <- list()
  
  #Microdata results
  mcmc.indod.res$micro <- merge(as.data.table(prop.table(table(data.ind))),
                                as.data.table(prop.table(table(mcmc.indod[,!c("otract","dtract")]))),
                                by = colnames(micro[,!"N"]))
  setnames(mcmc.indod.res$micro, c("N.x","N.y"), c("PUMS","MCMC"))
  #Marginal results
  vn <- names(mcmc.indod[,!c("otract","dtract")])
  #Tally up frequencies
  mcmc.indod.res$marg <- rbindlist(lapply(vn, function(n) setNames(mcmc.indod[,.N, by = c('otract',n) ], c("otract","Variable","MCMC"))))
  mcmc.indod.res$marg[ , MCMC := MCMC/sum(MCMC)]
  #Merging with marginals to compare
  mcmc.indod.res$marg <- merge(marg, mcmc.indod.res$marg, by = c('otract','Variable'), all.x = T)
  mcmc.indod.res$marg[is.na(MCMC), MCMC := 0]
  
  #OD
  odi <- as.data.table(prop.table(table(mcmc.indod[,.(otract,dtract,INDUS)])))
  lodes <- setNames(as.data.table(conds$ODI)[otract!="00000000000" & value>0,],c("otract","dtract","INDUS","LODES"))
  odi <- merge(odi,lodes, by = c("otract","dtract","INDUS"))
  mcmc.indod.res$ODI <- odi
  rm(odi,vn)
  
  # Plot results ####
  #Microdata plot
  ggplot(data = mcmc.indod.res$micro) + geom_point(aes(x=MCMC, y=PUMS)) + theme_bw() + coord_fixed()
  
  #Marginal plot
  ggplot(data=mcmc.indod.res$marg[ , lapply(.SD,sum), .SDcols = c("MCMC","Marginals"), by = Variable]) + 
    geom_point(aes(x=Marginals, y=MCMC, color=Variable)) +
    theme_bw() + coord_fixed() + theme(legend.position="none")  
  
  #Marginal plot w/ tracts
  ggplot(data=mcmc.indod.res$marg) + 
    geom_point(aes(x=Marginals, y=MCMC, color=Var), alpha = 0.1) +
    theme_bw() + coord_fixed() #+ theme(legend.position="none")  
  
  
  mcmc.indod.res$OD <- mcmc.indod.res$ODI[,lapply(.SD,sum), .SDcols = c("N","LODES"), by = .(otract,dtract)]
  mcmc.indod.res$DI <- mcmc.indod.res$ODI[,lapply(.SD,sum), .SDcols = c("N","LODES"), by = .(otract,INDUS)]
  mcmc.indod.res$OI <- mcmc.indod.res$ODI[,lapply(.SD,sum), .SDcols = c("N","LODES"), by = .(dtract,INDUS)]
  
  #ODI plot
  ggplot() + 
    geom_point(data=mcmc.indod.res$OD[sample(mcmc.indod.res$OD[,.I],10000,replace = F), ], aes(x=N, y=LODES, color="OD"), alpha=0.5) +
    geom_point(data=mcmc.indod.res$DI, aes(x=N, y=LODES, color="DI"), alpha=0.5) +
    geom_point(data=mcmc.indod.res$OI, aes(x=N, y=LODES, color="OI"), alpha=0.5) +
    theme_bw() + coord_fixed()
  
  #Output
  save(mcmc.indod, mcmc.indod.res, file = "./db/d3.2_MCMCindod.RData")
  #cleanup
  rm(list = ls()[grepl("mcmc.",ls())])
}

#### Gibbs sample Individuals ####
if(!file.exists("./db/d3.2_MCMCind.RData")) {
  # Gibbs individual ####
  data.ind = microdata$pums[ ,.(AGE,GEND,RELATE,OCC,INDUS)]
  #Conditional data
  pjoint.int = prop.table(table(data.ind[ , lapply(.SD,as.integer)]))
  #vars = c(dimnames(prop.table(table(data.ind))), LODES$ODIvars[c('otract','dtract')])
  vars = dimnames(prop.table(table(data.ind)))
  
  # Initial values, we know the two isolated optima!!
  starts <- rbind(as.data.table(table(data.ind))[OCC=="NOOCC" & INDUS == "NOIND",][which.max(N),!"N"],
                  as.data.table(table(data.ind)[,,,dimnames(table(data.ind))$OCC[-1],dimnames(table(data.ind))$INDUS[-1]])[which.max(N),!"N"])
  
  #conditionals
  conds = list(pjoint.int)
  condex = list(c(1,2,3,4,5))
  targets = list(c(1,2,3,4,5))
  
  #Generate individuals with OD
  mcmc.ind <- fun_gibbs(N=2500000, burnin=10000, thn=1, starts, vars, targets, condex, conds, lang='C')[,!"draw"]
  
  # Check results ####
  mcmc.ind.res <- list()
  #Check microdata results
  mcmc.ind.res$micro = merge(as.data.table(prop.table(table(mcmc.ind[,.(AGE,GEND,RELATE,OCC,INDUS)]))),
                             as.data.table(prop.table(table(data.ind))), by = colnames(data.ind), all = T, suffixes = c("Gibbs","Data"))
  mcmc.ind.res$micro = mcmc.ind.res$micro[ , lapply(.SD, function(x) ifelse(is.na(x),0,x))]
  ggplot(data = mcmc.ind.res$micro) + geom_point(aes(x=NGibbs, y=NData)) + theme_bw() + coord_fixed()
  
  #Check marginal results
  mcmc.ind.res$marg <- list("AGE" = apply(margins$indagesex, c(1,3), sum),
                            "GEND" = apply(margins$indagesex, c(1,2), sum),
                            "RELATE" = as.matrix(margins$indrelate),
                            "OCC"  = margins$indocc,
                            "INDUS"  = margins$indindus)
  mcmc.ind.res$marg <- rbindlist(lapply(names(mcmc.ind.res$marg), function(n) {
    setNames(as.data.table(colSums(mcmc.ind.res$marg[[n]]), keep.rownames = T),c("Variable","Marginals"))
  }))
  #aggregate data
  res.mcmc.ind <- rbindlist(lapply(names(mcmc.ind), function(n) {
    setNames(as.data.table(table(mcmc.ind[,.(AGE,GEND,RELATE,OCC,INDUS)]))[ , sum(N), by = n], c("Variable","MCMC"))
  }))
  
  #Merging together and melting
  mcmc.ind.res$marg <- merge(res.mcmc.ind, mcmc.ind.res$marg, by = 'Variable', all = T)
  mcmc.ind.res$marg[ , lapply(.SD,sum), .SDcols = c("MCMC","Marginals"), by = Variable]
  mcmc.ind.res$marg[ , (c("MCMC","Marginals")) := lapply(.SD,function(x) x/sum(x)), .SDcols = c("MCMC","Marginals")]
  
  ggplot(data=mcmc.ind.res$marg[ , lapply(.SD,sum), .SDcols = c("MCMC","Marginals"), by = Variable]) + 
    geom_point(aes(x=Marginals, y=MCMC, color=Variable)) +
    theme_bw() + coord_fixed() + theme(legend.position="none")
  
  rm(list = setdiff(ls(),c("mcmc.ind","mcmc.ind.res","mcmc.odi","mcmc.odi.res",
                           "margins","microdata","tables","tracts","fun_gibbs",
                           "fun_raketract.singledim","fun_rakebulk.singledim","fun_raketract.multidim",
                           "matkey","RcppGibbs")))
  
  rm(list = ls()[grepl("mcmc.",ls())])
  
  #Output
  save(mcmc.ind, mcmc.ind.res, file = "./db/d3.2_MCMCind.RData")
}

#### Gibbs sample Households ####
if(!file.exists("./db/d3.2_MCMChh.RData")) {
  # Gibbs household ####
  data.hh = microdata$pums[!duplicated(SERIALNO),.(INCOME,RESTY,HHVEH,HHSIZ)]
  
  #Conditional data
  pjoint.int = prop.table(table(data.hh[ , lapply(.SD,as.integer)]))
  
  # Initial values, we know the two isolated optima!!
  starts = as.data.table(table(data.hh))[which.max(N),!"N"]
  
  #conditionals
  conds = list(pjoint.int)
  condex = list(c(1,2,3,4))
  targets = list(c(1,2,3,4))
  vars = dimnames(prop.table(table(data.hh)))

  #Generate individuals with OD
  mcmc.hh <- fun_gibbs(N=2500000, burnin=1000, thn=1, starts, vars, targets, condex, conds, lang='C')[,!"draw"]
  
  # Check results ####
  #Check microdata results
  mcmc.hh.res <- list()
  mcmc.hh.res$micro = merge(as.data.table(prop.table(table(mcmc.hh[,.(INCOME,RESTY,HHVEH,HHSIZ)]))),
                            as.data.table(prop.table(table(data.hh))), by = colnames(data.hh), all = T, suffixes = c("Gibbs","Data"))
  mcmc.hh.res$micro = mcmc.hh.res$micro[ , lapply(.SD, function(x) ifelse(is.na(x),0,x))]
  ggplot(data = mcmc.hh.res$micro) + geom_point(aes(x=NGibbs, y=NData)) + theme_bw() + coord_fixed()
  
  
  #Check marginal results
  mcmc.hh.res$marg <- list("INCOME" = as.matrix(margins$hhinc),
                           "RESTY" = margins$hhdwell,
                           "HHVEH" = margins$hhveh,
                           "HHSIZ" = rowSums(margins$hhvehsize,dims=2))
  
  mcmc.hh.res$marg <- rbindlist(lapply(names(mcmc.hh.res$marg), function(n) {
    setNames(as.data.table(colSums(mcmc.hh.res$marg[[n]]), keep.rownames = T),c("Variable","Marginals"))
  }))
  # mcmc.hh.res$marg <- rbindlist(lapply(names(mcmc.hh.res$marg), function(n) {
  #   setNames(as.data.table(mcmc.hh.res$marg[[n]][tr,], keep.rownames = T),c("Variable","Marginals"))
  # }))
  
  #aggregate data
  res.mcmc.hh <- rbindlist(lapply(names(mcmc.hh), function(n) {
    setNames(as.data.table(table(mcmc.hh[,.(INCOME,RESTY,HHVEH,HHSIZ)]))[ , sum(N), by = n], c("Variable","MCMC"))
  }))
  
  #Merging together and melting
  mcmc.hh.res$marg <- merge(res.mcmc.hh, mcmc.hh.res$marg, by = 'Variable', all = T)
  mcmc.hh.res$marg[ , lapply(.SD,sum), .SDcols = c("MCMC","Marginals"), by = Variable]
  mcmc.hh.res$marg[ , (c("MCMC","Marginals")) := lapply(.SD,function(x) x/sum(x)), .SDcols = c("MCMC","Marginals")]
  
  ggplot(data=mcmc.hh.res$marg) +#[ , lapply(.SD,sum), .SDcols = c("MCMC","Marginals"), by = Variable]) + 
    geom_point(aes(x=Marginals, y=MCMC, color=Variable)) +
    theme_bw() + coord_fixed() #+ theme(legend.position="none")
  
  #Output
  save(mcmc.hh, mcmc.hh.res, file = "./db/d3.2_MCMChh.RData")
  
  rm(list = ls()[grepl("mcmc.",ls())])
}

#### Raking the results to fit marginals ####
if(!file.exists("./db/d3.2_MCMCrake.RData")) {
  if(!("mcmc.hh" %in% ls())) {load("./db/d3.2_MCMChh.RData")}
  if(!("mcmc.ind" %in% ls())) {load("./db/d3.2_MCMCind.RData")}

  
  #Setup totals and data ####
  freq.hh <- as.data.table(table(mcmc.hh[,.(INCOME,RESTY,HHVEH,HHSIZ)]))
  freq.ind <- as.data.table(table(mcmc.ind[,.(AGE,GEND,RELATE,OCC,INDUS)]))
  
                     
  #Setup ind marginals
  marg.ind <- list("AGE" = apply(margins$indagesex, c(1,3), sum),
                   "GEND" = apply(margins$indagesex, c(1,2), sum),
                   "RELATE" = as.matrix(margins$indrelate),
                   "OCC"  = margins$indocc,
                   "INDUS"  = margins$indindus)
  
  nm <- names(marg.ind)
  marg.ind <- lapply(names(marg.ind), function(n) as.data.table(setNames(melt(marg.ind[[n]]),c("otract",n,"Freq"))) )
  names(marg.ind) <- nm
  
  #Setup hh marginals
  marg.hh <- list("INCOME" = as.matrix(margins$hhinc),
                  "RESTY"  = margins$hhdwell,
                  "HHVEH"  = apply(margins$hhvehsize, c(1,3), sum),
                  "HHSIZ" = apply(margins$hhvehsize, c(1,2), sum))
  
  nm <- names(marg.hh)
  marg.hh <- lapply(names(marg.hh), function(n) as.data.table(setNames(melt(marg.hh[[n]]),c("otract",n,"Freq"))) )
  names(marg.hh) <- nm
  
  
  ##Individuals ####
  #Set up parallel cores
  if(.Platform$OS.type == "unix") {
    cl <- makeCluster(detectCores()-1, "FORK") 
  } else {
    cl <- makeCluster(detectCores(logical = F)-1, "PSOCK")
    clusterEvalQ(cl,c(library(data.table),library(survey)))
    clusterExport(cl,
                  list("fun_raketract.singledim","marg.ind","freq.ind"),
                  envir=environment())
  }
  
  #Raking
  mcmc.ind.raked <- pblapply(tracts, function(tr) fun_raketract.singledim(tr, freq.ind, marg.ind, bnd=c(0.1,8), trm=c(0.01,1000)), cl=cl)
  names(mcmc.ind.raked) <- tracts
  stopCluster(cl)
  
  #Check results
  ind.raked.res <- rbindlist(mcmc.ind.raked)
  ind.raked.res <- lapply(names(marg.ind), function(n) {
    out <- ind.raked.res[,lapply(.SD,sum), .SDcols=c("weights.raked","weights.sim"), by = c('otract',n)]
    return(setNames(out, c("otract","Variable","Raked","MCMC")))
  })
  ind.raked.res <- rbindlist(ind.raked.res)
  #Marginals
  ind.raked.res.marg <- as.data.table(rbindlist(lapply(marg.ind, function(x) setNames(x, c("otract","Variable","Marginals")))))
  ind.raked.res.marg[ , otract := as.character(otract)]
  #Merging together and melting
  mcmc.ind.raked.res <- merge(ind.raked.res, ind.raked.res.marg, by = c('otract','Variable'), all = T)
  mcmc.ind.raked.res[ , (c("Raked","MCMC","Marginals")) := lapply(.SD, function(x) ifelse(is.na(x),0,x)), .SDcols = c("Raked","MCMC","Marginals")]
  
  ggplot(data = mcmc.ind.raked.res) + geom_point(aes(x=Marginals, y=Raked, color=Variable), alpha=0.5) + theme_bw() + coord_fixed()
  ggplot(data = mcmc.ind.raked.res[,lapply(.SD,sum), by=Variable, .SDcols = c("Raked","MCMC","Marginals")]) + 
    geom_point(aes(x=MCMC, y=Marginals, color=Variable), alpha=0.5) + theme_bw() + coord_fixed()
  
  
  ## Households ####
  #Set up parallel cores
  if(.Platform$OS.type == "unix") {
    cl <- makeCluster(detectCores()-1, "FORK") 
  } else {
    cl <- makeCluster(detectCores(logical = F)-1, "PSOCK")
    clusterEvalQ(cl,c(library(data.table),library(survey)))
    clusterExport(cl,
                  list("fun_raketract.singledim","marg.hh","freq.hh"),
                  envir=environment())
  }
  #Raking on at a time
  mcmc.hh.raked <- pblapply(tracts, function(tr) fun_raketract.singledim(tr, freq.hh, marg.hh, bnd=c(0.2,4), trm=c(0.01,100)), cl=cl)
  names(mcmc.hh.raked) <- tracts
  stopCluster(cl)
  
  # Checking the raked results #
  #Households
  hh.raked.res <- rbindlist(mcmc.hh.raked)
  hh.raked.res <- lapply(names(marg.hh), function(n) {
    out <- hh.raked.res[,lapply(.SD,sum), .SDcols=c("weights.raked","weights.sim"), by = c('otract',n)]
    return(setNames(out, c("otract","Variable","Raked","MCMC")))
  })
  hh.raked.res <- rbindlist(hh.raked.res)
  #Marginals
  hh.raked.marg <- as.data.table(rbindlist(lapply(marg.hh, function(x) setNames(x, c("otract","Variable","Marginals")))))
  hh.raked.marg[ , otract := as.character(otract)]
  #Merging together and melting
  mcmc.hh.raked.res <- merge(hh.raked.res, hh.raked.marg, by = c('otract','Variable'), all = T)
  mcmc.hh.raked.res[ , (c("Raked","MCMC","Marginals")) := lapply(.SD, function(x) ifelse(is.na(x),0,x)), .SDcols = c("Raked","MCMC","Marginals")]

  ggplot(data = mcmc.hh.raked.res) + geom_point(aes(x=Marginals, y=Raked, color=Variable), alpha=0.5) + theme_bw() + coord_fixed()
  ggplot(data = mcmc.hh.raked.res[,lapply(.SD,sum), by=Variable, .SDcols = c("Raked","MCMC","Marginals")]) + 
    geom_point(aes(x=Raked, y=MCMC, color=Variable), alpha=0.5) + theme_bw() + coord_fixed()
  
  
  #Saving the output
  save(mcmc.hh.raked, mcmc.ind.raked, mcmc.hh.raked.res, mcmc.ind.raked.res, file = "./db/d3.2_MCMCrake.RData")
  rm(list = ls()[grepl("mcmc.",ls())])
}

#### Raking the OD results to fit marginals ####
if(!file.exists("./db/d3.2_MCMCindodrake.RData")) {
  if(!("mcmc.ind" %in% ls())) {load("./db/d3.2_MCMCind.RData")}
  if(!("mcmc.indod" %in% ls())) {load("./db/d3.2_MCMCindod.RData")}
  if(!("mcmc.odi" %in% ls())) {load("./db/d3.2_MCMClodes.RData")}
  
  #grab OD for tracts we need
  dtracts <- apply(mcmc.odi[-1,,], 1:2, sum)
  #Get conditional probability
  dtracts <- sweep(dtracts, 1, rowSums(dtracts), "/")
  #Get frequency
  dtracts <- rowSums(margins$indindus) * dtracts
  
  #Setup ind marginals
  marg.indod <- list("AGE" = apply(margins$indagesex, c(1,3), sum),
                     "GEND" = apply(margins$indagesex, c(1,2), sum),
                     "RELATE" = as.matrix(margins$indrelate),
                     "OCC"  = margins$indocc,
                     "INDUS"  = margins$indindus,
                     "dtract" = dtracts)
  
  nm <- names(marg.indod)
  marg.indod <- lapply(names(marg.indod), function(n) as.data.table(setNames(melt(marg.indod[[n]]),c("otract",n,"Freq"))) )
  names(marg.indod) <- nm
  
  marg.indod$dtract[ , dtract := as.character(dtract)]
  marg.indod$dtract[dtract == "0", dtract := "00000000000"]
  marg.indod$dtract[ , dtract := factor(dtract, levels = colnames(dtracts))]
  
  
  ##Individuals with OD ####
  #Set up parallel cores
  if(.Platform$OS.type == "unix") {
    cl <- makeCluster(detectCores()-1, "FORK") 
  } else {
    cl <- makeCluster(detectCores(logical = F)-1, "PSOCK")
    clusterEvalQ(cl,c(library(data.table),library(survey)))
    clusterExport(cl,
                  list("fun_raketract.singledim","marg.indod","mcmc.indod"),
                  envir=environment())
  }
  
  #Raking
  mcmc.indod.raked <- pblapply(tracts, function(tr) {
    freq.dat <- as.data.table(table(mcmc.indod[otract==tr,.(AGE,GEND,RELATE,OCC,INDUS,dtract)]))
    marg.dat <- marg.indod
    #Pulling the marginal
    nm <- names(marg.dat)
    levs.og <- lapply(nm, function(n) levels(marg.dat[[n]][[n]]))
    names(levs.og) <- nm
    
    #Pulling the marginal & data for current tract
    marg.dat <- lapply(marg.dat, function(x) x[otract == tr, !"otract"])
    if('otract' %in% colnames(freq.dat)) freq.dat <- freq.dat[otract == tr, !"otract"]
    
    #Removing margin zeros from data
    zeros <- lapply(marg.dat, function(x) x[Freq==0, ])
    for(n in names(zeros)) freq.dat[['N']][ freq.dat[[n]] %in% zeros[[n]][[n]] ] <- 0
    freq.dat <- freq.dat[N>0,]
    
    #Removing the zeros from the margins
    marg.dat <- lapply(marg.dat, function(x) x[Freq!=0, ])
    
    #Relevel the margin factors
    for(n in names(marg.dat)) {
      levs = levels(marg.dat[[n]][[n]])[levels(marg.dat[[n]][[n]]) %in% marg.dat[[n]][[n]]]
      marg.dat[[n]][[n]] <- factor(marg.dat[[n]][[n]], levels = levs)
    }
    #Matching factors in data to margin factors
    for(n in names(marg.dat)) freq.dat[[n]] <- factor(freq.dat[[n]], levels = levels(marg.dat[[n]][[n]]))
    
    #Scaling the weights
    total <- sum(marg.dat[[1]]$Freq)
    colnames(freq.dat)[which("N" == colnames(freq.dat))] <- "weights.sim"
    freq.dat[ , weights.sim := total*weights.sim/sum(weights.sim)]
    
    #Check if only 1 factor level remaining, aka only one household type
    singles <- sapply(marg.dat, function(x) ifelse(nrow(x[Freq>0,])<=1, x[[1]],NA))
    
    if(sum(is.na(singles)) > 0 ) {
      #Removing the singles temporarily
      freq.dat.stripped <- freq.dat[ , c(names(singles[is.na(singles)]),"weights.sim"), with=F]
      #Setting up the raking data
      design.dat <- svydesign(ids = ~1, weights=~weights.sim, data=freq.dat.stripped)
      #Formula
      form <- as.formula(paste("~",paste(names(singles[is.na(singles)]),collapse = "+"),sep=""))
      #Population marginals
      pop <- c(`(Intercept)` = total, unlist(lapply(names(singles[is.na(singles)]), function(n) {
        setNames(marg.dat[[n]]$Freq[-1], paste(n,marg.dat[[n]][[n]][-1],sep=""))
      })))
      raked.dat = 'character'
      #Raking the data, strictest standards
      
      raked.dat <- suppressWarnings(
        tryCatch(
          calibrate(design.dat, formula = form,
                    population = pop, calfun = "raking",
                    bounds = c(0,quantile(freq.dat$weights.sim, probs=0.90)),
                    trim = quantile(freq.dat$weights.sim, probs=c(0.01,0.99)), #c(0.0001,10), 
                    #epsilon=1,
                    sparse=T,
                    force=T),
          error = function(e) "error"))
      #Putting adding the weights back in
      if(all(class(raked.dat) == "character")) {
        freq.dat.raked <- cbind("otract" = tr, freq.dat, "weights.raked" = freq.dat$weights.sim)
      } else {
        freq.dat.raked <- cbind("otract" = tr, freq.dat, "weights.raked" = weights(raked.dat))
      }
    } else {
      #No change possible, just return the original, with eps = 0
      freq.dat.raked <- cbind("otract" = tr, freq.dat, "weights.raked" = freq.dat$weights.sim)
    }
    
    #Convert back to original factors
    for(n in nm) freq.dat.raked[[n]] <- as.character(freq.dat.raked[[n]])
    for(n in names(marg.dat)) freq.dat.raked[[n]] <- factor(freq.dat.raked[[n]], levels = levs.og[[n]])
    #Out
    return(freq.dat.raked)
  
  }, cl=cl)
  names(mcmc.indod.raked) <- tracts
  stopCluster(cl)
  
  #Check results
  indod.raked.res <- rbindlist(mcmc.indod.raked)
  indod.raked.res <- lapply(names(marg.indod), function(n) {
    out <- indod.raked.res[,lapply(.SD,sum), .SDcols=c("weights.raked","weights.sim"), by = c('otract',n)]
    out$Category <- n
    return(setNames(out, c("otract","Variable","Raked","MCMC","Category")))
  })
  indod.raked.res <- rbindlist(indod.raked.res)
  #Marginals
  indod.raked.res.marg <- as.data.table(rbindlist(lapply(marg.indod, function(x) setNames(x, c("otract","Variable","Marginals")))))
  indod.raked.res.marg[ , otract := as.character(otract)]
  #Merging together and melting
  mcmc.indod.raked.res <- merge(indod.raked.res, indod.raked.res.marg, by = c('otract','Variable'), all = T)
  mcmc.indod.raked.res[ , (c("Raked","MCMC","Marginals")) := lapply(.SD, function(x) ifelse(is.na(x),0,x)), .SDcols = c("Raked","MCMC","Marginals")]
  
  mcmc.indod.raked.res <- mcmc.indod.raked.res[Raked>0 & MCMC>0 & Marginals>0,]
  
  ggplot(data = mcmc.indod.raked.res) + 
    geom_point(aes(x=Marginals, y=Raked), alpha=0.1) + theme_bw() + coord_fixed()
  
  ggplot(data = mcmc.indod.raked.res[,lapply(.SD,sum), by=Variable, .SDcols = c("Raked","MCMC","Marginals")]) + 
    geom_point(aes(x=MCMC, y=Marginals), alpha=0.5) + theme_bw() + coord_fixed()
  
  #Saving the output
  save(mcmc.indod.raked, mcmc.indod.raked.res, file = "./db/d3.2_MCMCindodrake.RData")
  rm(list = ls()[grepl("mcmc.",ls())])
}


#### Cleanup ####
rm(list=ls())
