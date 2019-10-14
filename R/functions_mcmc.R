
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



