##
# Population Synthesis for Boston Metropolitan Region, by Nicholas Marc Fournier
# Last updated January 2017.
#
# NOTICE:  All information, intellectual and technical concepts contained herein is,
# and remains the property of Nicholas Marc Fournier. Dissemination of this 
# information or reproduction of this material is strictly forbidden unless
# prior written permission is obtained from Nicholas Marc Fournier
##

#### Loading tables & settings ####
# sourceCpp("./cpp_fastsample.cpp")
source("./s1.0_packages.R", echo=F)

#### Setting up functions ####
int_trs <- function(x){ # For generalisation purpose, x becomes a vector
  xv <- as.vector(x) # allows trs to work on matrices
  xint <- floor(xv) # integer part of the weight
  r <- xv - xint # decimal part of the weight
  def <- round(sum(r)) # the deficit population # the weights be 'topped up' (+ 1 applied)
  #topup <- sample(length(x), size = def, prob = r, replace = F)
  topup <- RcppSample(length(x), size = def, prob = r)
  xint[topup] <- xint[topup] + 1
  dim(xint) <- dim(x)
  dimnames(xint) <- dimnames(x)
  xint
}
int_expand_array <- function(x){
  # Transform the array into a dataframe
  count_data <- as.data.frame.table(x)
  # Store the indices of categories for the final population
  indices <- rep(1:nrow(count_data), count_data$Freq)
  # Create the final individuals
  ind_data <- count_data[indices,]
  data.table(ind_data)
}
joint_samp <- function(jfit, dest, pums, tr) {
  #Finding best fit random joint sample
  X <- jfit$X
  X <- X[X>0]
  
  dest[dest<0] <- 0
  #print(tr)
  #Draw samples
  sampvec <- sample(names(X), sum(X), T, X)
  
  #Generating joint sample by merging
  jointsamp = merge(data.table(SERIALNO = sampvec, HHID = paste(tr,1:length(sampvec),sep="_"), otract=tr), pums, by = 'SERIALNO')
  
  #Sample destinations
  if(length(dim(dest))>3) {
    #When dest is integrated
    jointsamp[ , dtract := apply(jointsamp, 1, function(samp) {
      n = names(dimnames(dest))[which(names(dimnames(dest))!="dtract")]
      dees = dest[samp[n[1]], samp[n[2]], samp[n[3]], samp[n[4]], samp[n[5]], ]
      sample(names(dees), 1, T, dees)
    })]
  } else {
    #When dest is not integrated
    jointsamp[ , dtract := apply(jointsamp, 1, function(x) sample(names(dest[,x['INDUS']]), 1, T, dest[,x['INDUS']]))]
  }
  
  #Reorder and cleanup
  return(jointsamp[ , c("HHID",colnames(pums[,!"SERIALNO"]),"otract","dtract"),with=F])
}
joint_samp_replace <- function(jfit, di, pums, tr) {
  #Finding best fit random joint sample
  X <- jfit$X
  X <- X[X>0]
  #print(tr)
  #Draw samples
  sampvec <- sample(names(X), sum(X), T, X)
  
  # #Empty table
  # jointsamp <- data.table(SERIALNO = as.character(), HHID = as.character(), otract=as.character(), pums[0,-1], dtract=as.character())
  
  #Generating joint sample by merging
  jointsamp = merge(data.table(SERIALNO = sampvec, HHID = paste(tr,1:length(sampvec),sep="_"), otract=tr), pums, by = 'SERIALNO')
  
  #Sample destinations
  jointsamp[ , dtract := apply(jointsamp_new, 1, function(x) {
    if(sum(di[,x['INDUS']]) > 0) sample(names(di[,x['INDUS']]), 1, T, di[,x['INDUS']])
    else as.character(NA)
  }) ]
  
  # badhh <- unique(jointsamp_new[is.na(dtract),HHID])
  # jointsamp <- rbind(jointsamp, jointsamp_new[!(HHID %in% badhh),])
  # remainder = sum(X) - length(unique(jointsamp$HHID))
  
  #Reorder and cleanup
  return(jointsamp[ , c("HHID",colnames(pums[,!"SERIALNO"]),"otract","dtract"),with=F])
}
joint_samp_old <- function(ind, hh, jfit, tr) {
  #Finding best fit random joint sample
  X <- jfit$X
  X <- X[X>0]
  #Grabbing b vector
  b <- jfit$errors$b
  #Pulling N-samples, keeping best fit
  rmse_min <- 100
  for(i in 1:10){
    samp_i <- sample(x = names(X), size = sum(X), replace = T, prob = X)
    rmse <- sqrt(sum((rowSums(A[ , samp_i])-b)^2)/length(b))
    #qplot(x = b, y = rowSums(A[ , samp_i]), geom='point')
    if(rmse < rmse_min) rmse_min <- rmse ; samp <- samp_i
  }
  #Flatting the IPF weights for individuals
  indflat <- as.data.table(as.data.frame.table(ind))
  indflat <- merge(indflat, types$indtypes[,!"Freq"], by = indvars)
  indflat <- indflat[Freq>0,]
  indflat <- indflat[ , lapply(.SD, as.character)]
  #Flattening and expanding the joint sample
  jointflat <- lapply(1:length(samp), function(x) {
    hh <- types$pumstypes[SERIALNO==samp[x],]
    hh[ , HHID := paste(tr,x,sep="_")]})
  #Adding unique household ID
  names(jointflat) <- paste(tr,1:length(samp), sep="_")
  jointflat <- lapply(names(jointflat), function(x) jointflat[[x]][ , HHID := x])
  jointflat <- rbindlist(jointflat)
  jointflat[ , dtract := as.character()]
  #person types to cycle through
  pertypes <- types[['indtypes']]$PERTYPE
  #Must in BOTH the joint sample and the flat IPF fit
  pertypes <- pertypes[pertypes %in% jointflat$PERTYPE & pertypes %in% indflat$PERTYPE]
  #Assigning destination for persons in chunks of person types
  jointflat <- lapply(pertypes, function(p) {
    #Matched person type subset
    persub <- indflat[PERTYPE == p, ]
    #Subsetting for person types
    jointsub <- jointflat[PERTYPE == p, ]
    #Sample size
    N <- nrow(jointflat[PERTYPE == p, ])
    #Assigning destination using MC sample
    jointsub[, dtract := sample(x = persub$dtract, size = N, replace = T, prob = persub$Freq)]
  })
  jointflat <- rbindlist(jointflat)
  #Merging in the other attributes
  jointflat <- merge(jointflat, types$hhtypes[,!"Freq"], by = "HHTYPE", all.y = F)
  jointflat <- merge(jointflat, types$indtypes[,!"Freq"], by = "PERTYPE", all.y = F)
  #
  jointflat[ , otract := tr]
  return(jointflat)
}

#Load data
load("./db/d2.3_microdata.RData")
load("./db/d2.5_margins.RData")

#Setup generic data inputs
cols = c("AGE","GEND","RELATE","OCC","INDUS","INCOME","RESTY","HHSIZ","HHVEH")
pums = microdata$pums[,c("SERIALNO",cols),with=F]
#Major variables
indvars  <- c("AGE","GEND","RELATE","OCC","INDUS")
hhvars <- c("INCOME","RESTY","HHVEH","HHSIZ")

#### Generate IPF with OD ####
if(!file.exists("./db/d5.1_IPFODpop.RData")) {
  load("./db/d4.1_IPFODjointfit.RData")
  # Loading database
  dbipfod.ind <- dbInit("./db/d3.1_dbipfindod", "RDS")
  
  #Set up parallel cores
  if(.Platform$OS.type == "unix") {
    cl <- makeCluster(detectCores()-1, "FORK") 
  } else {
    cl <- makeCluster(detectCores(logical = F)-1, "PSOCK")
    clusterEvalQ(cl,c(library(data.table),library(filehash)))
    clusterExport(cl,
                  list("ipfod.jointfit","dbipfod.ind","tracts","joint_samp","pums"),
                  envir=environment())
  }
  #Running generation in parallel
  ipfod.pop = pblapply(tracts, function(tr) joint_samp(ipfod.jointfit[[tr]], dbipfod.ind[[tr]], pums, tr), cl = cl)
  stopCluster(cl)
  
  #Combining
  ipfod.pop <- rbindlist(ipfod.pop)
  
  #Saving & cleanup
  save(ipfod.pop, file = "./db/d5.1_IPFODpop.RData")
}

#### Generate IPF ####
if(!file.exists("./db/d5.1_IPFpop.RData")) {
  load("./db/d4.1_IPFjointfit.RData")
  load("./db/d3.1_IPFlodes.RData")
  
  #Set up parallel cores
  if(.Platform$OS.type == "unix") {
    cl <- makeCluster(detectCores()-1, "FORK") 
  } else {
    cl <- makeCluster(detectCores(logical = F)-2, "PSOCK")
    clusterEvalQ(cl,c(library(data.table)))
    clusterExport(cl,
                  list("ipf.jointfit","ipf.lodes","tracts","joint_samp","pums"),
                  envir=environment())
  }
  #Running generation in parallel
  ipf.pop = pblapply(tracts, function(tr) joint_samp(ipf.jointfit[[tr]], ipf.lodes$p.hat[tr,,], pums, tr), cl = cl)
  stopCluster(cl)
  
  #Combining
  ipf.pop <- rbindlist(ipf.pop)
  
  #Saving & cleanup
  save(ipf.pop, file = "./db/d5.1_IPFpop.RData")
}

#### Generate MCMC with OD ####
if(!file.exists("./db/d5.2_MCMCODpop.RData")) {
  load("./db/d4.2_MCMCODjointfit.RData")
  load("./db/d3.2_MCMCindod.RData")
  #load("./db/d3.2_MCMClodes.RData")
  #load("./db/d3.1_IPFlodes.RData")
  
  #Loading stored results
  if(!dir.exists("./db/d5.2_dbMCMCindod")){ dbCreate("./db/d5.2_dbMCMCindod", "RDS") }
  dbmcmc.indod <- dbInit("./db/d5.2_dbMCMCindod", "RDS")
  
  
  #Set up parallel cores
  if(.Platform$OS.type == "unix") {
    cl <- makeCluster(detectCores()-1, "FORK") 
  } else {
    cl <- makeCluster(detectCores(logical = F)-1, "PSOCK")
    clusterEvalQ(cl,c(library(data.table),library(filehash)))
    clusterExport(cl,
                  list("margins","mcmc.indod","dbmcmc.indod"),
                  envir=environment())
  }
  
  #Creating integrated OD
  dummy <- pblapply(tracts, function(tr) {
    od = prop.table(table(mcmc.indod[otract==tr,!"otract"])) * sum(margins$indagesex[tr,,])
    dbInsert(dbmcmc.indod, tr, od)
    return(tr)
  }, cl = cl)
  stopCluster(cl)
  
  # mcmc.indod <- pblapply(tracts, function(tr) 
  #   prop.table(table(mcmc.indod[otract==tr,!"otract"])) * sum(margins$indagesex[tr,,]))
  # names(mcmc.indod) <- tracts
  
  #Set up parallel cores
  if(.Platform$OS.type == "unix") {
    cl <- makeCluster(detectCores()-1, "FORK") 
  } else {
    cl <- makeCluster(detectCores(logical = F)-2, "PSOCK")
    clusterEvalQ(cl,c(library(data.table),library(filehash)))
    clusterExport(cl,
                  list("mcmc.odjointfit","dbmcmc.indod","tracts","joint_samp","pums"),
                  envir=environment())
  }
  #Running generation in parallel
  mcmc.odpop = pblapply(tracts, function(tr) joint_samp(mcmc.odjointfit[[tr]], dbmcmc.indod[[tr]], pums, tr), cl = cl)
  stopCluster(cl)
  
  #Combining
  mcmc.odpop <- rbindlist(mcmc.odpop)
  
  #Saving & cleanup
  save(mcmc.odpop, file = "./db/d5.2_MCMCODpop.RData")
}

#### Generate raked MCMC with OD ####
if(!file.exists("./db/d5.2_MCMCODrakedpop.RData")) {
  load("./db/d4.2_MCMCODrakejointfit.RData")
  load("./db/d3.2_MCMCindodrake.RData")
  #load("./db/d3.2_MCMClodes.RData")
  #load("./db/d3.1_IPFlodes.RData")
  
  #Loading stored results
  if(!dir.exists("./db/d5.2_dbMCMCindodrake")){ dbCreate("./db/d5.2_dbMCMCindodrake", "RDS") }
  dbmcmc.indod.raked <- dbInit("./db/d5.2_dbMCMCindodrake", "RDS")
  
  #
  fm.ind <- as.formula(paste0("weights.raked~",paste(c(indvars,"dtract"), collapse = '+')))
  
  #Set up parallel cores
  if(.Platform$OS.type == "unix") {
    cl <- makeCluster(detectCores()-1, "FORK") 
  } else {
    cl <- makeCluster(detectCores(logical = F)-1, "PSOCK")
    clusterEvalQ(cl,c(library(data.table),library(filehash)))
    clusterExport(cl,
                  list("fm.ind","margins","mcmc.indod.raked","dbmcmc.indod.raked"),
                  envir=environment())
  }
  
  #Creating integrated OD
  dummy <- pblapply(tracts, function(tr) {
    od = prop.table(xtabs(fm.ind, data=mcmc.indod.raked[[tr]])) * sum(margins$indagesex[tr,,])
    dbInsert(dbmcmc.indod.raked, tr, od)
    return(tr)
  }, cl = cl)
  stopCluster(cl)
  
  # mcmc.indod.raked <- pblapply(tracts, function(tr) 
  #   prop.table(table(mcmc.indod.raked[otract==tr,!"otract"])) * sum(margins$indagesex[tr,,]))
  # names(mcmc.indod.raked) <- tracts
  
  #Set up parallel cores
  if(.Platform$OS.type == "unix") {
    cl <- makeCluster(detectCores()-1, "FORK") 
  } else {
    cl <- makeCluster(detectCores(logical = F)-2, "PSOCK")
    clusterEvalQ(cl,c(library(data.table),library(filehash)))
    clusterExport(cl,
                  list("mcmc.odrakejointfit","dbmcmc.indod.raked","tracts","joint_samp","pums"),
                  envir=environment())
  }
  #Running generation in parallel
  mcmc.odrakepop = pblapply(tracts, function(tr) {
    print(tr)
    joint_samp(mcmc.odrakejointfit[[tr]], dbmcmc.indod.raked[[tr]], pums, tr)
    })#, cl = cl)
  stopCluster(cl)
  
  #Combining
  mcmc.odrakepop <- rbindlist(mcmc.odrakepop)
  
  #Saving & cleanup
  save(mcmc.odpop, file = "./db/d5.2_MCMCODrakepop.RData")
}

#### Generate MCMC ####
if(!file.exists("./db/d5.2_MCMCpop.RData")) {
  load("./db/d4.2_MCMCjointfit.RData")
  load("./db/d3.2_MCMClodes.RData")
  #load("./db/d3.1_IPFlodes.RData")
  
  #Set up parallel cores
  if(.Platform$OS.type == "unix") {
    cl <- makeCluster(detectCores()-1, "FORK") 
  } else {
    cl <- makeCluster(detectCores(logical = F)-1, "PSOCK")
    clusterEvalQ(cl,c(library(data.table)))
    clusterExport(cl,
                  list("mcmc.jointfit","mcmc.odi","tracts","joint_samp","pums"),
                  envir=environment())
  }
  #Running generation in parallel
  mcmc.pop = pblapply(tracts, function(tr) joint_samp(mcmc.jointfit[[tr]], mcmc.odi[tr,,], pums, tr), cl = cl)
  stopCluster(cl)
  
  #Combining
  mcmc.pop <- rbindlist(mcmc.pop)
  
  #Saving & cleanup
  save(mcmc.pop, file = "./db/d5.2_MCMCpop.RData")
}

#### Generate BN ####
if(!file.exists("./db/d5.3_BNpop.RData")) {
  load("./db/d4.3_BNjointfit.RData")
  load("./db/d3.2_MCMClodes.RData")
  
  #Set up parallel cores
  if(.Platform$OS.type == "unix") {
    cl <- makeCluster(detectCores()-1, "FORK") 
  } else {
    cl <- makeCluster(detectCores(logical = F)-2, "PSOCK")
    clusterEvalQ(cl,c(library(data.table)))
    clusterExport(cl,
                  list("bn.jointfit","mcmc.odi","tracts","joint_samp","pums"),
                  envir=environment())
  }
  #Running generation in parallel
  bn.pop = pblapply(tracts, function(tr) joint_samp(bn.jointfit[[tr]], mcmc.odi[tr,,], pums, tr), cl = cl)
  stopCluster(cl)
  
  #Combining
  bn.pop <- rbindlist(bn.pop)
  
  #Saving & cleanup
  save(bn.pop, file = "./db/d5.3_BNpop.RData")
}

#### Generate BN with OD ####
if(!file.exists("./db/d5.3_BNODpop.RData")) {
  load("./db/d4.3_BNODjointfit.RData")
  load("./db/d3.3_BNindod.RData")
  #load("./db/d3.2_MCMClodes.RData")
  
  #Loading stored results
  if(!dir.exists("./db/d5.2_dbBNindod")){ dbCreate("./db/d5.2_dbBNindod", "RDS") }
  dbbn.indod <- dbInit("./db/d5.2_dbBNindod", "RDS")
  
  #Set up parallel cores
  if(.Platform$OS.type == "unix") {
    cl <- makeCluster(detectCores()-1, "FORK") 
  } else {
    cl <- makeCluster(detectCores(logical = F)-1, "PSOCK")
    clusterEvalQ(cl,c(library(data.table),library(filehash)))
    clusterExport(cl,
                  list("margins","bn.indod","dbbn.indod"),
                  envir=environment())
  }
  
  #Creating integrated OD
  dummy <- pblapply(tracts, function(tr) {
    od = prop.table(table(bn.indod[otract==tr,!"otract"])) * sum(margins$indagesex[tr,,])
    dbInsert(dbbn.indod, tr, od)
    return(tr)
  }, cl = cl)
  stopCluster(cl)
  
  
  #Set up parallel cores
  if(.Platform$OS.type == "unix") {
    cl <- makeCluster(detectCores()-1, "FORK") 
  } else {
    cl <- makeCluster(detectCores(logical = F)-2, "PSOCK")
    clusterEvalQ(cl,c(library(data.table),library(filehash)))
    clusterExport(cl,
                  list("bn.odjointfit","dbbn.indod","tracts","joint_samp","pums"),
                  envir=environment())
  }
  #Running generation in parallel
  bn.odpop = pblapply(tracts, function(tr) joint_samp(bn.odjointfit[[tr]], dbbn.indod[[tr]], pums, tr), cl = cl)
  stopCluster(cl)
  
  #Combining
  bn.odpop <- rbindlist(bn.odpop)
  
  #Saving & cleanup
  save(bn.odpop, file = "./db/d5.3_BNODpop.RData")
}


#### Cleanup ####

#rm(list=setdiff(ls(),c("clock","runtime","csvwriteout")))
