#### Loading tables & settings ####
source("./s1.0_packages.R", echo=F)

sourceCpp("./cpp_IPU.cpp")
sourceCpp("./cpp_fastsample.cpp")
#sourceCpp("./cpp_dupes.cpp")

#### Setting up functions ####
format_fulltypes <- function(ind, hh, indvars, hhvars){
  # #Loading pums data
  # Getting HH and Person types from IPF
  if(is.list(hh)){
    print("Aggregating person matrices")
    #Summing into single matrix
    indtypes <- Reduce("+", ind)
    hhtypes <- Reduce("+", hh)
  } else {
    print("Using existing total weights")
    indtypes <- indtots
    hhtypes <- hhtots
  }
  print("Formatting sparse matrix from Northeast US PUMS")
  pumstypes <- microdata[["pums"]][,c("SERIALNO", indvars, hhvars), with=F]
  pumstypes[ , SERIALNO := as.character(SERIALNO)]
  #expand IPF results into hh/person type frequencies
  indtypes <- as.data.table(indtypes)
  hhtypes <- as.data.table(hhtypes)
  #labeling types
  hhtypes[,HHTYPE := paste("HHTYPE",1:nrow(hhtypes),sep="")]
  indtypes[,PERTYPE := paste("PERTYPE",1:nrow(indtypes),sep="")]
  #Merging labels to person/household types in pums pop
  pumstypes <- merge(pumstypes, hhtypes[,c(hhvars,"HHTYPE"), with=F], by=hhvars, all.x = T)
  pumstypes <- merge(pumstypes, indtypes[,c(indvars,"PERTYPE"), with=F], by=indvars, all.x = T)
  #Compact DT
  pumstypes <- pumstypes[,.(SERIALNO, PERTYPE, HHTYPE)]
  #Getting frequencies
  pertypefreq <- pumstypes[, .N, by = .(SERIALNO,PERTYPE)]
  hhtypefreq <- pumstypes[!duplicated(SERIALNO), .N, by = .(SERIALNO,HHTYPE)]
  #Making uniform names and stacking results for sparse matrix
  colnames(pertypefreq) <- c("SERIALNO","TYPE","N")
  colnames(hhtypefreq) <- c("SERIALNO","TYPE","N")
  #Stacking long
  typefreq <- rbind(hhtypefreq, pertypefreq)
  
  #Releveling the factor levels
  typefreq[ , TYPE := as.factor(TYPE)]
  typefreq[ , SERIALNO := as.factor(SERIALNO)]
  
  #Sending to wide sparse matrix
  sparse <- sparseMatrix(i = as.integer(typefreq$TYPE),
                         j = as.integer(typefreq$SERIALNO),
                         x = typefreq$N,
                         dimnames = list(type = levels(typefreq$TYPE), serial = levels(typefreq$SERIALNO)))
  
  bigsparse <- list(hhtypes = hhtypes,
                    indtypes = indtypes,
                    pumstypes = pumstypes,
                    sparse = sparse)
  return(bigsparse)
}
format_bvector <- function(indcons, hhcons, sparsetypes, typevars){
  #expand IPF results into hh/person type frequencies
  indcons <- as.data.table(indcons)
  hhcons <- as.data.table(hhcons)
  #labeling types
  indcons <- merge(indcons, sparsetypes$indtypes[,c(indvars,"PERTYPE"),with=F], by=indvars, all=T)
  hhcons <- merge(hhcons, sparsetypes$hhtypes[,c(hhvars,"HHTYPE"),with=F], by=hhvars, all=T)
  #creating constraint vector
  if("value" %in% colnames(indcons)) {
    cons <- c(hhcons$value, indcons$value) 
  } else {
    cons <- c(hhcons$N, indcons$N)
  }
  names(cons) <- c(hhcons$HHTYPE, indcons$PERTYPE)
  #Removing corresponding zero columns from constraints
  cons <- cons[typevars]
  return(cons)
}
nnld_fit <- function(A, b, int = F) {
  print("Running least deviation fitting")
  time <- Sys.time()
  #Format for LAD
  Ac <- rbind(cbind(A, Diagonal(nrow(A), -1)), cbind(-A, Diagonal(nrow(A), -1)))
  if(int==T) b <- int_trs(b) ; names(b) <- rownames(A)
  #Objective coefficients
  f <- rep(0,nrow(A))
  f[which(grepl("HHTYPE",rownames(A)))] <- 1
  f[which(grepl("PERTYPE",rownames(A)))] <- 1
  f <- c(rep(0,ncol(A)), f)
  
  #Prepare data structures for the problem object. Number of columns and rows:
  prob <- initProbCLP()
  setObjDirCLP(prob, 1)
  #Matrix size
  nc <- ncol(Ac)
  nr <- nrow(Ac)
  #The constraint matrix is passed in column major order format. all indices start with 0!
  ia <- Ac@i
  #Column indices.
  ja <- Ac@p
  #Non-zero elements.
  ar <- Ac@x
  #Lower bounds for the variables (columns).
  clb <- rep(0, ncol(Ac))
  #Right hand side (row upper bounds for the rows).
  rub <- c(b,-b)
  #Objective coefficients.
  obj <- f
  #Load problem data into the problem object.
  loadProblemCLP(prob, nc, nr, ia, ja, ar, clb, NULL, obj, NULL, rub)
  #Solve the problem using the simplex algorithm.
  solveInitialCLP(prob)
  #Retrieve the value of the objective function after optimization.
  X <- getColPrimCLP(prob)[1:ncol(A)]
  names(X) <- colnames(A)
  #Free memory, allocated to the problem object.
  delProbCLP(prob)
  
  #getting results
  bhat <- as.vector(A %*% X)
  res <- data.table(varname = names(b),
                    bhat,
                    b,
                    vartype = gsub('[0-9]+', '', names(b)),
                    esttype = "Least-Deviation")
  res[ , sqerror := (b-bhat)^2]
  RMSE <- sqrt(sum(res$sqerror)/nrow(res))
  RMSN <- RMSE/mean(res$b)
  qplot(data=res, x=b, y=bhat) + theme_classic()
  time <- difftime(Sys.time(), time, units = 'secs')
  print(paste("Completed in", time, "seconds"))
  return(list(X = X, errors = res, RMSE = RMSE, RMSN = RMSN, time = time))
}
nnls_fit <- function(A, b) {
  print("Running least-squares fitting")
  time <- Sys.time()
  #Run NNLS
  LS <- nnls(A, b)
  
  #Retrieve the value of the objective function after optimization.
  X <- LS$x
  names(X) <- colnames(A)
  
  #getting results
  bhat <- as.vector(A %*% X)
  res <- data.table(varname = names(b), bhat, b, vartype = gsub('[0-9]+', '', names(b)))
  res[ , sqerror := (b-bhat)^2]
  res[ , esttype := "Least-Squares"]
  RMSE <- sqrt(sum(res$sqerror)/nrow(res))
  RMSN <- RMSE/mean(res$b)
  
  time <- difftime(Sys.time(), time, units = 'secs')
  print(paste("Completed in", time, "seconds"))
  return(list(X = X, errors = res, RMSE = RMSE, RMSN = RMSN, time = time))
}
glmnet_fit <- function(A, b, pumstypes, lambda = 0) {
  #print("Running least-squares fitting")
  time <- Sys.time()
  #Run NNLS with glmnet
  glmnetmod <- glmnet(A, b, alpha=1, lambda = 0, lower.limits=0, intercept=F)#, thresh = 1e-10)
  
  #Retrieve the value of the objective function after optimization.
  X <- coef(glmnetmod)[,1][-1]
  #names(X) <- colnames(A)
  
  #getting results
  bhat <- as.vector(A %*% X)
  res <- data.table(varname = names(b), bhat, b, vartype = gsub('[0-9]+', '', names(b)))
  res[ , sqerror := (b-bhat)^2]
  res[ , esttype := "Least-Squares (GDML)"]
  RMSE <- sqrt(sum(res$sqerror)/nrow(res))
  RMSN <- RMSE/mean(res$b)
  plot(bhat~b, res)
  plot(bhat~b, res[!(bhat>0 & b==0), ])
  
  #Finding the few erroneously weighted samples and removing
  zeros = res[bhat>0 & b==0,varname]
  #zeros = pumstypes[PERTYPE %in% zeros | HHTYPE %in% zeros, SERIALNO]
  #zeros = pumstypes[HHTYPE %in% zeros, SERIALNO]
  zeros = pumstypes[PERTYPE %in% zeros, SERIALNO]
  X <- X[!(names(X) %in% zeros)]
  #Remove 0
  X <- X[X>0]
  
  time <- difftime(Sys.time(), time, units = 'secs')
  #print(paste("Completed in", time, "seconds"))
  return(list(X = X, errors = res, RMSE = RMSE, RMSN = RMSN, time = time))
}
ipu_fit <- function(A, b, verb = T, c = 1e-5, cc = 1, maxit=100) {
  print("Running IPU fitting")
  time <- Sys.time()
  #base 0 hh type col index
  hh_idx <- which(grepl("HHTYPE", names(b)))-1
  
  #setting 0's to 0.01
  b[b==0] <- 1e-5
  
  #Running IPU
  X <- IPU_cpp_sparse(t(A), b, hh_idx, corncrit=cc, crit=c, maxit=maxit, print=verb)
  names(X) <- colnames(A)
  
  #If results in NA because too small then set to 0
  X[!is.finite(X)]<-0
  
  #getting results
  bhat <- as.vector(A %*% X)
  res <- data.table(varname = names(b), bhat, b, vartype = gsub('[0-9]+', '', names(b)))
  res[ , sqerror := (b-bhat)^2]
  res[ , esttype := "IPU"]
  #
  RMSE <- sqrt(sum(res$sqerror)/nrow(res))
  RMSN <- RMSE/mean(res$b)
  #
  time <- difftime(Sys.time(), time, units = 'secs')
  print(paste("Completed in", time, "seconds"))
  return(list(X = X, errors = res, RMSE = RMSE, RMSN = RMSN, time = time))
}
joint_fit <- function(tracts, ind, hh, sparsetypes, indvars, hhvars, par = F){
  
  #Extract A matrix
  A <- sparsetypes$sparse
  #Set up parallel cores
  if(par == T) {
    if(.Platform$OS.type == "unix") {
      cl <- makeCluster(detectCores()-1, "FORK") 
    } else {
      cl <- makeCluster(detectCores(logical = F)-1, "PSOCK")
      clusterEvalQ(cl,c(library(data.table),library(glmnet)))
      clusterExport(cl,
                    list("sparsetypes","A","ind","hh",
                         "indvars","hhvars","glmnet_fit","format_bvector"),
                    envir=environment())
    }
  } else {
    cl = NULL
  }
  #Running re-weighting
  print("Running joint fitting")
  jointfit <- pblapply(tracts, function(tr) {
    b <- format_bvector(ind[[tr]], hh[[tr]], sparsetypes, rownames(A))
    res <- glmnet_fit(A, b, sparsetypes$pumstypes)
    #res <- nnld_fit(A, b, sparsetypes$pumstypes)
    return(res)
  }, cl = cl)
  names(jointfit) <- tracts
  #Cleanup
  if(par==T) stopCluster(cl)
  gc()
  #Out
  return(jointfit)
}

#Load data
load("./db/d2.3_microdata.RData")
load("./db/d2.5_margins.RData")

#Major variables
indvars  <- c("AGE","GEND","RELATE","OCC","INDUS")
hhvars <- c("INCOME","RESTY","HHVEH","HHSIZ")
  
#### Fast Re-weighting of IPF w/ OD ####
if(!file.exists("./db/d4.1_IPFODjointfit.RData")) {
  # Loading database
  dbipfod.ind <- dbInit("./db/d3.1_dbipfindod", "RDS")
  if(!("ipf.hh" %in% ls())) {load("./db/d3.1_IPFindhh.RData")}
  rm("ipf.ind")
  
  dims = which(names(dimnames(dbipfod.ind[[tracts[1]]])) %in% indvars)
  #Separating out the individual from OD
  if(.Platform$OS.type == "unix") {
    cl <- makeCluster(detectCores()-1, "FORK") 
  } else {
    cl <- makeCluster(detectCores(logical = F)-1, "PSOCK")
    clusterEvalQ(cl,c(library(data.table),library(filehash)))
    clusterExport(cl,
                  list("dbipfod.ind","tracts","dims","indvars"),
                  envir=environment())
  }
  
  ipfod.ind = pblapply(tracts, function(tr) apply(dbipfod.ind[[tr]], dims, sum), cl = cl)
  names(ipfod.ind) <- tracts
  stopCluster(cl)
  gc()
  
  #Getting types
  if(!file.exists('./db/d4.1_ipfodsparse.RData')) {
    ipfod.sparse <- format_fulltypes(ipf.indod, ipf.hh, indvars, hhvars)
    save(ipfod.sparse, file = './db/d4.1_ipfodsparse.RData')
    gc()
  } else {
    load('./db/d4.1_ipfodsparse.RData')
  }
  
  #Joint fit
  ipfod.jointfit <- joint_fit(tracts, ipfod.ind, ipf.hh, ipfod.sparse, indvars, hhvars, par=T)
  
  # Check results
  ipfod.fitplots <- pblapply(ipfod.jointfit[1:6], function(x) {
    qplot(data=x$errors, x = b, y = bhat, geom='point') + theme_classic() + coord_fixed()
  })
  grid.arrange(grobs=ipfod.fitplots)
  
  # Fit results
  ipfod.errors <- rbindlist(pblapply(names(ipfod.jointfit), function(x) {
    tmp <- ipfod.jointfit[[x]]$errors
    tmp[ , otract := x]
  }))
  
  #Saving
  save(ipfod.jointfit, file = "./db/d4.1_IPFODjointfit.RData")
  
  #Cleanup
  rm(list=grepl("ipfod",ls()))
} 

#### Fast Re-weighting of IPF w/o OD ####
if(!file.exists("./db/d4.1_IPFjointfit.RData")) {
  # Loading IPF database
  if(all(!(c("ipf.hh","ipf.ind") %in% ls()))) {load("./db/d3.1_IPFindhh.RData")}
  
  #Getting types
  if(!file.exists('./db/d4.1_ipfsparse.RData')) {
    ipf.sparse <- format_fulltypes(ipf.ind, ipf.hh, indvars, hhvars)
    save(ipf.sparse, file = './db/d4.1_ipfsparse.RData')
    gc()
  } else {
    load('./db/d4.1_ipfsparse.RData')
  }
  
  #Joint fit
  ipf.jointfit <- joint_fit(tracts, ipf.ind, ipf.hh, ipf.sparse, indvars, hhvars, par=T)
  
  # Check results 
  ipf.fitplots <- pblapply(ipf.jointfit[1:6], function(x) {
    qplot(data=x$errors, x = b, y = bhat, geom='point') + theme_classic() + coord_fixed()
  })
  grid.arrange(grobs=ipf.fitplots)
  
  # Fit results
  ipf.errors <- rbindlist(pblapply(names(ipf.jointfit), function(x) {
    tmp <- ipf.jointfit[[x]]$errors
    tmp[ , otract := x]
  }))
  
  #Saving
  save(ipf.jointfit, file = "./db/d4.1_IPFjointfit.RData")
  
  #Cleanup
  rm(ipf.ind, ipf.hh, ipf.hhtots, ipf.indtots, cl, types)
} 

#### Fast Re-weighting of MCMC ####
if(!file.exists("./db/d4.2_MCMCjointfit.RData")) {
  # Loading database
  if(all(!(c("mcmc.hh.raked","mcmc.ind.raked") %in% ls()))) {load("./db/d3.2_MCMCrake.RData")}
  
  #Extracting the raked weights
  fm.hh <- as.formula(paste0("weights.raked~",paste(hhvars, collapse = '+')))
  fm.ind <- as.formula(paste0("weights.raked~",paste(indvars, collapse = '+')))
  
  mcmc.hh <- pblapply(mcmc.hh.raked, function(x) xtabs(fm.hh, data=x))
  mcmc.ind <- pblapply(mcmc.ind.raked, function(x) xtabs(fm.ind, data=x))
  rm(fm.hh,fm.ind)
  
  #Getting types
  if(!file.exists('./db/d4.2_MCMCsparse.RData')) {
    mcmc.sparse <- format_fulltypes(mcmc.ind, mcmc.hh, indvars, hhvars)
    save(mcmc.sparse, file= './db/d4.2_MCMCsparse.RData')
    gc()
  } else {
    load('./db/d4.2_MCMCsparse.RData')
  }
  
  #Running joint fit
  mcmc.jointfit <- joint_fit(tracts, mcmc.ind, mcmc.hh, mcmc.sparse, indvars, hhvars, par=T)
  
  
  # Check results 
  mcmc.fitplots <- pblapply(mcmc.jointfit[1:6], function(x) {
    qplot(data=x$errors, x = b, y = bhat, geom='point') + theme_classic() + coord_fixed()
  })
  grid.arrange(grobs=mcmc.fitplots)
  
  # Fit results
  mcmc.errors <- rbindlist(pblapply(names(mcmc.jointfit), function(x) {
    tmp <- mcmc.jointfit[[x]]$errors
    tmp[ , otract := x]
  }))
  
  #Saving
  save(mcmc.jointfit, file = "./db/d4.2_MCMCjointfit.RData")
  
  #Cleanup
  rm(list = ls()[grepl("mcmc",ls())])
  
}

#### Fast Re-weighting of MCMC w/ OD ####
if(!file.exists("./db/d4.2_MCMCODjointfit.RData")) {
  # Loading database
  if(all(!("mcmc.hh.raked" %in% ls()))) {load("./db/d3.2_MCMCrake.RData")}
  if(all(!("mcmc.indod" %in% ls()))) {load("./db/d3.2_MCMCindod.RData")}
  
  #Extracting the raked weights
  fm.hh <- as.formula(paste0("weights.raked~",paste(hhvars, collapse = '+')))
  
  mcmc.hh <- pblapply(mcmc.hh.raked, function(x) xtabs(fm.hh, data=x))
  mcmc.indod <- pblapply(tracts, function(tr) 
    prop.table(table(mcmc.indod[otract==tr,!c("otract","dtract")])) * sum(margins$indagesex[tr,,]))
  names(mcmc.indod) <- tracts
  rm(fm.hh)
  
  #Getting types
  if(!file.exists('./db/d4.2_MCMCODsparse.RData')) {
    mcmc.odsparse <- format_fulltypes(mcmc.indod, mcmc.hh, indvars, hhvars)
    save(mcmc.odsparse, file= './db/d4.2_MCMCODsparse.RData')
    gc()
  } else {
    load('./db/d4.2_MCMCODsparse.RData')
  }
  
  #Running joint fit
  mcmc.odjointfit <- joint_fit(tracts, mcmc.indod, mcmc.hh, mcmc.odsparse, indvars, hhvars, par=T)
  
  # Check results 
  mcmc.odfitplots <- pblapply(mcmc.odjointfit[1:6], function(x) {
    qplot(data=x$errors, x = b, y = bhat, geom='point') + theme_classic() + coord_fixed()
  })
  grid.arrange(grobs=mcmc.odfitplots)
  
  # Fit results
  mcmc.errors <- rbindlist(pblapply(names(mcmc.odjointfit), function(x) {
    tmp <- mcmc.odjointfit[[x]]$errors
    tmp[ , otract := x]
  }))
  
  #Saving
  save(mcmc.odjointfit, file = "./db/d4.2_MCMCODjointfit.RData")
  
  #Cleanup
  rm(list = ls()[grepl("mcmc.",ls())])
  
}

#### Fast Re-weighting of MCMC w/ OD ####
if(!file.exists("./db/d4.2_MCMCODjointfit.RData")) {
  # Loading database
  if(all(!("mcmc.hh.raked" %in% ls()))) {load("./db/d3.2_MCMCrake.RData")}
  if(all(!("mcmc.indod.raked" %in% ls()))) {load("./db/d3.2_MCMCindodrake.RData")}
  
  #Extracting the raked weights
  fm.hh <- as.formula(paste0("weights.raked~",paste(hhvars, collapse = '+')))
  fm.ind <- as.formula(paste0("weights.raked~",paste(indvars, collapse = '+')))
  
  mcmc.hh <- pblapply(mcmc.hh.raked, function(x) xtabs(fm.hh, data=x))
  mcmc.indod <- pblapply(mcmc.indod.raked, function(x) xtabs(fm.ind, data=x))
  rm(fm.hh, fm.ind)
  
  #Getting types
  if(!file.exists('./db/d4.2_MCMCODrakesparse.RData')) {
    mcmc.odsparse <- format_fulltypes(mcmc.indod, mcmc.hh, indvars, hhvars)
    save(mcmc.odsparse, file = './db/d4.2_MCMCODrakesparse.RData')
    gc()
  } else {
    load('./db/d4.2_MCMCODrakesparse.RData')
  }
  
  #Running joint fit
  mcmc.odrakejointfit <- joint_fit(tracts, mcmc.indod, mcmc.hh, mcmc.odsparse, indvars, hhvars, par=T)
  
  # Check results 
  mcmc.odfitplots <- pblapply(mcmc.odrakejointfit[1:6], function(x) {
    qplot(data=x$errors, x = b, y = bhat, geom='point') + theme_classic() + coord_fixed()
  })
  grid.arrange(grobs=mcmc.odfitplots)
  
  # Fit results
  mcmc.errors <- rbindlist(pblapply(names(mcmc.odrakejointfit), function(x) {
    tmp <- mcmc.odrakejointfit[[x]]$errors
    tmp[ , otract := x]
  }))
  
  #Saving
  save(mcmc.odrakejointfit, file = "./db/d4.2_MCMCODrakejointfit.RData")
  
  #Cleanup
  rm(list = ls()[grepl("mcmc.",ls())])
  
}

#### Fast Re-weighting of BN ####
if(!file.exists("./db/d4.3_BNjointfit.RData")) {
  # Loading database
  if(all(!(c("bn.hh.raked","bn.ind.raked") %in% ls()))) {load("./db/d3.3_BNrake.RData")}

  #Extracting the raked weights
  fm.hh <- as.formula(paste0("weights.raked~",paste(hhvars, collapse = '+')))
  fm.ind <- as.formula(paste0("weights.raked~",paste(indvars, collapse = '+')))
  
  bn.hh <- pblapply(bn.hh.raked, function(x) xtabs(fm.hh, data=x))
  bn.ind <- pblapply(bn.ind.raked, function(x) xtabs(fm.ind, data=x))
  rm(fm.hh,fm.ind)

  #Getting types
  if(!file.exists('./db/d4.3_BNsparse.RData')) {
    bn.sparse <- format_fulltypes(bn.ind, bn.hh, indvars, hhvars)
    save(bn.sparse, file= './db/d4.3_BNsparse.RData')
    gc()
  } else {
    load('./db/d4.3_BNsparse.RData')
  }
  
  #Running joint fit
  bn.jointfit <- joint_fit(tracts,bn.ind, bn.hh, bn.sparse, indvars, hhvars, par=T)
  
  # Check results 
  bn.fitplots <- pblapply(bn.jointfit[1:6], function(x) {
    qplot(data=x$errors, x = b, y = bhat, geom='point') + theme_classic() + coord_fixed()
  })
  grid.arrange(grobs=bn.fitplots)
  
  # Fit results
  bn.errors <- rbindlist(pblapply(names(bn.jointfit), function(x) {
    tmp <- bn.jointfit[[x]]$errors
    tmp[ , otract := x]
  }))
  
  #Saving
  save(bn.jointfit, file = "./db/d4.3_BNjointfit.RData")
  
  #Cleanup
  rm(list = ls()[grepl("bn",ls())])
  
}

#### Fast Re-weighting of BN w/ OD ####
if(!file.exists("./db/d4.3_BNODjointfit.RData")) {
  
  # Loading database
  if(all(!("bn.hh.raked" %in% ls()))) {load("./db/d3.3_BNrake.RData")}
  if(all(!("bn.indod" %in% ls()))) {load("./db/d3.3_BNindod.RData")}
  
  #Extracting the raked weights
  fm.hh <- as.formula(paste0("weights.raked~",paste(hhvars, collapse = '+')))
  
  bn.hh <- pblapply(bn.hh.raked, function(x) xtabs(fm.hh, data=x))
  bn.indod <- pblapply(tracts, function(tr) 
    prop.table(table(bn.indod[otract==tr,!c("otract","dtract")])) * sum(margins$indagesex[tr,,]))
  names(bn.indod) <- tracts
  rm(fm.hh)

  #Getting types
  if(!file.exists('./db/d4.3_BNODsparse.RData')) {
    bn.odsparse <- format_fulltypes(bn.indod, bn.hh, indvars, hhvars)
    save(bn.odsparse, file= './db/d4.3_BNODsparse.RData')
    gc()
  } else {
    load('./db/d4.3_BNODsparse.RData')
  }
  
  #Running joint fit
  bn.odjointfit <- joint_fit(tracts,bn.indod, bn.hh, bn.odsparse, indvars, hhvars, par=T)
  
  # Check results 
  bn.odfitplots <- pblapply(bn.odjointfit[1:6], function(x) {
    qplot(data=x$errors, x = b, y = bhat, geom='point') + theme_classic() + coord_fixed()
  })
  grid.arrange(grobs=bn.odfitplots)
  
  # Fit results
  bn.oderrors <- rbindlist(pblapply(names(bn.odjointfit), function(x) {
    tmp <- bn.odjointfit[[x]]$errors
    tmp[ , otract := x]
  }))
  
  #Saving
  save(bn.odjointfit, file = "./db/d4.3_BNODjointfit.RData")
  
  #Cleanup
  rm(list = ls()[grepl("bn",ls())])
  
}

#### Fast Re-weighting of raked BN w/ OD ####
if(!file.exists("./db/d4.3_BNODrakedjointfit.RData")) {
  
  # Loading database
  if(all(!("bn.hh.raked" %in% ls()))) {load("./db/d3.3_BNrake.RData")}
  if(all(!("bn.indod.raked" %in% ls()))) {load("./db/d3.3_BNindodrake.RData")}
  
  #Extracting the raked weights
  fm.hh <- as.formula(paste0("weights.raked~",paste(hhvars, collapse = '+')))
  fm.ind <- as.formula(paste0("weights.raked~",paste(indvars, collapse = '+')))
  
  bn.hh <- pblapply(bn.hh.raked, function(x) xtabs(fm.hh, data=x))
  bn.indod <- pblapply(bn.indod.raked, function(x) xtabs(fm.ind, data=x))
  
  #Getting types
  if(!file.exists('./db/d4.3_BNODrakesparse.RData')) {
    bn.odsparse <- format_fulltypes(bn.indod, bn.hh, indvars, hhvars)
    save(bn.odsparse, file = './db/d4.3_BNODrakesparse.RData')
    gc()
  } else {
    load('./db/d4.3_BNODrakesparse.RData')
  }
  
  #Running joint fit
  bn.odjointfit <- joint_fit(tracts,bn.indod, bn.hh, bn.odsparse, indvars, hhvars, par=T)
  
  # Check results 
  bn.odfitplots <- pblapply(bn.odjointfit[1:6], function(x) {
    qplot(data=x$errors, x = b, y = bhat, geom='point') + theme_classic() + coord_fixed()
  })
  grid.arrange(grobs=bn.odfitplots)
  
  # Fit results
  bn.oderrors <- rbindlist(pblapply(names(bn.odjointfit), function(x) {
    tmp <- bn.odjointfit[[x]]$errors
    tmp[ , otract := x]
  }))
  
  #Saving
  save(bn.odjointfit, file = "./db/d4.3_BNODrakejointfit.RData")
  
  #Cleanup
  rm(list = ls()[grepl("bn",ls())])
  
}

#### Cleanup ####
rm(list=ls())
