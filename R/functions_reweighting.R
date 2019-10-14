
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


