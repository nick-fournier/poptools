
#### Load packages
source("./s1.0_packages.R", echo=F)

#Fonts
library(extrafont)
#library(clpAPI)
library(Cairo)

if(.Platform$OS.type == "unix") face = "Times" else face = "Times New Roman"


#Save directory
dir <- "../pop-synth-paper/figures"

#### Checking each generation ####
if(!file.exists("./db/d6.0_genresults.RData")) {
  load("./db/d2.3_microdata.RData")
  load("./db/d2.4_lodes.RData")
  load("./db/d3.1_IPFlodes.RData")
  load("./db/d2.5_margins.RData")
  
  load("./db/d5.1_IPFpop.RData")
  load("./db/d5.1_IPFODpop.RData")
  load("./db/d5.2_MCMCpop.RData")
  load("./db/d5.2_MCMCODpop.RData")
  load("./db/d5.3_BNpop.RData")
  load("./db/d5.3_BNODpop.RData")
  
  indcols = c("AGE","GEND","RELATE","OCC","INDUS")
  hhcols = c("INCOME","RESTY","HHVEH","HHSIZ")
  
  pums = microdata$pums[,c("SERIALNO",c(indcols,hhcols)),with=F]

  check_results <- function(pop, pums, odi, indcols, hhcols) {
    #All columns combined
    cols = c(indcols,hhcols)
    #### Marginals ####
    print("Formatting marginals")
    #Check marginal results
    marg <- list("AGE" = apply(margins$indagesex, c(1,3), sum),
                 "GEND" = apply(margins$indagesex, c(1,2), sum),
                 "RELATE" = as.matrix(margins$indrelate),
                 "OCC"  = margins$indocc,
                 "INDUS"  = margins$indindus,
                 "INCOME" = as.matrix(margins$hhinc),
                 "RESTY" = margins$hhdwell,
                 "HHVEH" = margins$hhveh,
                 "HHSIZ" = rowSums(margins$hhvehsize,dims=2))
    
    #Aggregating population
    agg.hh <- rbindlist(lapply(hhcols, function(n) {
      cbind('Variable'=n,setNames(pop[!duplicated(HHID),.N,by=c(n,'otract')],c("var","otract","Synthetic")))
    }))
    agg.ind <- rbindlist(lapply(indcols, function(n) {
      cbind('Variable'=n,setNames(pop[,.N,by=c(n,'otract')],c("var","otract","Synthetic")))
    }))
    
    #Aggregating marginals
    marg <- rbindlist(lapply(marg, function(n) setNames(melt(n),c("otract","var","Marginals"))))
    marg[, otract := as.character(otract)]
    #Merging
    marg <- merge(rbind(agg.hh,agg.ind),marg, by=c("otract","var"))
    
    #### Cells ####
    print("Formatting microdata cells")
    #IND alone
    cells.ind<- merge(as.data.table(prop.table(table(pums[,indcols,with=F]))),
                      as.data.table(prop.table(table(pop[,indcols,with=F]))),
                      by = indcols)
    setnames(cells.ind,c("N.x","N.y"),c("PUMS","Synthetic"))
    #HH alone
    cells.hh <- merge(as.data.table(prop.table(table(pums[!duplicated(SERIALNO),hhcols,with=F]))),
                      as.data.table(prop.table(table(pop[!duplicated(HHID),hhcols,with=F]))),
                      by = hhcols)
    setnames(cells.hh,c("N.x","N.y"),c("PUMS","Synthetic"))
    #joint
    joint.cells <- merge(as.data.table(prop.table(table(pums[,!"SERIALNO"])))[N>0,],
                         as.data.table(prop.table(table(pop[,cols,with=F]))), by = cols, all.x=T)
    setnames(joint.cells,c("N.x","N.y"),c("PUMS","Synthetic"))
    
    #### OD ####
    print("Formatting origin-destinations")
    #Getting proportions
    pop <- pop[ , lapply(.SD, as.character)]
    res.od = as.data.table(table(pop[INDUS!="NOIND",.(otract,dtract)]))
    res.oi = as.data.table(table(pop[INDUS!="NOIND",.(otract,INDUS)]))
    res.di = as.data.table(table(pop[INDUS!="NOIND",.(dtract,INDUS)]))
    
    lodes.od = setNames(melt(as.data.table(lodes$mtxOD,keep.rownames = T), id.vars='rn'),
                        c("otract","dtract","LODES"))[,  dtract := as.character(dtract)]
    lodes.oi = setNames(melt(as.data.table(lodes$mtxOI,keep.rownames = T), id.vars='rn'),
                        c("otract","INDUS","LODES"))
    lodes.di = setNames(melt(as.data.table(lodes$mtxDI,keep.rownames = T), id.vars='rn'),
                        c("dtract","INDUS","LODES"))
    
    #merging proportions to gibbs results
    res.od = merge(res.od[otract!="00000000000" & dtract!="00000000000",],
                   lodes.od[otract!="00000000000" & dtract!="00000000000",], by=c("otract","dtract"), all=T)
    res.oi = merge(res.oi[otract!="00000000000",], lodes.oi[otract!="00000000000",], by=c("otract","INDUS"), all=T)
    res.di = merge(res.di[dtract!="00000000000",], lodes.di[dtract!="00000000000",], by=c("dtract","INDUS"), all=T)
    
    #Checking for zeros
    res.od = res.od[ , (c('N','LODES')) := lapply(.SD, function(x) ifelse(is.na(x),0,x)), .SDcols = c('N','LODES')]
    res.oi = res.oi[ , (c('N','LODES')) := lapply(.SD, function(x) ifelse(is.na(x),0,x)), .SDcols = c('N','LODES')]
    res.di = res.di[ , (c('N','LODES')) := lapply(.SD, function(x) ifelse(is.na(x),0,x)), .SDcols = c('N','LODES')]
    
    
    "hh cells" = cells.hh[,.(PUMS,Synthetic)][ , sqrt(sum((PUMS-Synthetic)^2)/.N)/mean(PUMS)]
    "ind cells" = cells.ind[,.(PUMS,Synthetic)][ , sqrt(sum((PUMS-Synthetic)^2)/.N)/mean(PUMS)]
    
    
    #### SRMSE ####
    print("Finalizing")
    SRMSE <- c("Marginal" = marg[,sqrt(sum((Marginals - Synthetic)^2)/.N)/mean(Marginals)],
               "Multilevel joint cells" = joint.cells[ , sqrt(sum((PUMS-Synthetic)^2)/.N)/mean(PUMS)],
               "Independent joint cells" = rbind(cells.ind[,.(PUMS,Synthetic)],cells.hh[,.(PUMS,Synthetic)])[ , sqrt(sum((PUMS-Synthetic)^2)/.N)/mean(PUMS)],
               "hh cells" = cells.hh[,.(PUMS,Synthetic)][ , sqrt(sum((PUMS-Synthetic)^2)/.N)/mean(PUMS)],
               "ind cells" = cells.ind[,.(PUMS,Synthetic)][ , sqrt(sum((PUMS-Synthetic)^2)/.N)/mean(PUMS)],
               "OI" = res.oi[,sqrt(sum((LODES-N)^2)/.N)/mean(LODES)],
               "DI" = res.di[,sqrt(sum((LODES-N)^2)/.N)/mean(LODES)],
               "OD" = res.od[,sqrt(sum((LODES-N)^2)/.N)/mean(LODES)],
               "ODI" = rbindlist(list(res.od[,.(N,LODES)],res.oi[,.(N,LODES)],res.di[,.(N,LODES)]))[,sqrt(sum((LODES-N)^2)/.N)/mean(LODES)])
    SRMSE <- as.data.table(SRMSE, keep.rownames = T)
    
    #Output
    output <- list('marginals' = marg,
                   'cells' = rbind(cells.ind[,.(PUMS,Synthetic)],cells.hh[,.(PUMS,Synthetic)]),
                   'jointcells' = joint.cells,
                   'odi' = list('od' = res.od, 'oi' = res.oi, 'di' = res.di),
                   'SRMSE' = SRMSE)
    return(output)
  }
  
  results.ipfod = check_results(ipfod.pop, pums, ipf.lodes, indcols, hhcols)
  results.ipf   = check_results(ipf.pop, pums, lodes, indcols, hhcols)
  
  results.mcmc  = check_results(mcmc.pop, pums, lodes, indcols, hhcols)
  results.mcmcod  = check_results(mcmc.odpop, pums, lodes, indcols, hhcols)
  
  results.bn    = check_results(bn.pop, pums, lodes, indcols, hhcols)
  results.bnod    = check_results(bn.odpop, pums, lodes, indcols, hhcols)
  
  results <- list('IPF' = results.ipf,
                  'IPFOD' = results.ipfod,
                  'MCMC' = results.mcmc, 
                  'MCMCOD' = results.mcmcod, 
                  'BN' = results.bn, 
                  'BNOD' = results.bnod)
  
  save(results, file = "./db/d6.0_genresults.RData")
} else {
  load("./db/d6.0_genresults.RData")
  
  #Errors
  SRMSE <- rbindlist(lapply(results,
                            function(x) setNames(data.table(t(x$SRMSE$SRMSE)),x$SRMSE$rn)))
  SRMSE <- cbind("Pop" = names(results), SRMSE)
  SRMSE <- t(SRMSE[,-1])
  colnames(SRMSE) <- names(results)
  SRMSE
}


#### Compare methods ####
if(!file.exists("./db/d6.0_compare_all.RData")){
  sourceCpp("./cpp_IPU.cpp")
  format_fulltypes <- function(ind, hh, indvars, hhvars){
    # Getting HH and Person types from IPF
    if(is.list(hh)){
      print("Aggregating person matrices")
      #Summing into single matrix
      indtypes <- Reduce("+", ind)
      hhtypes <- Reduce("+", hh)
    } else {
      print("Using existing total weights")
      indtypes <- ind
      hhtypes <- hh
    }
    print("Formatting sparse matrix from PUMS")
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
  format_bvector <- function(indcons, hhcons, sparsetypes, typevars, indvars, hhvars){
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
  glmnet_fit <- function(A, b, pumstypes, l = 0, a = 1, tol=1e-10, maxit=1e5) {
    #print("Running least-squares fitting")
    time <- Sys.time()
    #Run NNLS with glmnet
    glmnetmod <- glmnet(A, b, alpha = a, lambda = l, lower.limits=0, intercept=F, thresh = tol, maxit=maxit)
    
    #Retrieve the value of the objective function after optimization.
    X <- coef(glmnetmod)[,1][-1]
    #names(X) <- colnames(A)
    
    #getting results
    bhat <- as.vector(A %*% X)
    res <- data.table(varname = names(b), bhat, b, vartype = gsub('[0-9]+', '', names(b)))
    res[ , sqerror := (b-bhat)^2]
    res[ , esttype := "CCD"]
    RMSE <- sqrt(sum(res$sqerror)/nrow(res))
    RMSN <- RMSE/mean(res$b)
    #plot(bhat~b, res)
    #plot(bhat~b, res[!(bhat>0 & b==0), ])
    
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
  
  # Formatting full constraints and matrix
  if(!file.exists("./db/d6.0_sparse.RData")) {
    #Load general data
    load("./db/d2.3_microdata.RData")
    load('./db/d3.1_IPFindhh.RData')
    indcols = c("AGE","GEND","RELATE","OCC","INDUS")
    hhcols = c("INCOME","RESTY","HHSIZ","HHVEH")
    pums = microdata$pums[,c("SERIALNO",c(indcols,hhcols)),with=F]
    
    types <- format_fulltypes(ipf.indtots, ipf.hhtots, indcols, hhcols)
    
    # sparse <- types$sparse
    # 
    # chnk <- split(sample(1:ncol(sparse),ncol(sparse)), ceiling(seq_along(1:ncol(sparse))/50000))
    # 
    # if(.Platform$OS.type == "unix") {
    #   cl <- makeCluster(detectCores()-1, "FORK")
    # } else {
    #   cl <- makeCluster(detectCores(logical = F)-1, "PSOCK")
    #   clusterEvalQ(cl,c(library(Matrix)))
    #   clusterExport(cl,
    #                 list("sparse"),
    #                 envir=environment())
    # }
    # 
    # NoDup <- pblapply(chnk, function(x) t(mgcv::uniquecombs(t(sparse[,x]))), cl = cl)
    # stopCluster(cl)
    # 
    # sn <- unlist(lapply(NoDup, function(x) dimnames(x)$serial))
    # 
    # sparse <- sparse[ , dimnames(sparse)$serial %in% sn]
    # 
    # types$sparse <- sparse
    
    save(types, file = "./db/d6.0_sparse.RData")
  } else {
    load("./db/d6.0_sparse.RData")
    A <- types$sparse
    b <- format_bvector(ipf.indtots, ipf.hhtots, types, rownames(A), indcols, hhcols)
  }

  
  #####Comparison for whole region
  if(!file.exists("./db/d6.0_compare_ipu.RData")) {
    res_ipu <- ipu_fit(A, b, verb=T, c=5e-4, cc=1, maxit=100)
    save(res_ipu, file = "./db/d6.0_compare_ipu.RData")
  } else {
    load("./db/d6.0_compare_ipu.RData")
  }
  
  # if(!file.exists("./db/d6.0_compare_nnls.RData")) {
  #   res_nnls <- nnls_fit(A, b)
  #   save(res_nnls, file = "./db/d6.0_compare_nnls.RData")
  # } else {
  #   load("./db/d6.0_compare_nnls.RData")
  # }
  
  if(!file.exists("./db/d6.0_compare_nnld.RData")) {
    res_nnld <- nnld_fit(A, b)
    save(res_nnld, file = "./db/d6.0_compare_nnld.RData")
  } else {
    load("./db/d6.0_compare_nnld.RData")
  }
  
  if(!file.exists("./db/d6.0_compare_ccd.RData")) {
    res_ccd <- glmnet_fit(A, b, types$pumstypes, a=1, l=0, tol=1e-14)
    save(res_ccd, file = "./db/d6.0_compare_ccd.RData")
  } else {
    load("./db/d6.0_compare_ccd.RData")
  }
  
  if(!file.exists("./db/d6.0_compare_lasso.RData")) {
    res_lasso <- glmnet_fit(A, b, types$pumstypes, a=1, l=1)
    save(res_lasso, file = "./db/d6.0_compare_lasso.RData")
  } else {
    load("./db/d6.0_compare_lasso.RData")
  }
  if(!file.exists("./db/d6.0_compare_ridge.RData")) {
    res_ridge <- glmnet_fit(A, b, types$pumstypes, a=0, l=1, tol=1e-12, maxit=1e7)
    res_ridge$errors$esttype <- "NNLS"
    save(res_ridge, file = "./db/d6.0_compare_ridge.RData")
  } else {
    load("./db/d6.0_compare_ridge.RData")
  }
  
  compare <- c(list("ipu"=res_ipu,
                    "nnls"=res_ridge,
                    "nlad"=res_nnld,
                    "ccd"=res_ccd))
  
  ####Combining
  # compare <- c(list("nnls"=res_nnls,"nnld"=res_nnld,"ipu"=res_ipu,"ccd"=res_ccd))
  
  #Creating performance table
  perform <- lapply(compare, function(res) data.table("Method"=res$errors$esttype[1],
                                                      "RMSN"=res$RMSN,
                                                      "RMSE"=res$RMSE,
                                                      "MAPE"=res$MAPE,
                                                      "Time"=res$time))
  perform <- do.call(rbind, perform)
  perform[ , Scale := Time/as.numeric(min(Time))]
  perform
  
  save(compare, perform, file= "./db/d6.0_compare_all.RData")
} else{
  load("./db/d6.0_compare_all.RData")
}


#### PLOT COMPARE METHODS ####
plots.compare <- lapply(compare, function(res) {
  dat = res$errors
  dat$bhat = dat$bhat/1000
  dat$b = dat$b/1000
  
  RMSN =  format(dat[ , sqrt(sum((bhat - b)^2)/.N)/mean(b)], digits=3, scientific=T)
  slope = round(coef(lm(bhat~b, data=dat))[2],5)
  R = round(summary(lm(bhat~b, data=dat))$r.squared,5)
  
  ggplot(data=dat, aes(x=b, y=bhat)) + 
    geom_smooth(method="lm", color="gray", size=0.5, fullrange=T) +
    geom_point(shape=1, size=1.5, color="black") +
    scale_y_continuous('Estimated IPF types (thousands)', breaks = seq(0,1000,50), label = scales::comma) +
    scale_x_continuous('Target IPF types (thousands)', breaks = seq(0,1000,50), label = scales::comma) +
    #annotate("text", x=0.7*max(dat$bhat), y=50, label = paste("Slope=",slope), family=face) +
    #annotate("text", x=0.7*max(dat$bhat), y=25, label = paste("R^{2}~'='~",R), parse=T, family=face) +
    annotate("text", x=0.7*max(dat$bhat), y=0, label = paste("RMSN =",RMSN), family=face) +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + coord_fixed() +
    theme(legend.position = "none", text=element_text(family=face, size=14))})
grid.arrange(grobs=plots.compare, ncol = 2)

#save
for(n in names(compare)) {
  ggsave(paste0(dir,"/method_",n,".png"), type="cairo", plot = plots.compare[[n]], width = 3, height = 3, units = "in", dpi=300)  
}

#Zoom
plots.comparezoom <- lapply(compare, function(res) {
  ggplot(data=res$errors, aes(x=b, y=bhat)) +
    geom_smooth(method="lm", color="gray", size=0.5, fullrange=T) +
    geom_point(size=1.5, shape=1) +
    scale_y_continuous('Estimated IPF types', limits = c(0,5), breaks = seq(0,5,1)) + 
    scale_x_continuous('Target IPF types', limits = c(0,5), breaks = seq(0,5,1)) +
    #annotate("text", x=2, y=0.50, label = paste("Slope=",round(coef(lm(bhat~b, data=res$errors))[2],5)), family=face) +
    #annotate("text", x=2, y=0.25, label = paste("R^{2}~'='~",round(summary(lm(bhat~b, res$errors))$r.squared,5)), parse=T, family=face) +
    annotate("text", x=3, y=0, label = paste("RMSN =",format(res$RMSN,digits=3,scientific=T)), family=face) +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + coord_fixed() +
    theme(legend.position = "none", text=element_text(family=face, size=14))})
grid.arrange(grobs=plots.comparezoom, ncol = 2)

#save
for(n in names(compare)) {
  ggsave(paste0(dir,"/methodzoom_",n,".png"), type="cairo", plot = plots.comparezoom[[n]], width = 3, height = 3, units = "in", dpi=300)
}



#### Plot final gens ####
keepers <- c("IPF","IPFOD")

####Marginal totals
plots.margtot <- lapply(results, function(res) {
  
  dat = res$marginals[ , lapply(.SD, sum), by = c("var","Variable"), .SDcols = c("Synthetic","Marginals")]
  RMSN = round(dat[ , sqrt(sum((Marginals - Synthetic)^2)/.N)/mean(Marginals)], 4)
  
  ggplot(data=dat, aes(x=Marginals/1000000, y=Synthetic/1000000, color=Variable)) +
    geom_smooth(method="lm", color="gray", size=0.5, fullrange=T) +
    geom_point(size=1.5, shape=1) + theme_classic() + coord_fixed() + 
    # annotate("text", x=0.8*max(dat$Marginals/1000), y=0.8,
    #          label = paste("Slope=",round(coef(lm(Synthetic~Marginals, data=dat))[2],4)), family=face) +
    # annotate("text", x=0.8*max(dat$Marginals/1000), y=0.4,
    #          label = paste("R^{2}~'='~",round(summary(lm(Synthetic~Marginals, dat))$r.squared,4)), parse=T, family=face) +
    annotate("text", x=0.8*max(dat$Marginals/1000000), y=0, 
             label = paste("RMSN =",RMSN), family=face) +
    scale_x_continuous("Marginal totals (millions)", breaks = seq(0, 3, by = 0.5), limits = c(0,3)) +
    scale_y_continuous("Synthetic totals (millions)", breaks = seq(0, 3, by = 0.5), limits = c(0,3)) +
    scale_color_manual("Variable",
                       breaks = c("HHVEH","HHSIZ","RESTY","INCOME","INDUS","OCC","RELATE","GEND","AGE"),
                       labels = c("Household vehicles","Household size","Dwelling type","Income","Industry","Occupation","Relationship","Gender","Age"),
                       values = c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666','#e31a1c')) +
    theme(text=element_text(family=face,size=12),
          axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "pt")),
          legend.title = element_blank(),
          legend.background = element_rect(fill="transparent"),
          legend.key.height=unit(0.65,"line"),
          legend.position = c(.3,0.7))
})
plots.margtot$IPF

#save
for(n in keepers) {
  ggsave(paste0(dir,"/margtot_",n,".png"), type="cairo", plot = plots.margtot[[n]], width = 3, height = 3, units = "in", dpi=300)  
}

####Marginal tracts
plots.margtracts <- lapply(results, function(res) {
         RMSN = round(res$marginals[ , sqrt(sum((Marginals - Synthetic)^2)/.N)/mean(Marginals)], 4)
         
         ggplot(data=res$marginals, aes(x=Marginals/1000, y=Synthetic/1000)) +
           geom_smooth(method="lm", color="gray", size=0.5, fullrange=T) +
           geom_point(size=0.001) + theme_classic() + coord_fixed() + 
           # annotate("text", x=0.8*max(res$marginals$Marginals/1000), y=0.8,
           #          label = paste("Slope=",round(coef(lm(Synthetic~Marginals, data=res$marginals))[2],4)), family=face) +
           # annotate("text", x=0.8*max(res$marginals$Marginals/1000), y=0.4,
           #          label = paste("R^{2}~'='~",round(summary(lm(Synthetic~Marginals, res$marginals))$r.squared,4)), parse=T, family=face) +
           annotate("text", x=0.8*max(res$marginals$Marginals/1000), y=0, 
                    label = paste("RMSN =",RMSN), family=face) +
           scale_x_continuous("Marginal totals (thousands)", breaks = seq(0, 6, by = 1), label = scales::comma, limits = c(0,6)) +
           scale_y_continuous("Synthetic totals (thousands)", breaks = seq(0, 6, by = 1), label = scales::comma, limits = c(0,6)) +
         theme(text=element_text(family=face,size=12), 
               axis.title.y = element_text(margin = unit(c(0, 14, 0, 0), "pt"))
               )
       })
#plots.margtracts$IPF

#save
for(n in keepers) {
  ggsave(paste0(dir,"/margtracts_",n,".png"), type="cairo", plot = plots.margtracts[[n]], width = 3, height = 3, units = "in", dpi=300)  
}


####Cells
plots.cells <- lapply(results,
                     function(res) {
                       RMSN = round(res$cells[ , sqrt(sum((PUMS - Synthetic)^2)/.N)/mean(PUMS)], 4)
                       
                       ggplot(data=res$cells, aes(x=PUMS, y=Synthetic)) +
                         geom_smooth(method="lm", color="gray", size=0.5, fullrange=T) +
                         geom_point(size=0.001) + theme_classic() + coord_fixed() + 
                         # annotate("text", x=0.8*max(res$cells$PUMS), y=0.01,
                         #          label = paste("Slope=",round(coef(lm(Synthetic~PUMS, data=res$cells))[2],4)), family=face) +
                         # annotate("text", x=0.8*max(res$cells$PUMS), y=0.005,
                         #          label = paste("R^{2}~'='~",round(summary(lm(Synthetic~PUMS, res$cells))$r.squared,4)), parse=T, family=face) +
                         annotate("text", x=0.7*max(res$cells$PUMS), y=0, 
                                  label = paste("RMSN =",RMSN), family=face) +
                         scale_x_continuous("PUMS proportions", breaks = seq(0, 1, by = 1e-2), limits = c(0,0.07)) +
                         scale_y_continuous("Synthetic proportions", breaks = seq(0, 1, by = 1e-2), limits = c(0,0.07)) +
                         theme(text=element_text(family=face,size=12), 
                               #axis.title.y = element_text(margin = unit(c(0, 10, 0, 0), "pt"))
                               )
                     })
plots.cells$IPF

#save
for(n in keepers) {
  ggsave(paste0(dir,"/cells_",n,".png"), type="cairo", plot = plots.cells[[n]], width = 3, height = 3, units = "in", dpi=300)  
}

####OD
plots.odi <- lapply(results,
                   function(res) {
                     #Origin-destination
                     RMSN = lapply(res$odi, function(x) round(x[ , sqrt(sum((LODES - N)^2)/.N)/mean(N)], 4))
                     slope = lapply(res$odi, function(x) round(coef(lm(LODES~N, data=x))[2],4))
                     R = lapply(res$odi, function(x) round(summary(lm(LODES~N, data=x))$r.squared,4))
                     
                     plot.OD <- ggplot(data=res$odi$od,aes(x=LODES, y=N)) +
                       geom_smooth(method="lm", color="gray", size=0.5, fullrange=T) +
                       geom_point(size=0.1) +
                       # annotate("text", x=0.7*max(res$odi$od$LODES), y=1000, label = paste("Slope=",slope$od), family=face) +
                       # annotate("text", x=0.7*max(res$odi$od$LODES), y=500, label = paste("R^{2}~'='~",R$od), parse=T, family=face) +
                       annotate("text", x=0.7*max(res$odi$od$LODES), y=0, label = paste("RMSN =",RMSN$od), family = face) + 
                       scale_x_continuous("LODES", labels = scales::comma, breaks = seq(0, max(res$odi$od$LODES), 100), limits=c(0,max(res$odi$od$LODES))) + 
                       scale_y_continuous("Synthetic", labels = scales::comma, breaks = seq(0, max(res$odi$od$LODES), 100), limits=c(0,max(res$odi$od$LODES))) +
                       theme_classic() + coord_fixed() + theme(text=element_text(family=face, size=12),
                                                               axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "pt")))
                     #plot.OD
                     #Destination by industry
                     plot.DI <- ggplot(data=res$odi$di, aes(x=LODES, y=N)) +#, color=INDUS)) + 
                       geom_smooth(method="lm", color="gray", size=0.5, fullrange=T) +
                       geom_point(size=0.5) +
                       # annotate("text", x=0.7*max(res$odi$di$LODES), y=5000, label = paste("Slope=",slope$di), family=face) +
                       # annotate("text", x=0.7*max(res$odi$di$LODES), y=2500, label = paste("R^{2}~'='~",R$di), parse=T, family=face) +
                       annotate("text", x=0.7*max(res$odi$di$LODES), y=0, label = paste("RMSN =",RMSN$di), family = face) + 
                       scale_x_continuous("LODES", labels = scales::comma, breaks = seq(0, max(res$odi$di$LODES), 10000)) + 
                       scale_y_continuous("Synthetic", labels = scales::comma, breaks = seq(0, max(res$odi$di$LODES),10000)) +
                       theme_classic() + coord_fixed() + theme(text=element_text(family=face, size=12))
                     #plot.DI
                     #Origin by industry
                     plot.OI <- ggplot(data=res$odi$oi, aes(x=LODES, y=N)) + #, color = INDUS)) + 
                       geom_smooth(method="lm", color="gray", size=0.5, fullrange=T) +
                       geom_point(size=0.5) +
                       #annotate("text", x=0.7*max(res$odi$oi$LODES), y=1000, label = paste("Slope=",slope$oi), family=face) +
                       #annotate("text", x=0.7*max(res$odi$oi$LODES), y=500, label = paste("R^{2}~'='~",R$oi), parse=T, family=face) +
                       annotate("text", x=0.7*max(res$odi$oi$LODES), y=0, label = paste("RMSN =",RMSN$oi), family = face) + 
                       scale_x_continuous("LODES", labels = scales::comma, breaks = seq(0, max(res$odi$oi$LODES), 500)) + 
                       scale_y_continuous("Synthetic", labels = scales::comma, breaks = seq(0, max(res$odi$oi$LODES), 500)) +
                       theme_classic() + coord_fixed() + theme(text=element_text(family=face, size=12),
                                                               axis.title.y = element_text(margin = unit(c(0, 8, 0, 0), "pt")))
                     #plot.OI
                     
                     return(list("OD"=plot.OD,"DI"=plot.DI,"OI"=plot.OI))
                     })

plots.odi$IPFOD$OI

#save
for(n in keepers) {
  for(od in names(plots.odi[[n]])){ 
    ggsave(paste0(dir,"/odi_",n,"_",od,".png"), type="cairo", plot = plots.odi[[n]][[od]], width = 3, height = 3, units = "in", dpi=300)
  }
}


