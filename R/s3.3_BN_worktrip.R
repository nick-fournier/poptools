##
# Population Synthesis for Boston Metropolitan Region, by Nicholas Marc Fournier
# Last updated January 2017.
#
# NOTICE:  All information, intellectual and technical concepts contained herein is,
# and remains the property of Nicholas Marc Fournier. Dissemination of this 
# information or reproduction of this material is strictly forbidden unless
# prior written permission is obtained from Nicholas Marc Fournier
##

#### Setting up functions ####
source("./s1.0_packages.R", echo=F)
sourceCpp("cpp_gibbs.cpp")

fun_raketract.singledim <- function(tr, freq.dat, marg.dat, bnd=c(0,Inf), trm=NULL) {
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
    #Raking the data, strictest standards
    raked.dat <- suppressWarnings(
      tryCatch(calibrate(design.dat, formula = form, population = pop, calfun = "raking", bounds = bnd, trim = trm),
               error = function(e) "error"))
    #Try again with no lower boundaries
    if(all(raked.dat == "error")) {
      print("Relaxed boundaries")
      raked.dat <- suppressWarnings(
        tryCatch(calibrate(design.dat, formula = form, population = pop, calfun = "linear", bounds = c(0,1000), force=T),
                 error = function(e) "error"))
    }
    #Putting the weights back into the original data
    freq.dat.raked <- cbind("otract" = tr, freq.dat, "weights.raked" = weights(raked.dat))
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
fun_plotnetwork <- function(structure, ht = "400px"){
  nodes.uniq <- unique(c(structure$arcs[,1], structure$arcs[,2]))
  nodes <- data.frame(id = nodes.uniq,
                      label = nodes.uniq,
                      color = "darkturquoise",
                      font.face = "serif",
                      shadow = TRUE,
                      size = 10)
  
  edges <- data.frame(from = structure$arcs[,1],
                      to = structure$arcs[,2],
                      arrows = "to",
                      smooth = TRUE,
                      shadow = TRUE,
                      color = "black")
  
  return(visNetwork(nodes, edges, height = ht, width = "100%"))
}


#Load general data
set.seed(12345)
load("./db/d2.3_microdata.RData")
load("./db/d2.5_margins.RData")

#### Bayesian Network for Households & Individuals ####
if(!file.exists("./db/d3.3_BNindhh.RData")) {
  #Data
  data.ind <- microdata$pums[,.(AGE,GEND,RELATE,OCC,INDUS)]
  data.hh <- microdata$pums[!duplicated(SERIALNO),.(INCOME,RESTY,HHVEH,HHSIZ)]

  # Individuals ####
  #Set up root graph guess
  dag.ind <- empty.graph(names(data.ind))
  modelstring(dag.ind) <- "[AGE][GEND][RELATE|AGE:GEND][OCC|AGE:GEND:INDUS][INDUS|AGE:GEND]"
  #Learning the structure
  dag.ind.tabu <- tabu(data.ind, start = dag.ind, tabu = 100)
  fun_plotnetwork(dag.ind.tabu)
  #Estimating the probabilities
  bn.ind.fit <- bn.fit(dag.ind.tabu, data.ind)
  #Simulate full population
  bn.ind = as.data.table(rbn(bn.ind.fit, sum(margins$indagesex), dag.ind.tabu))

  # Households ####
  #Set up root graph guess
  dag.hh.start <- empty.graph(names(data.hh))
  modelstring(dag.hh.start) <- "[HHSIZ][INCOME|HHSIZ][RESTY|HHSIZ:INCOME][HHVEH|INCOME:HHSIZ]"
  #Learning the structure
  dag.hh.tabu <- tabu(data.hh, start = dag.hh.start, tabu = 100)
  fun_plotnetwork(dag.hh.tabu)
  
  modelstring(dag.hh.tabu)
  fun_plotnetwork(model2network("[Income][Household size|Income][Household vehicles|Income:Household size][Dwelling Type|Income:Household vehicles:Household size]"))
  
  #Estimating the probabilities
  bn.hh.fit <- bn.fit(dag.hh.tabu, data.hh)
  #Simulate full population
  bn.hh = as.data.table(rbn(bn.hh.fit, sum(margins$hhinc), dag.hh.tabu))
  
  # Check results ####
  bn.ind.res <- merge(as.data.table(prop.table(table(bn.ind))),
                      as.data.table(prop.table(table(data.ind))),by=c("AGE","GEND","RELATE","INDUS","OCC"))
  
  bn.hh.res <- merge(as.data.table(prop.table(table(bn.hh))),
                      as.data.table(prop.table(table(data.hh))),by=c("INCOME","RESTY","HHVEH","HHSIZ"))
  #
  ggplot(data=bn.ind.res, aes(x=N.x, y=N.y)) + geom_point() + coord_fixed() + theme_classic()
  ggplot(data=bn.hh.res, aes(x=N.x, y=N.y)) + geom_point() + coord_fixed() + theme_classic()
  
  #Output
  save(bn.ind, bn.hh, bn.ind.res, bn.hh.res, file = "./db/d3.3_BNindhh.RData")
  #
  rm(list = setdiff(ls(),c("bn.hh","bn.ind","bn.ind.fit","bn.hh.fit","bn.ind.res","bn.hh.res",
                           "margins","microdata","tables","tracts",
                           "fun_raketract.singledim","fun_plotnetwork")))
}

#### Bayesian Network for individuals with OD ####
if(!file.exists("./db/d3.3_BNindod.RData")) {
  load("./db/d3.3_BNindhh.RData")
  load("./db/d2.4_lodes.RData")
  load("./db/d3.2_MCMClodes.RData")
  
  #Data
  data.ind <- microdata$pums[,.(AGE,GEND,RELATE,OCC,INDUS)]
  
  # Individuals ####
  #Set up root graph guess
  dag.ind.guess <- empty.graph(names(data.ind))
  modelstring(dag.ind.guess) <- "[AGE][GEND][RELATE|AGE:GEND][OCC|AGE:GEND:INDUS][INDUS|AGE:GEND]"
  
  fun_plotnetwork(dag.ind.guess)
  
  #Learning the structure
  dag.ind.tabu <- tabu(data.ind, start = dag.ind.guess, tabu = 100)
  fun_plotnetwork(dag.ind.tabu)
  #Estimating the probabilities
  bn.ind.fit <- bn.fit(dag.ind.tabu, data.ind)
  #Simulate full population
  bn.ind.sim = as.data.table(rbn(bn.ind.fit, sum(margins$indagesex), dag.ind.tabu))
  
  # Check BN results ####
  bn.ind.res <- merge(as.data.table(prop.table(table(bn.ind.sim))),
                      as.data.table(prop.table(table(bn.ind.sim))),by=c("AGE","GEND","RELATE","INDUS","OCC"))
  
  #
  ggplot(data=bn.ind.res, aes(x=N.x, y=N.y)) + geom_point() + coord_fixed() + theme_classic()
  
  
  # Setup LODES conditionals ####
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
  mtxOD <- mtxOD[c("00000000000",tracts),c("00000000000",tracts)]
  
  #Get conditional probabilities
  INDUS.otract = aperm(mtxOI[c("00000000000",tracts),], c(2,1))
  INDUS.otract = sweep(INDUS.otract, 2, colSums(INDUS.otract), "/")

  INDUS.dtract = aperm(mtxDI[c("00000000000",tracts),], c(2,1))
  INDUS.dtract = sweep(INDUS.dtract, 2, colSums(INDUS.dtract), "/")
  
  INDUS.dtract.otract = aperm(mcmc.odi[c("00000000000",tracts),c("00000000000",tracts),],c(3,2,1))
  INDUS.dtract.otract = sweep(INDUS.dtract.otract, 2:3, colSums(INDUS.dtract.otract), "/")
  
  dtract.otract = aperm(mtxOD[c("00000000000",tracts),c("00000000000",tracts)], c(2,1))
  dtract.otract = sweep(dtract.otract, 2, colSums(dtract.otract), "/")
  
  #Adding into conditional
  conds <- list("GEND" = aperm(cpt(GEND~INDUS,data.ind),c(2,1)),
                "OCC" = aperm(cpt(OCC~GEND+INDUS,data.ind),c(3,1,2)),
                "AGE" = aperm(cpt(AGE~GEND+OCC+INDUS,data.ind),c(4,1,2,3)),
                "RELATE" = aperm(cpt(RELATE~AGE+GEND+OCC,data.ind),c(4,1,2,3)),
                # "AGE" = aperm(cpt(AGE~GEND+OCC,data.ind),c(3,1,2)),
                # "RELATE" = aperm(cpt(RELATE~AGE+GEND,data.ind),c(3,1,2)),
                "otract" = as.table(rowSums(mtxOI[c("00000000000",tracts),])),
                "INDUS" = INDUS.dtract.otract,
                #"INDUS" = INDUS.dtract,
                "dtract" = dtract.otract)
  
  #If there are any NaN, fill in with very very very tiny number to get it to run.
  conds$INDUS[!is.finite(conds$INDUS)] <- .Machine$double.xmin
  conds$INDUS <- sweep(conds$INDUS, 2:3, colSums(conds$INDUS), "/")
  
  conds$AGE <- sweep(conds$AGE, 2:4, colSums(conds$AGE), "/")
  conds$AGE[!is.finite(conds$AGE)] <- .Machine$double.xmin
  conds$AGE <- sweep(conds$AGE, 2:4, colSums(conds$AGE), "/")
  
  conds$RELATE <- sweep(conds$RELATE, 2:4, colSums(conds$RELATE), "/")
  conds$RELATE[!is.finite(conds$RELATE)] <- .Machine$double.xmin
  conds$RELATE <- sweep(conds$RELATE, 2:4, colSums(conds$RELATE), "/")

  rm(INDUS.dtract.otract, dtract.otract, INDUS.dtract, INDUS.otract)

  #Custom fit based on modifed Tabu search
  fun_plotnetwork(dag.ind.tabu)
  modelstring(dag.ind.tabu)
  
  dag.indod <- model2network("[GEND|INDUS][OCC|GEND:INDUS][AGE|GEND:OCC:INDUS][RELATE|AGE:GEND:OCC][otract][INDUS|dtract:otract][dtract|otract]")
  #dag.indod <- model2network("[GEND|INDUS][OCC|GEND:INDUS][AGE|GEND:OCC][RELATE|AGE:GEND][otract][INDUS|dtract:otract][dtract|otract]")
  #dag.indod <- model2network("[GEND|INDUS][OCC|GEND:INDUS][AGE|GEND:OCC][RELATE|AGE:GEND][otract][INDUS|dtract][dtract|otract]")
  fun_plotnetwork(dag.indod)
  
  fun_plotnetwork(model2network("[Gender|Industry][Occupation|Gender:Industry][Age|Gender:Occupation:Industry][Relationship|Age:Gender:Occupation][Origin][Industry|Destination:Origin][Destination|Origin]"))
  
  
  #
  bn.indod.fit <- custom.fit(x=dag.indod, dist=conds)
  
  #Simulate
  bn.indod <- as.data.table(rbn(bn.indod.fit, sum(margins$indagesex), dag.indod))
  
  #Reorder
  bn.indod <- bn.indod[ , c(names(data.ind),"otract","dtract"), with=F]
  
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
  
  # Check results ####
  #results list
  bn.indod.res <- list()
  
  #Microdata results
  bn.indod.res$micro <- merge(as.data.table(prop.table(table(data.ind))),
                              as.data.table(prop.table(table(bn.indod[,!c("otract","dtract")]))),
                              by = colnames(micro[,!"N"]))
  setnames(bn.indod.res$micro, c("N.x","N.y"), c("PUMS","MCMC"))
  #Marginal results
  vn <- names(bn.indod[,!c("otract","dtract")])
  #Tally up frequencies
  bn.indod.res$marg <- rbindlist(lapply(vn, function(n) setNames(bn.indod[,.N, by = c('otract',n) ], c("otract","Variable","MCMC"))))
  bn.indod.res$marg[ , MCMC := MCMC/sum(MCMC)]
  #Merging with marginals to compare
  bn.indod.res$marg <- merge(marg, bn.indod.res$marg, by = c('otract','Variable'), all.x = T)
  bn.indod.res$marg[is.na(MCMC), MCMC := 0]
  
  #OD
  odi <- as.data.table(prop.table(table(bn.indod[,.(otract,dtract,INDUS)])))
  lodes <- setNames(as.data.table(mcmc.odi)[otract!="00000000000" & value>0,],c("otract","dtract","INDUS","LODES"))
  odi <- merge(odi,lodes, by = c("otract","dtract","INDUS"))
  bn.indod.res$ODI <- odi
  rm(odi,vn)
  
  # Plot results ####
  #Microdata plot
  ggplot(data = bn.indod.res$micro) + geom_point(aes(x=MCMC, y=PUMS)) + theme_bw() + coord_fixed()
  
  #Marginal plot
  ggplot(data=bn.indod.res$marg[ , lapply(.SD,sum), .SDcols = c("MCMC","Marginals"), by = Variable]) + 
    geom_point(aes(x=Marginals, y=MCMC, color=Variable)) +
    theme_bw() + coord_fixed() + theme(legend.position="none")  
  
  #Marginal plot w/ tracts
  ggplot(data=bn.indod.res$marg) + 
    geom_point(aes(x=Marginals, y=MCMC, color=Var), alpha = 0.1) +
    theme_bw() + coord_fixed() #+ theme(legend.position="none")

  bn.indod.res$OD <- bn.indod.res$ODI[,lapply(.SD,sum), .SDcols = c("N","LODES"), by = .(otract,dtract)]
  bn.indod.res$DI <- bn.indod.res$ODI[,lapply(.SD,sum), .SDcols = c("N","LODES"), by = .(otract,INDUS)]
  bn.indod.res$OI <- bn.indod.res$ODI[,lapply(.SD,sum), .SDcols = c("N","LODES"), by = .(dtract,INDUS)]
  
  #ODI plot
  ggplot() + 
    geom_point(data=bn.indod.res$OD[sample(bn.indod.res$OD[,.I],10000,replace = F), ], aes(x=N, y=LODES, color="OD"), alpha=0.1) +
    geom_point(data=bn.indod.res$DI, aes(x=N, y=LODES, color="DI"), alpha=0.1) +
    geom_point(data=bn.indod.res$OI, aes(x=N, y=LODES, color="OI"), alpha=0.1) +
    #xlim(c(0,1e-5)) + ylim(c(0,1e-5)) +
    theme_bw() + coord_fixed()
  
  
  #Output
  save(bn.indod, bn.indod.res, file = "./db/d3.3_BNindod.RData")
  #cleanup
  rm(list = ls()[grepl("dag.|bn.",ls())])
  
}

#### Raking the results to fit marginals ####
if(!file.exists("./db/d3.3_BNrake.RData")) {
  if( all(!(c("bn.hh","bn.ind") %in% ls())) ) {load("./db/d3.3_BNindhh.RData")}
  
  #Setup totals and data
  freq.hh <- as.data.table(table(bn.hh[,.(INCOME,RESTY,HHVEH,HHSIZ)]))
  freq.ind <- as.data.table(table(bn.ind[,.(AGE,GEND,RELATE,OCC,INDUS)]))
  
  #Setup marginals
  marg.hh <- list("INCOME" = as.matrix(margins$hhinc),
                  "RESTY"  = margins$hhdwell,
                  "HHVEH"  = apply(margins$hhvehsize, c(1,3), sum),
                  "HHSIZ" = apply(margins$hhvehsize, c(1,2), sum))
  
  nm <- names(marg.hh)
  marg.hh <- lapply(names(marg.hh), function(n) as.data.table(setNames(melt(marg.hh[[n]]),c("otract",n,"Freq"))) )
  names(marg.hh) <- nm
  
  marg.ind <- list("AGE" = apply(margins$indagesex, c(1,3), sum),
                   "GEND" = apply(margins$indagesex, c(1,2), sum),
                   "RELATE" = as.matrix(margins$indrelate),
                   "OCC"  = margins$indocc,
                   "INDUS"  = margins$indindus)
  nm <- names(marg.ind)
  marg.ind <- lapply(names(marg.ind), function(n)  as.data.table(setNames(melt(marg.ind[[n]]),c("otract",n,"Freq"))) )
  names(marg.ind) <- nm
  
  # Households ####
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
  bn.hh.raked <- pblapply(tracts, function(tr) fun_raketract.singledim(tr, freq.hh, marg.hh, bnd=c(0.2,4), trm=c(0.01,100)), cl=cl)
  names(bn.hh.raked) <- tracts
  stopCluster(cl)
  
  #Checking results
  hh.raked.res <- rbindlist(bn.hh.raked)
  hh.raked.res <- lapply(names(marg.hh), function(n) {
    out <- hh.raked.res[,lapply(.SD,sum), .SDcols=c("weights.raked","weights.sim"), by = c('otract',n)]
    return(setNames(out, c("otract","Variable","Raked","BN")))
  })
  hh.raked.res <- rbindlist(hh.raked.res)
  #Marginals
  hh.raked.marg <- as.data.table(rbindlist(lapply(marg.hh, function(x) setNames(x, c("otract","Variable","Marginals")))))
  hh.raked.marg[ , otract := as.character(otract)]
  #Merging together and melting
  bn.hh.raked.res <- merge(hh.raked.res, hh.raked.marg, by = c('otract','Variable'), all = T)
  bn.hh.raked.res[ , (c("Raked","BN","Marginals")) := lapply(.SD, function(x) ifelse(is.na(x),0,x)), .SDcols = c("Raked","BN","Marginals")]
  
  ggplot(data = bn.hh.raked.res) + geom_point(aes(x=Marginals, y=Raked, color=Variable), alpha=0.5) + theme_bw() + coord_fixed()
  ggplot(data = bn.hh.raked.res[,lapply(.SD,sum), by=Variable, .SDcols = c("Raked","BN","Marginals")]) + 
    geom_point(aes(x=BN, y=Raked, color=Variable), alpha=0.5) + theme_bw() + coord_fixed()  
  
  # Individuals ####
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
  
  bn.ind.raked <- pblapply(tracts, function(tr) fun_raketract.singledim(tr, freq.ind, marg.ind, bnd=c(0.1,8), trm=c(0.01,1000)), cl=cl)
  names(bn.ind.raked) <- tracts
  stopCluster(cl)
  
  # Checking results
  ind.raked.res <- rbindlist(bn.ind.raked)
  ind.raked.res <- lapply(names(marg.ind), function(n) {
    out <- ind.raked.res[,lapply(.SD,sum), .SDcols=c("weights.raked","weights.sim"), by = c('otract',n)]
    return(setNames(out, c("otract","Variable","Raked","BN")))
  })
  ind.raked.res <- rbindlist(ind.raked.res)
  #Marginals
  ind.raked.res.marg <- as.data.table(rbindlist(lapply(marg.ind, function(x) setNames(x, c("otract","Variable","Marginals")))))
  ind.raked.res.marg[ , otract := as.character(otract)]
  #Merging together and melting
  bn.ind.raked.res <- merge(ind.raked.res, ind.raked.res.marg, by = c('otract','Variable'), all = T)
  bn.ind.raked.res[ , (c("Raked","BN","Marginals")) := lapply(.SD, function(x) ifelse(is.na(x),0,x)), .SDcols = c("Raked","BN","Marginals")]
  
  ggplot(data = bn.ind.raked.res) + geom_point(aes(x=Marginals, y=Raked, color=Variable), alpha=0.5) + theme_bw() + coord_fixed()
  ggplot(data = bn.ind.raked.res[,lapply(.SD,sum), by=Variable, .SDcols = c("Raked","BN","Marginals")]) + 
    geom_point(aes(x=BN, y=Raked, color=Variable), alpha=0.5) + theme_bw() + coord_fixed()  
  

  #Saving the output
  save(bn.hh.raked, bn.ind.raked, bn.hh.raked.res, bn.ind.raked.res, file = "./db/d3.3_BNrake.RData")
  #
  rm(list = setdiff(ls(),c("bn.hh.raked", "bn.ind.raked", "bn.hh.raked.res", "bn.ind.raked.res",
                           "bn.hh","bn.ind","bn.ind.fit","bn.hh.fit","bn.ind.res","bn.hh.res",
                           "margins","microdata","tables","tracts",
                           "fun_raketract.singledim","fun_plotnetwork")))
}

#### Raking the OD results to fit marginals ####
if(!file.exists("./db/d3.3_BNindodrake.RData")) {
  if(!("bn.indod" %in% ls())) {load("./db/d3.3_BNindod.RData")}
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
                  list("fun_raketract.singledim","marg.indod","bn.indod"),
                  envir=environment())
  }
  
  #Raking
  bn.indod.raked <- pblapply(tracts, function(tr) {
    freq.dat <- as.data.table(table(bn.indod[otract==tr,.(AGE,GEND,RELATE,OCC,INDUS,dtract)]))
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
  names(bn.indod.raked) <- tracts
  stopCluster(cl)
  
  #Check results
  indod.raked.res <- rbindlist(bn.indod.raked)
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
  bn.indod.raked.res <- merge(indod.raked.res, indod.raked.res.marg, by = c('otract','Variable'), all = T)
  bn.indod.raked.res[ , (c("Raked","MCMC","Marginals")) := lapply(.SD, function(x) ifelse(is.na(x),0,x)), .SDcols = c("Raked","MCMC","Marginals")]
  
  bn.indod.raked.res <- bn.indod.raked.res[Raked>0 & MCMC>0 & Marginals>0,]
  
  ggplot(data = bn.indod.raked.res) + 
    geom_point(aes(x=Marginals, y=Raked), alpha=0.1) + theme_bw() + coord_fixed()
  
  ggplot(data = bn.indod.raked.res[,lapply(.SD,sum), by=Variable, .SDcols = c("Raked","MCMC","Marginals")]) + 
    geom_point(aes(x=MCMC, y=Marginals), alpha=0.5) + theme_bw() + coord_fixed()
  
  #Saving the output
  save(bn.indod.raked, bn.indod.raked.res, file = "./db/d3.3_BNindodrake.RData")
  rm(list = ls()[grepl("bn.",ls())])
}




#### Cleanup ####
rm(list=ls())
