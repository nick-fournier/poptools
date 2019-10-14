#### Setting up functions ####
seeder <- function(x, tr){ #this function creates a matrix for each tract, allows for parallel loops
  #puma the tract is in
  puma <- tract2puma[TRACTID==tr,PUMA5CE]
  pumas <- tract2puma[, PUMA5CE]
  #puma frequencies
  tab.puma <- table(x[PUMA==puma, !"PUMA", with=F])
  #Regional frequencies
  tab.gba <- table(x[PUMA %in% pumas, !"PUMA", with=F])
  #Regional frequencies
  tab.region <- table(x[, !"PUMA", with=F])
  
  #Store in new table starting with pumas
  tab <- tab.puma
  
  #Find missing, fill from GBA
  #First proportionally scale the matrix up by the amount being added
  fill <- tab.gba[tab==0] * sum(tab) / sum(tab.gba)
  tab <- tab * (1 + sum(fill) / sum(tab))
  #Second add the replacement fill probabilities
  tab[tab==0] <- tab.gba[tab==0] * sum(tab) / sum(tab.gba)
  
  #Find missing, fill from region
  #First proportionally scale the matrix up by the amount being added
  fill <- tab.region[tab==0] * sum(tab) / sum(tab.region)
  tab <- tab * (1 + sum(fill) / sum(tab))
  #Second add the replacement fill probabilities
  tab[tab==0] <- tab.region[tab==0] * sum(tab) / sum(tab.region)
  
  #Proportionally scale the filled seed back down to puma total
  tab <- tab * sum(tab.puma) / sum(tab)
  
  #fill in remaining 0's with very very very small values, just in case.
  #tablocal[which(tablocal==0)] <- 1e-5
  #convert to array with dimnames
  array(tab, dim = dim(tab), dimnames = dimnames(tab))
}
format_pums <- function(indvars, hhvars) {
  pums <- microdata[["pums"]][,c("SERIALNO", indvars, hhvars), with=F]
  #add puma
  pumas <- tables[["pumsind"]][,c("SERIALNO","PUMA")]
  pumas <- pumas[!duplicated(SERIALNO),]
  pums <- merge(pums, pumas, by="SERIALNO")
  # #renumber hhid
  #pums <- merge(pums, data.table(HHID = 1:length(unique(pums$SERIALNO)), SERIALNO = unique(pums$SERIALNO)), by = "SERIALNO")
  #separating
  pums_ind <- pums[, c("PUMA",indvars), with = F]
  pums_hh <- pums[!duplicated(SERIALNO), c("PUMA",hhvars), with = F]
  return(list(joint = pums, ind = pums_ind, hh = pums_hh))
}
ipf_worktripdist <- function(crit=1e-8, reset = F){
  # Total O-D IPF
  if(!file.exists("./db/d3.1_IPFlodes.RData") | reset == T) {
    print("Running IPF for work trips")
    #Data setup LODES ####
    #lodes <- setuplodes()
    load("./db/d2.4_lodes.RData")
    
    LODES=lodes$LODES
    mtxOD=lodes$mtxOD
    mtxDI=lodes$mtxDI
    mtxOI=lodes$mtxOI
    
    #Adding in non-work from marginal census table
    margOI <- apply(margins$indindusocc, c(1,3), sum)
    dimnames(margOI) <- list(otract = dimnames(margOI)[[1]], INDUS = dimnames(margOI)[[2]])
    
    #Sorting row order for consistency!!!
    margOI <- margOI[tracts,]
    mtxOI <- mtxOI[tracts,]       
    mtxDI <- mtxDI[tracts,]
    mtxOD <- mtxOD[tracts,tracts]
    
    ODnames <- dimnames(mtxOD)
    OInames <- dimnames(mtxOI)
    DInames <- dimnames(mtxDI)
    
    #Checking totals
    data.table("Indus" = colnames(mtxDI),
               "DI" = colSums(mtxDI),
               "OI" = colSums(mtxOI),
               "marginOI" = colSums(margOI)[colnames(mtxOI)],
               "Ratio" = colSums(mtxOI)/colSums(mtxDI),
               "MargRatio" = colSums(mtxOI)/colSums(margOI)[colnames(mtxOI)])
    c("DI" = sum(mtxDI), "OD" = sum(mtxOD), "OI" = sum(mtxOI), "margOD" = sum(colSums(margOI)[-1]))
    
    #Plot checks ####
    plot(x=rowSums(mtxOI), y=rowSums(mtxOD))
    plot(x=rowSums(mtxDI), y=colSums(mtxOD))
    plot(x=colSums(mtxOI), y=colSums(mtxDI))
    
    #Filling 0's with 1's in marginals so that they can be fitted if necessary
    mtxOI[mtxOI==0] <- 1
    mtxDI[mtxDI==0] <- 1
    #Adjusting OD to match marg totals ####
    mtxOD <- mtxOD*rowSums(margOI[,-which(colnames(margOI)=='NOIND')])/rowSums(mtxOD)
    mtxOI[!is.finite(mtxOI)] <- 0
    #Adjusting to match OD totals
    mtxOI <- mtxOI*rowSums(mtxOD)/rowSums(mtxOI)
    mtxOI[!is.finite(mtxOI)] <- 0
    #Adjusting to match DI totals
    mtxDI <- mtxDI*colSums(mtxOD)/rowSums(mtxDI)
    mtxDI[!is.finite(mtxDI)] <- 0
    
    #Adding in NOIND ####
    mtxOI <- cbind(mtxOI, NOIND=margOI[,'NOIND'])
    mtxDI <- cbind(mtxDI, NOIND=margOI[,'NOIND'])
    for(tr in tracts) mtxOD[tr,tr] <- margOI[tr, 'NOIND'] + mtxOD[tr,tr]
    
    #putting names back
    names(dimnames(mtxOI)) <- names(OInames)
    names(dimnames(mtxDI)) <- names(DInames)
    names(dimnames(mtxOD)) <- names(ODnames)
    
    #recheck
    data.table("Indus" = colnames(mtxDI),
               "DI" = colSums(mtxDI),
               "OI" = colSums(mtxOI),
               "marginOI" = colSums(margOI)[colnames(mtxOI)],
               "Ratio" = colSums(mtxOI)/colSums(mtxDI),
               "MargRatio" = colSums(mtxOI)/colSums(margOI)[colnames(mtxOI)])
    c("DI" = sum(mtxDI), "OD" = sum(mtxOD), "OI" = sum(mtxOI), "margOD" = sum(colSums(margOI)))#<- this needs to match
    all_equal(sum(mtxDI),sum(mtxOD),sum(mtxOI), sum(colSums(margOI)))
    
    #Plot checks ####
    plot(x=rowSums(mtxOI), y=rowSums(mtxOD))
    plot(x=rowSums(mtxDI), y=colSums(mtxOD))
    plot(x=colSums(mtxOI), y=colSums(mtxDI))
    
    #
    mtxOI <- mtxOI[tracts, dimnames(margOI)[[2]]]
    mtxDI <- mtxDI[tracts, dimnames(margOI)[[2]]]
    #Setting up uniform array for seed ####
    seed <- array(data = 1,
                  dim = c(dim(mtxOD), ncol(mtxOI)),
                  dimnames = c(dimnames(mtxOD), list(INDUS = dimnames(mtxDI)$INDUS)))
    #Applying weights to seed
    seed <- aperm(aaply(seed, 3, '*', mtxOD), c(2,3,1))
    #Fill zeros
    seed[seed==0] <- 1e-10
    
    #Filling zeros along all but diagonal for no industry, must stay in home tract if no work!
    seed[,,'NOIND'] <- seed[,,'NOIND'] * diag(1, length(tracts))
    #Scale totals back up for unity
    seed[,,'NOIND'] <- seed[,,'NOIND'] * sum(mtxOD) / sum(seed[,,'NOIND'])
    
    #Scale by industry
    scale <- (colSums(mtxDI)/sum(mtxOD) + colSums(mtxDI)/sum(mtxOD))/2
    #Scale along industry dimension
    for(i in names(scale)) seed[,,i] <- seed[,,i]*scale[i]
    rm(i)
    
    #Setting up target data (marginals)
    target <- list(mtxOI, mtxDI, mtxOD)
    #Setting up dimensional references between marginals and seed array
    descript <- list(
      c(which(names(dimnames(seed))=='otract'),which(names(dimnames(seed))=='INDUS')),
      c(which(names(dimnames(seed))=='dtract'),which(names(dimnames(seed))=='INDUS')),
      c(which(names(dimnames(seed))=='otract'),which(names(dimnames(seed))=='dtract')))
    
    #Running IPF ####
    ipf.lodes <- Ipfp(seed, descript, target, print = T, tol = crit)
    #check errors
    print(ipf.lodes$error.margins)
    print("Saving total flows data")
    save(ipf.lodes, file = "./db/d3.1_IPFlodes.RData")
    gc()
  } else {
    print("Loading existing total flows data")
    load("./db/d3.1_IPFlodes.RData")
    print(ipf.lodes$error.margins)
  }
  return(ipf.lodes$x.hat)
}
ipf_indwork <- function(x, tr, verb=F, crit=1e-5){
  #seed
  seed <-  seeder(x, tr)
  names <- c(names(dimnames(seed)),'dtract')
  #Proportion of trips to scale by
  workprop <- ipf.lodes[tr,,] / sum(ipf.lodes[tr,,])
  #Duplicating seed matrix for every possible destination into a list
  seedlist <- lapply(tracts, function(i) seed)
  names(seedlist) <- tracts
  #binding seed list into array along the new dimension
  seed <- abind(seedlist, along = length(dim(seed))+1, new.names = tracts, use.dnns = T)
  names(dimnames(seed))[which(names(dimnames(seed)) == "")] <- "dtract"
  
  #Scaling seed to dest-industry matrix 
  for(j in tracts) { for(i in dimnames(seed)$INDUS) { seed[,,,,i,j] <- seed[,,,,i,j] * workprop[j,i] }}
  rm(i,j)

  #pointers
  descript<-list(which(names=="RELATE"), #relate, no hfam
                 #which(names=="SCHOL"), #school enrollment
                 c(which(names=="OCC"),which(names=="INDUS")), #industry occ
                 #which(names=="INDUS"),
                 c(which(names=="GEND"),which(names=="AGE")), #agesex
                 c(which(names=="dtract"),which(names=="INDUS"))) 
  #marginals
  target <- list(unlist(margins[["indrelate"]][tr,]), #relate no hfam
                 #unlist(margins[["indschool"]][tr,,]), #school
                 margins[["indindusocc"]][tr,,], #industry occ
                 #unlist(margins[["indindus"]][tr,]),
                 #unlist(colSums(ipf.lodes[tr,,])),
                 margins[["indagesex"]][tr,,], #agesex
                 ipf.lodes[tr,,]) #destination and industry
  
  ipfres <- Ipfp(seed, descript, target, print = verb, iter = 1000, tol = crit)
  print(ipfres$error.margins)
  #weights <- int_trs(ipfres$x.hat)
  weights <- ipfres$x.hat
  return(weights)
}
ipf_ind <- function(x, tr, verb=F, crit=1e-5){
  #seed
  seed <-  seeder(x, tr)
  names <- names(dimnames(seed))
  #pointers
  descript<-list(which(names=="RELATE"), #relate, no hfam
                 c(which(names=="OCC"),which(names=="INDUS")), #industry occ
                 c(which(names=="GEND"),which(names=="AGE"))) #agesex
  #marginals
  target <- list(unlist(margins[["indrelate"]][tr,]), #relate no hfam
                 margins[["indindusocc"]][tr,,], #industry occ
                 margins[["indagesex"]][tr,,])
  
  ipfres <- Ipfp(seed, descript, target, print = F, iter = 1000, tol=1e-10)
  weights <- ipfres$x.hat
  #print(ipfres$error.margins)
  return(weights)
  
}
ipf_totind <- function(x){
  tab <- table(x[,!"PUMA",with=F])
  seed <- array(tab, dim = dim(tab), dimnames = dimnames(tab))
  names <- names(dimnames(seed))
  
  descript<-list(which(names=="RELATE"), #relate, no hfam
                 c(which(names=="OCC"),which(names=="INDUS")), #industry occ
                 c(which(names=="GEND"),which(names=="AGE"))) #agesex
  
  target <- list(unlist(colSums(margins[["indrelate"]])), #relate no hfam
                 colSums(margins[["indindusocc"]]), #industry occ
                 colSums(margins[["indagesex"]])) #agesex
  
  ipfres <- Ipfp(seed, descript, target, print = F, iter = 1000, tol = 1e-10)
  print(ipfres$error.margins)
  weights <- ipfres$x.hat
  return(weights)
}
ipf_hh <- function(x, tr){
  seed <-  seeder(x, tr)
  names <- names(dimnames(seed))
  descript<-list(which(names=="INCOME"), #inc
                 which(names=="RESTY"), #hh type
                 c(which(names=="HHSIZ"),which(names=="HHVEH"))) #vehicles size
  
  target <- list(unlist(margins[["hhinc"]][tr,]),
                 unlist(margins[["hhdwell"]][tr,]),
                 margins[["hhvehsize"]][tr,,])
  
  ipfres <- Ipfp(seed, descript, target, print = F, iter = 1000, tol=1e-10)
  weights <- ipfres$x.hat
  #print(ipfres$error.margins)
  return(weights)
}
ipf_tothh <- function(x){
  tab <- table(x[,!"PUMA",with=F])
  seed <- array(tab, dim = dim(tab), dimnames = dimnames(tab))
  names <- names(dimnames(seed))
  #seed[seed==0] <- 1e-5
  descript<-list(which(names=="INCOME"), #inc
                 which(names=="RESTY"),
                 c(which(names=="HHSIZ"),which(names=="HHVEH"))) #hh vehicles & hhsize size
  
  target <- list(unlist(colSums(margins[["hhinc"]])),
                 unlist(colSums(margins[["hhdwell"]])),
                 colSums(margins[["hhvehsize"]]))
  
  ipfres <- Ipfp(seed, descript, target, print = F, iter = 1000, tol=1e-10)
  print(ipfres$error.margins)
  weights <- ipfres$x.hat
  return(weights)
}

load("./db/d2.3_microdata.RData")
load("./db/d2.5_margins.RData")

#### Setting data ####
#formatting PUMS for IPF
indvars=c("GEND","AGE","RELATE","OCC","INDUS")#,"HOURS","TTIME", "SCHOL"
hhvars=c("INCOME","RESTY","HHSIZ","HHVEH") #RACE
pums <- format_pums(indvars, hhvars)

#Setup tracts
tract2puma <- data.table(tables[['tract2puma']])
tract2puma[, TRACTID := paste(STATEFP, COUNTYFP, TRACTCE, sep="")]
tract2puma <- tract2puma[STATEFP=="25",]

#### Work trip distribution IPF ####
ipf.lodes <- ipf_worktripdist(crit=1e-2, reset = F)

#### Person & Household IPF ####
ipf.indtot <- ipf_totind(pums[['ind']])
ipf.hhtot <- ipf_tothh(pums[['hh']])

#Zone-by-zone IPF without work destination
if(!file.exists('./db/d3.1_IPFindhh.RData')) {
  print("Running IPF for households")
  ipf.hh <- pblapply(tracts, function(tr) ipf_hh(pums[['hh']], tr))
  print("Running IPF for individuals without destination")
  ipf.ind <- pblapply(tracts, function(tr) ipf_ind(pums[['ind']], tr))
  #Add names
  names(ipf.hh) <- tracts
  names(ipf.ind) <- tracts
  #saving
  save(ipf.hh, ipf.hhtots, ipf.ind, ipf.indtots, file='./db/d3.1_IPFindhh.RData') #indpoplist #
} else {
  load('./db/d3.1_IPFindhh.RData')
}

#Zone-by-zone IPF WITH work destination
#Loading stored results
if(!dir.exists("./db/d3.1_dbIPFindwork")){ dbCreate("./db/d3.1_dbIPFindwork", "RDS") }
dbindworkpoplist <- dbInit("./db/d3.1_dbIPFindwork", "RDS")

#Set up parallel cores
if(.Platform$OS.type == "unix") {
  cl <- makeCluster(detectCores()-1, "FORK") 
} else {
  cl <- makeCluster(detectCores(logical = F)-1, "PSOCK")
  clusterEvalQ(cl,c(library(filehash),library(data.table),library(abind),library(mipfp)))
  clusterExport(cl,
                list("dbindworkpoplist","ipf_indwork","seeder",
                     "tract2puma","pums","ipf.lodes","tracts","margins"),
                envir=environment())
}

#running IPF and storing straight to hard drive
print(paste("Running IPF for individuals", "Tracts remaining:", sum(!(tracts %in% dbList(dbindworkpoplist)))))
#dummy variable is an arbitary export to get the lapply function to operate without consuming memory
if(sum(!(tracts %in% dbList(dbindworkpoplist)))>0){
  dummy <- pblapply(tracts[!(tracts %in% dbList(dbindworkpoplist))], function(tr) {
    if(!(tr %in% dbList(dbindworkpoplist))) dbInsert(dbindworkpoplist, tr, ipf_indwork(pums[['ind']], tr, verb=T, crit=1e-8))
    return(tr)
  }, cl = cl)
}
stopCluster(cl)

#### Cleanup ####
rm(list=ls())


