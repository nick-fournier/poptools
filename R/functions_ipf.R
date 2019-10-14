

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
setuplodes <- function() {
  #Setting up lodes data
  
  #loading data ####
  print("Loading LODES data")
  #list of lodes files
  lslodes <- list(
    'OD' = fread("./rdata/LODES/ma_od_main_JT00_2015.csv", sep = ",", integer64 = "character"),
    'dtotal' = fread("./rdata/LODES/ma_wac_S000_JT00_2015.csv", sep = ",", integer64 = "character"),
    'ototal' = fread("./rdata/LODES/ma_rac_S000_JT00_2015.csv", sep = ",", integer64 = "character"))
  
  print("Loading conversion tables")
  #List of support files
  lsconv <- list(
    header = fread("./rdata/LODES/wacheaders.csv", sep = ",", integer64 = "character"),
    BLK2TAZ = fread("./rdata/spatial/BLOCKS2TAZ.csv", stringsAsFactors = F, sep = ",", integer64 = "character"),
    naicskey = fread("./rdata/LODES/naicskey.csv", sep = ","))
  
  #Setting up lodes ####
  print("Setting up LODES data")
  #Loading NAICS ID key
  naicskey <- lsconv[["naicskey"]]
  naicskey <- naicskey[!duplicated(naicskey[,-1]),-1] #clean up key for use
  rownames(naicskey) <- NULL
  #loading workplace tables
  headers <- lsconv[["header"]]
  headers <- headers[which(grepl("CNS|000",headers$Variable)),c("Variable")]
  headers <- headers[[1]]
  
  #trimming columns and aggregating worker counts per tract ####
  LODES <- lapply(lslodes, function(x) { #subsets columns, aggregates blocks to tracts, consolidates to INDUS variables
    #OD pairs
    if(all(c("h_geocode","w_geocode") %in% names(x))) {
      print("Setting up O-D table")
      x <- x[,c("h_geocode","w_geocode","S000","SA01","SA02","SA03")]
      #convert block to tracts
      print("Aggregating OD matrix into tracts")
      #trimming geocode to tract level
      x[, otract := substr(h_geocode,0,11)]
      x[, dtract := substr(w_geocode,0,11)]
      x <- x[,!c("h_geocode","w_geocode"), with=F]
      #label destinations outside work tracts as "00000000000"
      x[!(dtract %in% tracts), dtract := "00000000000"]
      x[!(otract %in% tracts), otract := "00000000000"]
      #convert to number
      x <- x[, lapply(.SD, as.numeric), by = c("otract","dtract")]
      #Aggregating
      x <- x[,lapply(.SD,sum),by=c("otract","dtract")]
      #renumbering
      rownames(x) <- NULL
      colnames(x)[which(colnames(x) %in% c("S000","SA01","SA02","SA03"))] <- c('total','under29','30to54','55over')
    }
    #work D-zones
    if("w_geocode" %in% names(x) & !("h_geocode" %in% names(x))) {
      print("Aggregating destinations into tracts")
      x <- x[,c("w_geocode",headers), with=F]
      #convert block to tract
      x[, dtract := substr(w_geocode,0,11)]
      x <- x[,-which(colnames(x) %in% c("w_geocode","h_geocode")), with = F]
      #label destinations outside work tracts as "00000000000"
      x[!(dtract %in% tracts), dtract := "00000000000"]
      #convert to number
      x <- x[, lapply(.SD, as.numeric), by = "dtract"]
      #Aggregate over tracts's
      x <- x[,lapply(.SD,sum),by="dtract"]
      #renumbering
      rownames(x) <- NULL
      x2 <- as.data.table(sapply(unique(naicskey$INDUS), function(h)
        if(sum(h == naicskey$INDUS)>1) 
        { rowSums(x[, naicskey[INDUS == h, VAR], with = F]) }
        else
        { x[[naicskey[INDUS == h, VAR]]] }
      ))
      #putting back
      x <- data.table(dtract = x$dtract, x2, stringsAsFactors = F)
      rownames(x) <- NULL
    }
    #home O-zones
    if("h_geocode" %in% names(x) & !("w_geocode" %in% names(x))) {
      print("Aggregating origins into tracts")
      x <- x[,c("h_geocode",headers), with=F]
      #convert block to tract
      x[, otract := substr(h_geocode,0,11)]
      x <- x[,-which(colnames(x) %in% c("w_geocode","h_geocode")), with = F]
      #label origins outside home tracts as "00000000000"
      x[!(otract %in% tracts), otract := "00000000000"]
      #convert to number
      x <- x[, lapply(.SD, as.numeric), by = "otract"]
      #Aggregate over tract's
      x <- x[,lapply(.SD,sum), by="otract"]
      #renumbering
      rownames(x) <- NULL
      x2 <- as.data.table(sapply(unique(naicskey$INDUS), function(h)
        if(sum(h == naicskey$INDUS)>1) 
        { rowSums(x[, naicskey[naicskey$INDUS == h, VAR], with = F]) }
        else
        { x[[naicskey[naicskey$INDUS == h, VAR]]] }
      ))
      #putting back
      x <- data.table(otract = x$otract, x2, stringsAsFactors = F)
      rownames(x) <- NULL
    }
    return(x)
  })
  
  #setting up totals matrices ####
  print("Setting up totals tables")
  mtxDI <- as.matrix(LODES[['dtotal']][,-1], dimnames = list(LODES[['dtotal']]$dtract, colnames(LODES[['dtotal']])[-1]))
  rownames(mtxDI) <- LODES[['dtotal']]$dtract
  mtxOI <- as.matrix(LODES[['ototal']][,-1], dimnames = list(LODES[['ototal']]$otract, colnames(LODES[['ototal']])[-1]))
  rownames(mtxOI) <- LODES[['ototal']]$otract
  mtxOD <- acast(LODES[['OD']], otract~dtract, value.var = "total", fill = 0)
  
  dimnames(mtxDI) <- list(dtract = dimnames(mtxDI)[[1]], INDUS = dimnames(mtxDI)[[2]])
  dimnames(mtxOI) <- list(otract = dimnames(mtxOI)[[1]], INDUS = dimnames(mtxOI)[[2]])
  dimnames(mtxOD) <- list(otract = dimnames(mtxOD)[[1]], dtract = dimnames(mtxOD)[[2]])
  return(list(LODES=LODES, mtxOD=mtxOD, mtxDI=mtxDI, mtxOI=mtxOI))
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
    lodes <- setuplodes()
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


