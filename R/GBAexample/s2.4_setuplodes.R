

load("./data/microdata.RData")
load("./data/marginal_ind.RData")
load("./data/marginal_hh.RData")

if(!exists("tracts")) tracts <- marginals_hh[[1]][['tracts']]

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

# #Save for example data
# lodes_od = as.data.table(LODES$OD)
# lodes_di = as.data.table(LODES$dtotal)
# lodes_oi = as.data.table(LODES$ototal)
# fwrite(lodes_od, file = "./data/lodes_od.csv")
# fwrite(lodes_di, file = "./data/lodes_di.csv")
# fwrite(lodes_oi, file = "./data/lodes_oi.csv")
# save(lodes_od, file = "./data/lodes_od.RData")
# save(lodes_di, file = "./data/lodes_di.RData")
# save(lodes_oi, file = "./data/lodes_oi.RData")


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

#Match lodes to census marginals?
margOI <- margins$indindus
names(dimnames(margOI)) <- c('otract','INDUS')# <- list(otract = dimnames(margOI)[[1]], INDUS = dimnames(margOI)[[2]])

#Sorting row order for consistency!!!
mtxOI <- mtxOI[c("00000000000",tracts),]
mtxDI <- mtxDI[c("00000000000",tracts),]
mtxOD <- mtxOD[c("00000000000",tracts),c("00000000000",tracts)]
# #We only care out outside region as destination, not for origin.
# mtxOD["00000000000",] <- 0
# mtxOI["00000000000",] <- 0

#store name for later
ODnames <- dimnames(mtxOD)
OInames <- dimnames(mtxOI)
DInames <- dimnames(mtxDI)

#Plot checks
#OD vs OI
qplot(data=data.table(OI=rowSums(mtxOI), OD=rowSums(mtxOD))[-1,], x=OD,y=OI) + theme_classic() + coord_fixed()
#OD vs DI
qplot(data=data.table(DI=rowSums(mtxDI), OD=colSums(mtxOD))[-1,], x=OD,y=DI) + theme_classic() + coord_fixed()
#OI vs DI (industry)
qplot(data=data.table(DI=colSums(mtxDI), OI=colSums(mtxOI)), x=OI,y=DI) + theme_classic() + coord_fixed()
#OI vs margOI
qplot(data=data.table(OI=rowSums(mtxOI[tracts,]), margOI=rowSums(margOI[tracts,-1])), x=OI,y=margOI) + theme_classic() + coord_fixed()

lodes <- list("mtxOI" = mtxOI, "mtxOD" = mtxOD, "mtxDI" = mtxDI)

save(lodes, file = "./db/d2.4_lodes.RData")
rm(list=ls())
