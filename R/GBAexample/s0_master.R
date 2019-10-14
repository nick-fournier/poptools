#### Create data base folder if not already
# if(!dir.exists("./db")){ dir.create("./db") }

#### Installing necessary packages ####
print("Loading packages")
system.time(source("./s1.0_packages.R", echo=F))

#### Loading raw data in ####
print("Loading raw data")
if(!file.exists("./db/d2.1_tables.RData")) {
  system.time(source("./s2.1_loadtables.R", echo=F))
  save(tables, file = "./db/d2.1_tables.RData")
}

#### Formatting the microdata ####
print("Loading micro data")
if(!file.exists("./db/d2.3_microdata.RData")) {
  # Checking for census tract consistency
  #counties to ignore 
  # MA,25,001,Barnstable County,H1
  # MA,25,003,Berkshire County,H4 <<<
  # MA,25,005,Bristol County,H1
  # MA,25,007,Dukes County,H1 <<<
  # MA,25,009,Essex County,H4
  # MA,25,011,Franklin County,H4 <<<
  # MA,25,013,Hampden County,H4 <<<
  # MA,25,015,Hampshire County,H4 <<<
  # MA,25,017,Middlesex County,H4
  # MA,25,019,Nantucket County,H4 <<<
  # MA,25,021,Norfolk County,H1
  # MA,25,023,Plymouth County,H1
  # MA,25,025,Suffolk County,H4
  # MA,25,027,Worcester County,H4
  #Tracts used by CTPS
  tracts <- fread("./rdata/spatial/BLOCKS2TAZ.csv", colClasses = "character", stringsAsFactors = F)
  tracts <- tracts[,substr(GEOID10,0,11)]
  tracts <- unique(tracts)
  
  # making list of empty tracts to avoid divide by 0 error
  tmp <- rbind(with(tables[['indagesex']], data.table(GEO.id2, HD01_VD01)[-1,]),
               with(tables[['hhinc']], data.table(GEO.id2, HD01_VD01)[-1,]))
  tmp[ , HD01_VD01 := as.numeric(HD01_VD01)]
  tmp <- tmp[ , min(HD01_VD01), by = GEO.id2]
  emptytracts <- tmp[V1==0, GEO.id2]
  emptytracts <- c(emptytracts, c("25025980101", "25025980700")) #Adding mistmatched tracts found manually in OD tables.
  # emptytracts <- c(emptytracts, c("25001014900", "25017354800", "25021417300", "25025020303",
  #                                 "25025981502", "25025980300", "25025981100", "25027732902"))
  emptytracts<-unique(emptytracts)
  tracts<-tracts[!(tracts %in% emptytracts)]
  rm(tmp)
  
  #Setup PUMS data
  print("Setting up PUMS microdata")
  system.time(source("./s2.3_pumsmicroformat.R", echo=F))
  #Setup MTS data
  print("Setting up MTS microdata")
  system.time(source("./s2.4_mtsmicrosample.R", echo=F))
  #Saving
  save(microdata, tracts, file = "./db/d2.3_microdata.RData")
}

#### Set up LODES ####
print("Loading LODES data")
if(!file.exists("./db/d2.4_setuplodes.RData")) {
  source("./s2.4_setuplodes.R", echo=F)
}

#### Set up marginals ####
print("Loading marginal data")
if(!file.exists("./db/d2.5_margins.RData")) {
  margins <- list()
  
  #Setting up margin tables for individuals
  print("Setting up individual margins")
  source("./s2.5_indmargins.R", echo=F)
  
  #Setting up margin tables for households
  print("Setting up household margins")
  source("./s2.6_hhmargins.R", echo=F)
  
  #Storing margins
  print("Saving margins")
  save(margins, file = "./db/d2.5_margins.RData")
}

#### IPF ####
if( !file.exists("./db/d3.1_IPFindhh.RData") & !all(tracts %in% list.files("./db/d3.1_dbIPFindwork"))) {
  system.time(source("./s3.1_IPF_worktrip.R", echo=F))
}

#### MCMC ####
system.time(source("./s3.2_MCMC_worktrip.R", echo=F))

#### BN ####
system.time(source("./s3.3_BN_worktrip.R", echo=F))


#
system.time(source("./s4.1_NLAD_IPF.R", echo=F))




# #Write to csv
# csvwriteout <- function() {
#   library(data.table)
#   print("Writing data to CSV:")
#   print("Loading vehicle data")
#   load("./db/d7.0_vehpop.RData")
#   print("Loading population data")
#   load("./db/d9_finalpop.RData")
#   #cleaning up individuals
#   print("Checking for errors")
#   indpop$postcode_id <- indpop$taz_id
#   indpop$taz_id <- NULL
#   err <- F
#   #check if any null
#   for(i in 1:length(indpop)) { if(any(is.na(indpop[[i]]))){
#     print(paste("NA's in ", colnames(indpop)[i]))
#     err <- T
#     } }
#   for(i in 1:length(hhpop)) { if(any(is.na(hhpop[[i]]))){
#     print(paste("NA's in ", colnames(indpop)[i]))
#     err <- T
#   } }
#   if(err==F) {
#     #writing
#     print("No errors found, Writing data")
#     fwrite(indpop, paste("./output/CSV/individual_", format(Sys.Date(), format="%m-%d-%Y"), ".csv", sep = ""))
#     fwrite(hhpop, paste("./output/CSV/household_", format(Sys.Date(), format="%m-%d-%Y"), ".csv", sep = ""))
#     fwrite(vehpop, paste("./output/CSV/vehicles_", format(Sys.Date(), format="%m-%d-%Y"), ".csv", sep = ""))
#     print("Done saving data")
#   } else {
#     print("Errors found, breaking from program.")
#   }
# }
# 
# system.time(csvwriteout())

#Refresh
source("./cleanup.R", echo = F)
