


#### Create a local data base folder on hard drive for caching
if(!dir.exists("./db")){ dir.create("./db") }

#### Installing necessary packages ####
print("Loading packages")
system.time(source("./s1.0_packages.R", echo=F))

#### Loading raw data in ####
print("Loading raw data")
system.time(source("./s2.1_loadtables.R", echo=F))

#### Formatting the microdata ####
print("Loading micro data")
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

#### Set up LODES ####
print("Loading LODES data")
source("./s2.4_setuplodes.R", echo=F)


#### Set up marginals ####
print("Loading marginal data")
margins <- list()

#Setting up margin tables for individuals
print("Setting up individual margins")
source("./s2.5_indmargins.R", echo=F)

#Setting up margin tables for households
print("Setting up household margins")
source("./s2.6_hhmargins.R", echo=F)

#### Synthesizing separate populations  ####
system.time(source("./s3.1_IPF_worktrip.R", echo=F))
system.time(source("./s3.2_MCMC_worktrip.R", echo=F))
system.time(source("./s3.3_BN_worktrip.R", echo=F))

#### Running joint fitting ####
system.time(source("./s4.0_JointFit.R", echo=F))

#### Sampling the final population ####
system.time(source("./s5.0_Generation", echo=F))

