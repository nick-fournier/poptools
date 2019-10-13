##
# Population Synthesis for Boston Metropolitan Region, by Nicholas Marc Fournier
# Last updated January 2017.
#
# NOTICE:  All information, intellectual and technical concepts contained herein is,
# and remains the property of Nicholas Marc Fournier. Dissemination of this 
# information or reproduction of this material is strictly forbidden unless
# prior written permission is obtained from Nicholas Marc Fournier
##


# Create data base folder 
if(!dir.exists("./db")){ dir.create("./db") }
#Loads packages
source("./s1.0_packages.R", echo=F)

#Loading tables
source("./s2.1_loadtables.R", echo=F)
print("Saving")
save(tables, file = "./db/d2.1_tables.RData")

# Checking for census tract consistency ####
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
# emptytracts <- c(emptytracts, c("25001014900", "25017354800", "25021417300",
#                                 "25025020303", "25025981502", "25025980300",
#                                 "25025981100", "25027732902"))
emptytracts<-unique(emptytracts)
tracts<-tracts[!(tracts %in% emptytracts)]

#Set up PUMS data ####
print("Setting up PUMS microdata")
source("./s2.3_pumsmicroformat.R", echo=F)
microdata <- pumsformatting(tables[["pumshh"]], tables[["pumsind"]])

#Set up MTS data ####
print("Setting up MTS microdata")
source("./s2.4_mtsmicrosample.R", echo=F)
print("Saving microdata")
save(microdata, pwork, tracts, file = "./db/d2.3_microdata.RData")

# #Set up PUMS US data ####
# if(!file.exists("./db/d2.4_microdataus.RData")) {
#   if(!file.exists("./db/d2.2_tablesus.RData")) {
#     source("./s2.2_loaduspums.R", echo=F)
#     print("Saving US PUMS")
#     save(tablesus, file = "./db/d2.2_tablesus.RData")
#   } else {
#     load("./db/d2.2_tablesus.RData")
#     print("s1.4 - Existing tables loaded")
#   }
#   microdata.us <- pumsformatting(tablesus[["uspumshh"]], tablesus[["uspumsind"]])
#   save(microdata.us, file = "./db/d2.4_microdataus.RData")
# }

#Set up marginals ####
margins<-list()

#Setting up margin tables for individuals
print("Setting up individual margins")
source("./s2.5_indmargins.R", echo=F)

#Setting up margin tables for households
print("Setting up household margins")
source("./s2.6_hhmargins.R", echo=F)

#Storing margins
print("Saving margins")
save(margins, tracts, microdata, file = "./db/d2.5_margins.RData")

rm(list=setdiff(ls(),c("clock","runtime","csvwriteout")))

# Microdata tables complete