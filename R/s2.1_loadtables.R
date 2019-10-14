#### Loading Tables ####
# Microsample
tables<-list()
tables[["mtshh"]]  <- read.csv("./rdata/MTS/HH.csv", stringsAsFactors=F)
tables[["mtsind"]] <- read.csv("./rdata/MTS/PER.csv", stringsAsFactors=F)
tables[["mtsplace"]] <- read.csv("./rdata/MTS/PLACE.csv", stringsAsFactors=F)
tables[["mtsveh"]] <- read.csv("./rdata/MTS/VEH.csv", stringsAsFactors=F)
#pums
hhcols <- c("SERIALNO","NP","VEH","BLD","MRGX","HINCP","HHT","HUPAC","NR","NOC")
indcols <- c("SERIALNO","SPORDER","SEX","AGEP","RELP","ESR","WKHP","JWTR","JWMNP","INDP",
             "OCCP10","OCCP12","SOCP10","SOCP12","PWGTP","RAC1P","HISP","DRIVESP","SCHG",
             "SCH","PINCP","PUMA00","PUMA10","POWPUMA00","POWPUMA10")

mahh <- fread("./rdata/PUMS/ss15hma.csv", stringsAsFactors=F, sep = ",", integer64 = "character",
                            na.strings = c("N.A.","N.A.//"), data.table = T, select = hhcols)
maind <- fread("./rdata/PUMS/ss15pma.csv",stringsAsFactors=F, sep = ",", integer64 = "character",
               na.strings = c("N.A.","N.A.//"), data.table = T, select = indcols)

rihh <- fread("./rdata/PUMS/ss15hri.csv", stringsAsFactors=F, sep = ",", integer64 = "character",
                            na.strings = c("N.A.","N.A.//"), data.table = T, select = hhcols)
riind <- fread("./rdata/PUMS/ss15pri.csv",stringsAsFactors=F, sep = ",", integer64 = "character",
               na.strings = c("N.A.","N.A.//"), data.table = T, select = indcols)

cthh <- fread("./rdata/PUMS/ss15hct.csv", stringsAsFactors=F, sep = ",", integer64 = "character",
                            na.strings = c("N.A.","N.A.//"), data.table = T, select = hhcols)
ctind <- fread("./rdata/PUMS/ss15pct.csv",stringsAsFactors=F, sep = ",", integer64 = "character",
               na.strings = c("N.A.","N.A.//"), data.table = T, select = indcols)

nhhh <- fread("./rdata/PUMS/ss15hnh.csv", stringsAsFactors=F, sep = ",", integer64 = "character",
                            na.strings = c("N.A.","N.A.//"), data.table = T, select = hhcols)
nhind <- fread("./rdata/PUMS/ss15pnh.csv",stringsAsFactors=F, sep = ",", integer64 = "character",
               na.strings = c("N.A.","N.A.//"), data.table = T, select = indcols)

meind <- fread("./rdata/PUMS/ss15pme.csv",stringsAsFactors=F, sep = ",", integer64 = "character",
                             na.strings = c("N.A.","N.A.//"), data.table = T, select = indcols)
mehh <- fread("./rdata/PUMS/ss15hme.csv", stringsAsFactors=F, sep = ",", integer64 = "character",
              na.strings = c("N.A.","N.A.//"), data.table = T, select = hhcols)

vtind <- fread("./rdata/PUMS/ss15pvt.csv",stringsAsFactors=F, sep = ",", integer64 = "character",
               na.strings = c("N.A.","N.A.//"), data.table = T, select = indcols)
vthh <- fread("./rdata/PUMS/ss15hvt.csv", stringsAsFactors=F, sep = ",", integer64 = "character",
              na.strings = c("N.A.","N.A.//"), data.table = T, select = hhcols)

nyind <- fread("./rdata/PUMS/ss15pny.csv",stringsAsFactors=F, sep = ",", integer64 = "character",
               na.strings = c("N.A.","N.A.//"), data.table = T, select = indcols)
nyhh <- fread("./rdata/PUMS/ss15hny.csv", stringsAsFactors=F, sep = ",", integer64 = "character",
              na.strings = c("N.A.","N.A.//"), data.table = T, select = hhcols)

pumshh <- rbindlist(list(mahh, rihh, cthh, nhhh, mehh, vthh, nyhh))
pumsind <- rbindlist(list(maind, riind, ctind, nhind, meind, vtind, nyind))
rm(cthh, ctind, mahh, maind, nhhh,nhind, rihh, riind, mehh, meind, vthh, vtind, nyhh, nyind)
# pumshh <- mahh
# pumsind <- maind
# rm(mahh, maind)

pumsind[PUMA00==-9, PUMA00 :=  NA]
pumsind[PUMA10==-9, PUMA10 :=  NA]
pumsind[ is.na(PUMA10), PUMA := PUMA00]
pumsind[ is.na(PUMA00), PUMA := PUMA10]
pumsind[ , PUMA := sprintf("%05d", PUMA)]

pumsind[POWPUMA00==-9, POWPUMA00 :=  NA]
pumsind[POWPUMA10==-9, POWPUMA10 :=  NA]
pumsind[ is.na(POWPUMA10), POWPUMA := POWPUMA00]
pumsind[ is.na(POWPUMA00), POWPUMA := POWPUMA10]
pumsind[ , POWPUMA := ifelse(is.na(POWPUMA),NA,sprintf("%05d", POWPUMA))]

sum(is.na(pumsind$PUMA))

tables[["pumshh"]] <- pumshh
tables[["pumsind"]] <- pumsind

tables[["pumsocc"]] <- read.csv("./rdata/PUMS/occupationkey.csv", stringsAsFactors=F)
tables[["pumsindus"]] <- read.csv("./rdata/PUMS/industrykey.csv", stringsAsFactors=F)

# Marginal Census Tables

print("Loading marginal data for individuals")
tables[["indagesex"]]   <- read.csv("./rdata/marginals/ACS_15_5YR_B01001_with_ann.csv", stringsAsFactors=F) # <<<
tables[["indrelate"]]   <- read.csv("./rdata/marginals/ACS_15_5YR_B09019_with_ann.csv", stringsAsFactors=F) # <<<
tables[["indindusocc"]] <- read.csv("./rdata/marginals/ACS_15_5YR_C24050_with_ann.csv", stringsAsFactors=F) # <<<
tables[["hhinc"]]       <- read.csv("./rdata/marginals/ACS_15_5YR_B19001_with_ann.csv", stringsAsFactors=F) # <<<
tables[["hhdwellsize"]] <- read.csv("./rdata/marginals/ACS_15_5YR_B25124_with_ann.csv", stringsAsFactors=F) # <<<
tables[["hhvehsize"]]   <- read.csv("./rdata/marginals/ACS_15_5YR_B08201_with_ann.csv", stringsAsFactors=F) # <<<
tables[["hhraceA"]]     <- read.csv("./rdata/marginals/ACS_15_5YR_B11001A_with_ann.csv", stringsAsFactors=F) # <<<
tables[["hhraceB"]]     <- read.csv("./rdata/marginals/ACS_15_5YR_B11001B_with_ann.csv", stringsAsFactors=F) # <<<
tables[["hhraceC"]]     <- read.csv("./rdata/marginals/ACS_15_5YR_B11001C_with_ann.csv", stringsAsFactors=F) # <<<
tables[["hhraceD"]]     <- read.csv("./rdata/marginals/ACS_15_5YR_B11001D_with_ann.csv", stringsAsFactors=F) # <<<
tables[["hhraceE"]]     <- read.csv("./rdata/marginals/ACS_15_5YR_B11001E_with_ann.csv", stringsAsFactors=F) # <<<
tables[["hhraceF"]]     <- read.csv("./rdata/marginals/ACS_15_5YR_B11001F_with_ann.csv", stringsAsFactors=F) # <<<
tables[["hhraceG"]]     <- read.csv("./rdata/marginals/ACS_15_5YR_B11001G_with_ann.csv", stringsAsFactors=F) # <<<
tables[["hhraceH"]]     <- read.csv("./rdata/marginals/ACS_15_5YR_B11001H_with_ann.csv", stringsAsFactors=F) # <<<


# Individuals
# tables[["indagesex"]]   <- read.csv("./rdata/marginals/DEC_10_SF1_SF1DP1_with_ann.csv", stringsAsFactors=F)  # <<<
# tables[["indrelate"]]   <- read.csv("./rdata/marginals/DEC_10_SF1_P29_with_ann.csv", stringsAsFactors=F) # <<<
# tables[["indworkagesex"]]  <- read.csv("./rdata/marginals/ACS_10_5YR_B23001_with_ann.csv", stringsAsFactors=F)
# tables[["indoccsex"]] <- read.csv("./rdata/marginals/ACS_10_5YR_C24010_with_ann.csv", stringsAsFactors=F)
# tables[["indindusocc"]] <- read.csv("./rdata/marginals/ACS_10_5YR_C24050_with_ann.csv", stringsAsFactors=F) # <<<
# 
# tables[["indindussex"]]  <- read.csv("./rdata/marginals/ACS_10_5YR_C24030_with_ann.csv", stringsAsFactors=F)
# tables[["indindusmode"]]  <- read.csv("./rdata/marginals/ACS_10_5YR_B08126_with_ann.csv", stringsAsFactors=F)
# 
# tables[["indschoolagesex"]]<- read.csv("./rdata/marginals/ACS_10_5YR_B14003_with_ann.csv", stringsAsFactors=F)#
# tables[["indschool"]]<- read.csv("./rdata/marginals/ACS_10_5YR_B14001_with_ann.csv", stringsAsFactors=F)
# tables[["indmodeage"]]  <- read.csv("./rdata/marginals/ACS_10_5YR_B08101_with_ann.csv", stringsAsFactors=F)
# tables[["indmodesex"]]  <- read.csv("./rdata/marginals/ACS_10_5YR_B08006_with_ann.csv", stringsAsFactors=F)
# tables[["indmodettime"]]<- read.csv("./rdata/marginals/ACS_10_5YR_C08134_with_ann.csv", stringsAsFactors=F)#
# tables[["indttimesex"]] <- read.csv("./rdata/marginals/ACS_10_5YR_B08012_with_ann.csv", stringsAsFactors=F)
# tables[["indttime"]]    <- read.csv("./rdata/marginals/ACS_10_5YR_B08303_with_ann.csv", stringsAsFactors=F)
# tables[["indhrssex"]]   <- read.csv("./rdata/marginals/ACS_15_5YR_B23022_with_ann.csv", stringsAsFactors=F) 
# tables[["indhrssex2"]]  <- read.csv("./rdata/marginals/ACS_15_5YR_B23026_with_ann.csv", stringsAsFactors=F)
# tables[["indparttime"]] <- read.csv("./rdata/marginals/ACS_10_5YR_B17004_with_ann.csv", stringsAsFactors=F)
# 
# # Households
# print("Loading marginal data for households")
# tables[["hhinc"]]       <- read.csv("./rdata/marginals/ACS_10_5YR_B19001_with_ann.csv", stringsAsFactors=F) # <<<
# tables[["hhworksize"]]  <- read.csv("./rdata/marginals/ACS_10_5YR_B08202_with_ann.csv", stringsAsFactors=F) 
# tables[["hhdwellsize"]] <- read.csv("./rdata/marginals/ACS_10_5YR_B25124_with_ann.csv", stringsAsFactors=F) # <<<
# tables[["hhfam"]]       <- read.csv("./rdata/marginals/DEC_10_SF1_QTP11_with_ann.csv", stringsAsFactors=F)
# tables[["hhnonrel"]]    <- read.csv("./rdata/marginals/DEC_10_SF1_P27_with_ann.csv", stringsAsFactors=F)
# tables[["hhownage"]]    <- read.csv("./rdata/marginals/ACS_10_5YR_B25007_with_ann.csv", stringsAsFactors=F)
# tables[["hhvehsize"]]   <- read.csv("./rdata/marginals/ACS_10_5YR_B08201_with_ann.csv", stringsAsFactors=F) # <<<
# tables[["hhraceA"]]     <- read.csv("./rdata/marginals/ACS_10_5YR_B11001A_with_ann.csv", stringsAsFactors=F) # <<<
# tables[["hhraceB"]]     <- read.csv("./rdata/marginals/ACS_10_5YR_B11001B_with_ann.csv", stringsAsFactors=F) # <<<
# tables[["hhraceC"]]     <- read.csv("./rdata/marginals/ACS_10_5YR_B11001C_with_ann.csv", stringsAsFactors=F) # <<<
# tables[["hhraceD"]]     <- read.csv("./rdata/marginals/ACS_10_5YR_B11001D_with_ann.csv", stringsAsFactors=F) # <<<
# tables[["hhraceE"]]     <- read.csv("./rdata/marginals/ACS_10_5YR_B11001E_with_ann.csv", stringsAsFactors=F) # <<<
# tables[["hhraceF"]]     <- read.csv("./rdata/marginals/ACS_10_5YR_B11001F_with_ann.csv", stringsAsFactors=F) # <<<
# tables[["hhraceG"]]     <- read.csv("./rdata/marginals/ACS_10_5YR_B11001G_with_ann.csv", stringsAsFactors=F) # <<<
# tables[["hhraceH"]]     <- read.csv("./rdata/marginals/ACS_10_5YR_B11001H_with_ann.csv", stringsAsFactors=F) # <<<
# 
# tables[["hhsize"]]      <- read.csv("./rdata/marginals/DEC_10_SF1_QTH2_with_ann.csv", stringsAsFactors=F)
# Household size for proportional assignment
# 
# #Vehicles
tables[["aggvehs"]]     <- read.csv("./rdata/marginals/ACS_15_5YR_B25046_with_ann.csv", stringsAsFactors=F)

# Workplace
tables[["workplace"]]   <- read.csv("./rdata/marginals/BP_2010_00A1_with_ann.csv", stringsAsFactors=F)
#schools
tables[["enroll"]]   <- read.csv("./rdata/spatial/enrollment.csv", stringsAsFactors=F, numeral="no.loss")
tables[["taz_xy"]]   <- read.csv("./rdata/spatial/taz_xy.csv", stringsAsFactors=F, numeral="no.loss")
tables[["schools_xy"]]   <- read.csv("./rdata/spatial/schools_tazxy.csv", stringsAsFactors=F, numeral="no.loss")
tables[["colleges_xy"]]   <- read.csv("./rdata/spatial/colleges_tazxy.csv", stringsAsFactors=F, numeral="no.loss")
tables[["tract2puma"]]   <- read.csv("./rdata/spatial/2010_Census_Tract_to_2010_PUMA.txt",stringsAsFactors=F, colClasses = "character")

rm(hhcols, indcols, pumshh, pumsind)
#checking if tracts in tables match
all(unique(tables[["indagesex"]]$GEO.id2)==unique(tables[["hhinc"]]$GEO.id2))
all(unique(tables[["hhdwellsize"]]$GEO.id2)==unique(tables[["hhinc"]]$GEO.id2))
all(unique(tables[["hhdwellsize"]]$GEO.id2)==unique(tables[["hhworksize"]]$GEO.id2))
all(unique(tables[["indrelate"]]$GEO.id2)==unique(tables[["hhworksize"]]$GEO.id2))
all(unique(tables[["indrelate"]]$GEO.id2)==unique(tables[["hhnonrel"]]$GEO.id2))

