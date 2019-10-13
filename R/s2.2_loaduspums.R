##
# Population Synthesis for Boston Metropolitan Region, by Nicholas Marc Fournier
# Last updated January 2017.
#
# NOTICE:  All information, intellectual and technical concepts contained herein is,
# and remains the property of Nicholas Marc Fournier. Dissemination of this 
# information or reproduction of this material is strictly forbidden unless
# prior written permission is obtained from Nicholas Marc Fournier
##

#### Loading Tables ####
# Microsample
print("Loading micro data for US")
#pums
hhcols <- c("SERIALNO","NP","VEH","BLD","MRGX","HINCP","HHT","HUPAC","NR","NOC")
indcols <- c("SERIALNO","SPORDER","SEX","AGEP","RELP","ESR","WKHP","JWTR","JWMNP","INDP",
             "OCCP10","OCCP12","SOCP10","SOCP12","PWGTP","RAC1P","HISP","DRIVESP","SCHG",
             "SCH","PINCP","PUMA00","PUMA10","POWPUMA00","POWPUMA10")

pumsusind <- fread("./rdata/PUMS/ss15pus.csv", stringsAsFactors=F, sep = ",", integer64 = "character",
                            na.strings = c("N.A.","N.A.//"), data.table = T, select = indcols)
pumsushh <- fread("./rdata/PUMS/ss15hus.csv", stringsAsFactors=F, sep = ",", integer64 = "character",
                            na.strings = c("N.A.","N.A.//"), data.table = T, select = hhcols)

pumsusind[PUMA00==-9, PUMA00 :=  NA]
pumsusind[PUMA10==-9, PUMA10 :=  NA]
pumsusind[ is.na(PUMA10), PUMA := PUMA00]
pumsusind[ is.na(PUMA00), PUMA := PUMA10]
pumsusind[ , PUMA := sprintf("%05d", PUMA)]

pumsusind[POWPUMA00==-9, POWPUMA00 :=  NA]
pumsusind[POWPUMA10==-9, POWPUMA10 :=  NA]
pumsusind[ is.na(POWPUMA10), POWPUMA := POWPUMA00]
pumsusind[ is.na(POWPUMA00), POWPUMA := POWPUMA10]
pumsusind[ , POWPUMA := ifelse(is.na(POWPUMA),NA,sprintf("%05d", POWPUMA))]
sum(is.na(pumsusind$PUMA))


tablesus <- list()
tablesus[["uspumshh"]] <- pumsushh
tablesus[["uspumsind"]] <- pumsusind

rm(hhcols, indcols, pumsushh, pumsusind)

