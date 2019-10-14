hh <- tables[["pumshh"]]
ind <- tables[["pumsind"]]

#### Households ####
print("Setting up household PUMS")
hh <- hh[, .(SERIALNO,NP,VEH,BLD,MRGX,HINCP,HHT,HUPAC,NR)]
#Removing NA's
hh <- hh[!is.na(HINCP),]
#Formatting HH to match MTS & IPF
hh[ , OWN := ifelse(MRGX == 3 | is.na(MRGX), "RENT", "OWN")]
hh[ , HHSIZ := paste("HHSIZ", ifelse(NP > 4, 4, NP), sep = "")]
hh[ , HHVEH := paste("HHVEH", ifelse(VEH > 4, 4, VEH), sep = "")]

#hh units type
hh[ , RESTY := ifelse(BLD == 2 | BLD == 3, "UNITS1", NA)] # 1	Single Family Detached Dwelling
hh[ BLD == 4 | BLD == 5, RESTY := "UNITS2TO4"]            # 2	Building with 2-4 Units
hh[ BLD == 6 | BLD == 7, RESTY := "UNITS5TO19"]           # 3	Building with 5-19 Units
hh[ BLD == 8 | BLD == 9, RESTY := "UNITS20ORMORE"]        # 4	Building with 20 or More Units
hh[ BLD == 1 | BLD == 10, RESTY := "UNITSOTHER"]          # 4	Building with 20 or More Units
# Income1 <$15,000, # Income2 $15,000-$24,999, Income3 $25,000-$34,999, Income4 $35,000-$49,999,
# Income5 $50,000-$74,999, Income6 $75,000-$99,999, Income7 $100,000-$149,999, Income8 $150,000<
hh[ , INCOME := ifelse(HINCP < 15000, "INCOME1", NA)]
hh[ HINCP >= 15000 & HINCP < 25000, INCOME := "INCOME2", INCOME]
hh[ HINCP >= 25000 & HINCP < 35000, INCOME := "INCOME3", INCOME]
hh[ HINCP >= 35000 & HINCP < 50000, INCOME := "INCOME4", INCOME]
hh[ HINCP >= 50000 & HINCP < 75000, INCOME := "INCOME5", INCOME]
hh[ HINCP >= 75000 & HINCP < 100000, INCOME := "INCOME6", INCOME]
hh[ HINCP >= 100000 & HINCP < 150000, INCOME := "INCOME7", INCOME]
hh[ HINCP >= 150000, INCOME := "INCOME8", INCOME]

# #Personal income
# hh[ , PINCOME := ifelse(PINCP < 15000, "INCOME1", NA)]
# hh[ HINCP >= 15000 & HINCP < 25000, PINCOME := "INCOME2", PINCOME]
# hh[ HINCP >= 25000 & HINCP < 35000, PINCOME := "INCOME3", PINCOME]
# hh[ HINCP >= 35000 & HINCP < 50000, PINCOME := "INCOME4", PINCOME]
# hh[ HINCP >= 50000 & HINCP < 75000, PINCOME := "INCOME5", PINCOME]
# hh[ HINCP >= 75000 & HINCP < 100000, PINCOME := "INCOME6", PINCOME]
# hh[ HINCP >= 100000 & HINCP < 150000, PINCOME := "INCOME7", PINCOME]
# hh[ HINCP >= 150000, PINCOME := "INCOME8", PINCOME]

#hgend
hh <- merge(hh, ind[SPORDER == 1, .(SERIALNO,SEX)], by = "SERIALNO")
hh[ , HHOLDER := ifelse(SEX == 1, "MALEHEAD", "FEMALEHEAD")]

#hh race
hh <- merge(hh, ind[SPORDER == 1, .(SERIALNO,RAC1P)], by = "SERIALNO")
hh <- merge(hh, ind[SPORDER == 1, .(SERIALNO,HISP)], by = "SERIALNO")
hh[ , RACE := ifelse(RAC1P == 1, "WHITE", NA)]
hh[ RAC1P == 2, RACE := "BLACK"]
hh[ RAC1P == 3, RACE := "NATIVE"]
hh[ RAC1P == 4, RACE := "NATIVE"]
hh[ RAC1P == 5, RACE := "NATIVE"]
hh[ RAC1P == 6, RACE := "ASIAN"]
hh[ RAC1P == 7, RACE := "PACIFIC"]
hh[ RAC1P == 8, RACE := "OTHER"]
hh[ RAC1P == 9, RACE := "MULTI"]
hh[ HISP != 1, RACE := "HISPLAT"]

#hhworkers
hh <- merge(hh, data.frame(table(ind[ !is.na(WKHP), SERIALNO])),
            by.x = "SERIALNO", by.y = "Var1", all.x = T)
hh[ is.na(Freq), Freq := 0]
hh[ Freq > 3, Freq := 3]
hh[ , HHWRK := paste("HHWRK", Freq, sep = "")]

#hfam
hh[ , HFAM := ifelse(hh$HHT < 4 , "FAM", "NONFAM")]
#married?
hh[ , HMATE := ifelse(hh$HHT == 1, "MARRIED", "SINGLE")]
#Alone?
hh[ , HALONE := ifelse(hh$HHT == 4 | hh$HHT == 6, "ALONE", "NOTALONE")]
#kids?
hh[ , HKIDS := ifelse(HUPAC %in% 1:3, "KIDS", "NOKIDS")]
#nonrels
hh[ , HNONREL := ifelse(NR == 1, "NONREL", "NONONREL")]
#HH fam type category
hh[ , HFAMCAT := ifelse(HHT %in% c(4,6), "NONFAMALONE", NA)]
hh[ HHT %in% c(5,7), HFAMCAT := "NONFAMNOTALONE"]
hh[ HHT %in% 2:3 & HUPAC %in% 1:3, HFAMCAT := "FAMSINGLEKIDS"]
hh[ HHT %in% 2:3 & HUPAC == 4, HFAMCAT := "FAMSINGLENOKIDS"]
hh[ HHT == 1 & HUPAC %in% 1:3, HFAMCAT := "FAMMARRIEDKIDS"]
hh[ HHT == 1 & HUPAC == 4, HFAMCAT := "FAMMARRIEDNOKIDS"]
##
hh[ , HFAMCOMP := ifelse(HHT == 1, "FAMMARRIED", NA)]
hh[ , HFAMCOMP := ifelse(HHT %in% c(2,3), "FAMSINGLE", NA)]
hh[ HFAMCAT == "NONFAMALONE" , HFAMCOMP := "NONFAMALONE"]
hh[ HFAMCAT == "NONFAMNOTALONE" , HFAMCOMP := "NONFAMNOTALONE"]


#### Individuals ####
print("Setting up individuals PUMS")
#Subsetting columns we need
ind <- ind[ , .(SERIALNO,SPORDER,SEX,AGEP,RELP,ESR,WKHP,JWTR,JWMNP,INDP,
                OCCP10,OCCP12,SOCP10,SOCP12,PWGTP,RAC1P,HISP,DRIVESP,SCHG,PINCP)]
#ind <- ind[ind$SERIALNO %in% hh$SERIALNO,]
#gender
ind[ , GEND := ifelse(SEX == 1, "MALE", "FEMALE")]
#ind race
ind[ , RACE := ifelse(RAC1P == 1, "WHITE", NA)]
ind[ RAC1P == 2, RACE := "BLACK"]
ind[ RAC1P == 3, RACE := "NATIVE"]
ind[ RAC1P == 4, RACE := "NATIVE"]
ind[ RAC1P == 5, RACE := "NATIVE"]
ind[ RAC1P == 6, RACE := "ASIAN"]
ind[ RAC1P == 7, RACE := "PACIFIC"]
ind[ RAC1P == 8, RACE := "OTHER"]
ind[ RAC1P == 9, RACE := "MULTI"]
ind[ HISP != 1, RACE := "HISPLAT"]

#age
ind[ , AGE := ifelse(AGEP < 10, "AGE0TO9", NA)]
ind[ AGEP >= 10 & AGEP < 15, AGE := "AGE10TO14"]
ind[ AGEP >= 15 & AGEP < 20, AGE := "AGE15TO19"]
ind[ AGEP >= 20 & AGEP < 25, AGE := "AGE20TO24"]
ind[ AGEP >= 25 & AGEP < 45, AGE := "AGE25TO44"]
ind[ AGEP >= 45 & AGEP < 55, AGE := "AGE45TO54"]
ind[ AGEP >= 55 & AGEP < 65, AGE := "AGE55TO64"]
ind[ AGEP >= 65, AGE := "AGE65UP"]
#relate
ind[ , RELATE := ifelse(RELP == 0, "HEAD", NA)]
ind[ RELP == 1, RELATE := "SPOUSE"]
ind[ RELP %in% c(2,3,4,7,9,14), RELATE := "CHILD"]
ind[ RELP %in% c(5,6,8,10), RELATE := "RELATIVE"]
ind[ is.na(RELATE), RELATE := "NONRELATIVE"]
#Hours
ind[ is.na(WKHP), WKHP := 0]
ind[ , HOURS := ifelse(WKHP == 0, "HOURS0", NA)]
ind[WKHP < 35 & WKHP >= 1, HOURS := "HOURS1TO34"]
ind[WKHP >= 35, HOURS := "HOURS35UP"]
#Works
ind[ , WORKS := ifelse(HOURS != "HOURS0", "EMP", "UNEMP")]
#Work mode
ind[ is.na(JWTR), JWTR := 0]
ind[ is.na(DRIVESP), DRIVESP := 0]
ind[ , WMODE := ifelse(WORKS == "UNEMP", "NONWORK", NA)]
ind[ JWTR == 0 & WORKS == "EMP" , WMODE := "WORKSHOME"]
ind[ JWTR == 1, WMODE := "DRIVE"]
ind[ JWTR == 10, WMODE := "WALK"]
ind[ JWTR == 11, WMODE := "WORKSHOME"]
ind[ JWTR > 1 & JWTR < 7, WMODE := "TRANSIT"]
ind[ JWTR >= 7 & JWTR < 10 | JWTR == 12, WMODE := "OTHERMODE"]
ind[ DRIVESP > 1, WMODE := "RIDEPOOL"]
#School enroll
ind[ is.na(ind$SCHG), SCHG := 0]
ind[ , SCHOL := ifelse(SCHG == 0, "NOTENROLLED", "ENROLLED")]
#Hfam
ind <- merge(ind, hh[, .(SERIALNO,HFAM)], by = "SERIALNO", all = T)
#Travel time
ind[ is.na(JWMNP), JWMNP := 0]
ind[ , TTIME := ifelse(JWMNP == 0, "MINS0", NA)]
ind[ JWMNP > 0   & JWMNP < 15, TTIME := "MINS1TO14"]
ind[ JWMNP >= 15 & JWMNP < 35, TTIME := "MINS15TO34"]
ind[ JWMNP >= 35 & JWMNP < 45, TTIME := "MINS35TO44"]
ind[ JWMNP >= 45 & JWMNP < 60, TTIME := "MINS45TO59"]
ind[ JWMNP >= 60, TTIME := "MINS60UP"]
#occupation
occkey <- data.table(tables[["pumsocc"]])
ind[ , OCCP := OCCP12]
ind[ is.na(OCCP12), OCCP := OCCP10[is.na(OCCP12)]]
ind[ is.na(OCCP), OCCP := 0]
ind[ , OCC := NA]
for(i in 1:nrow(occkey)) {
  ind[ , OCC := ifelse(OCCP >= occkey[i,START] & OCCP <= occkey[i,END], occkey[i,OCC], OCC)]
}

#industry
induskey <- data.table(tables[["pumsindus"]])
ind[is.na(INDP), INDP := 0]
ind[ , INDUS := NA]
for(i in 1:nrow(induskey)) {
  ind[ , INDUS := ifelse(INDP >= induskey[i,START] & INDP <= induskey[i,END], induskey[i,INDUS], INDUS)]
}
#Retired/unemployed but reported industry/occ
ind[HOURS=="HOURS0", OCC := "NOOCC"]
ind[HOURS=="HOURS0", INDUS := "NOIND"]


#### Factor order ####
print("Cleaning up")
hh[ , HHSIZ := factor(HHSIZ, levels = c("HHSIZ1","HHSIZ2","HHSIZ3","HHSIZ4"))]
hh[ , HHVEH := factor(HHVEH, levels = c("HHVEH0","HHVEH1","HHVEH2","HHVEH3","HHVEH4"))]
hh[ , RACE := factor(RACE, levels = c("ASIAN","BLACK","HISPLAT","MULTI","NATIVE","OTHER","PACIFIC","WHITE"))]
hh[ , RESTY := factor(RESTY, levels = c("UNITS1","UNITS2TO4","UNITS5TO19","UNITS20ORMORE","UNITSOTHER"))]
hh[ , OWN := factor(OWN, levels = c("OWN","RENT"))]
hh[ , INCOME := factor(INCOME, levels = c("INCOME1","INCOME2","INCOME3","INCOME4","INCOME5","INCOME6","INCOME7","INCOME8"))]
hh[ , HHWRK := factor(HHWRK, levels = c("HHWRK0","HHWRK1","HHWRK2","HHWRK3"))]
hh[ , HFAM := factor(HFAM, levels = c("FAM","NONFAM"))]
hh[ , HHOLDER := factor(HHOLDER, levels = c("MALEHEAD","FEMALEHEAD"))]
hh[ , HMATE := factor(HMATE, levels = c("MARRIED","SINGLE"))]
hh[ , HKIDS := factor(HKIDS, levels = c("KIDS","NOKIDS"))]
hh[ , HALONE := factor(HALONE, levels = c("ALONE","NOTALONE"))]
hh[ , HNONREL := factor(HNONREL, levels = c("NONREL","NONONREL"))]
hh[ , HFAMCAT := factor(HFAMCAT, levels = c("NONFAMNOTALONE","NONFAMALONE","FAMSINGLEKIDS","FAMSINGLENOKIDS","FAMMARRIEDKIDS","FAMMARRIEDNOKIDS"))]
hh[ , HFAMCOMP := factor(HFAMCOMP, levels = c("NONFAMNOTALONE","NONFAMALONE","FAMSINGLE","FAMMARRIED"))]

ind[ , GEND := factor(GEND, levels = c("MALE","FEMALE"))]
ind[ , AGE := factor(AGE, levels = c("AGE0TO9","AGE10TO14","AGE15TO19","AGE20TO24","AGE25TO44","AGE45TO54","AGE55TO64","AGE65UP"))]
ind[ , RELATE := factor(RELATE, levels = c("HEAD","SPOUSE","CHILD","RELATIVE","NONRELATIVE"))]
ind[ , WORKS := factor(WORKS, levels = c("EMP", "UNEMP"))]
ind[ , HOURS := factor(HOURS, levels = c("HOURS0","HOURS1TO34","HOURS35UP"))]
ind[ , WMODE := factor(WMODE, levels = c("DRIVE","RIDEPOOL","TRANSIT","WALK","WORKSHOME","OTHERMODE","NONWORK"))]
ind[ , SCHOL := factor(SCHOL, levels = c("ENROLLED","NOTENROLLED"))]
ind[ , RACE := factor(RACE, levels = c("ASIAN","BLACK","HISPLAT","MULTI","NATIVE","OTHER","PACIFIC","WHITE"))]
ind[ , HFAM := factor(HFAM, levels = c("FAM","NONFAM"))]
ind[ , OCC := factor(OCC, levels = occkey$OCC[-which(duplicated(occkey$OCC, fromLast = T))])]
ind[ , INDUS := factor(INDUS, levels = induskey$INDUS[-which(duplicated(induskey$INDUS, fromLast = T))])]
ind[ , TTIME := factor(TTIME, levels =  c("MINS0","MINS1TO14","MINS15TO34","MINS35TO44","MINS45TO59","MINS60UP"))]

#### Store ####
#trimming columns
microdata <- list()
microdata[["pums"]] <- merge(ind[,.(SERIALNO,SPORDER,GEND,AGE,RELATE,WORKS,HOURS,WMODE,SCHOL,OCC,INDUS,TTIME,PINCP)],
                         hh[,.(SERIALNO,HHSIZ,HHVEH,RACE,RESTY,INCOME,OWN,HHWRK,HFAM,HHOLDER,HMATE,HKIDS,HALONE,HNONREL,HFAMCAT,HFAMCOMP)],
                         by = "SERIALNO")

microdata[["pums_hh"]] <- microdata[["pums"]][ ,.(HHSIZ,HHVEH,RACE,RESTY,INCOME,OWN,HHWRK,HFAM,HHOLDER,HMATE,HKIDS,HALONE,HNONREL,HFAMCAT,HFAMCOMP)]
microdata[["pums_ind"]] <- microdata[["pums"]][ ,.(GEND,AGE,RELATE,WORKS,HOURS,WMODE,SCHOL,HFAM,OCC,INDUS,TTIME)]

rm(i, occkey, induskey, ind, hh)
# Microdata tables complete
