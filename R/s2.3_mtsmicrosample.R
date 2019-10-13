##
# Population Synthesis for Boston Metropolitan Region, by Nicholas Marc Fournier
# Last updated January 2017.
#
# NOTICE:  All information, intellectual and technical concepts contained herein is,
# and remains the property of Nicholas Marc Fournier. Dissemination of this 
# information or reproduction of this material is strictly forbidden unless
# prior written permission is obtained from Nicholas Marc Fournier
##

##### Ininitalize Individual Variables for IPF #####
hh <- tables[["mtshh"]]
ind  <- tables[["mtsind"]]
place <- tables[["mtsplace"]]
#
##### Creating list of all tracts from Census ####
#List of usable tracts, must match between HH and Inds



##### Individuals #####
  ##### Employment Status
  ind$WORKS<-ifelse(ind$WORKS==1, "EMP", ind$WORKS)
  ind$WORKS<-ifelse(ind$WORKS==2, "UNEMP", ind$WORKS)
  ind$WORKS<-ifelse(ind$AGE<=15, "UNEMP", ind$WORKS)
  ind$WORKS<-ifelse(ind$AGE==99 & ind$AGEB==1, "UNEMP", ind$WORKS)

  ##### School status
  ind$SCHOL <- ifelse(is.na(ind$SCHOL), 0, ind$SCHOL)
  ind$SCHOL <- ifelse(ind$ENROL==1, "ENROLLED","NOTENROLLED")
  
  # ind$SCHOL <- ifelse(!(ind$SCHOL %in% c(1:8)), "NONSTUDENT", ind$SCHOL)
  # ind$SCHOL <- ifelse(ind$SCHOL > 0 & ind$SCHOL < 5, "PTO12", ind$SCHOL)
  # ind$SCHOL <- ifelse(ind$SCHOL==7, "UGRAD", ind$SCHOL)
  # ind$SCHOL <- ifelse(ind$SCHOL==8|ind$SCHOL==5|ind$SCHOL==6, "GRADPROF", ind$SCHOL)
  
  ##### Work hours
  ind$HOURS<-ifelse(ind$WORKS=="UNEMP", 0, ind$HOURS)
  ind$WDAYS<-ifelse(ind$WORKS=="UNEMP", 0, ind$WDAYS)
  ind$WDAYS<-ifelse(ind$HOURS==24 & ind$WDAYS > 3, 3, ind$WDAYS)
  ind$WDAYS<-ifelse(ind$HOURS>24 & ind$WDAYS > 7, 5, ind$WDAYS)
  ind$HOURS<-ifelse(ind$HOURS>24 & ind$WDAYS > 7, 8, ind$HOURS)
  dayhours<-table(subset(ind, ind$WDAYS<8 & ind$WDAYS>0 & ind$HOURS < 25, select=c("WDAYS","HOURS")))
  dayhours<-dayhours/sum(dayhours)
  
  #sampling in work hours per day
  for(i in unique(ind$WDAYS[ind$HOURS>24 & ind$WDAYS<8])) {
    days<-as.character(i)
    sel <- which(ind$WDAYS[ind$HOURS>24] == i)
    ind$HOURS[ind$HOURS>24][sel] <- as.numeric(sample(names(dayhours[days,]), length(sel), T, prob = dayhours[days,]))
  }
  
  #sampling in work days per week
  for(i in unique(ind$HOURS[ind$WDAYS>7 & ind$HOURS<98])) {
    hrs<-as.character(i)
    sel <- which(ind$HOURS[ind$WDAYS==8|ind$WDAYS==9] == i)
    ind$WDAYS[ind$WDAYS==8|ind$WDAYS==9][sel] <- as.numeric(sample(1:7, length(sel), T, prob = dayhours[,hrs]))
  }
  rm(hrs, sel, days)
  
  #weekly working hours
  ind$WHOURS<-ind$HOURS * ind$WDAYS
  ind$HOURS<-ifelse(ind$WHOURS == 0, "HOURS0", ind$HOURS)
  ind$HOURS<-ifelse(ind$WHOURS > 0 & ind$WHOURS < 35, "HOURS1TO34", ind$HOURS)
  ind$HOURS<-ifelse(ind$WHOURS >= 35, "HOURS35UP", ind$HOURS)
  
  ##### Mode to Work
  # Need to go thru WMODE==97
  ind$WMODE<-ifelse(ind$WMODE==0, "WORKSHOME", ind$WMODE)
  ind$WMODE<-ifelse(ind$WMODE==1, "WALK", ind$WMODE)
  #ind$WMODE<-ifelse(ind$WMODE==2, "BICYCLE", ind$WMODE)
  ind$WMODE<-ifelse(ind$WMODE==3, "DRIVE", ind$WMODE)
  ind$WMODE<-ifelse(ind$WMODE==4 | ind$WMODE==9, "RIDEPOOL", ind$WMODE)
  ind$WMODE<-ifelse(ind$WMODE==5, "TRANSIT", ind$WMODE)
  ind$WMODE<-ifelse(ind$WMODE==6|ind$WMODE==7|ind$WMODE==8|ind$WMODE==97|ind$WMODE==98|ind$WMODE==99|ind$WMODE==2, "OTHERMODE", ind$WMODE)
  ind$WMODE<-ifelse(is.na(ind$WMODE) & ind$AGE<16, "NONWORK", ind$WMODE)
  ind$WMODE<-ifelse(is.na(ind$WMODE) & ind$WORKS=="UNEMP", "NONWORK", ind$WMODE)
  ind<-ind[!is.na(ind$WMODE),]
  
  ##### Age
  ind$AGE<-ifelse(ind$AGE<9, "AGE0TO9", ind$AGE)
  ind$AGE<-ifelse(ind$AGE>=10 & ind$AGE<15, "AGE10TO14", ind$AGE)
  ind$AGE<-ifelse(ind$AGE>=15 & ind$AGE<20, "AGE15TO19", ind$AGE)
  ind$AGE<-ifelse(ind$AGE>=20 & ind$AGE<25, "AGE20TO24", ind$AGE)
  ind$AGE<-ifelse(ind$AGE>=25 & ind$AGE<45, "AGE25TO44", ind$AGE)
  ind$AGE<-ifelse(ind$AGE>=45 & ind$AGE<55, "AGE45TO54", ind$AGE)
  ind$AGE<-ifelse(ind$AGE>=55 & ind$AGE<65, "AGE55TO64", ind$AGE)
  ind$AGE<-ifelse(ind$AGE>=65 & ind$AGE<99, "AGE65UP", ind$AGE)
  ind<-ind[!ind$AGE==99,] #ditching 98"s and 99"s for now
  
  ##### Gender
  ind$GEND<-ifelse(ind$GEND==1, "MALE", ind$GEND)
  ind$GEND<-ifelse(ind$GEND==2, "FEMALE", ind$GEND)
  ind<-ind[ind$GEND!=9,] #ditching 9"s for now
  
  ##### Ind Household status
  ind$RELATE<-ifelse(ind$RELATE==1, "HEAD", ind$RELATE)
  ind$RELATE<-ifelse(ind$RELATE==2, "SPOUSE", ind$RELATE)
  ind$RELATE<-ifelse(ind$RELATE==3, "CHILD", ind$RELATE)
  ind$RELATE<-ifelse(ind$RELATE==4 | ind$RELATE==5, "RELATIVE", ind$RELATE)
  ind$RELATE<-ifelse(ind$RELATE==6, "NONRELATIVE", ind$RELATE)
  ind<-ind[!ind$RELATE==9,] #ditching 98"s and 99"s for now
  
  ##### Ind/HH home and work Tract/block ID
  #home tract
  hh$tractID<-paste(sprintf("%02d",hh$HSTATE10),
                    sprintf("%03d",hh$HCOUNTY10),
                    sprintf("%06d",hh$HTRACT10),sep = "")
  #home block
  hh$hblock<-paste(sprintf("%02d",hh$HSTATE10),
                    sprintf("%03d",hh$HCOUNTY10),
                    sprintf("%06d",hh$HTRACT10),
                    sprintf("%04d",hh$HBLOCK10),sep = "")
  #work block
  ind$wblock<-paste(sprintf("%02d",ind$WSTATE10),
                    sprintf("%03d",ind$WCOUNTY10),
                    sprintf("%06d",ind$WTRACT10),
                    sprintf("%04d",ind$WBLOCK10),sep = "")
  ind$wblock[is.na(ind$WBLOCK10)] <- NA
  
  #send home tract to individuals
  ind<-merge(ind, hh[, c("SAMPN", "tractID")], by="SAMPN")
  ind<-merge(ind, hh[, c("SAMPN", "hblock")], by="SAMPN")

##### Households ####
  ##### Dwelling unit type  #$RESTY is hometype
  hh<-hh[hh$RESTY!=8 & hh$RESTY!=9,]
  hh$RESTY<-ifelse(hh$RESTY==1, "UNITS1", hh$RESTY)   # 1	Single Family Detached Dwelling
  hh$RESTY<-ifelse(hh$RESTY==2, "UNITS2TO4", hh$RESTY)     # 2	Building with 2-4 Units
  hh$RESTY<-ifelse(hh$RESTY==3, "UNITS5TO19", hh$RESTY)   # 3	Building with 5-19 Units
  hh$RESTY<-ifelse(hh$RESTY==4, "UNITS20ORMORE", hh$RESTY)   # 4	Building with 20 or More Units
  hh$RESTY<-ifelse(hh$RESTY==7, "UNITSOTHER", hh$RESTY)   # 4	Building with 20 or More Units
  ##### Rent or own home
  hh<-hh[!hh$OWN>=7,]
  hh$OWN<-ifelse(hh$OWN==1, "OWN", hh$OWN)   # 1	Own or mortage
  hh$OWN<-ifelse(hh$OWN==2, "RENT", hh$OWN)  # 2 rents
  
  ##### Household Size
  hh<-hh[!hh$HHSIZ>=98,]
  hh$HHSIZ<-ifelse(hh$HHSIZ>=4, 4, hh$HHSIZ)   # Household with 4 or more people
  hh$HHSIZ<-ifelse(hh$HHSIZ==1, "HHSIZ1", hh$HHSIZ)
  hh$HHSIZ<-ifelse(hh$HHSIZ==2, "HHSIZ2", hh$HHSIZ)
  hh$HHSIZ<-ifelse(hh$HHSIZ==3, "HHSIZ3", hh$HHSIZ)
  hh$HHSIZ<-ifelse(hh$HHSIZ==4, "HHSIZ4", hh$HHSIZ)
  
  ##### Household Vehicles
  hh<-hh[!hh$HHVEH>=98,]
  hh$HHVEH<-ifelse(hh$HHVEH>=4, 4, hh$HHVEH)   # Household with 4 or more people
  hh$HHVEH<-ifelse(hh$HHVEH==0, "HHVEH0", hh$HHVEH)
  hh$HHVEH<-ifelse(hh$HHVEH==1, "HHVEH1", hh$HHVEH)
  hh$HHVEH<-ifelse(hh$HHVEH==2, "HHVEH2", hh$HHVEH)
  hh$HHVEH<-ifelse(hh$HHVEH==3, "HHVEH3", hh$HHVEH)
  hh$HHVEH<-ifelse(hh$HHVEH==4, "HHVEH4", hh$HHVEH)
  
  ##### Household Income
  hh<-hh[hh$INCOME!=99,]
  hh$INCOME<-ifelse(hh$INCOME==1, "INCOME1", hh$INCOME)
  hh$INCOME<-ifelse(hh$INCOME==2, "INCOME2", hh$INCOME)
  hh$INCOME<-ifelse(hh$INCOME==3, "INCOME3", hh$INCOME)
  hh$INCOME<-ifelse(hh$INCOME==4, "INCOME4", hh$INCOME)
  hh$INCOME<-ifelse(hh$INCOME==5, "INCOME5", hh$INCOME)
  hh$INCOME<-ifelse(hh$INCOME==6, "INCOME6", hh$INCOME)
  hh$INCOME<-ifelse(hh$INCOME==7, "INCOME7", hh$INCOME)
  hh$INCOME<-ifelse(hh$INCOME==8, "INCOME8", hh$INCOME)
  
  ##### Household Race
  #Race/Ethnicity: 1-Asian, 2-Black/African American, 3-Hispanic/Latino, 4-Native American/American Indian, 5-Pacific Islander, 6-White
  hh$RACE<-ifelse(hh$RACE==1, "WHITE", hh$RACE)
  hh$RACE<-ifelse(hh$RACE==2, "BLACK", hh$RACE)
  hh$RACE<-ifelse(hh$RACE==3, "NATIVE", hh$RACE)
  hh$RACE<-ifelse(hh$RACE==4, "ASIAN", hh$RACE)
  hh$RACE<-ifelse(hh$RACE==5, "PACIFIC", hh$RACE)
  hh$RACE<-ifelse(hh$RACE==6 | hh$RACE==9, "OTHER", hh$RACE)
  hh$RACE<-ifelse(hh$RACE==7, "MULTI", hh$RACE)
  hh$RACE<-ifelse(hh$HISP==1, "HISPLAT", hh$RACE)
  
  ##### Household Worker
  hh<-hh[hh$HHWRK!=99,]
  #storing hh size >4 proportions for later
  pmtx <- table(data.frame(tractID = tracts, HHWRK = sample(3:6, length(tracts), T)))*0 #creating empty table of proper dimensions
  pwork <- subset(hh, select=c("HHWRK")) #subsetting hhwrkers
  pwork <- table(pwork)[4:7] / sum(table(pwork)[4:7]) #calculating frequencies
  pwork <- as.matrix(pwork) #converting to matrix
  for(i in 1:length(tracts)){pmtx[i,] <- t(pwork)} #filling in proportion for all tracts...
  pwork<-pmtx
  rm(pmtx)
  
  hh$HHWRK<-ifelse(hh$HHWRK>=3, 3, hh$HHWRK)          # Household with 3 or more workers
  hh$HHWRK<-ifelse(hh$HHWRK==0, "HHWRK0", hh$HHWRK)   # Household with 0
  hh$HHWRK<-ifelse(hh$HHWRK==1, "HHWRK1", hh$HHWRK)   # Household with 1 or more workers
  hh$HHWRK<-ifelse(hh$HHWRK==2, "HHWRK2", hh$HHWRK)   # Household with 2 or more workers
  hh$HHWRK<-ifelse(hh$HHWRK==3, "HHWRK3", hh$HHWRK)   # Household with 3 or more workers
  
  ##### tallying individual variables to household  ####
  #merging inds and hh to determine family structure
  indiehh<-merge(ind[, c("SAMPN","PERNO","tractID","RELATE","AGE","GEND","HOURS","WMODE","WORKS","SCHOL")],
                 hh[, c("SAMPN", "tractID", "RESTY", "HHSIZ", "INCOME", "OWN", "HHWRK", "HHVEH", "RACE")],
                 by=c("SAMPN", "tractID"), all.x=T)
  indiehh<-indiehh[which(!is.na(indiehh$HHSIZ)),]
  
  #talling up genders, member types, hh sizes
  tempgend<-as.data.frame.matrix(table(subset(indiehh, select=c("SAMPN","GEND"))))
  temprela<-as.data.frame.matrix(table(subset(indiehh, select=c("SAMPN","RELATE"))))
  temphgend<-subset(indiehh, RELATE=="HEAD", select=c("SAMPN","GEND","AGE"))     #head gender & age
  tempsize<-subset(hh, select=c("SAMPN", "HHSIZ"))
  colnames(temphgend)<-c("SAMPN","HGEND","HAGE")
  
  #Merging tallies into one
  tallies<-merge(temprela,tempgend,by="row.names")
  colnames(tallies)[which(colnames(tallies)=="Row.names")]<-"SAMPN"
  tallies<-merge(tempsize,tallies,by="SAMPN")
  tallies<-merge(tallies,temphgend,by="SAMPN")
  rm(tempsize,temprela,tempgend,temphgend)
  
  #Naming tally types for factors
  tallies$HFAM    <-ifelse(tallies$CHILD>0 | tallies$SPOUSE>0 | tallies$RELATIVE>0 ,"FAM","NONFAM") #htype
  tallies$HHOLDER <-paste(tallies$HGEND,"HEAD",sep="")
  tallies$HMATE   <-ifelse(tallies$SPOUSE>0, "MARRIED","SINGLE")
  tallies$HKIDS   <-ifelse(tallies$CHILD, "KIDS", "NOKIDS")
  tallies$HALONE  <-ifelse(tallies$HHSIZ!="HHSIZ1", "NOTALONE", "ALONE")
  tallies$HNONREL <-ifelse(tallies$NONRELATIVE>0, "NONREL","NONONREL")
  tallies<-tallies[, c("SAMPN", "HFAM", "HHOLDER", "HMATE", "HKIDS", "HALONE", "HNONREL")]
  
  ##### merging tallies and trip durations ####
  hh<-merge(hh,tallies,by="SAMPN")
  ind<-merge(ind,tallies[,c("SAMPN","HFAM")],by="SAMPN")
  
  #Getting travel times from place table
  place$ptract<-paste(sprintf("%02d", place$STATE10),sprintf("%03d", place$COUNTY10),sprintf("%06d", place$TRACT10),sep = "")
  place <- subset(place, place$TPURP == 3, select = c("SAMPN","PERNO","TRPDUR","ptract"))
  place <- place[!duplicated(place),]
  place <- place[!is.na(place),]
  place <- aggregate(TRPDUR ~., data=place, mean)
  
  #adding trip time
  ind<-merge(ind, place, by=c("SAMPN","PERNO"), all.x = T)
  ind$TRPDUR[ind$WORKS=="UNEMP"] <- 0
  ind$TRPDUR[ind$WMODE=="WORKSHOME"] <- 0
  
  #subsetting times
  times <- subset(ind, !is.na(ind$TRPDUR) & !is.na(ind$ptract), select = c("tractID", "TRPDUR"))
  
  #table of possibilities
  times <- times[which(times$tractID %in% tracts),]
  times$hcnty <- substr(times$tractID,0,5)
  times <- subset(times, !is.na(times), select=c("hcnty", "TRPDUR"))
  times <- aggregate(TRPDUR~., data=times, mean)
  
  #formatting times
  ttime <- function(x) {
    x$TTIME <- ifelse(x$TRPDUR == 0, "MINS0", x$TTIME)
    x$TTIME <- ifelse(x$TRPDUR > 0   & x$TRPDUR < 15 , "MINS1TO14", x$TTIME)
    x$TTIME <- ifelse(x$TRPDUR >= 15 & x$TRPDUR < 35, "MINS15TO34", x$TTIME)
  x$TTIME <- ifelse(x$TRPDUR >= 35 & x$TRPDUR < 45, "MINS35TO44", x$TTIME)
  x$TTIME <- ifelse(x$TRPDUR >= 45 & x$TRPDUR < 60, "MINS45TO59", x$TTIME)
  x$TTIME <- ifelse(x$TRPDUR >= 60, "MINS60UP", x$TTIME)

  return(x)
  }
ind$TTIME <- NA
times$TTIME <- NA
ind <- ttime(ind)
times <- ttime(times)

#fillin gmissing travel times at the county level?
for(i in 1:length(times$hcnty)) {
  ind$TTIME[which(is.na(ind$TRPDUR) & substr(ind$tractID,0,5) == times$hcnty[i])] <- times$TTIME[i]
}
ind$TTIME[is.na(ind$TTIME)] <- sample(times$TTIME, length(ind$TTIME[is.na(ind$TTIME)]), T)
rm(times, i, ttime, place)

##### Cleaning up data #####

#reorder
mts_ind <- ind[, c("tractID","RELATE","AGE","GEND","HOURS","WMODE","WORKS","SCHOL","HFAM","TTIME")]
mts_hh <- hh[, c("tractID", "RESTY", "HHSIZ", "INCOME", "OWN", "HHWRK", "HHVEH", "RACE","HFAM","HHOLDER","HMATE","HKIDS","HALONE","HNONREL")]

#### Factor order ####
print("Cleaning up")
mts_hh$HHSIZ <- factor(mts_hh$HHSIZ, levels = c("HHSIZ1","HHSIZ2","HHSIZ3","HHSIZ4"))
mts_hh$HHVEH <- factor(mts_hh$HHVEH, levels = c("HHVEH0","HHVEH1","HHVEH2","HHVEH3","HHVEH4"))
mts_hh$RACE <- factor(mts_hh$RACE, levels = c("ASIAN","BLACK","HISPLAT","MULTI","NATIVE","OTHER","PACIFIC","WHITE") )
mts_hh$RESTY <- factor(mts_hh$RESTY, levels = c("UNITS1","UNITS2TO4","UNITS5TO19","UNITS20ORMORE","UNITSOTHER"))
mts_hh$OWN <- factor(mts_hh$OWN, levels = c("OWN","RENT"))
mts_hh$INCOME <- factor(mts_hh$INCOME, levels = c("INCOME1","INCOME2","INCOME3","INCOME4","INCOME5","INCOME6","INCOME7","INCOME8"))
mts_hh$HHWRK <- factor(mts_hh$HHWRK, levels = c("HHWRK0","HHWRK1","HHWRK2","HHWRK3"))
mts_hh$HFAM <- factor(mts_hh$HFAM, levels = c("FAM","NONFAM"))
mts_hh$HHOLDER <- factor(mts_hh$HHOLDER, levels = c("MALEHEAD","FEMALEHEAD"))
mts_hh$HMATE <- factor(mts_hh$HMATE, levels = c("MARRIED","SINGLE"))
mts_hh$HKIDS <- factor(mts_hh$HKIDS, levels = c("KIDS","NOKIDS"))
mts_hh$HALONE <- factor(mts_hh$HALONE, levels = c("ALONE","NOTALONE"))
mts_hh$HNONREL <- factor(mts_hh$HNONREL, levels = c("NONREL","NONONREL"))

mts_ind$GEND <- factor(mts_ind$GEND, levels = c("MALE","FEMALE"))
mts_ind$AGE <- factor(mts_ind$AGE, levels = c("AGE0TO9","AGE10TO14","AGE15TO19","AGE20TO24","AGE25TO44","AGE45TO54","AGE55TO64","AGE65UP"))
mts_ind$RELATE <- factor(mts_ind$RELATE, levels = c("HEAD","SPOUSE","CHILD","RELATIVE","NONRELATIVE"))
mts_ind$WORKS <- factor(mts_ind$WORKS, levels = c("EMP", "UNEMP"))
mts_ind$HOURS <- factor(mts_ind$HOURS, levels = c("HOURS0","HOURS1TO34","HOURS35UP"))
mts_ind$WMODE <- factor(mts_ind$WMODE, levels = c("DRIVE","RIDEPOOL","TRANSIT","WALK","WORKSHOME","OTHERMODE","NONWORK"))
mts_ind$SCHOL <- factor(mts_ind$SCHOL, levels = c("ENROLLED","NOTENROLLED"))
mts_ind$HFAM <- factor(mts_ind$HFAM, levels = c("FAM","NONFAM"))
mts_ind$TTIME <- factor(mts_ind$TTIME, levels =  c("MINS0","MINS1TO14","MINS15TO34","MINS35TO44","MINS45TO59","MINS60UP"))

#adding dummy column for jobs
occkey<-tables[["pumsocc"]]$OCC
induskey<-tables[["pumsindus"]]$INDUS
mts_ind$OCC <- factor(NA, levels = occkey[-which(duplicated(occkey, fromLast = T))])
mts_ind$INDUS <- factor(NA, levels = induskey[-which(duplicated(induskey, fromLast = T))])

#Merged MTS sample for later
mts <- merge(ind[,c("SAMPN","wblock","RELATE","AGE","GEND","HOURS","WMODE","WORKS","SCHOL","TTIME")],
             hh[,c("SAMPN","tractID","hblock","RESTY","HHSIZ","INCOME","OWN","HHWRK","HHVEH","RACE","HFAM","HHOLDER","HMATE","HKIDS","HALONE","HNONREL")],
             by = "SAMPN")

mts$HHSIZ <- factor(mts$HHSIZ, levels = c("HHSIZ1","HHSIZ2","HHSIZ3","HHSIZ4"))
mts$HHVEH <- factor(mts$HHVEH, levels = c("HHVEH0","HHVEH1","HHVEH2","HHVEH3","HHVEH4"))
mts$RACE <- factor(mts$RACE, levels = c("ASIAN","BLACK","HISPLAT","MULTI","NATIVE","OTHER","PACIFIC","WHITE") )
mts$RESTY <- factor(mts$RESTY, levels = c("UNITS1","UNITS2TO4","UNITS5TO19","UNITS20ORMORE","UNITSOTHER"))
mts$OWN <- factor(mts$OWN, levels = c("OWN","RENT"))
mts$INCOME <- factor(mts$INCOME, levels = c("INCOME1","INCOME2","INCOME3","INCOME4","INCOME5","INCOME6","INCOME7","INCOME8"))
mts$HHWRK <- factor(mts$HHWRK, levels = c("HHWRK0","HHWRK1","HHWRK2","HHWRK3"))
mts$HFAM <- factor(mts$HFAM, levels = c("FAM","NONFAM"))
mts$HHOLDER <- factor(mts$HHOLDER, levels = c("MALEHEAD","FEMALEHEAD"))
mts$HMATE <- factor(mts$HMATE, levels = c("MARRIED","SINGLE"))
mts$HKIDS <- factor(mts$HKIDS, levels = c("KIDS","NOKIDS"))
mts$HALONE <- factor(mts$HALONE, levels = c("ALONE","NOTALONE"))
mts$HNONREL <- factor(mts$HNONREL, levels = c("NONREL","NONONREL"))

mts$GEND <- factor(mts$GEND, levels = c("MALE","FEMALE"))
mts$AGE <- factor(mts$AGE, levels = c("AGE0TO9","AGE10TO14","AGE15TO19","AGE20TO24","AGE25TO44","AGE45TO54","AGE55TO64","AGE65UP"))
mts$RELATE <- factor(mts$RELATE, levels = c("HEAD","SPOUSE","CHILD","RELATIVE","NONRELATIVE"))
mts$WORKS <- factor(mts$WORKS, levels = c("EMP", "UNEMP"))
mts$HOURS <- factor(mts$HOURS, levels = c("HOURS0","HOURS1TO34","HOURS35UP"))
mts$WMODE <- factor(mts$WMODE, levels = c("DRIVE","RIDEPOOL","TRANSIT","WALK","WORKSHOME","OTHER","NONWORK"))
mts$SCHOL <- factor(mts$SCHOL, levels = c("ENROLLED","NOTENROLLED"))
mts$HFAM <- factor(mts$HFAM, levels = c("FAM","NONFAM"))
mts$TTIME <- factor(mts$TTIME, levels = c("MINS0","MINS1TO14","MINS15TO34","MINS35TO44","MINS45TO59","MINS60UP"))
#Removing MTS sample number
#mts_hh<-mts_hh[,-which(names(mts_hh)=="SAMPN")]
#mts_ind<-mts_ind[,-which(names(mts_ind)=="SAMPN")]
mts<-mts[,-which(names(mts)=="SAMPN")]

#### updating to remove erroneous tracts ####
mts_ind <- mts_ind[mts_ind$tractID %in% tracts,]
mts_hh  <- mts_hh[mts_hh$tractID %in% tracts,]
mts <- mts[substr(mts$wblock,0,11) %in% c(NA,tracts),]

pwork <- pwork[tracts,]

#checking if all match
length(unique(mts_ind$tractID))==length(unique(mts_hh$tractID))

#putting mts micro data into microdata file####
microdata[["mts_hh"]] <- mts_hh
microdata[["mts_ind"]] <- mts_ind
microdata[["mts"]] <- mts
rm(hh, ind, indiehh, mts, mts_ind, mts_hh, tallies, dayhours, emptytracts, pwork, induskey, occkey)
# Microdata tables complete
