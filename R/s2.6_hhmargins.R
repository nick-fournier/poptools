# Setting up margins for households
pums_hh <- microdata[["pums_hh"]]
pums_ind <- microdata[["pums_ind"]]

dimsize <- function(tmp){
  tmp <- table(tmp)
  array(tmp, dim = dim(tmp), dimnames = dimnames(tmp))
}

dimshh <- dimsize(pums_hh)
dimsind <- dimsize(pums_ind)
rm(dimsize)

##### HH Income #####
print("Setting up hh income")
con_hhinc<-tables[["hhinc"]]
#removing excess header rows
colnames(con_hhinc)<-con_hhinc[1,]
con_hhinc<-con_hhinc[-1,]
rownames(con_hhinc)<-con_hhinc$Id2

#removing excess columns and column text
con_hhinc<-con_hhinc[,-c(1,2,3,4,which(grepl("Margin of Error",colnames(con_hhinc))))]
colnames(con_hhinc)<-gsub("Estimate; Total: - ", "", colnames(con_hhinc))

#converting char to numeric
con_hhinc[,]<-apply(con_hhinc[,], 2, as.numeric)

#subsetting only by matching mts tracts
con_hhinc<-con_hhinc[tracts,]

# Grouping columns to match MTS income groups
# Income1 <$15,000, # Income2 $15,000-$24,999, Income3 $25,000-$34,999, Income4 $35,000-$49,999,
# Income5 $50,000-$74,999, Income6 $75,000-$99,999, Income7 $100,000-$149,999, Income8 $150,000<
con_hhinc$INCOME1<-rowSums(con_hhinc[,which(grepl("Less than", colnames(con_hhinc))):which(grepl("\\$14,999", colnames(con_hhinc)))])
con_hhinc$INCOME2<-rowSums(con_hhinc[,which(grepl("\\$15,000", colnames(con_hhinc))):which(grepl("\\$24,999", colnames(con_hhinc)))])
con_hhinc$INCOME3<-rowSums(con_hhinc[,which(grepl("\\$25,000", colnames(con_hhinc))):which(grepl("\\$34,999", colnames(con_hhinc)))])
con_hhinc$INCOME4<-rowSums(con_hhinc[,which(grepl("\\$35,000", colnames(con_hhinc))):which(grepl("\\$49,999", colnames(con_hhinc)))])
con_hhinc$INCOME5<-rowSums(con_hhinc[,which(grepl("\\$50,000", colnames(con_hhinc))):which(grepl("\\$74,999", colnames(con_hhinc)))])
con_hhinc$INCOME6<-con_hhinc[,which(grepl("\\$75,000 to \\$99,999", colnames(con_hhinc)))]
con_hhinc$INCOME7<-rowSums(con_hhinc[,which(grepl("\\$100,000", colnames(con_hhinc))):which(grepl("\\$149,999", colnames(con_hhinc)))])
con_hhinc$INCOME8<-rowSums(con_hhinc[,which(grepl("150,000", colnames(con_hhinc))):which(grepl("200,000", colnames(con_hhinc)))])
con_hhinc<-con_hhinc[,which(names(con_hhinc)=="INCOME1"):which(names(con_hhinc)=="INCOME8")]

##### HH dwell & Size ###############
print("Setting up hh size and building type")
con_hhdwellsize<-tables[["hhdwellsize"]]
#removing excess header rows
colnames(con_hhdwellsize)<-con_hhdwellsize[1,]
con_hhdwellsize<-con_hhdwellsize[-1,]
rownames(con_hhdwellsize)<-con_hhdwellsize$Id2

#removing excess columns & column text
con_hhdwellsize<-con_hhdwellsize[,-c(1,2,3, which(grepl("Margin of Error",colnames(con_hhdwellsize))))]
colnames(con_hhdwellsize)<-gsub("Estimate;| |:|,|", "", colnames(con_hhdwellsize))
colnames(con_hhdwellsize)<-gsub("Owneroccupied-", "OWNHHSIZ", colnames(con_hhdwellsize))
colnames(con_hhdwellsize)<-gsub("Renteroccupied-", "RENTHHSIZ", colnames(con_hhdwellsize))
colnames(con_hhdwellsize)<-gsub("personhousehold|detachedorattached|or-more", "", colnames(con_hhdwellsize))
colnames(con_hhdwellsize)<-gsub("--", "UNITS", colnames(con_hhdwellsize))
colnames(con_hhdwellsize)<-gsub("MobilehomeboatRVvanetc.", "OTHER", colnames(con_hhdwellsize))
colnames(con_hhdwellsize)<-toupper(colnames(con_hhdwellsize))

con_hhdwellsize <- con_hhdwellsize[,-which(colnames(con_hhdwellsize) %in% c("TOTAL","OWNEROCCUPIED","RENTEROCCUPIED"))]
con_hhdwellsize <- con_hhdwellsize[,-which(grepl("-",colnames(con_hhdwellsize)))]

#converting char to numeric
con_hhdwellsize[,] <- apply(con_hhdwellsize[,], 2, as.numeric)

#subsetting only by matching mts tracts
con_hhdwellsize<-con_hhdwellsize[tracts,]

#Grouping columns to match mts
for(size in c(dimnames(dimshh)$HHSIZ, "HHSIZ5")){
  restysum <- con_hhdwellsize[,paste("OWN",size,"UNITS20TO49",sep="")] + con_hhdwellsize[,paste("OWN",size,"UNITS50ORMORE",sep="")]
  con_hhdwellsize[,paste("OWN",size,"UNITS20ORMORE",sep="")] <- restysum
  
  restysum <- con_hhdwellsize[,paste("RENT",size,"UNITS20TO49",sep="")] + con_hhdwellsize[,paste("RENT",size,"UNITS50ORMORE",sep="")]
  con_hhdwellsize[,paste("RENT",size,"UNITS20ORMORE",sep="")] <- restysum
}
con_hhdwellsize <- con_hhdwellsize[,-which(grepl("50ORMORE|20TO49", colnames(con_hhdwellsize)))]

#Capping HHSIZ > 4
#own
con_hhdwellsize$OWNHHSIZ4UNITS1        <- con_hhdwellsize$OWNHHSIZ4UNITS1 + con_hhdwellsize$OWNHHSIZ5UNITS1
con_hhdwellsize$OWNHHSIZ4UNITS2TO4     <- con_hhdwellsize$OWNHHSIZ4UNITS2TO4 + con_hhdwellsize$OWNHHSIZ5UNITS2TO4
con_hhdwellsize$OWNHHSIZ4UNITS5TO19    <- con_hhdwellsize$OWNHHSIZ4UNITS5TO19 + con_hhdwellsize$OWNHHSIZ5UNITS5TO19
con_hhdwellsize$OWNHHSIZ4UNITS20ORMORE <- con_hhdwellsize$OWNHHSIZ4UNITS20ORMORE + con_hhdwellsize$OWNHHSIZ5UNITS20ORMORE
con_hhdwellsize$OWNHHSIZ4UNITSOTHER    <- con_hhdwellsize$OWNHHSIZ4UNITSOTHER + con_hhdwellsize$OWNHHSIZ5UNITSOTHER
#rent
con_hhdwellsize$RENTHHSIZ4UNITS1        <- con_hhdwellsize$RENTHHSIZ4UNITS1 + con_hhdwellsize$RENTHHSIZ5UNITS1
con_hhdwellsize$RENTHHSIZ4UNITS2TO4     <- con_hhdwellsize$RENTHHSIZ4UNITS2TO4 + con_hhdwellsize$RENTHHSIZ5UNITS2TO4
con_hhdwellsize$RENTHHSIZ4UNITS5TO19    <- con_hhdwellsize$RENTHHSIZ4UNITS5TO19 + con_hhdwellsize$RENTHHSIZ5UNITS5TO19
con_hhdwellsize$RENTHHSIZ4UNITS20ORMORE <- con_hhdwellsize$RENTHHSIZ4UNITS20ORMORE + con_hhdwellsize$RENTHHSIZ5UNITS20ORMORE
con_hhdwellsize$RENTHHSIZ4UNITSOTHER    <- con_hhdwellsize$RENTHHSIZ4UNITSOTHER + con_hhdwellsize$RENTHHSIZ5UNITSOTHER
#deleting old columns and reordering them
con_hhdwellsize <- con_hhdwellsize[,-which(grepl("HHSIZ5", colnames(con_hhdwellsize)))]

# Transform into an 3D-array
names <- c(list(tracts),dimnames(dimshh)[c("OWN","HHSIZ","RESTY")])
con_hhdwellsize2 <- array(NA, dim=c(length(tracts),length(dimnames(dimshh)$OWN),
                                    length(dimnames(dimshh)$HHSIZ),
                                    length(dimnames(dimshh)$RESTY)), dimnames = names)

for(tract in rownames(con_hhdwellsize)){
  for (own in dimnames(con_hhdwellsize2)$OWN){
    for (size in dimnames(con_hhdwellsize2)$HHSIZ){
      for (resty in dimnames(con_hhdwellsize2)$RESTY){
        con_hhdwellsize2[tract,own,size,resty] <- con_hhdwellsize[tract,paste(own,size,resty,sep="")]
}}}}
rm(tract, own, size, resty, names, con_hhdwellsize)

con_hhowndwellsize2 <- con_hhdwellsize2
con_hhdwellsize2 <- colSums(aperm(con_hhdwellsize2, c(2,1,3,4)), dims=1)

##### HH vehicles & Size ###############
print("Setting up hh vehicles and size")
con_hhvehsize<-tables[["hhvehsize"]]
#removing excess header rows
colnames(con_hhvehsize)<-con_hhvehsize[1,]
con_hhvehsize<-con_hhvehsize[-1,]
rownames(con_hhvehsize)<-con_hhvehsize$Id2

#removing excess columns & column text
con_hhvehsize<-con_hhvehsize[,-which(grepl("Margin of Error",colnames(con_hhvehsize)))]
con_hhvehsize<-con_hhvehsize[,-c(1:9)]
colnames(con_hhvehsize)<-gsub(" ", "", colnames(con_hhvehsize))
colnames(con_hhvehsize)<-gsub(":", "", colnames(con_hhvehsize))
colnames(con_hhvehsize)<-gsub("Estimate;", "", colnames(con_hhvehsize))
colnames(con_hhvehsize)<-gsub("Total-", "", colnames(con_hhvehsize))
colnames(con_hhvehsize)<-gsub("1-personhousehold-", "HHSIZ1", colnames(con_hhvehsize))
colnames(con_hhvehsize)<-gsub("2-personhousehold-", "HHSIZ2", colnames(con_hhvehsize))
colnames(con_hhvehsize)<-gsub("3-personhousehold-", "HHSIZ3", colnames(con_hhvehsize))
colnames(con_hhvehsize)<-gsub("4-or-more-personhousehold-", "HHSIZ4", colnames(con_hhvehsize))
colnames(con_hhvehsize)<-gsub("Novehicleavailable", "HHVEH0", colnames(con_hhvehsize))
colnames(con_hhvehsize)<-gsub("1vehicleavailable", "HHVEH1", colnames(con_hhvehsize))
colnames(con_hhvehsize)<-gsub("2vehiclesavailable", "HHVEH2", colnames(con_hhvehsize))
colnames(con_hhvehsize)<-gsub("3vehiclesavailable", "HHVEH3", colnames(con_hhvehsize))
colnames(con_hhvehsize)<-gsub("4ormorevehiclesavailable", "HHVEH4", colnames(con_hhvehsize))
con_hhvehsize<-con_hhvehsize[,-which(grepl("personhousehold",colnames(con_hhvehsize)))]

#converting char to numeric
con_hhvehsize[,] <- apply(con_hhvehsize[,], 2, as.numeric)

#subsetting only by matching mts tracts
con_hhvehsize<-con_hhvehsize[tracts,]

# Transform into an 3D-array
names <- c(list(tracts),dimnames(dimshh)[c("HHSIZ","HHVEH")])
con_hhvehsize2 <- array(NA, dim=c(length(tracts),4,5), dimnames = names)

for(tract in rownames(con_hhvehsize)){
  for (size in dimnames(con_hhvehsize2)$HHSIZ){
    for (veh in dimnames(con_hhvehsize2)$HHVEH){
      con_hhvehsize2[tract,size,veh] <- con_hhvehsize[tract,paste(size,veh,sep="")]
    }}}
rm(tract,size,veh,names, con_hhvehsize)

##### HH Race #####
print("Setting up hh race")

rtab <- c("A","B","C","D","E","F","G","H")
race <- c("WHITEHISP", "BLACK","NATIVE","ASIAN","PACIFIC","OTHER","MULTI","WHITE")
con_hhrace <- foreach(i = 1:length(race)) %do% {
  hhrace<-tables[[paste("hhrace", rtab[i], sep="")]]
  #removing excess header rows
  colnames(hhrace)<-hhrace[1,]
  hhrace<-hhrace[-1,]
  rownames(hhrace)<-hhrace$Id2
  
  #removing excess columns and column text
  hhrace<-hhrace[,-c(1,2,3,4,which(grepl("Margin of Error",colnames(hhrace))))]
  hhrace<-hhrace[,which(grepl("households:\\b",colnames(hhrace)))]
  colnames(hhrace)<-gsub("Estimate; ", "", colnames(hhrace))
  colnames(hhrace)<-gsub(" ", "", colnames(hhrace))
  colnames(hhrace)<-gsub(":", "", colnames(hhrace))
  colnames(hhrace)<-gsub("ilyhouseholds", "", colnames(hhrace))
  colnames(hhrace)<-toupper(colnames(hhrace))
  colnames(hhrace)<-paste(colnames(hhrace), race[i], sep="")
  
  return(hhrace)
}
con_hhrace <- do.call(cbind, con_hhrace)
rm(i, race, rtab, hhrace)
#converting char to numeric
con_hhrace[,]<-apply(con_hhrace[,], 2, as.numeric)

#subsetting only by matching mts tracts
con_hhrace<-con_hhrace[tracts,]
#subtracting white from white+hispanic
con_hhrace$FAMWHITEHISP <- con_hhrace$FAMWHITEHISP - con_hhrace$FAMWHITE
con_hhrace$NONFAMWHITEHISP <- con_hhrace$NONFAMWHITEHISP - con_hhrace$NONFAMWHITE
#renaming
colnames(con_hhrace)[which(colnames(con_hhrace) == c("FAMWHITEHISP","NONFAMWHITEHISP"))] <- c("FAMHISPLAT", "NONFAMHISPLAT")

# Transform into an 3D-array
names <- c(list(tracts),c(dimnames(dimshh)["HFAM"], dimnames(dimshh)["RACE"]))
con_hhrace2 <- array(NA, dim=c(length(tracts),length(dimnames(dimshh)$HFAM),length(dimnames(dimshh)$RACE)), dimnames = names)

for(tract in rownames(con_hhrace)){
  for (fam in dimnames(con_hhrace2)$HFAM){
    for (race in dimnames(con_hhrace2)$RACE){
      con_hhrace2[tract,fam,race] <- con_hhrace[tract,paste(fam,race,sep="")]
    }}}
rm(tract, fam, race, names, con_hhrace)


#Collapsing along HFAM dimensions
con_hhracecomp <- aaply(con_hhrace2, c(1,3), sum)


##### check the totals per zone ####
table(rowSums(con_hhinc)==rowSums(con_hhracecomp))
table(rowSums(con_hhinc)==rowSums(con_hhvehsize2))
table(rowSums(con_hhinc)==rowSums(con_hhdwellsize2))

#check totals
sum(con_hhinc)
sum(con_hhdwellsize2)
sum(con_hhvehsize2)
sum(con_hhracecomp)

#### Putting into neat list form ####
margins[["hhinc"]] <- con_hhinc
margins[["hhracecomp"]] <- con_hhracecomp
margins[["hhdwellsize"]] <- con_hhdwellsize2
margins[["hhvehsize"]] <- con_hhvehsize2
margins[["hhdwell"]] <- apply(margins[['hhdwellsize']], c(1,3), sum)
margins[["hhveh"]] <- apply(margins[['hhvehsize']], c(1,3), sum)

for(i in names(margins)){
  if(any(margins[[i]]<0)){
    print(i)
  }
}

#cleanup ####

rm(list=setdiff(ls(),c("clock","runtime","csvwriteout",
                       "tables","microdata","tracts","margins","settings")))

print("HH Margins done")
