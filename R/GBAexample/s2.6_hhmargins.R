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
marg_hhinc<-tables[["hhinc"]]
#removing excess header rows
colnames(marg_hhinc)<-marg_hhinc[1,]
marg_hhinc<-marg_hhinc[-1,]
rownames(marg_hhinc)<-marg_hhinc$Id2

#removing excess columns and column text
marg_hhinc<-marg_hhinc[,-c(1,2,3,4,which(grepl("Margin of Error",colnames(marg_hhinc))))]
colnames(marg_hhinc)<-gsub("Estimate; Total: - ", "", colnames(marg_hhinc))

#converting char to numeric
marg_hhinc[,]<-apply(marg_hhinc[,], 2, as.numeric)

#subsetting only by matching mts tracts
marg_hhinc<-marg_hhinc[tracts,]

# Grouping columns to match MTS income groups
# Income1 <$15,000, # Income2 $15,000-$24,999, Income3 $25,000-$34,999, Income4 $35,000-$49,999,
# Income5 $50,000-$74,999, Income6 $75,000-$99,999, Income7 $100,000-$149,999, Income8 $150,000<
marg_hhinc$INCOME1<-rowSums(marg_hhinc[,which(grepl("Less than", colnames(marg_hhinc))):which(grepl("\\$14,999", colnames(marg_hhinc)))])
marg_hhinc$INCOME2<-rowSums(marg_hhinc[,which(grepl("\\$15,000", colnames(marg_hhinc))):which(grepl("\\$24,999", colnames(marg_hhinc)))])
marg_hhinc$INCOME3<-rowSums(marg_hhinc[,which(grepl("\\$25,000", colnames(marg_hhinc))):which(grepl("\\$34,999", colnames(marg_hhinc)))])
marg_hhinc$INCOME4<-rowSums(marg_hhinc[,which(grepl("\\$35,000", colnames(marg_hhinc))):which(grepl("\\$49,999", colnames(marg_hhinc)))])
marg_hhinc$INCOME5<-rowSums(marg_hhinc[,which(grepl("\\$50,000", colnames(marg_hhinc))):which(grepl("\\$74,999", colnames(marg_hhinc)))])
marg_hhinc$INCOME6<-marg_hhinc[,which(grepl("\\$75,000 to \\$99,999", colnames(marg_hhinc)))]
marg_hhinc$INCOME7<-rowSums(marg_hhinc[,which(grepl("\\$100,000", colnames(marg_hhinc))):which(grepl("\\$149,999", colnames(marg_hhinc)))])
marg_hhinc$INCOME8<-rowSums(marg_hhinc[,which(grepl("150,000", colnames(marg_hhinc))):which(grepl("200,000", colnames(marg_hhinc)))])
marg_hhinc<-marg_hhinc[,which(names(marg_hhinc)=="INCOME1"):which(names(marg_hhinc)=="INCOME8")]

##### HH dwell & Size ###############
print("Setting up hh size and building type")
marg_hhdwellsize<-tables[["hhdwellsize"]]
#removing excess header rows
colnames(marg_hhdwellsize)<-marg_hhdwellsize[1,]
marg_hhdwellsize<-marg_hhdwellsize[-1,]
rownames(marg_hhdwellsize)<-marg_hhdwellsize$Id2

#removing excess columns & column text
marg_hhdwellsize<-marg_hhdwellsize[,-c(1,2,3, which(grepl("Margin of Error",colnames(marg_hhdwellsize))))]
colnames(marg_hhdwellsize)<-gsub("Estimate;| |:|,|", "", colnames(marg_hhdwellsize))
colnames(marg_hhdwellsize)<-gsub("Owneroccupied-", "OWN_HHSIZ", colnames(marg_hhdwellsize))
colnames(marg_hhdwellsize)<-gsub("Renteroccupied-", "RENT_HHSIZ", colnames(marg_hhdwellsize))
colnames(marg_hhdwellsize)<-gsub("personhousehold|detachedorattached|or-more", "", colnames(marg_hhdwellsize))
colnames(marg_hhdwellsize)<-gsub("--", "_UNITS", colnames(marg_hhdwellsize))
colnames(marg_hhdwellsize)<-gsub("MobilehomeboatRVvanetc.", "OTHER", colnames(marg_hhdwellsize))
colnames(marg_hhdwellsize)<-toupper(colnames(marg_hhdwellsize))

marg_hhdwellsize <- marg_hhdwellsize[,-which(colnames(marg_hhdwellsize) %in% c("TOTAL","OWNEROCCUPIED","RENTEROCCUPIED"))]
marg_hhdwellsize <- marg_hhdwellsize[,-which(grepl("-",colnames(marg_hhdwellsize)))]

#converting char to numeric
marg_hhdwellsize[,] <- apply(marg_hhdwellsize[,], 2, as.numeric)

#subsetting only by matching mts tracts
marg_hhdwellsize<-marg_hhdwellsize[tracts,]

#Grouping columns to match mts
for(size in c(dimnames(dimshh)$HHSIZ, "HHSIZ5")){
  restysum <- marg_hhdwellsize[,paste("OWN",size,"UNITS20TO49",sep="_")] + marg_hhdwellsize[,paste("OWN",size,"UNITS50ORMORE",sep="_")]
  marg_hhdwellsize[,paste("OWN",size,"UNITS20ORMORE",sep="_")] <- restysum
  
  restysum <- marg_hhdwellsize[,paste("RENT",size,"UNITS20TO49",sep="_")] + marg_hhdwellsize[,paste("RENT",size,"UNITS50ORMORE",sep="_")]
  marg_hhdwellsize[,paste("RENT",size,"UNITS20ORMORE",sep="_")] <- restysum
}
marg_hhdwellsize <- marg_hhdwellsize[,-which(grepl("50ORMORE|20TO49", colnames(marg_hhdwellsize)))]

#Capping HHSIZ > 4
#own
marg_hhdwellsize$OWN_HHSIZ4_UNITS1        <- marg_hhdwellsize$OWN_HHSIZ4_UNITS1 + marg_hhdwellsize$OWN_HHSIZ5_UNITS1
marg_hhdwellsize$OWN_HHSIZ4_UNITS2TO4     <- marg_hhdwellsize$OWN_HHSIZ4_UNITS2TO4 + marg_hhdwellsize$OWN_HHSIZ5_UNITS2TO4
marg_hhdwellsize$OWN_HHSIZ4_UNITS5TO19    <- marg_hhdwellsize$OWN_HHSIZ4_UNITS5TO19 + marg_hhdwellsize$OWN_HHSIZ5_UNITS5TO19
marg_hhdwellsize$OWN_HHSIZ4_UNITS20ORMORE <- marg_hhdwellsize$OWN_HHSIZ4_UNITS20ORMORE + marg_hhdwellsize$OWN_HHSIZ5_UNITS20ORMORE
marg_hhdwellsize$OWN_HHSIZ4_UNITSOTHER    <- marg_hhdwellsize$OWN_HHSIZ4_UNITSOTHER + marg_hhdwellsize$OWN_HHSIZ5_UNITSOTHER
#rent
marg_hhdwellsize$RENT_HHSIZ4_UNITS1        <- marg_hhdwellsize$RENT_HHSIZ4_UNITS1 + marg_hhdwellsize$RENT_HHSIZ5_UNITS1
marg_hhdwellsize$RENT_HHSIZ4_UNITS2TO4     <- marg_hhdwellsize$RENT_HHSIZ4_UNITS2TO4 + marg_hhdwellsize$RENT_HHSIZ5_UNITS2TO4
marg_hhdwellsize$RENT_HHSIZ4_UNITS5TO19    <- marg_hhdwellsize$RENT_HHSIZ4_UNITS5TO19 + marg_hhdwellsize$RENT_HHSIZ5_UNITS5TO19
marg_hhdwellsize$RENT_HHSIZ4_UNITS20ORMORE <- marg_hhdwellsize$RENT_HHSIZ4_UNITS20ORMORE + marg_hhdwellsize$RENT_HHSIZ5_UNITS20ORMORE
marg_hhdwellsize$RENT_HHSIZ4_UNITSOTHER    <- marg_hhdwellsize$RENT_HHSIZ4_UNITSOTHER + marg_hhdwellsize$RENT_HHSIZ5_UNITSOTHER
#deleting old columns and reordering them
marg_hhdwellsize <- marg_hhdwellsize[,-which(grepl("HHSIZ5", colnames(marg_hhdwellsize)))]

# Transform into an 3D-array
names <- c(list(tracts),dimnames(dimshh)[c("OWN","HHSIZ","RESTY")])
marg_hhdwellsize2 <- array(NA, dim=c(length(tracts),length(dimnames(dimshh)$OWN),
                                    length(dimnames(dimshh)$HHSIZ),
                                    length(dimnames(dimshh)$RESTY)), dimnames = names)

for(tract in rownames(marg_hhdwellsize)){
  for (own in dimnames(marg_hhdwellsize2)$OWN){
    for (size in dimnames(marg_hhdwellsize2)$HHSIZ){
      for (resty in dimnames(marg_hhdwellsize2)$RESTY){
        marg_hhdwellsize2[tract,own,size,resty] <- marg_hhdwellsize[tract,paste(own,size,resty,sep="_")]
}}}}
rm(tract, own, size, resty, names)

marg_hhowndwellsize2 <- marg_hhdwellsize2
marg_hhdwellsize2 <- colSums(aperm(marg_hhdwellsize2, c(2,1,3,4)), dims=1)

##### HH vehicles & Size ###############
print("Setting up hh vehicles and size")
marg_hhvehsize<-tables[["hhvehsize"]]
#removing excess header rows
colnames(marg_hhvehsize)<-marg_hhvehsize[1,]
marg_hhvehsize<-marg_hhvehsize[-1,]
rownames(marg_hhvehsize)<-marg_hhvehsize$Id2

#removing excess columns & column text
marg_hhvehsize<-marg_hhvehsize[,-which(grepl("Margin of Error",colnames(marg_hhvehsize)))]
marg_hhvehsize<-marg_hhvehsize[,-c(1:9)]
colnames(marg_hhvehsize)<-gsub(" ", "", colnames(marg_hhvehsize))
colnames(marg_hhvehsize)<-gsub(":", "", colnames(marg_hhvehsize))
colnames(marg_hhvehsize)<-gsub("Estimate;", "", colnames(marg_hhvehsize))
colnames(marg_hhvehsize)<-gsub("Total-", "", colnames(marg_hhvehsize))
colnames(marg_hhvehsize)<-gsub("1-personhousehold-", "HHSIZ1_", colnames(marg_hhvehsize))
colnames(marg_hhvehsize)<-gsub("2-personhousehold-", "HHSIZ2_", colnames(marg_hhvehsize))
colnames(marg_hhvehsize)<-gsub("3-personhousehold-", "HHSIZ3_", colnames(marg_hhvehsize))
colnames(marg_hhvehsize)<-gsub("4-or-more-personhousehold-", "HHSIZ4_", colnames(marg_hhvehsize))
colnames(marg_hhvehsize)<-gsub("Novehicleavailable", "HHVEH0", colnames(marg_hhvehsize))
colnames(marg_hhvehsize)<-gsub("1vehicleavailable", "HHVEH1", colnames(marg_hhvehsize))
colnames(marg_hhvehsize)<-gsub("2vehiclesavailable", "HHVEH2", colnames(marg_hhvehsize))
colnames(marg_hhvehsize)<-gsub("3vehiclesavailable", "HHVEH3", colnames(marg_hhvehsize))
colnames(marg_hhvehsize)<-gsub("4ormorevehiclesavailable", "HHVEH4", colnames(marg_hhvehsize))
marg_hhvehsize<-marg_hhvehsize[,-which(grepl("personhousehold",colnames(marg_hhvehsize)))]

#converting char to numeric
marg_hhvehsize[,] <- apply(marg_hhvehsize[,], 2, as.numeric)

#subsetting only by matching mts tracts
marg_hhvehsize<-marg_hhvehsize[tracts,]

# Transform into an 3D-array
names <- c(list(tracts),dimnames(dimshh)[c("HHSIZ","HHVEH")])
marg_hhvehsize2 <- array(NA, dim=c(length(tracts),4,5), dimnames = names)

for(tract in rownames(marg_hhvehsize)){
  for (size in dimnames(marg_hhvehsize2)$HHSIZ){
    for (veh in dimnames(marg_hhvehsize2)$HHVEH){
      marg_hhvehsize2[tract,size,veh] <- marg_hhvehsize[tract,paste(size,veh,sep="_")]
    }}}
rm(tract,size,veh,names)

##### check the totals per zone ####
table(rowSums(marg_hhinc)==rowSums(marg_hhvehsize2))
table(rowSums(marg_hhinc)==rowSums(marg_hhdwellsize2))

#check totals
sum(marg_hhinc)
sum(marg_hhdwellsize2)
sum(marg_hhvehsize2)

#### Putting into neat list form ####
margins[["hhinc"]] <- marg_hhinc
margins[["hhdwellsize"]] <- marg_hhdwellsize2
margins[["hhvehsize"]] <- marg_hhvehsize2
margins[["hhdwell"]] <- apply(margins[['hhdwellsize']], c(1,3), sum)
margins[["hhveh"]] <- apply(margins[['hhvehsize']], c(1,3), sum)

for(i in names(margins)){
  if(any(margins[[i]]<0)){
    print(i)
  }
}



marginal_hhinc = as.data.table(margins[["hhinc"]],keep.rownames = T)
marginal_hhvehsize = as.data.table(margins[["hhvehsize"]],keep.rownames = T)
marginal_hhdwell = as.data.table(margins[["hhdwell"]],keep.rownames = T)

setnames(marginal_hhdwell, "rn","tract")
setnames(marginal_hhinc, "rn","tract")
setnames(marginal_hhvehsize, "V3","tract")
fwrite(marginal_hhdwell, file = "./data/marginal_hh_dwellingtype.csv")
fwrite(marginal_hhinc, file = "./data/marginal_hh_income.csv")
fwrite(marginal_hhvehsize, file = "./data/marginal_hh_vehiclesandsize.csv")

names(dimnames(marg_hhvehsize2)) <- c("tracts","HHSIZ","HHVEH")

marginals_hh <- list("inc"=marg_hhinc, "dwell" = marginal_hhdwell, "vehsize" = marg_hhvehsize2)
save(marginals_hh, file = "./data/marginal_hh.RData")


print("HH Margins done")
