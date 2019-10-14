# Setting up margins for individuals
pums_hh <- microdata[["pums_hh"]]
pums_ind <- microdata[["pums_ind"]]

dimsize <- function(tmp){
  tmp <- table(tmp)
  array(tmp, dim = dim(tmp), dimnames = dimnames(tmp))
}

dimshh <- dimsize(pums_hh)
dimsind <- dimsize(pums_ind)
rm(dimsize)

##### Ind Age & Sex ################
print("Setting up age & sex tables")
con_indagesex<-tables[["indagesex"]]
#removing excess header rows
colnames(con_indagesex)<-con_indagesex[1,]
con_indagesex<-con_indagesex[-1,]
rownames(con_indagesex)<-con_indagesex$Id2

#removing excess columns and column text
con_indagesex<-con_indagesex[,-c(4, which(grepl("Id|Id2|Geography|Margin",colnames(con_indagesex))))]
colnames(con_indagesex)<-gsub("Estimate;| |:|-|,", "", colnames(con_indagesex))
# con_indagesex<-con_indagesex[,-c(1,2,3,which(grepl("Percent",colnames(con_indagesex))))]
# con_indagesex<-con_indagesex[,which(grepl("SEX AND AGE",colnames(con_indagesex)))]
# con_indagesex<-con_indagesex[,which(grepl("Male",colnames(con_indagesex))|grepl("Female",colnames(con_indagesex)))]
# con_indagesex<-con_indagesex[,-which(grepl("Median",colnames(con_indagesex)))]
#saving this for proportioning later
# agesex16<-con_indagesex[,which(grepl("16 years and over",colnames(con_indagesex)))]
# agesex16[,]<-apply(agesex16[,], 2, as.numeric)
# colnames(agesex16) <- c("MALE", "FEMALE")
# agesex16<-agesex16[tracts,]
#
# con_indagesex<-con_indagesex[,-which(grepl("yearsandover",colnames(con_indagesex)) & !grepl("85",colnames(con_indagesex)))]
# con_indagesex<-con_indagesex[,which(grepl("population -",colnames(con_indagesex)))]
# colnames(con_indagesex)<-gsub("Number; SEX AND AGE| |-|population|years", "", colnames(con_indagesex))
# colnames(con_indagesex)<-gsub("Under", "0TO", colnames(con_indagesex))
colnames(con_indagesex)<-gsub("Under", "0TO", colnames(con_indagesex))
colnames(con_indagesex)<-gsub("years", "", colnames(con_indagesex))
colnames(con_indagesex)<-gsub("and", "TO", colnames(con_indagesex))
colnames(con_indagesex)<-gsub("TOover", "ANDOVER", colnames(con_indagesex))
colnames(con_indagesex)<-toupper(colnames(con_indagesex))

#subsetting only by matching mts tracts
con_indagesex<-con_indagesex[tracts,]

#converting to numeric
con_indagesex[,]<-apply(con_indagesex[,], 2, as.numeric)

#grouping by ages
for(i in dimnames(dimsind)$GEND) {
  con_indagesex[,paste(i,"AGE0TO9",sep="")]<-rowSums(con_indagesex[, paste(i, c("0TO5","5TO9"),sep="")])
  con_indagesex[,paste(i,"AGE10TO14",sep="")]<-con_indagesex[, paste(i, c("10TO14"),sep="")]
  con_indagesex[,paste(i,"AGE15TO19",sep="")]<-rowSums(con_indagesex[, paste(i, c("15TO17","18TO19"),sep="")])
  #con_indagesex[,paste(i,"AGE15TO19",sep="")]<-con_indagesex[, paste(i, c("15TO19"),sep="")]
  #con_indagesex[,paste(i,"AGE20TO24",sep="")]<-con_indagesex[, paste(i, c("20TO24"),sep="")]
  con_indagesex[,paste(i,"AGE20TO24",sep="")]<-rowSums(con_indagesex[, paste(i, c("20","21","22TO24"),sep="")])
  con_indagesex[,paste(i,"AGE25TO44",sep="")]<-rowSums(con_indagesex[, paste(i, c("25TO29","30TO34","35TO39","40TO44"),sep="")])
  con_indagesex[,paste(i,"AGE45TO54",sep="")]<-rowSums(con_indagesex[, paste(i, c("45TO49","50TO54"),sep="")])
  con_indagesex[,paste(i,"AGE55TO64",sep="")]<-rowSums(con_indagesex[, paste(i, c("55TO59","60TO61","62TO64"),sep="")])
  con_indagesex[,paste(i,"AGE65UP",sep="")]<-rowSums(con_indagesex[, paste(i, c("65TO66","67TO69","70TO74","75TO79","80TO84","85ANDOVER"),sep="")])
}

con_indagesex<-con_indagesex[,which(grepl("AGE", colnames(con_indagesex)))]

# Transform age-sex into an 3D-array
names <- c(list(tracts),dimnames(dimsind)[c("GEND","AGE")])
con_indagesex2 <- array(NA, dim=c(length(tracts),length(dimnames(dimsind)$GEND),length(dimnames(dimsind)$AGE)), dimnames = names)

for(tract in rownames(con_indagesex)){
  for (sex in dimnames(con_indagesex2)$GEND){
    for (age in dimnames(con_indagesex2)$AGE){
      con_indagesex2[tract,sex,age] <- con_indagesex[tract,paste(sex,age,sep="")]
    }}}

##### Ind Relationship ################
print("Setting relationship tables")
#removing excess header rows
con_indrelate<-tables[["indrelate"]]
colnames(con_indrelate)<-con_indrelate[1,]
con_indrelate<-con_indrelate[-1,]
rownames(con_indrelate)<-con_indrelate$Id2
#subsetting only by matching mts tracts
con_indrelate<-con_indrelate[tracts,]

#removing excess columns and column text
con_indrelate<-con_indrelate[,-c(which(grepl("Id|Id2|Geography|Margin",colnames(con_indrelate))))]
colnames(con_indrelate)<-gsub("Estimate;| |:|-|,", "", colnames(con_indrelate))
#con_indrelate<-con_indrelate[,-c(1:3)]
colnames(con_indrelate)<-toupper(colnames(con_indrelate))
#colnames(con_indrelate)<-gsub(" ", "", colnames(con_indrelate))
#colnames(con_indrelate)<-gsub(":", "", colnames(con_indrelate))
#colnames(con_indrelate)<-gsub("-", "", colnames(con_indrelate))
colnames(con_indrelate)<-gsub("INHOUSEHOLDSIN", "", colnames(con_indrelate))
con_indrelate<-con_indrelate[ , -which(colnames(con_indrelate) %in%
                                         c("TOTAL","INHOUSEHOLDS","FAMILYHOUSEHOLDS","NONFAMILYHOUSEHOLDS",
                                           "FAMILYHOUSEHOLDSHOUSEHOLDERMALE","FAMILYHOUSEHOLDSHOUSEHOLDERFEMALE",
                                           "NONFAMILYHOUSEHOLDSHOUSEHOLDERMALE","NONFAMILYHOUSEHOLDSHOUSEHOLDERFEMALE"))]
con_indrelate<-con_indrelate[ , -which(grepl("BIOLOGICAL|ADOPTED|FOSTER|STEPCHILD", colnames(con_indrelate)))]
con_indrelate<-con_indrelate[ , -which(grepl("ROOMER|HOUSEMATE|UNMARRIED|FOSTER|OTHERNONRELATIVE|LIVINGALONE", colnames(con_indrelate)))]
colnames(con_indrelate)<-gsub("NONFAMILYHOUSEHOLDS|FAMILYHOUSEHOLDS", "", colnames(con_indrelate))

#converting to numeric and removing bad text
con_indrelate[,]<-data.frame(apply(con_indrelate[,], 2, function(y) as.numeric(gsub("*\\(.*?\\)", "", y))))

#Grouping
con_indrelate$HEAD <- rowSums(con_indrelate[ , which(grepl("HOUSEHOLDER",colnames(con_indrelate)))])
con_indrelate$NONRELATIVE <- rowSums(con_indrelate[ , grepl("NONRELATIVES|INGROUPQUARTERS",colnames(con_indrelate))])
con_indrelate$CHILD <- rowSums(con_indrelate[ , which(grepl("CHILD",colnames(con_indrelate)))])
con_indrelate$RELATIVE <- rowSums(con_indrelate[ , which(grepl("SISTER|PARENT|SON|OTHERRELATIVES",colnames(con_indrelate)))])

#consolidating
con_indrelate <- con_indrelate[ , c("HEAD","SPOUSE","CHILD","RELATIVE","NONRELATIVE")]

##### Ind occupation & Industry #####
print("Setting up industry & occupation tables")
con_indindusocc<-tables[["indindusocc"]]
#removing excess header rows
colnames(con_indindusocc)<-con_indindusocc[1,]
con_indindusocc<-con_indindusocc[-1,]
rownames(con_indindusocc)<-con_indindusocc$Id2
#cleanup
con_indindusocc<-con_indindusocc[,-c(4, which(grepl("Id|Id2|Geography|Margin",colnames(con_indindusocc))))]
colnames(con_indindusocc)<-gsub("Estimate; Total| |:|-|,", "", colnames(con_indindusocc))
#Major occ categories
colnames(con_indindusocc)<-gsub("managementbusinessscienceandartsoccupations", "MGMTBIZSCIART", colnames(con_indindusocc), ignore.case = T)
colnames(con_indindusocc)<-gsub("serviceoccupations", "SERVICE", colnames(con_indindusocc), ignore.case = T)
colnames(con_indindusocc)<-gsub("salesandofficeoccupations", "SALESOFFICEADMIN", colnames(con_indindusocc), ignore.case = T)
colnames(con_indindusocc)<-gsub("naturalresourcesconstructionandmaintenanceoccupations", "NATRESCONSTMAINT", colnames(con_indindusocc), ignore.case = T)
colnames(con_indindusocc)<-gsub("productiontransportationandmaterialmovingoccupations", "PRODTRANS", colnames(con_indindusocc), ignore.case = T)
#major indus categories
colnames(con_indindusocc)<-gsub("Agricultureforestryfishingandhuntingandmining", "NATRES", colnames(con_indindusocc), ignore.case = T)
colnames(con_indindusocc)<-gsub("Construction", "CONST", colnames(con_indindusocc), ignore.case = T)
colnames(con_indindusocc)<-gsub("Manufacturing", "MFG", colnames(con_indindusocc), ignore.case = T)
colnames(con_indindusocc)<-gsub("Wholesaletrade", "WHLTRADE", colnames(con_indindusocc), ignore.case = T)
colnames(con_indindusocc)<-gsub("Retailtrade", "RETTRADE", colnames(con_indindusocc), ignore.case = T)
colnames(con_indindusocc)<-gsub("Transportationandwarehousingandutilities", "TRANSUTILS", colnames(con_indindusocc), ignore.case = T)
colnames(con_indindusocc)<-gsub("Information", "INFO", colnames(con_indindusocc), ignore.case = T)
colnames(con_indindusocc)<-gsub("Financeandinsuranceandrealestateandrentalandleasing", "FINESTATE", colnames(con_indindusocc), ignore.case = T)
colnames(con_indindusocc)<-gsub("Professionalscientificandmanagementandadministrativeandwastemanagementservices", "PROFSCIMGMT", colnames(con_indindusocc), ignore.case = T)
colnames(con_indindusocc)<-gsub("Educationalservicesandhealthcareandsocialassistance", "EDUCSOC", colnames(con_indindusocc), ignore.case = T)
colnames(con_indindusocc)<-gsub("Artsentertainmentandrecreationandaccommodationandfoodservices", "ARTSACCOM", colnames(con_indindusocc), ignore.case = T)
colnames(con_indindusocc)<-gsub("Otherservicesexceptpublicadministration", "OTHER", colnames(con_indindusocc), ignore.case = T)
colnames(con_indindusocc)<-gsub("Publicadministration", "PUBADMIN", colnames(con_indindusocc), ignore.case = T)
#subsetting only by matching mts tracts
con_indindusocc<-con_indindusocc[tracts,]
#converting to numeric
con_indindusocc[,]<-apply(con_indindusocc[,], 2, as.numeric)

# #Combining duplicates
# dupes <- data.table(table(gsub('[[:digit:]]+', '', colnames(con_indindusocc))))
# dupes <- dupes[N>1,V1]
# for(i in dupes) con_indindusocc[, i] <- rowSums(con_indindusocc[, grepl(i,colnames(con_indindusocc))])
# con_indindusocc <- con_indindusocc[ , -which(grepl('[[:digit:]]+', colnames(con_indindusocc)))]

#trim columns and separate into dedicated tables
con_indindus <- con_indindusocc[,c("NATRES","TRANSUTILS","CONST","MFG","WHLTRADE","RETTRADE","INFO","FINESTATE","PROFSCIMGMT","EDUCSOC","ARTSACCOM","OTHER","PUBADMIN")]
con_indocc <- con_indindusocc[,c("MGMTBIZSCIART","SALESOFFICEADMIN","NATRESCONSTMAINT","PRODTRANS","SERVICE")]
con_indindusocc <- con_indindusocc[,-which(colnames(con_indindusocc) %in%
                                             c("MGMTBIZSCIART","SALESOFFICEADMIN","NATRESCONSTMAINT","PRODTRANS","SERVICE",
                                               "NATRES","TRANSUTILS","CONST","MFG","WHLTRADE","RETTRADE","INFO","FINESTATE","PROFSCIMGMT","EDUCSOC","ARTSACCOM","OTHER","PUBADMIN"))]


#placeholder columns, known to be 0, impossible
con_indindusocc[,paste("NOOCC",colnames(con_indindus),sep="")] <- 0
con_indindusocc[,paste(colnames(con_indocc),"NOIND",sep="")] <- 0
#Adding Unemp
con_indindus$NOIND <- rowSums(con_indagesex2) - rowSums(con_indindus)
con_indocc$NOOCC <- rowSums(con_indagesex2) - rowSums(con_indocc)
con_indindusocc$NOOCCNOIND <- rowSums(con_indagesex2) - rowSums(con_indindusocc)

# Transform into an 3D-array
names <- c(list(tracts),dimnames(dimsind)["OCC"], dimnames(dimsind)["INDUS"])
con_indindusocc2 <- array(NA, dim=c(length(tracts),length(dimnames(dimsind)$OCC),length(dimnames(dimsind)$INDUS)), dimnames = names)
for(tract in rownames(con_indindusocc)){
  for (occ in dimnames(con_indindusocc2)$OCC){
    for (indus in dimnames(con_indindusocc2)$INDUS){
      con_indindusocc2[tract,occ,indus] <- con_indindusocc[tract,paste(occ,indus,sep="")]
    }}}

##### Row sum check ####
table(rowSums(con_indagesex2)==rowSums(con_indindusocc)) 
table(rowSums(con_indagesex2)==rowSums(con_indrelate))

##### Observe the global total ####
sum(con_indagesex2)
sum(con_indrelate)
sum(con_indindusocc2)

##### Putting into neat list form ####
margins[["indagesex"]] <- con_indagesex2
margins[["indrelate"]] <- con_indrelate
margins[["indindusocc"]] <- con_indindusocc2
margins[["indindus"]] <- apply(con_indindusocc2, c(1,3), sum)
margins[["indocc"]] <- apply(con_indindusocc2, c(1,2), sum)



for(i in names(margins)){
  if(any(margins[[i]]<0)){
    print(i)
  }
}

# cleanup ####
rm(list=setdiff(ls(),c("clock","runtime","csvwriteout",
                       "tables","microdata","tracts","margins","settings")))

print("Ind margins done")
