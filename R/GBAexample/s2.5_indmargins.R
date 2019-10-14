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
marg_indagesex<-tables[["indagesex"]]
#removing excess header rows
colnames(marg_indagesex)<-marg_indagesex[1,]
marg_indagesex<-marg_indagesex[-1,]
rownames(marg_indagesex)<-marg_indagesex$Id2

#removing excess columns and column text
marg_indagesex<-marg_indagesex[,-c(4, which(grepl("Id|Id2|Geography|Margin",colnames(marg_indagesex))))]
colnames(marg_indagesex)<-gsub("Estimate;| |:|-|,", "", colnames(marg_indagesex))
# marg_indagesex<-marg_indagesex[,-c(1,2,3,which(grepl("Percent",colnames(marg_indagesex))))]
# marg_indagesex<-marg_indagesex[,which(grepl("SEX AND AGE",colnames(marg_indagesex)))]
# marg_indagesex<-marg_indagesex[,which(grepl("Male",colnames(marg_indagesex))|grepl("Female",colnames(marg_indagesex)))]
# marg_indagesex<-marg_indagesex[,-which(grepl("Median",colnames(marg_indagesex)))]
#saving this for proportioning later
# agesex16<-marg_indagesex[,which(grepl("16 years and over",colnames(marg_indagesex)))]
# agesex16[,]<-apply(agesex16[,], 2, as.numeric)
# colnames(agesex16) <- c("MALE", "FEMALE")
# agesex16<-agesex16[tracts,]
#
# marg_indagesex<-marg_indagesex[,-which(grepl("yearsandover",colnames(marg_indagesex)) & !grepl("85",colnames(marg_indagesex)))]
# marg_indagesex<-marg_indagesex[,which(grepl("population -",colnames(marg_indagesex)))]
# colnames(marg_indagesex)<-gsub("Number; SEX AND AGE| |-|population|years", "", colnames(marg_indagesex))
# colnames(marg_indagesex)<-gsub("Under", "0TO", colnames(marg_indagesex))
colnames(marg_indagesex)<-gsub("Under", "0TO", colnames(marg_indagesex))
colnames(marg_indagesex)<-gsub("years", "", colnames(marg_indagesex))
colnames(marg_indagesex)<-gsub("and", "TO", colnames(marg_indagesex))
colnames(marg_indagesex)<-gsub("TOover", "ANDOVER", colnames(marg_indagesex))
colnames(marg_indagesex)<-toupper(colnames(marg_indagesex))

#subsetting only by matching mts tracts
marg_indagesex<-marg_indagesex[tracts,]

#converting to numeric
marg_indagesex[,]<-apply(marg_indagesex[,], 2, as.numeric)

#grouping by ages
for(i in dimnames(dimsind)$GEND) {
  marg_indagesex[,paste(i,"_AGE0TO9",sep="")]<-rowSums(marg_indagesex[, paste(i, c("0TO5","5TO9"),sep="")])
  marg_indagesex[,paste(i,"_AGE10TO14",sep="")]<-marg_indagesex[, paste(i, c("10TO14"),sep="")]
  marg_indagesex[,paste(i,"_AGE15TO19",sep="")]<-rowSums(marg_indagesex[, paste(i, c("15TO17","18TO19"),sep="")])
  #marg_indagesex[,paste(i,"A_GE15TO19",sep="")]<-marg_indagesex[, paste(i, c("15TO19"),sep="")]
  #marg_indagesex[,paste(i,"_AGE20TO24",sep="")]<-marg_indagesex[, paste(i, c("20TO24"),sep="")]
  marg_indagesex[,paste(i,"_AGE20TO24",sep="")]<-rowSums(marg_indagesex[, paste(i, c("20","21","22TO24"),sep="")])
  marg_indagesex[,paste(i,"_AGE25TO44",sep="")]<-rowSums(marg_indagesex[, paste(i, c("25TO29","30TO34","35TO39","40TO44"),sep="")])
  marg_indagesex[,paste(i,"_AGE45TO54",sep="")]<-rowSums(marg_indagesex[, paste(i, c("45TO49","50TO54"),sep="")])
  marg_indagesex[,paste(i,"_AGE55TO64",sep="")]<-rowSums(marg_indagesex[, paste(i, c("55TO59","60TO61","62TO64"),sep="")])
  marg_indagesex[,paste(i,"_AGE65UP",sep="")]<-rowSums(marg_indagesex[, paste(i, c("65TO66","67TO69","70TO74","75TO79","80TO84","85ANDOVER"),sep="")])
}

marg_indagesex<-marg_indagesex[,which(grepl("AGE", colnames(marg_indagesex)))]

# Transform age-sex into an 3D-array
names <- c(list(tracts),dimnames(dimsind)[c("GEND","AGE")])
marg_indagesex2 <- array(NA, dim=c(length(tracts),length(dimnames(dimsind)$GEND),length(dimnames(dimsind)$AGE)), dimnames = names)

for(tract in rownames(marg_indagesex)){
  for (sex in dimnames(marg_indagesex2)$GEND){
    for (age in dimnames(marg_indagesex2)$AGE){
      marg_indagesex2[tract,sex,age] <- marg_indagesex[tract,paste(sex,age,sep="_")]
    }}}

##### Ind Relationship ################
print("Setting relationship tables")
#removing excess header rows
marg_indrelate<-tables[["indrelate"]]
colnames(marg_indrelate)<-marg_indrelate[1,]
marg_indrelate<-marg_indrelate[-1,]
rownames(marg_indrelate)<-marg_indrelate$Id2
#subsetting only by matching mts tracts
marg_indrelate<-marg_indrelate[tracts,]

#removing excess columns and column text
marg_indrelate<-marg_indrelate[,-c(which(grepl("Id|Id2|Geography|Margin",colnames(marg_indrelate))))]
colnames(marg_indrelate)<-gsub("Estimate;| |:|-|,", "", colnames(marg_indrelate))
#marg_indrelate<-marg_indrelate[,-c(1:3)]
colnames(marg_indrelate)<-toupper(colnames(marg_indrelate))
#colnames(marg_indrelate)<-gsub(" ", "", colnames(marg_indrelate))
#colnames(marg_indrelate)<-gsub(":", "", colnames(marg_indrelate))
#colnames(marg_indrelate)<-gsub("-", "", colnames(marg_indrelate))
colnames(marg_indrelate)<-gsub("INHOUSEHOLDSIN", "", colnames(marg_indrelate))
marg_indrelate<-marg_indrelate[ , -which(colnames(marg_indrelate) %in%
                                         c("TOTAL","INHOUSEHOLDS","FAMILYHOUSEHOLDS","NONFAMILYHOUSEHOLDS",
                                           "FAMILYHOUSEHOLDSHOUSEHOLDERMALE","FAMILYHOUSEHOLDSHOUSEHOLDERFEMALE",
                                           "NONFAMILYHOUSEHOLDSHOUSEHOLDERMALE","NONFAMILYHOUSEHOLDSHOUSEHOLDERFEMALE"))]
marg_indrelate<-marg_indrelate[ , -which(grepl("BIOLOGICAL|ADOPTED|FOSTER|STEPCHILD", colnames(marg_indrelate)))]
marg_indrelate<-marg_indrelate[ , -which(grepl("ROOMER|HOUSEMATE|UNMARRIED|FOSTER|OTHERNONRELATIVE|LIVINGALONE", colnames(marg_indrelate)))]
colnames(marg_indrelate)<-gsub("NONFAMILYHOUSEHOLDS|FAMILYHOUSEHOLDS", "", colnames(marg_indrelate))

#converting to numeric and removing bad text
marg_indrelate[,]<-data.frame(apply(marg_indrelate[,], 2, function(y) as.numeric(gsub("*\\(.*?\\)", "", y))))

#Grouping
marg_indrelate$HEAD <- rowSums(marg_indrelate[ , which(grepl("HOUSEHOLDER",colnames(marg_indrelate)))])
marg_indrelate$NONRELATIVE <- rowSums(marg_indrelate[ , grepl("NONRELATIVES|INGROUPQUARTERS",colnames(marg_indrelate))])
marg_indrelate$CHILD <- rowSums(marg_indrelate[ , which(grepl("CHILD",colnames(marg_indrelate)))])
marg_indrelate$RELATIVE <- rowSums(marg_indrelate[ , which(grepl("SISTER|PARENT|SON|OTHERRELATIVES",colnames(marg_indrelate)))])

#consolidating
marg_indrelate <- marg_indrelate[ , c("HEAD","SPOUSE","CHILD","RELATIVE","NONRELATIVE")]

##### Ind occupation & Industry #####
print("Setting up industry & occupation tables")
marg_indindusocc<-tables[["indindusocc"]]
#removing excess header rows
colnames(marg_indindusocc)<-marg_indindusocc[1,]
marg_indindusocc<-marg_indindusocc[-1,]
rownames(marg_indindusocc)<-marg_indindusocc$Id2
#cleanup
marg_indindusocc<-marg_indindusocc[,-c(4, which(grepl("Id|Id2|Geography|Margin",colnames(marg_indindusocc))))]
colnames(marg_indindusocc)<-gsub("Estimate; Total| |:|-|,", "", colnames(marg_indindusocc))
#Major occ categories
colnames(marg_indindusocc)<-gsub("managementbusinessscienceandartsoccupations", "MGMTBIZSCIART_", colnames(marg_indindusocc), ignore.case = T)
colnames(marg_indindusocc)<-gsub("serviceoccupations", "SERVICE_", colnames(marg_indindusocc), ignore.case = T)
colnames(marg_indindusocc)<-gsub("salesandofficeoccupations", "SALESOFFICEADMIN_", colnames(marg_indindusocc), ignore.case = T)
colnames(marg_indindusocc)<-gsub("naturalresourcesconstructionandmaintenanceoccupations", "NATRESCONSTMAINT_", colnames(marg_indindusocc), ignore.case = T)
colnames(marg_indindusocc)<-gsub("productiontransportationandmaterialmovingoccupations", "PRODTRANS_", colnames(marg_indindusocc), ignore.case = T)
#major indus categories
colnames(marg_indindusocc)<-gsub("Agricultureforestryfishingandhuntingandmining", "NATRES", colnames(marg_indindusocc), ignore.case = T)
colnames(marg_indindusocc)<-gsub("Construction", "CONST", colnames(marg_indindusocc), ignore.case = T)
colnames(marg_indindusocc)<-gsub("Manufacturing", "MFG", colnames(marg_indindusocc), ignore.case = T)
colnames(marg_indindusocc)<-gsub("Wholesaletrade", "WHLTRADE", colnames(marg_indindusocc), ignore.case = T)
colnames(marg_indindusocc)<-gsub("Retailtrade", "RETTRADE", colnames(marg_indindusocc), ignore.case = T)
colnames(marg_indindusocc)<-gsub("Transportationandwarehousingandutilities", "TRANSUTILS", colnames(marg_indindusocc), ignore.case = T)
colnames(marg_indindusocc)<-gsub("Information", "INFO", colnames(marg_indindusocc), ignore.case = T)
colnames(marg_indindusocc)<-gsub("Financeandinsuranceandrealestateandrentalandleasing", "FINESTATE", colnames(marg_indindusocc), ignore.case = T)
colnames(marg_indindusocc)<-gsub("Professionalscientificandmanagementandadministrativeandwastemanagementservices", "PROFSCIMGMT", colnames(marg_indindusocc), ignore.case = T)
colnames(marg_indindusocc)<-gsub("Educationalservicesandhealthcareandsocialassistance", "EDUCSOC", colnames(marg_indindusocc), ignore.case = T)
colnames(marg_indindusocc)<-gsub("Artsentertainmentandrecreationandaccommodationandfoodservices", "ARTSACCOM", colnames(marg_indindusocc), ignore.case = T)
colnames(marg_indindusocc)<-gsub("Otherservicesexceptpublicadministration", "OTHER", colnames(marg_indindusocc), ignore.case = T)
colnames(marg_indindusocc)<-gsub("Publicadministration", "PUBADMIN", colnames(marg_indindusocc), ignore.case = T)
#subsetting only by matching mts tracts
marg_indindusocc<-marg_indindusocc[tracts,]
#converting to numeric
marg_indindusocc[,]<-apply(marg_indindusocc[,], 2, as.numeric)

# #Combining duplicates
# dupes <- data.table(table(gsub('[[:digit:]]+', '', colnames(marg_indindusocc))))
# dupes <- dupes[N>1,V1]
# for(i in dupes) marg_indindusocc[, i] <- rowSums(marg_indindusocc[, grepl(i,colnames(marg_indindusocc))])
# marg_indindusocc <- marg_indindusocc[ , -which(grepl('[[:digit:]]+', colnames(marg_indindusocc)))]

#trim columns and separate into dedicated tables
marg_indindus <- marg_indindusocc[,c("NATRES","TRANSUTILS","CONST","MFG","WHLTRADE","RETTRADE","INFO","FINESTATE","PROFSCIMGMT","EDUCSOC","ARTSACCOM","OTHER","PUBADMIN")]
marg_indocc <- marg_indindusocc[,c("MGMTBIZSCIART_","SALESOFFICEADMIN_","NATRESCONSTMAINT_","PRODTRANS_","SERVICE_")]
colnames(marg_indocc) <- gsub("_","",colnames(marg_indocc))

marg_indindusocc <- marg_indindusocc[,-which(colnames(marg_indindusocc) %in%
                                             c("MGMTBIZSCIART_","SALESOFFICEADMIN_","NATRESCONSTMAINT_","PRODTRANS_","SERVICE_",
                                               "NATRES","TRANSUTILS","CONST","MFG","WHLTRADE","RETTRADE","INFO","FINESTATE","PROFSCIMGMT","EDUCSOC","ARTSACCOM","OTHER","PUBADMIN"))]


#placeholder columns, known to be 0, impossible
marg_indindusocc[,paste("NOOCC_",colnames(marg_indindus),sep="")] <- 0
marg_indindusocc[,paste(colnames(marg_indocc),"_NOIND",sep="")] <- 0
#Adding Unemp
marg_indindus$NOIND <- rowSums(marg_indagesex2) - rowSums(marg_indindus)
marg_indocc$NOOCC <- rowSums(marg_indagesex2) - rowSums(marg_indocc)
marg_indindusocc$NOOCC_NOIND <- rowSums(marg_indagesex2) - rowSums(marg_indindusocc)

# Transform into an 3D-array
names <- c(list(tracts),dimnames(dimsind)["OCC"], dimnames(dimsind)["INDUS"])
marg_indindusocc2 <- array(NA, dim=c(length(tracts),length(dimnames(dimsind)$OCC),length(dimnames(dimsind)$INDUS)), dimnames = names)
for(tract in rownames(marg_indindusocc)){
  for (occ in dimnames(marg_indindusocc2)$OCC){
    for (indus in dimnames(marg_indindusocc2)$INDUS){
      marg_indindusocc2[tract,occ,indus] <- marg_indindusocc[tract,paste(occ,indus,sep="_")]
    }}}

##### Row sum check ####
table(rowSums(marg_indagesex2)==rowSums(marg_indindusocc)) 
table(rowSums(marg_indagesex2)==rowSums(marg_indrelate))

##### Observe the global total ####
sum(marg_indagesex2)
sum(marg_indrelate)
sum(marg_indindusocc2)

##### Putting into neat list form ####
margins[["indagesex"]] <- marg_indagesex2
margins[["indrelate"]] <- marg_indrelate
margins[["indindusocc"]] <- marg_indindusocc2
margins[["indindus"]] <- apply(marg_indindusocc2, c(1,3), sum)
margins[["indocc"]] <- apply(marg_indindusocc2, c(1,2), sum)



for(i in names(margins)){
  if(any(margins[[i]]<0)){
    print(i)
  }
}

# marginal_indagesex = as.data.table(marg_indagesex,keep.rownames = T)
# marginal_indrelate = as.data.table(marg_indrelate,keep.rownames = T)
# marginal_indindusocc = as.data.table(marg_indindusocc,keep.rownames = T)
# setnames(marginal_indagesex, "rn","tract")
# setnames(marginal_indrelate, "rn","tract")
# setnames(marginal_indindusocc, "rn","tract")
# fwrite(marginal_indindusocc, file = "./data/marginal_ind_relate.csv")
# fwrite(marginal_indindusocc, file = "./data/marginal_ind_agesex.csv")
# fwrite(marginal_indindusocc, file = "./data/marginal_ind_indusocc.csv")
# 
# names(dimnames(marg_indagesex2)) <- c("tracts","GEND","AGE")
# names(dimnames(marg_indindusocc2)) <- c("tracts","OCC","INDUS")
# 
# marginals_ind <- list("agesex"=marg_indagesex2, "indusocc" = marg_indindusocc2, "relate" = marginal_indrelate)
# save(marginals_ind, file = "./data/marginal_ind.RData")



print("Ind margins done")
