
int_trs <- function(x){ # For generalisation purpose, x becomes a vector
  xv <- as.vector(x) # allows trs to work on matrices
  xint <- floor(xv) # integer part of the weight
  r <- xv - xint # decimal part of the weight
  def <- round(sum(r)) # the deficit population # the weights be 'topped up' (+ 1 applied)
  #topup <- sample(length(x), size = def, prob = r, replace = F)
  topup <- RcppSample(length(x), size = def, prob = r)
  xint[topup] <- xint[topup] + 1
  dim(xint) <- dim(x)
  dimnames(xint) <- dimnames(x)
  xint
}
int_expand_array <- function(x){
  # Transform the array into a dataframe
  count_data <- as.data.frame.table(x)
  # Store the indices of categories for the final population
  indices <- rep(1:nrow(count_data), count_data$Freq)
  # Create the final individuals
  ind_data <- count_data[indices,]
  data.table(ind_data)
}
joint_samp <- function(jfit, dest, pums, tr) {
  #Finding best fit random joint sample
  X <- jfit$X
  X <- X[X>0]
  
  dest[dest<0] <- 0
  #print(tr)
  #Draw samples
  sampvec <- sample(names(X), sum(X), T, X)
  
  #Generating joint sample by merging
  jointsamp = merge(data.table(SERIALNO = sampvec, HHID = paste(tr,1:length(sampvec),sep="_"), otract=tr), pums, by = 'SERIALNO')
  
  #Sample destinations
  if(length(dim(dest))>3) {
    #When dest is integrated
    jointsamp[ , dtract := apply(jointsamp, 1, function(samp) {
      n = names(dimnames(dest))[which(names(dimnames(dest))!="dtract")]
      dees = dest[samp[n[1]], samp[n[2]], samp[n[3]], samp[n[4]], samp[n[5]], ]
      sample(names(dees), 1, T, dees)
    })]
  } else {
    #When dest is not integrated
    jointsamp[ , dtract := apply(jointsamp, 1, function(x) sample(names(dest[,x['INDUS']]), 1, T, dest[,x['INDUS']]))]
  }
  
  #Reorder and cleanup
  return(jointsamp[ , c("HHID",colnames(pums[,!"SERIALNO"]),"otract","dtract"),with=F])
}
joint_samp_replace <- function(jfit, di, pums, tr) {
  #Finding best fit random joint sample
  X <- jfit$X
  X <- X[X>0]
  #print(tr)
  #Draw samples
  sampvec <- sample(names(X), sum(X), T, X)
  
  # #Empty table
  # jointsamp <- data.table(SERIALNO = as.character(), HHID = as.character(), otract=as.character(), pums[0,-1], dtract=as.character())
  
  #Generating joint sample by merging
  jointsamp = merge(data.table(SERIALNO = sampvec, HHID = paste(tr,1:length(sampvec),sep="_"), otract=tr), pums, by = 'SERIALNO')
  
  #Sample destinations
  jointsamp[ , dtract := apply(jointsamp_new, 1, function(x) {
    if(sum(di[,x['INDUS']]) > 0) sample(names(di[,x['INDUS']]), 1, T, di[,x['INDUS']])
    else as.character(NA)
  }) ]
  
  # badhh <- unique(jointsamp_new[is.na(dtract),HHID])
  # jointsamp <- rbind(jointsamp, jointsamp_new[!(HHID %in% badhh),])
  # remainder = sum(X) - length(unique(jointsamp$HHID))
  
  #Reorder and cleanup
  return(jointsamp[ , c("HHID",colnames(pums[,!"SERIALNO"]),"otract","dtract"),with=F])
}
joint_samp_old <- function(ind, hh, jfit, tr) {
  #Finding best fit random joint sample
  X <- jfit$X
  X <- X[X>0]
  #Grabbing b vector
  b <- jfit$errors$b
  #Pulling N-samples, keeping best fit
  rmse_min <- 100
  for(i in 1:10){
    samp_i <- sample(x = names(X), size = sum(X), replace = T, prob = X)
    rmse <- sqrt(sum((rowSums(A[ , samp_i])-b)^2)/length(b))
    #qplot(x = b, y = rowSums(A[ , samp_i]), geom='point')
    if(rmse < rmse_min) rmse_min <- rmse ; samp <- samp_i
  }
  #Flatting the IPF weights for individuals
  indflat <- as.data.table(as.data.frame.table(ind))
  indflat <- merge(indflat, types$indtypes[,!"Freq"], by = indvars)
  indflat <- indflat[Freq>0,]
  indflat <- indflat[ , lapply(.SD, as.character)]
  #Flattening and expanding the joint sample
  jointflat <- lapply(1:length(samp), function(x) {
    hh <- types$pumstypes[SERIALNO==samp[x],]
    hh[ , HHID := paste(tr,x,sep="_")]})
  #Adding unique household ID
  names(jointflat) <- paste(tr,1:length(samp), sep="_")
  jointflat <- lapply(names(jointflat), function(x) jointflat[[x]][ , HHID := x])
  jointflat <- rbindlist(jointflat)
  jointflat[ , dtract := as.character()]
  #person types to cycle through
  pertypes <- types[['indtypes']]$PERTYPE
  #Must in BOTH the joint sample and the flat IPF fit
  pertypes <- pertypes[pertypes %in% jointflat$PERTYPE & pertypes %in% indflat$PERTYPE]
  #Assigning destination for persons in chunks of person types
  jointflat <- lapply(pertypes, function(p) {
    #Matched person type subset
    persub <- indflat[PERTYPE == p, ]
    #Subsetting for person types
    jointsub <- jointflat[PERTYPE == p, ]
    #Sample size
    N <- nrow(jointflat[PERTYPE == p, ])
    #Assigning destination using MC sample
    jointsub[, dtract := sample(x = persub$dtract, size = N, replace = T, prob = persub$Freq)]
  })
  jointflat <- rbindlist(jointflat)
  #Merging in the other attributes
  jointflat <- merge(jointflat, types$hhtypes[,!"Freq"], by = "HHTYPE", all.y = F)
  jointflat <- merge(jointflat, types$indtypes[,!"Freq"], by = "PERTYPE", all.y = F)
  #
  jointflat[ , otract := tr]
  return(jointflat)
}
