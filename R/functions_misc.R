

fun_raketract.singledim <- function(tr, freq.dat, marg.dat, bnd=c(0,Inf), trm=NULL, forcestrict=F, skip=F) {
  #Pulling the marginal
  nm <- names(marg.dat)
  levs.og <- lapply(nm, function(n) levels(marg.dat[[n]][[n]]))
  names(levs.og) <- nm
  
  #Pulling the marginal & data for current tract
  marg.dat <- lapply(marg.dat, function(x) x[otract == tr, !"otract"])
  if('otract' %in% colnames(freq.dat)) freq.dat <- freq.dat[otract == tr, !"otract"]
  
  #Removing margin zeros from data
  zeros <- lapply(marg.dat, function(x) x[Freq==0, ])
  for(n in names(zeros)) freq.dat[['N']][ freq.dat[[n]] %in% zeros[[n]][[n]] ] <- 0
  freq.dat <- freq.dat[N>0,]
  
  #Removing the zeros from the margins
  marg.dat <- lapply(marg.dat, function(x) x[Freq!=0, ])
  
  #Relevel the margin factors
  for(n in names(marg.dat)) {
    levs = levels(marg.dat[[n]][[n]])[levels(marg.dat[[n]][[n]]) %in% marg.dat[[n]][[n]]]
    marg.dat[[n]][[n]] <- factor(marg.dat[[n]][[n]], levels = levs)
  }
  #Matching factors in data to margin factors
  for(n in names(marg.dat)) freq.dat[[n]] <- factor(freq.dat[[n]], levels = levels(marg.dat[[n]][[n]]))
  
  #Scaling the weights
  total <- sum(marg.dat[[1]]$Freq)
  colnames(freq.dat)[which("N" == colnames(freq.dat))] <- "weights.sim"
  freq.dat[ , weights.sim := total*weights.sim/sum(weights.sim)]
  
  #Check if only 1 factor level remaining, aka only one household type
  singles <- sapply(marg.dat, function(x) ifelse(nrow(x[Freq>0,])<=1, x[[1]],NA))
  
  if(sum(is.na(singles)) > 0 ) {
    #Removing the singles temporarily
    freq.dat.stripped <- freq.dat[ , c(names(singles[is.na(singles)]),"weights.sim"), with=F]
    #Setting up the raking data
    design.dat <- svydesign(ids = ~1, weights=~weights.sim, data=freq.dat.stripped)
    #Formula
    form <- as.formula(paste("~",paste(names(singles[is.na(singles)]),collapse = "+"),sep=""))
    #Population marginals
    pop <- c(`(Intercept)` = total, unlist(lapply(names(singles[is.na(singles)]), function(n) {
      setNames(marg.dat[[n]]$Freq[-1], paste(n,marg.dat[[n]][[n]][-1],sep=""))
    })))
    raked.dat = 'character'
    #Raking the data, strictest standards
    if(skip==F) {
      raked.dat <- suppressWarnings(
        tryCatch(
          calibrate(design.dat, formula = form, population = pop, calfun = "raking", bounds = bnd, trim = trm, force=forcestrict),
          error = function(e) "error"))
    }
    #Try again with relaxed constraints
    if(all(class(raked.dat) == "character")) {
      print("Relaxed boundaries")
      raked.dat <- suppressWarnings(
        tryCatch(
          calibrate(design.dat, formula = form, population = pop, calfun = "linear", bounds = c(0,1000), trim=trm, eps=1, force=T),
          error = function(e) "error"))
    }
    
    #Putting adding the weights back in
    if(all(class(raked.dat) == "character")) {
      freq.dat.raked <- cbind("otract" = tr, freq.dat, "weights.raked" = freq.dat$weights.sim)
    } else {
      freq.dat.raked <- cbind("otract" = tr, freq.dat, "weights.raked" = weights(raked.dat))
    }
  } else {
    #No change possible, just return the original, with eps = 0
    freq.dat.raked <- cbind("otract" = tr, freq.dat, "weights.raked" = freq.dat$weights.sim)
  }
  
  #Convert back to original factors
  for(n in nm) freq.dat.raked[[n]] <- as.character(freq.dat.raked[[n]])
  for(n in names(marg.dat)) freq.dat.raked[[n]] <- factor(freq.dat.raked[[n]], levels = levs.og[[n]])
  #Out
  return(freq.dat.raked)
}
fun_raketract.singledim <- function(tr, freq.dat, marg.dat, bnd=c(0,Inf), trm=NULL) {
  #Pulling the marginal
  nm <- names(marg.dat)
  levs.og <- lapply(nm, function(n) levels(marg.dat[[n]][[n]]))
  names(levs.og) <- nm
  
  #Pulling the marginal & data for current tract
  marg.dat <- lapply(marg.dat, function(x) x[otract == tr, !"otract"])
  if('otract' %in% colnames(freq.dat)) freq.dat <- freq.dat[otract == tr, !"otract"]
  
  #Removing margin zeros from data
  zeros <- lapply(marg.dat, function(x) x[Freq==0, ])
  for(n in names(zeros)) freq.dat[['N']][ freq.dat[[n]] %in% zeros[[n]][[n]] ] <- 0
  freq.dat <- freq.dat[N>0,]
  
  #Removing the zeros from the margins
  marg.dat <- lapply(marg.dat, function(x) x[Freq!=0, ])
  
  #Relevel the margin factors
  for(n in names(marg.dat)) {
    levs = levels(marg.dat[[n]][[n]])[levels(marg.dat[[n]][[n]]) %in% marg.dat[[n]][[n]]]
    marg.dat[[n]][[n]] <- factor(marg.dat[[n]][[n]], levels = levs)
  }
  #Matching factors in data to margin factors
  for(n in names(marg.dat)) freq.dat[[n]] <- factor(freq.dat[[n]], levels = levels(marg.dat[[n]][[n]]))
  
  #Scaling the weights
  total <- sum(marg.dat[[1]]$Freq)
  colnames(freq.dat)[which("N" == colnames(freq.dat))] <- "weights.sim"
  freq.dat[ , weights.sim := total*weights.sim/sum(weights.sim)]
  
  #Check if only 1 factor level remaining, aka only one household type
  singles <- sapply(marg.dat, function(x) ifelse(nrow(x[Freq>0,])<=1, x[[1]],NA))
  
  if(sum(is.na(singles)) > 0 ) {
    #Removing the singles temporarily
    freq.dat.stripped <- freq.dat[ , c(names(singles[is.na(singles)]),"weights.sim"), with=F]
    #Setting up the raking data
    design.dat <- svydesign(ids = ~1, weights=~weights.sim, data=freq.dat.stripped)
    #Formula
    form <- as.formula(paste("~",paste(names(singles[is.na(singles)]),collapse = "+"),sep=""))
    #Population marginals
    pop <- c(`(Intercept)` = total, unlist(lapply(names(singles[is.na(singles)]), function(n) {
      setNames(marg.dat[[n]]$Freq[-1], paste(n,marg.dat[[n]][[n]][-1],sep=""))
    })))
    #Raking the data, strictest standards
    raked.dat <- suppressWarnings(
      tryCatch(calibrate(design.dat, formula = form, population = pop, calfun = "raking", bounds = bnd, trim = trm),
               error = function(e) "error"))
    #Try again with no lower boundaries
    if(all(raked.dat == "error")) {
      print("Relaxed boundaries")
      raked.dat <- suppressWarnings(
        tryCatch(calibrate(design.dat, formula = form, population = pop, calfun = "linear", bounds = c(0,1000), force=T),
                 error = function(e) "error"))
    }
    #Putting the weights back into the original data
    freq.dat.raked <- cbind("otract" = tr, freq.dat, "weights.raked" = weights(raked.dat))
  } else {
    #No change possible, just return the original, with eps = 0
    freq.dat.raked <- cbind("otract" = tr, freq.dat, "weights.raked" = freq.dat$weights.sim)
  }
  
  #Convert back to original factors
  for(n in nm) freq.dat.raked[[n]] <- as.character(freq.dat.raked[[n]])
  for(n in names(marg.dat)) freq.dat.raked[[n]] <- factor(freq.dat.raked[[n]], levels = levs.og[[n]])
  #Out
  return(freq.dat.raked)
}
fun_plotnetwork <- function(structure, ht = "400px"){
  nodes.uniq <- unique(c(structure$arcs[,1], structure$arcs[,2]))
  nodes <- data.frame(id = nodes.uniq,
                      label = nodes.uniq,
                      color = "darkturquoise",
                      font.face = "serif",
                      shadow = TRUE,
                      size = 10)
  
  edges <- data.frame(from = structure$arcs[,1],
                      to = structure$arcs[,2],
                      arrows = "to",
                      smooth = TRUE,
                      shadow = TRUE,
                      color = "black")
  
  return(visNetwork(nodes, edges, height = ht, width = "100%"))
}


