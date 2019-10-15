

#Downloads the latest or specified LODES data for home-work, home-industry, and work-industry
dl.lodes <- function(state=NULL, year=NULL, filedest=NULL) {
  
  baseurl <- "https://lehd.ces.census.gov/data/lodes/LODES7"
  
  #Basic checks
  if(is.null(state)) stop("Invalid state specified")
  if(is.null(filedest)) {
    filedest <- "./data/raw"
    print(paste("No file destination specified, using default:", filedest))
  }
  if(is.null(year)) {
    url <- paste(baseurl,"ma/od/", sep="/")
    html <- paste(readLines(url), collapse="\n")
    matched <- unlist(str_match_all(html, "<a href=\"(.*?)\""))
    matched <- matched[!grepl("<a href=",matched)]
    matched <- matched[grepl(paste(files,collapse = "|"),matched)]
    year <- max(as.numeric(gsub("ma_od_main_JT00_|.csv.gz","",matched)))
    print(paste("No year specified, defaulting to", year, "as most recent year available"))
  }
  #Ensure lower case
  state <- tolower(state)
  
  #File base
  files <- paste0(state,c("_od_main_JT00_","_rac_S000_JT00_","_wac_S000_JT00_"),year,".csv.gz")
  
  #downloading files
  for(lode in files) {
    fdir <- paste0(filedest, lode)
    if(grepl("od", lode)) 
      prefix <- "od" 
    if(grepl("wac", lode)) 
      prefix <- "wac" 
    if(grepl("rac", lode)) 
      prefix <- "rac"
  
    fdir <- paste(filedest, lode, sep = "/")
    url <- paste(baseurl,state,prefix,lode,sep="/")
    if(!file.exists(fdir)) download.file(url, fdir)
  }
}


#Downloads the latest or specified PUMS data
dl.pums <- function(state=NULL, year=NULL, filedest=NULL) {
  
  baseurl <- "https://www2.census.gov/programs-surveys/acs/data/pums"
  
  #Basic checks
  if(is.null(state)) stop("Invalid state specified")
  if(is.null(filedest)) {
    filedest <- "./data/raw"
    print(paste("No file destination specified, using default:", filedest))
  }
  if(is.null(year)) {
    html <- paste(readLines(baseurl), collapse="\n")
    matched <- unlist(str_match_all(html, "<a href=\"(.*?)\""))
    matched <- matched[!grepl("<a href=",matched)]
    matched <- as.numeric(gsub("/","",matched[grepl("\\d{4}/",matched)]))
    years <- sort(matched)
    
    #Check contents of latest folder, if empty check next most recent
    for(year in rev(years)) {
      url <- paste(baseurl,year,"5-Year/", sep="/")
      html <- tryCatch(paste(readLines(url), collapse="\n"),
               warning = function(w) NA,
               error = function(e) NA)
      if(!is.na(html)) break;
    }
    print(paste("No year specified, defaulting to", year, "as most recent year available"))
  }
  #Ensure lower case
  state <- tolower(state)
  
  #File base
  files <- paste0(c("csv_h","csv_p"),state,".zip")
  
  #downloading files
  for(pums in files) {
    fdir <- paste(filedest, pums, sep = "/")
    url <- paste(baseurl,year,"5-Year",pums, sep="/")
    if(!file.exists(fdir)) download.file(url, fdir)
  }
}



#Downloads the latest or specified PUMS data
dl.tables <- function(state=NULL, year=NULL, filedest=NULL) {
  
  baseurl <- "https://www2.census.gov/programs-surveys/acs/summary_file"
  
  
  #Abbreviations
  statebrevs <- structure(.Data = state.name, .Names = state.abb)
  #Ensure lower case
  state <- tolower(state)
  names(statebrevs) <- tolower(names(statebrevs))
  
  
  #Basic checks
  if(is.null(state)) stop("Invalid state specified")
  if(is.null(filedest)) {
    filedest <- "./data/raw"
    print(paste("No file destination specified, using default:", filedest))
  }
  if(is.null(year)) {
    html <- paste(readLines(baseurl), collapse="\n")
    matched <- unlist(str_match_all(html, "<a href=\"(.*?)\""))
    matched <- matched[!grepl("<a href=",matched)]
    matched <- as.numeric(gsub("/","",matched[grepl("\\d{4}/",matched)]))
    years <- sort(matched)
    
    #Check contents of latest folder, if empty check next most recent
    for(year in rev(years)) {
      url <- paste(baseurl,year,"data/", sep="/")
      html <- paste(readLines(url), collapse="\n")
      matched <- unlist(str_match_all(html, "<a href=\"(.*?)\""))
      if(any(grepl("5_year_by_state",matched))) break;
    }
    print(paste("No year specified, defaulting to", year, "as most recent year available"))
  }
  
  ## download lookup file
  url <- paste(baseurl,year,"documentation/user_tools/ACS_5yr_Seq_Table_Number_Lookup.txt", sep = "/")
  fdir <- paste(filedest,"ACS_5yr_Seq_Table_Number_Lookup.txt", sep = "/")
  if(!file.exists(fdir)) download.file(url, fdir)
  lookup <- fread(fdir)
  
  #determine which files
  tableid <- c("B01001","B08201","B09019","B19001","B25124","C24050")
  tableseq <- unique(lookup[`Table ID` %in% tableid, `Sequence Number`])
  seqfiles <- paste0("e",year,"5",state,sprintf("%04d",tableseq),"000.txt")
  headfiles <- paste0("xls_temp/seq",tableseq,".xlsx")
  
  
  ## download header file
  file <- paste0(year,"_5yr_Summary_FileTemplates.zip")
  url <- paste(baseurl,year,"data", file, sep = "/")
  fdir <- paste(filedest, file, sep = "/")
  if(!file.exists(fdir)) download.file(url, fdir)
  
  #unzip and load header files
  unzip(fdir, files = c(headfiles,paste0("xls_temp/",year,"_SFGeoFileTemplate.xls")), exdir = sub(".zip","",fdir))
  headerdat <- lapply(paste(gsub(".zip","",fdir), headfiles, sep = "/"), function(x) read.xlsx(x, 1))
  names(headerdat) <- tableseq
  
  
  ## download the summary files
  file <- paste0(statebrevs[state],"_Tracts_Block_Groups_Only.zip")
  url <- paste(baseurl,year,"data/5_year_by_state",file, sep = "/")
  fdir <- paste(filedest, file, sep = "/")
  if(!file.exists(fdir)) download.file(url, fdir)
  
  #unzip and load summary files
  unzip(fdir, files = c(seqfiles,paste0("g",year,"5",state,".csv")), exdir = gsub(".zip","",fdir))
  censusdat <- lapply(paste(gsub(".zip","",fdir), seqfiles, sep = "/"), fread)
  names(censusdat) <- tableseq
  
  
  #
  for(x in as.character(tableseq)) colnames(censusdat[[x]]) <- as.matrix(headerdat[[x]])[1,]
  
  
}






