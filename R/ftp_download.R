
library(stringr)
library(data.table)
library(xlsx)


#### Setup the state name and codes
#state code number
statecode <- fread("./data/raw/fips.txt")
statecode[ , STATE := tolower(statecode[,STATE])]

#Abbreviations
statebrevs <- data.table(STATE = state.name, ABBREV = state.abb)

#Ensure lower case
state <- tolower(state)
statebrevs <- statebrevs[ , lapply(.SD, tolower)]
statebrevs <- rbind(statebrevs, data.table(STATE = "district of columbia", ABBREV = "dc"))

#merge state abbreviation to state code
statecode <- merge(statebrevs, statecode, by = "STATE")
rm(statebrevs)


#### FTP Download functions ####
#Downloads the latest or specified LODES data for home-work, home-industry, and work-industry
dl.lodes <- function(state=NULL, year=NULL, filedest="./data/raw") {
  
  baseurl <- "https://lehd.ces.census.gov/data/lodes/LODES7"
  
  #Ensure lower case
  state <- tolower(state)
  
  #Basic checks
  if(is.null(state) | !(state %in% statecode$ABBREV)) 
    stop("Invalid state specified")
  if(!file.exists(filedest)) {
    filedest="./data/raw"
    print(paste("Invalid file destination, using default:", filedest))
  }

  #Checking the year
  url <- paste(baseurl,"ma/od/", sep="/")
  html <- paste(readLines(url), collapse="\n")
  matched <- unlist(str_match_all(html, "<a href=\"(.*?)\""))
  matched <- matched[!grepl("<a href=",matched)]
  matched <- matched[grepl(paste(files,collapse = "|"),matched)]
  yr <- max(as.numeric(gsub("ma_od_main_JT00_|.csv.gz","",matched)))

  if(yr > year) 
    print(paste0("Latest data year available is ", yr,", but ", year, " is OK."))
  if(yr < year) 
    stop(paste0("Data for year ", year, " is unavailable, latest data year is ", yr,"."))
  if(is.null(year)) {
    print(paste("No year specified, defaulting to", yr, "as most recent year available"))
    year = yr
  }
  
  #File base
  files <- paste0(state,c("_od_main_JT00_","_rac_S000_JT00_","_wac_S000_JT00_"),year,".csv.gz")
  
  #downloading files
  for(lode in files) {
    fdir <- paste(filedest, lode, sep = "/")
    if(grepl("od", lode)) prefix <- "od" 
    if(grepl("wac", lode)) prefix <- "wac" 
    if(grepl("rac", lode)) prefix <- "rac"
    url <- paste(baseurl,state,prefix,lode,sep="/")
    if(!file.exists(fdir)) download.file(url, fdir)
  }
  
  #Load LODES into R
  lodes <- lapply(paste(filedest,files,sep = "/"), fread)
  names(lodes) <- c("od","rac","wac")
  
  #Return but do not print
  invisible(lodes)
}

#Downloads the latest or specified PUMS data
dl.pums <- function(state=NULL, year=NULL, filedest= "./data/raw") {
  
  baseurl <- "https://www2.census.gov/programs-surveys/acs/data/pums"
  
  #Ensure lower case
  state <- tolower(state)
  
  #Basic checks
  if(is.null(state) | !(state %in% statecode$ABBREV)) 
    stop("Invalid state specified")
  if(!file.exists(filedest)) {
    filedest="./data/raw"
    print(paste("Invalid file destination, using default:", filedest))
  }
  
  #Checking the data year
  html <- paste(readLines(baseurl), collapse="\n")
  matched <- unlist(str_match_all(html, "<a href=\"(.*?)\""))
  matched <- matched[!grepl("<a href=",matched)]
  matched <- as.numeric(gsub("/","",matched[grepl("\\d{4}/",matched)]))
  years <- sort(matched)
  
  #Check contents of latest folder, if empty check next most recent
  for(yr in rev(years)) {
    url <- paste(baseurl,yr,"5-Year/", sep="/")
    html <- tryCatch(paste(readLines(url), collapse="\n"),
                     warning = function(w) NA,
                     error = function(e) NA)
    if(!is.na(html)) break;
  }
  if(yr > year) 
    print(paste0("Latest data year available is ", yr,", but ", year, " is OK."))
  if(yr < year) 
    stop(paste0("Data for year ", year, " is unavailable, latest data year is ", yr,"."))
  if(is.null(year)) {
    print(paste("No year specified, defaulting to", yr, "as most recent year available"))
    year = yr
  }
  

  #File base
  zips <- paste0(c("csv_h","csv_p"),state,".zip")
  
  #downloading files
  for(pums in zips) {
    fdir <- paste(filedest, pums, sep = "/")
    url <- paste(baseurl,year,"5-Year",pums, sep="/")
    if(!file.exists(fdir)) download.file(url, fdir)
    if(!file.exists(gsub(".zip","",fdir))) unzip(fdir, exdir = gsub(".zip","",fdir))
  }
  
  #finding the files extracted
  files <- list.files(paste(filedest, gsub(".zip","", zips), sep = "/"))
  
  #Reading PUMS data into R
  pums.data <- lapply(paste(filedest, gsub(".zip","", zips), files[grepl(".csv",files)], sep = "/"), fread)
  
  #Cleaning up the extracted zips
  for(pums in zips) unlink(paste(filedest, gsub(".zip","",pums), sep = "/"), recursive = T)
  
  #Return but do not print
  invisible(lodes)
}

#Downloads the latest or specified PUMS data
dl.tables <- function(state=NULL,
                      year=NULL,
                      tables.id = c("B01001","B08201","B09019","B19001","B25124","C24050"),
                      filedest="./data/raw") {
  
  baseurl <- "https://www2.census.gov/programs-surveys/acs/summary_file"
  
  #Basic checks
  if(is.null(state) | !(state %in% statecode$ABBREV)) 
    stop("Invalid state specified")
  if(!file.exists(filedest)) {
    filedest="./data/raw"
    print(paste("Invalid file destination, using default:", filedest))
  }
  
  #Checking the data year
  html <- paste(readLines(baseurl), collapse="\n")
  matched <- unlist(str_match_all(html, "<a href=\"(.*?)\""))
  matched <- matched[!grepl("<a href=",matched)]
  matched <- as.numeric(gsub("/","",matched[grepl("\\d{4}/",matched)]))
  years <- sort(matched)
  
  #Check contents of latest folder, if empty check next most recent
  for(yr in rev(years)) {
    url <- paste(baseurl,yr,"data/", sep="/")
    html <- paste(readLines(url), collapse="\n")
    matched <- unlist(str_match_all(html, "<a href=\"(.*?)\""))
    if(any(grepl("5_year_by_state",matched))) break;
  }
  if(yr > year) 
    print(paste0("Latest data year available is ", yr,", but ", year, " is OK."))
  if(yr < year) 
    stop(paste0("Data for year ", year, " is unavailable, latest data year is ", yr,"."))
  if(is.null(year)) {
    print(paste("No year specified, defaulting to", yr, "as most recent year available"))
    year = yr
  }
  
  ## download lookup file
  url <- paste(baseurl,year,"documentation/user_tools/ACS_5yr_Seq_Table_Number_Lookup.txt", sep = "/")
  fdir <- paste(filedest,"ACS_5yr_Seq_Table_Number_Lookup.txt", sep = "/")
  if(!file.exists(fdir)) download.file(url, fdir)
  lookup <- fread(fdir)
  
  #Check table IDs
  if( !all(tables.id %in% lookup$`Table ID`) ) stop("Invalid tables specified")
  
  #determine which files
  tables.seq <- unique(lookup[`Table ID` %in% tables.id, `Sequence Number`])
  tables.seqfiles <- paste0("e",year,"5",state,sprintf("%04d",tables.seq),"000.txt")
  tables.headfiles <- paste0("xls_temp/seq",tables.seq,".xlsx")
  
  
  ## download header file
  file <- paste0(year,"_5yr_Summary_FileTemplates.zip")
  url <- paste(baseurl,year,"data", file, sep = "/")
  fdir <- paste(filedest, file, sep = "/")
  if(!file.exists(fdir)) download.file(url, fdir)
  
  #unzip and load header files
  unzip(fdir, files = c(tables.headfiles,paste0("xls_temp/",year,"_SFGeoFileTemplate.xls")), exdir = sub(".zip","",fdir))
  tables.headers <- lapply(paste(gsub(".zip","",fdir), tables.headfiles, sep = "/"), function(x) as.character(as.matrix(read.xlsx(x, 1))[1,]))
  names(tables.headers) <- tables.seq
  
  #read geocode header file
  geocode.header <- read.xlsx(paste(gsub(".zip","",fdir), paste0("xls_temp/",year,"_SFGeoFileTemplate.xls"), sep = "/"),1)
  geocode.header <- structure(.Data = as.matrix(geocode.header)[1,], .Names = colnames(geocode.header))
  
  #cleanup extracted zip files
  unlink(gsub(".zip","",fdir), recursive = T)
  
  
  ## download the summary files
  file <- paste0(statecode[ABBREV==state, STATE],"_Tracts_Block_Groups_Only.zip")
  url <- paste(baseurl,year,"data/5_year_by_state",file, sep = "/")
  fdir <- paste(filedest, file, sep = "/")
  if(!file.exists(fdir)) download.file(url, fdir)
  
  #unzip and load summary files
  unzip(fdir, files = c(tables.seqfiles,paste0("g",year,"5",state,".csv")), exdir = gsub(".zip","",fdir))
  tables.data <- lapply(paste(gsub(".zip","",fdir), tables.seqfiles, sep = "/"), fread)
  names(tables.data) <- tables.id
  
  #Read the geocode locations data
  geocode.data <- fread(paste(gsub(".zip","",fdir), paste0("g",year,"5",state,".csv"), sep="/"))
  colnames(geocode.data) <- names(geocode.header)
  
  #cleanup extracted zip files
  unlink(gsub(".zip","",fdir), recursive = T)
  
  #Naming the census table data columns and adding GEOID
  for(x in as.character(tables.seq)) {
    colnames(tables.data[[x]]) <- tables.headers[[x]]
    tables.data[[x]] <- merge(geocode.data[ , .(LOGRECNO,GEOID,NAME,STATE,COUNTY,TRACT)], tables.data[[x]], by = "LOGRECNO")
  }
  
  #Return but do not print
  invisible(lodes)
}

#Download the latest shapefiles
dl.geodata <- function(state=NULL, year=NULL, filedest="./data/raw") {
  
  baseurl <- "https://www2.census.gov/geo/tiger"
 
  #Basic checks
  if(is.null(state) | !(state %in% statecode$ABBREV)) 
    stop("Invalid state specified")
  if(!file.exists(filedest)) {
    filedest="./data/raw"
    print(paste("Invalid file destination, using default:", filedest))
  }
  if(is.null(year)) {
    html <- paste(readLines(baseurl), collapse="\n")
    years <- unlist(str_match_all(html, "<a href=\"(.*?)\""))
    years <- years[!grepl("<a href=",years)]
    years <- as.numeric(gsub("GENZ|/","",years[grepl("GENZ\\d{4}/",years)]))
    years <- sort(years)
    year <- max(years)
    print(paste("No year specified, defaulting to", year, "as most recent year available"))
  }
  
  #
  file <- paste("cb",year,statecode[ABBREV==state,CODE],"bg_500k.zip", sep = "_")
  paste(baseurl,paste0("GENZ",year),
        paste("cb",year,statecode[ABBREV==state,CODE],"bg_500k.zip", sep = "_"), sep = "/")
  
  
}



