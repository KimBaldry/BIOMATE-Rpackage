#' @title Reformat profiling sensor (CTD) data to WHP-Exchange format.
#'
#' @author Kimberlee Baldry
#' @description BIOMATE reformats profiling sensor files and saves them as a file that is tagged by EXPOCODE, STATION and CAST. The resulting file follows WHP-Exchange format guidelines. Information is required from the user, stored within a user-entered metadata file particularly for batch runs, or provided via the BIOMATE submission app. The metadata file must be named "PROF_meta.csv". All PROF/CTD files must be singular per cast with the same extention, in one folder. Use the split_delim_file function to convert files to this format.
#' @note This function only supports well structured text delimited and NetCDF file formats.
#'
#' @return returns a successful run message and saves reformatted data to a "profiling_sensors" subdirectory in path_out.
#'
#' @param file_path The path to the processing metadata file, titled "PROF_meta.csv".
#' @param path_out The path where reformatted files are to be sent. A subdirectory "profiling_sensors" will be created in here, if one does not already exist.
#' @param userID An ID code for the user creating the files.
#' @param row_start The start row of the processing metadata table, for batch run worflows.
#' @param row_end The end row of the processing metadata table, for batch run worflows.
#' @param trace_file TRUE will print out each file name as it is reformatted, FALSE suppresses this. It helps with de-bugging and identifying formatting errors during development.
#'
#' @import dplyr
#' @import data.table
#' @import tools
#' @import ncdf4
#' @import stringi
#' @import tidyr
#' @import readr
#' @import geosphere
#' @import lubridate
#'
#' @export



PROF_to_WHPE = function(file_path, path_out,userID = "IMASUTASKB",row_start = 1,row_end = NA, trace_file = F){

  # this small function prevents the drop of midnight 00:00:00
  print.POSIXct2 <- function(x){format(x,"%Y-%m-%d %H:%M:%S %Z")}
  ### Output directory
  # if the output directory does not exist for this layer, create it
  out_dir = file.path(path_out,"profiling_sensors")
  if(!file.exists(out_dir))
  {dir.create(out_dir)}

  # headers that are in the metadata file
  out_header_data = c("STNNBR",	"CASTNO", "DATE","TIME_s",	"TIME_b",	"TIME_e",	"LATITUDE_s",	"LONGITUDE_s",	"LATITUDE_b",	"LONGITUDE_b",	"LATITUDE_e",	"LONGITUDE_e")

  ### read in metadata file
  # check that the file exists
  meta_path = file.path(file_path,"PROF_meta.csv")
  if(!file.exists(meta_path)){stop(paste("Error: There is no file named PROF_meta.csv located in the dirctory",file_path))}
  # read it in
  meta = fread(file = meta_path,header = T,strip.white = T,stringsAsFactors = F, keepLeadingZeros = T)
  heads = as.character(fread(meta_path, nrows = 1, header = F, stringsAsFactors = F))
  heads = heads[heads != "NA"]
  meta = meta[,..heads]
  #if(length(which(!is.na(meta[,ncol(meta)]))) == 0){meta = meta[,-ncol(meta)]}
  if(is.na(row_end)){row_end = nrow(meta)}

  ### reformat files ###
  # loop through entries in the meta file
  for(rw in row_start:row_end)
  {
    # relevent metadata for the row rw
    info = meta[rw,]
    ex = info$EXPOCODE
    if(info$header_sep == "colon"){info$header_sep = ":"}
    if(info$header_sep == "comma"){info$header_sep = ","}
    if(info$header_sep == "dash"){info$header_sep = "-"}
    if(info$header_sep == "equals"){info$header_sep = "="}
    if(info$header_sep == "space"){info$header_sep = " "}
    if(is.empty(info$header_sep)){info$header_sep = ""}


    # list ctd files
    ctd_files = list.files(path = info$path, pattern = info$extention,full.names = T)

    # list variables that have info in the ctd files
    vars = colnames(info)[-c(1:which(colnames(info) == "start_data_var"))]
    vars = vars[!is.na(info[,..vars])]
    vars = vars[!(info[,..vars] == "")]
    # variables where data is stored in file headers - if any
    if(any(grepl("header-",info[,..vars]))){header_vars = vars[grepl("header-",info[,..vars])]
    # variables
    data_vars = vars[!grepl("header-",info[,..vars])]}else{data_vars = vars}
    # information required to ID samples
    ID_info = info[,..out_header_data]
    if(any(grepl("header-",ID_info))){
      rm_head = colnames(ID_info)[!grepl("header-",ID_info)]
      ID_info = ID_info[,..rm_head]}


    # loop through each file
    for(fl in ctd_files){
      old_file = fl

      # data stored in headers
      if(exists("header_data")){rm(header_data)}
      header_data = data.frame(matrix("",nrow = 1, ncol = length(out_header_data)), stringsAsFactors = F)
      colnames(header_data)= out_header_data

      # text file reformating
      if(info$file_type == "text delim"){
      # open file and read lines


      # get number lines in the file
      f <- file( fl, open = "r" )
      lines = readLines( f, -1L)
      close(f)
      n_lines = length(lines)

      # read header lines
      f <- file( fl, open = "r" )
      n = 0
      while( TRUE ){
        line <- readLines( f, 1L )
        n = n+1
        if(is.empty(line)){next}
        if(n == n_lines){print("Error: Couldnt find header line. Headers must be mislabled in metadata?")
          break}

        header_in_line = unlist(lapply(data_vars,function(x){grepl(info[,..x],line, fixed = T)}))
        # if there are at least three data_vars that apear, break this is likeley the header line. print the header line.
        if(length(which(header_in_line == T) ) > 2 & substr(line,1,1) != "#"){
          # write to a file to check header line is picked and variables have been grabbed. This file is for debugging purposes
          # remove this later
          tf = file(file.path(out_dir,paste("Check_headers",Sys.Date(),".txt", sep = "")),open = "at")
          writeLines(paste(old_file),tf)
          writeLines(paste("The header line picked:",line),tf)
          writeLines(paste("The variables found:",data_vars[which(header_in_line ==T)]),tf)
          writeLines("All good?",tf)
          writeLines(" ",tf)
          header_line = line
          close(tf)
          break}
        if(length(line) == 0){break}

        ### get header data if it is in the header of the file
        if(exists("header_vars")){

        # for each header variable store the associated data
          for(hd in header_vars){
          head_pattern = sub("header-","",info[,..hd])
          if(grepl(head_pattern, line ,fixed = T)){
            # common header starters
              line = sub("[#]|[*]|[%]|[!]","",line)
              header_data[,hd] = trimws(sub(info$header_sep,"",sub(head_pattern,"",line, fixed = T)))
              if(is.logical(header_data[,hd])){header_data[,hd] = ""}}
          }

        }

      }
      close( f )


      # get number of the header line to skip
      n = grep(header_line,lines, fixed = T)
      if(length(n) == 0){ n = 1}
      # check for extra lines before the data table starts.
      # There should be some numeric data in the table, but not in the header or unit lines.
      b = n+1
      line = fread(fl,stringsAsFactors = F, skip = n, nrows = 1, header = F, keepLeadingZeros = T)
      while(!any(unlist(lapply(line,is.numeric)))){
        line = fread(fl,stringsAsFactors = F, skip = b, nrows = 1, header = F, keepLeadingZeros = T)
        b = b+1}

      ### Get data table
      # Use fread to read the data table - this will adapt if there is a units line or not.

      # catch the out-of-format lines here for debugging
      #tryCatch({data = as.data.frame(fread(fl,stringsAsFactors = F, skip = n, na.strings = info$missing_value,strip.white = T , header = F))}, warning=function(w) print(fl))

      # rectangular data shouldnt be read with fread
      if(info$delim == "rect"){
        data = read_table2(fl,col_names = F, skip = b-1, na = as.character(info$missing_value),col_types = cols())
        if(info$source == "MGDS"){
          ## this format has a manual insertion. The header lines arent available for all data columns and they break the rectangular format
          headers = as.character(fread(fl,stringsAsFactors = F, skip = n-1, nrows = 1, header = F, keepLeadingZeros = T, tz = info$TZ))
        }else{
         # get the header line
          headers = fread(fl, skip = n-1, nrows = 1, header = F, stringsAsFactors = F, keepLeadingZeros = T)
          headers = as.character(headers[1,])
        }


        }else{
      if(is.na(info$missing_value)){data = as.data.frame(fread(fl,stringsAsFactors = F, skip = b-1,strip.white = T, header = F, keepLeadingZeros = T, tz = info$TZ))
}else{
        data = as.data.frame(fread(fl,stringsAsFactors = F, skip = b-1, na.strings = as.character(info$missing_value),strip.white = T, header = F, keepLeadingZeros = T, tz = info$TZ))

      }
      # get the header line
      headers = as.character(fread(fl,stringsAsFactors = F, skip = n-1, nrows = 1, header = F))
      }

      headers = headers[headers != "NA"]
      colnames(data)[1:length(headers)] = headers

      # remove empty columns
      if(length(which(!is.na(data[,ncol(data)]))) == 0){data = data[,-ncol(data)]}
      # remove empty rows
      data = data[rowSums(matrix(unlist(lapply(as.matrix(data),is.empty)), ncol = ncol(data))) != ncol(data),]
      #data = data[rowSums(!is.na(data) | (data != "")) > 1,]

      # reassign missing value
      for(cl in 1:ncol(data)){
        data[which(data[,cl] == info$missing_value),cl] <- NA}

      ### get header info that is in the data table. This info should be repeated, so just grab the last entry.
      grabbed_vars = info[,..data_vars]
      idx = which(colnames(grabbed_vars) %in% colnames(header_data))
      grabbed_vars = grabbed_vars[,..idx]
      idx = which(grabbed_vars %in% colnames(data))
      grabbed_vars = grabbed_vars[,..idx]
      header_data[,colnames(grabbed_vars)] = data[nrow(data),as.character(grabbed_vars)]


      # 2 data variables to one header variable
      if(any(grepl("-",ID_info))){
        for(gb in which(grepl("-",ID_info))){
          ID_vars = unlist(strsplit(as.character(ID_info[,..gb]), "-"))
          header_data[,colnames(ID_info[,..gb])] = paste(data[1, ID_vars], collapse = "")
        }
      }

      }

      # specifically set to Marlin netcdf format - get header data
      if(info$file_type == "netcdf"){
        f = nc_open(fl)
        heads = out_header_data[(out_header_data %in% data_vars)]
        vars_to_extract = info[,..heads]
          for(vr in 1:ncol(vars_to_extract)){
            if(out_header_data[(out_header_data %in% data_vars)][vr] == "STNNBR"){header_data[,"STNNBR"]  = as.numeric(ncatt_get(f,varid = 0,attname = as.character(vars_to_extract[,..vr]))$value)}
            if(out_header_data[(out_header_data %in% data_vars)][vr] == "DATE"){header_data[,"DATE"]  = ncvar_get(f,as.character(vars_to_extract[,..vr]))[1]}
            if(out_header_data[(out_header_data %in% data_vars)][vr] == "TIME_s"){header_data[,"TIME_s"]  = ncvar_get(f,as.character(vars_to_extract[,..vr]))[1]}
            if(out_header_data[(out_header_data %in% data_vars)][vr] == "TIME_b"){header_data[,"TIME_b"]  = ncvar_get(f,as.character(vars_to_extract[,..vr]))[2]}
            if(out_header_data[(out_header_data %in% data_vars)][vr] == "TIME_e"){header_data[,"TIME_e"]  = ncvar_get(f,as.character(vars_to_extract[,..vr]))[3]}
            if(out_header_data[(out_header_data %in% data_vars)][vr] == "LATITUDE_s"){header_data[,"LATITUDE_s"]  = ncvar_get(f,as.character(vars_to_extract[,..vr]))[1]}
            if(out_header_data[(out_header_data %in% data_vars)][vr] == "LONGITUDE_s"){header_data[,"LONGITUDE_s"]  = ncvar_get(f,as.character(vars_to_extract[,..vr]))[1]}
            if(out_header_data[(out_header_data %in% data_vars)][vr] == "LATITUDE_b"){header_data[,"LATITUDE_b"]  = ncvar_get(f,as.character(vars_to_extract[,..vr]))[2]}
            if(out_header_data[(out_header_data %in% data_vars)][vr] == "LONGITUDE_b"){header_data[,"LONGITUDE_b"]  = ncvar_get(f,as.character(vars_to_extract[,..vr]))[2]}
            if(out_header_data[(out_header_data %in% data_vars)][vr] == "LATITUDE_e"){header_data[,"LATITUDE_e"]  = ncvar_get(f,as.character(vars_to_extract[,..vr]))[3]}
            if(out_header_data[(out_header_data %in% data_vars)][vr] == "LONGITUDE_e"){header_data[,"LONGITUDE_e"]  = ncvar_get(f,as.character(vars_to_extract[,..vr]))[3]}
          }
      }


      ### reformat date and times ###
      # date
      if(info$DATE_format != info$TIME_format){
      if(nchar(header_data$DATE)>12){header_data$DATE = trimws(substr(header_data$DATE,1,12))}}
      if(info$DATE_format == "sec_since_year_start"){
        Y = substr(ex,5,8)
        header_data$DATE = as.POSIXct(paste(Y,"01","01", sep = "-"), tz = info$TZ) + as.numeric(header_data$DATE)

        if(!is.empty(info$TIME_s))
        {date_time = as.POSIXct(paste(Y,"01","01", sep = "-"), tz = info$TZ) + as.numeric(header_data$TIME_s)
        if(attributes(date_time)$tzone != "UTC"){attributes(date_time)$tzone = "UTC"}
        header_data$TIME_s = sub(as.Date(date_time),"",print.POSIXct2(date_time))}
        if(!is.empty(info$TIME_b))
        {
          date_time = as.POSIXct(paste(Y,"01","01", sep = "-"), tz = info$TZ) + as.numeric(header_data$TIME_b)
          if(attributes(date_time)$tzone != "UTC"){attributes(date_time)$tzone = "UTC"}
          header_data$TIME_b = sub(as.Date(date_time),"",print.POSIXct2(date_time))}
        if(!is.empty(info$TIME_e))
        {date_time = as.POSIXct(paste(Y,"01","01", sep = "-"), tz = info$TZ) + as.numeric(header_data$TIME_e)
        if(attributes(date_time)$tzone != "UTC"){attributes(date_time)$tzone = "UTC"}
        header_data$TIME_e = sub(as.Date(date_time),"",print.POSIXct2(date_time))}

        }else{
          if(!is.POSIXct(header_data$DATE[1])){
      header_data$DATE = as.POSIXct(as.character(header_data$DATE),format = info$DATE_format,tz = info$TZ)}

      # times
      if(info$DATE_format != info$TIME_format){

        if(!is.empty(info$TIME_s))
        {dt = paste(header_data$DATE,header_data$TIME_s)
          date_time = as.POSIXct(dt,format = paste("%Y-%m-%d",info$TIME_format),tz = info$TZ)
          if(info$TIME_format == "%H"){date_time = date_time + as.numeric(header_data$TIME_s)%%1}
        if(attributes(date_time)$tzone != "UTC"){attributes(date_time)$tzone = "UTC"}
        header_data$TIME_s = sub(as.Date(date_time),"",print.POSIXct2(date_time))}

        if(!is.empty(info$TIME_b)){
          if(!is.empty(info$TIME_b_format)){
            if(info$TIME_b_format == "seconds"){header_data$TIME_b = as.POSIXct(paste(header_data$DATE,trimws(header_data$TIME_s)), format = "%Y-%m-%d %H:%M:%S", tz = info$TZ) +  as.numeric(header_data$TIME_b)
                header_data$TIME_b = sub(as.Date(header_data$TIME_b),"",header_data$TIME_b)}}else{
                  dt = paste(header_data$DATE,header_data$TIME_b)
                  date_time = as.POSIXct(dt,format = paste("%Y-%m-%d",info$TIME_format),tz = info$TZ)
                  if(info$TIME_format == "%H"){date_time = date_time + as.numeric(header_data$TIME_b)%%1}
                  if(attributes(date_time)$tzone != "UTC"){attributes(date_time)$tzone = "UTC"}
                  header_data$TIME_b = sub(as.Date(date_time),"",print.POSIXct2(date_time))}}

        if(!is.empty(info$TIME_e))
        {dt = paste(header_data$DATE,header_data$TIME_e)
          date_time = as.POSIXct(dt,format = paste("%Y-%m-%d",info$TIME_format),tz = info$TZ)
          if(info$TIME_format == "%H"){date_time = date_time + as.numeric(header_data$TIME_e)%%1}
        if(attributes(date_time)$tzone != "UTC"){attributes(date_time)$tzone = "UTC"}
        header_data$TIME_e = sub(as.Date(date_time),"",print.POSIXct2(date_time))}

      }else{
        if(!is.empty(info$TIME_s) )
        {if(!is.POSIXct(header_data$TIME_s[1])){
          date_time = as.POSIXct(as.character(header_data$TIME_s),format = info$DATE_format,tz = info$TZ)}
        if(attributes(date_time)$tzone != "UTC"){attributes(date_time)$tzone = "UTC"}
        header_data$TIME_s = sub(as.Date(date_time),"",print.POSIXct2(date_time))}
        if(!is.empty(info$TIME_b))
        {if(!is.POSIXct(header_data$TIME_b[1])){date_time = as.POSIXct(as.character(header_data$TIME_b),format = info$DATE_format,tz = info$TZ)}
        if(attributes(date_time)$tzone != "UTC"){attributes(date_time)$tzone = "UTC"}
        header_data$TIME_b = sub(as.Date(date_time),"",print.POSIXct2(date_time))}
        if(!is.empty(info$TIME_e))
        {if(!is.POSIXct(header_data$TIME_e[1])){date_time = as.POSIXct(as.character(header_data$TIME_e),format = info$DATE_format,tz = info$TZ)}
        if(attributes(date_time)$tzone != "UTC"){attributes(date_time)$tzone = "UTC"}
        header_data$TIME_e = sub(as.Date(date_time),"",print.POSIXct2(date_time))}
      }}

      if(attributes(header_data$DATE)$tzone != "UTC"){attributes(header_data$DATE)$tzone = "UTC"}
     header_data$DATE = as.Date(header_data$DATE)

      ### reformat positions ###
      if(info$LATITUDE_s == info$LONGITUDE_s & !is.empty(info$LATITUDE_s)){
        lat_lon = strsplit(header_data$LATITUDE_s, split = " ")
        idx = which(!is.empty(lat_lon[[1]], first.only = F))

        header_data$LATITUDE_s = lat_lon[[1]][idx[1]]
        header_data$LONGITUDE_s = lat_lon[[1]][idx[2]]
      }
      if(info$LATITUDE_b == info$LONGITUDE_b & !is.empty(info$LATITUDE_b)){
        lat_lon = strsplit(header_data$LATITUDE_b, split = " ")
        idx = which(!is.empty(lat_lon[[1]], first.only = F))

        header_data$LATITUDE_b = lat_lon[[1]][idx[1]]
        header_data$LONGITUDE_b = lat_lon[[1]][idx[2]]
      }
      if(info$LATITUDE_e == info$LONGITUDE_e & !is.empty(info$LATITUDE_e)){
        lat_lon = strsplit(header_data$LATITUDE_e, split = " ")
        idx = which(!is.empty(lat_lon[[1]], first.only = F))

        header_data$LATITUDE_e = lat_lon[[1]][idx[1]]
        header_data$LONGITUDE_e = lat_lon[[1]][idx[2]]
      }
      # if not in %deg
      # %deg %min %sec %pos format or %deg %min %pos format
        if(grepl("%pos",info$POSITION_format) & grepl("%min",info$POSITION_format)){
      # grab the delimiter
          delim = substr(sub("%pos","",sub("%deg","",sub("%min","",sub("%sec","",info$POSITION_format)))),1,1)

      ### start pos
      # split based on delimiter
          lat = gsub("[[:alpha:]]","",header_data$LATITUDE_s)
          lat = strsplit(lat,split = delim)
          lat = lat[[1]][which(!is.empty(lat[[1]],first.only = F))]
          lon = gsub("[[:alpha:]]","",header_data$LONGITUDE_s)
          lon = strsplit(lon,split = delim)
          lon = lon[[1]][which(!is.empty(lon[[1]],first.only = F))]
      # turn to numeric
          if(length(lat) == 2){
            lat = as.numeric(lat[1]) + as.numeric(lat[2])/60
            lon = as.numeric(lon[1]) + as.numeric(lon[2])/60
          }
          if(length(lat) == 3){
            lat = as.numeric(lat[1]) + as.numeric(lat[2])/60 + as.numeric(lat[3])/(60*60)
            lon = as.numeric(lon[1]) + as.numeric(lon[2])/60 + as.numeric(lon[3])/(60*60)
          }
      # turn to negative based on %pos
          if(!is.empty( header_data$LATITUDE_s)){


          if(gsub("[^[:alpha:]]+","", header_data$LATITUDE_s) == "S"){header_data$LATITUDE_s = -as.numeric(lat)}else{header_data$LATITUDE_s = as.numeric(lat)}
          if(gsub("[^[:alpha:]]+","", header_data$LONGITUDE_s) == "W"){header_data$LONGITUDE_s = -as.numeric(lon)}else{header_data$LONGITUDE_s = as.numeric(lon)}
          }
      ### bottom pos
      # split based on delimiter
          lat = gsub("[[:alpha:]]","",header_data$LATITUDE_b)
          lat = strsplit(lat,split = delim)
          lat = lat[[1]][which(!is.empty(lat[[1]],first.only = F))]
          lon = gsub("[[:alpha:]]","",header_data$LONGITUDE_b)
          lon = strsplit(lon,split = delim)
          lon = lon[[1]][which(!is.empty(lon[[1]],first.only = F))]
          # turn to numeric
          if(length(lat) == 2){
            lat = as.numeric(lat[1]) + as.numeric(lat[2])/60
            lon = as.numeric(lon[1]) + as.numeric(lon[2])/60
          }
          if(length(lat) == 3){
            lat = as.numeric(lat[1]) + as.numeric(lat[2])/60 + as.numeric(lat[3])/(60*60)
            lon = as.numeric(lon[1]) + as.numeric(lon[2])/60 + as.numeric(lon[3])/(60*60)
          }

          if(!is.empty( header_data$LATITUDE_b)){
      # turn to negative based on %pos


          if(gsub("[^[:alpha:]]+","", header_data$LATITUDE_b) == "S"){header_data$LATITUDE_b = -as.numeric(lat)}else{header_data$LATITUDE_b = as.numeric(lat)}
          if(gsub("[^[:alpha:]]+","", header_data$LONGITUDE_b) == "W"){header_data$LONGITUDE_b = -as.numeric(lon)}else{header_data$LONGITUDE_b = as.numeric(lon)}
          }
      ### finish pos
      # split based on delimiter
          lat = gsub("[[:alpha:]]","",header_data$LATITUDE_e)
          lat = strsplit(lat,split = delim)
          lat = lat[[1]][which(!is.empty(lat[[1]],first.only = F))]
          lon = gsub("[[:alpha:]]","",header_data$LONGITUDE_e)
          lon = strsplit(lon,split = delim)
          lon = lon[[1]][which(!is.empty(lon[[1]],first.only = F))]
          # turn to numeric
          if(length(lat) == 2){
            lat = as.numeric(lat[1]) + as.numeric(lat[2])/60
            lon = as.numeric(lon[1]) + as.numeric(lon[2])/60
          }
          if(length(lat) == 3){
            lat = as.numeric(lat[1]) + as.numeric(lat[2])/60 + as.numeric(lat[3])/(60*60)
            lon = as.numeric(lon[1]) + as.numeric(lon[2])/60 + as.numeric(lon[3])/(60*60)
          }

      # turn to negative based on %pos
          if(!is.empty(header_data$LATITUDE_e)){


          if(gsub("[^[:alpha:]]+","", header_data$LATITUDE_e) == "S"){header_data$LATITUDE_e = -as.numeric(lat)}else{header_data$LATITUDE_e = as.numeric(lat)}
          if(gsub("[^[:alpha:]]+","", header_data$LONGITUDE_e) == "W"){header_data$LONGITUDE_e = -as.numeric(lon)}else{header_data$LONGITUDE_e = as.numeric(lon)}
          }

      }



      # if there is no cast number, set it to 1
      if(is.empty(header_data$CASTNO)){header_data$CASTNO = 1}
      # make sure that the stnnbr is alphanumeric
      header_data$STNNBR = gsub("[^[:alnum:]]","",header_data$STNNBR)
      header_data$STNNBR = stri_replace_all_regex(header_data$STNNBR, "\\b0*(\\d+)\\b", "$1")
      header_data$CASTNO = gsub("[^[:alnum:]]","",header_data$CASTNO)
      header_data$CASTNO = stri_replace_all_regex(header_data$CASTNO, "\\b0*(\\d+)\\b", "$1")

      # # remove repeated character strings
      # # query STNNBR entries and find longest common (forward) substring. Taking the first and last means leading 0's in a count wont be removed
      # # This only deals with leading repitition. e.g TARA01, TARA02...TARA10 becomes 01,02,10
      # # Doesnt deal with internal or end repitition e.g PS01C2, PS02C1...PS10C1 becomes 01C2, 02C2, 10C1
      # sstr <- stri_sub(data2$STNNBR[1], 1, 1:nchar(data2$STNNBR[1]))
      # for(i in 2:nrow(data2)){
      #   #iterively compare strings
      #   sstr <- na.omit(stri_extract_all_coll(data2$STNNBR[i], sstr, simplify=TRUE))
      # }
      # if(length(sstr)>0){
      #   ## match the longest one
      #   lcs = sstr[which.max(nchar(sstr))]
      #   # double check and remove common substring
      #   if(all(grepl(lcs, data2$STNNBR))){data2$STNNBR = gsub(lcs, "", data2$STNNBR)}
      # }



      # if there is no station info, the station number will be assigned as the file number
      if(is.empty(info$STNNBR)){header_data$STNNBR = which(ctd_files == fl)}
      #### Write the new file
      # new file name and path
      new_file = paste(ex,"CTD",header_data$STNNBR,header_data$CASTNO,"ctd1.csv",sep = "_")
      new_file_path = file.path(out_dir,new_file)
      # if the file exists, delete it
      if(file.exists(new_file_path)){unlink(new_file_path)}

      # writing the file

      # file metadata
      subsource = source_info %>% filter(source == info$source)
      cite_tags = c(unlist(strsplit(info$citation,";" )), unlist(strsplit(subsource$citations,";" )))
      cite_tags = cite_tags[is.character(cite_tags)]

      fd <- file(new_file_path, open = "wt")
      writeLines(paste("CTD,",gsub("-","",Sys.Date()),userID, sep = ""), fd)
      writeLines("#semi-manual exchange", fd)
      writeLines(paste("#ORIGINAL_CTDFILE:", basename(old_file)), fd)
      writeLines(paste("#CTDFILE_MOD_DATE:",Sys.time(),"AEST"), fd)
      writeLines(paste("#SOURCED_FROM: ", info$source, "(",subsource$url,")", sep = ""), fd)
      if(is.empty(info$contact)){writeLines(paste("#DATASET_CONTACT:", info$PI),fd)}else{
        writeLines(paste("#DATASET_CONTACT: ", info$PI,"(",info$contact,")", sep = ""),fd)}
      writeLines(paste("#DOI/s:", paste(bib[unlist(strsplit(info$citation,";" ))]$doi, collapse = ",")), fd)
      writeLines(paste("#BIOMATE_CITE_TAGS:", paste(cite_tags, collapse = ",")), fd)

      dcite = format(bib[cite_tags], style = "text", .bibstyle = "BIOMATE")
      dcite = gsub(pattern = "\n", replacement = " ", dcite)
      writeLines(paste("#DATA_CITATION/S:", paste(dcite , collapse = "\n# and")), fd)
      if(!is.null(info$Notes)){writeLines(paste("#NOTE:", info$Notes), fd)}

      # CTD metadata
      writeLines(paste("NUMBER_HEADERS =", 18), fd)
      writeLines(paste("EXPOCODE =", ex), fd)
      writeLines(paste("SHIP =",platforms$`Platform Name`[match(substr(ex,1,4), platforms$`NODC code`)]),fd)
      writeLines(paste("STNNBR or EVENTNBR =", header_data$STNNBR), fd)
      writeLines(paste("CASTNO =", header_data$CASTNO), fd)
      writeLines(paste("CTDID =", paste(ex,"CTD",header_data$STNNBR,header_data$CASTNO, sep = "_")), fd)
      writeLines(paste("DATE =", trimws(header_data$DATE)), fd)
      writeLines(paste("TIMEZONE = UTC"), fd)
      writeLines(paste("CTD_START_TIME =", trimws(sub("UTC","",header_data$TIME_s))),fd)
      writeLines(paste("CTD_START_LATITUDE =", header_data$LATITUDE_s),fd)
      writeLines(paste("CTD_START_LONGITUDE =", header_data$LONGITUDE_s),fd)
      writeLines(paste("CTD_BOTTOM_TIME =", trimws(sub("UTC","",header_data$TIME_b))),fd)
      writeLines(paste("CTD_BOTTOM_LATITUDE =", header_data$LATITUDE_b),fd)
      writeLines(paste("CTD_BOTTOM_LONGITUDE =", header_data$LONGITUDE_b),fd)
      writeLines(paste("CTD_END_TIME =", trimws(sub("UTC","",header_data$TIME_e))),fd)
      writeLines(paste("CTD_END_LATITUDE =", header_data$LATITUDE_e),fd)
      writeLines(paste("CTD_END_LONGITUDE =", header_data$LONGITUDE_e),fd)
      writeLines("missing_value = -999", fd)

      if(info$file_type == "text delim"){
      grabbed_vars = info[,..data_vars]
      idx = which(!colnames(grabbed_vars) %in% colnames(header_data))
      grabbed_vars = grabbed_vars[,..idx]
      idx = which(grabbed_vars %in% colnames(data))
      grabbed_vars = grabbed_vars[,..idx]
      data2 = data[,as.character(grabbed_vars)]
      colnames(data2) = colnames(grabbed_vars)

      # conversions = colnames(info)[grepl("_c",colnames(info))]
      # if(length(conversions)>0){
      #   for(cn in conversions){
      #     if(!is.null(info[,..cn])){if(unlist(strsplit(cn,"_c")) %in% colnames(data2)){data2[,unlist(strsplit(cn,"_c"))] = data2[,unlist(strsplit(cn,"_c"))]*info[,..cn]}}
      #   }}

      }
      # netcdf - marlin format. Create new reformatted data frame
      if(info$file_type == "netcdf"){
        grabbed_vars = info[,..data_vars]
        idx = which(!colnames(grabbed_vars) %in% colnames(header_data))
        grabbed_vars = grabbed_vars[,..idx]
        idx = which(grabbed_vars %in% names(f[['var']]))
        grabbed_vars = grabbed_vars[,..idx]

        if(exists("data2", inherits = FALSE)){rm(data2)}
        for(vr in as.character(grabbed_vars)){
          if(!exists("data2", inherits = FALSE)){data2 = ncvar_get(f,vr)}else{data2 = cbind(data2, ncvar_get(f,vr))}
          }
        data2 = as.data.frame(data2)
        colnames(data2) = colnames(grabbed_vars)
        #conversions
        conversions = colnames(info)[grepl("_c",colnames(info))]
        if(length(conversions)>0){
          for(cn in conversions){
            if(!is.null(info[,..cn])){if(unlist(strsplit(cn,"_c")) %in% colnames(data2)){data2[,unlist(strsplit(cn,"_c"))] = data2[,unlist(strsplit(cn,"_c"))]*info[,..cn]}}
          }}
        nc_close(f)

        for(cl in 1:ncol(data2)){
          data2[which(data2[,cl] == info$missing_value),cl] <- NA}
      }
      data2 = data2[,colnames(data2) != "NA"]
      data2 = data2[,colSums(is.na(data2))<nrow(data2)]
      # downcast data only
      max_p = which.max(data2$CTDPRS)
      data2 = data2[1:max_p,]

      # write the data frame into the file
      writeLines(toString(colnames(data2)), fd)
      heads = paste(colnames(data2),"_u",sep = "")
      writeLines(toString(info[,..heads]), fd)
      data2[is.na(data2)] <- -999
      data2[data2 == "NA"] <- -999
      data2[data2 == ""] <- -999
      for(row in 1:nrow(data2))
      {writeLines(toString(data2[row,]),fd)}
      close(fd)
      rm(data2)
      if(exists("data", inherits = FALSE)){rm(data)}
      if(trace_file){print(fl)}

    }
    print(paste("row",rw,"EXPOCODE",ex,"successfully run"))
  }

}

