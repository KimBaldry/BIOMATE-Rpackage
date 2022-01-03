#' @title Reformat underway sensor data to WHP-Exchange format
#'
#' @author Kimberlee Baldry
#' @description BIOMATE reformats underway sensor files and saves them as a file that is tagged by EXPOCODE. The resulting file follows WHP-Exchange format guidelines. Information is required from the user, stored within a user-entered metadata file particularly for batch runs, or provided via the BIOMATE submission app. The metadata file must be named "UWY_meta.csv". Multiple, sequential underway sensor files can exist for each EXPOCODE, but one file cannot contain data from multiple EXPOCODES. Use the split_delim_file function to convert files if there are multiple voyages in a single file.
#' @note This function only supports well structured text delimited and NetCDF file formats.
#'
#' @return returns a successful run message and saves reformatted data to a "underway_sensors" subdirectory in path_out.
#'
#' @param file_path The path to the processing metadata file, titled "UWY_meta.csv".
#' @param path_out The path where reformatted files are to be sent. A subdirectory "underway_sensors" will be created in here, if one does not already exist.
#' @param userID An ID code for the user creating the files.
#' @param row_start The start row of the processing metadata table, for batch run worflows.
#' @param row_end The end row of the processing metadata table, for batch run worflows.
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


UWY_to_WHPE = function(file_path, path_out,userID = "IMASUTASKB",row_start = 1,row_end = NA, t_thresh = 6, d_thresh = 1000){

  ### output directory ###
  # if the output directory does not exist for this layer, create it
  out_dir = file.path(path_out,"underway_sensors")
  if(!file.exists(out_dir))
  {dir.create(out_dir)}

  ### metadata file ###
  # check that the file path exists
  meta_path = file.path(file_path,"UWY_meta.csv")
  if(!file.exists(meta_path)){stop(paste("Error: There is no file named UWY_meta.csv located in the dirctory",file_path))}
  # read it in
  meta = fread(file = meta_path,header = T,strip.white = T,stringsAsFactors = F)
  if(length(which(!is.na(meta[,ncol(meta)]))) == 0){meta = meta[,-ncol(meta)]}
  if(is.na(row_end)){row_end = nrow(meta)}

  # set up pigment names and the final data frame
  all_headers = c("DATE","TIME","LATITUDE","LONGITUDE")
  pig_names = colnames(meta)[-c(1:(which(colnames(meta) == "CTDSAL")-1))]
  pig_names = pig_names[pig_names != "Notes"] # dont require Notes
  all_headers = c(all_headers, pig_names)

  pig_id_data = c("DATE","TIME",	"LATITUDE",	"LONGITUDE")

  ### reformat files ###
  # loop through the entries in the metadata file
  for(rw in row_start:row_end)
  {
    # relevent metadata for the row rw
    info = meta[rw,]
    ex = info$EXPOCODE
    if(!is.empty(info$header_sep)){ # convert words to punctuation
      if(info$header_sep == "colon"){info$header_sep = ":"}
      if(info$header_sep == "comma"){info$header_sep = ","}
      if(info$header_sep == "dash"){info$header_sep = "-"}
      if(info$header_sep == "equals"){info$header_sep = "="}
      if(info$header_sep == "space"){info$header_sep = " "}}
    if(is.empty(info$header_sep)){info$header_sep = ""}

    # list uwy files
    uwy_files = list.files(path = info$path, pattern = info$extention,full.names = T)

    # list variables that have info in the files
    vars = colnames(info)[-c(1:which(colnames(info) == "start_data_var"))]
    vars = vars[!unlist(lapply(info[,..vars], is.empty))]
    # list variables where data is stored in file headers - if any
    if(any(grepl("header-",info[,..vars]))){header_vars = vars[grepl("header-",info[,..vars])]
    # variables
    data_vars = vars[!grepl("header-",info[,..vars])]}else{data_vars = vars}
    # information required to ID samples
    ID_info = info[,..pig_id_data]
    if(any(grepl("header-",ID_info))){
      rm_head = colnames(ID_info)[!grepl("header-",ID_info)]
      ID_info = ID_info[,..rm_head]}


    # list all of the files in the data stream

    files = list.files(path = info$path, pattern = info$extention,full.names = T,recursive = F)
    if(length(files) == 0){stop(paste("file_path in UWY_meta.csv is wrong:",info$path,"It does not exist or the file is not here"))}

    ## header data ##
    # loop through files to grab header data
    for(fl in files){
      # text files - at the moment this is the only option. No netcdf files were found in search.
      if(info$file_type == "text delim"){

        # get the number of headers that are specified in the meta file to check that enough have been found in the file
        names = data_vars[which(!unlist(lapply(info[,..data_vars],is.empty)))]
        n_headers = length(which(colnames(info[,..names])[!duplicated(as.character(info[,..names]))] %in% all_headers))

        # 2 data variables to one header variable - find less header matches
        if(any(grepl("-",ID_info))){
          n_headers = n_headers - length(which(grepl("-",ID_info)))
        }


        # get total number of lines in the file
        f <- file( fl, open = "r" )
        lines = readLines( f, -1L)
        close(f)
        n_lines = length(lines)

        # read the file to get the header line.
        # Identify the header by crosschecking the names recorded in the metadata information
        f <- file( fl, open = "r" )
        n = 0
        while( TRUE ){
          line <- readLines( f, 1L)
          n = n+1
          # if it is the end of the file?
          if(n == n_lines){
            if(exists("msg")){print(msg)}
            print("Error: Couldnt find header line. Headers must be mislabled in metadata?")
            break}
          if(length(line)== 0){next}
          # searching for data variables in headers
          dv = data_vars[which(!unlist(lapply(info[,..data_vars],is.empty)))]
          header_in_line = unlist(lapply(dv,function(x){grepl(info[,..x],line, fixed = T)}))

          if(length(which(header_in_line == T) ) < n_headers & length(which(header_in_line == T) ) > 2){missing = dv[!header_in_line]
          missing = missing[missing %in% pig_names[pig_names != "Notes" & !grepl("_u", pig_names)]]
          msg = paste("Header assignments might not be right. Missing underway headers in",fl,"for",ex, "are:", paste(missing, collapse = ","))
          }
          # if there are at least three data_vars that apear, break this is likely the header line
          if(length(which(header_in_line == T) ) >= n_headers  & !substr(line,1,1) %in% c("!","#")){
            # write to a file to check header line is picked and variables have been grabbed. This can be removed. Maybe put a switch on for debugging.
            tf = file(file.path(out_dir,paste("Check_headers",Sys.Date(),".txt", sep = "")),open = "at")
            writeLines(paste(fl),tf)
            writeLines(paste("The header line picked:",line),tf)
            writeLines(paste("The variables found:",data_vars[which(header_in_line ==T)]),tf)
            writeLines("All good?",tf)
            writeLines(" ",tf)
            close(tf)
            header_line = line

            break}

        }
        close( f )

        # get number of the header line to skip for data table
        n = grep(header_line,lines, fixed = T)
        if(length(n) == 0){ n = 1}
        # check for extra lines before the data table starts.
        # There should be some numeric data in the table, but not in the header or unit lines.
        b = n+1
        line = fread(fl,stringsAsFactors = F, skip = n, nrows = 1, header = F)
        while(!any(unlist(lapply(line,is.numeric)))){
          line = fread(fl,stringsAsFactors = F, skip = b, nrows = 1, header = F)
          b = b+1}

        ### Get data table ###
        # Use fread to read the data table - this will adapt if there is a units line or not.

        # catch the out-of-format lines here for debugging
        #tryCatch({data = as.data.frame(fread(fl,stringsAsFactors = F, skip = n, na.strings = info$missing_value,strip.white = T , header = F))}, warning=function(w) print(fl))

        # rectangular data shouldnt be read with fread
        if(info$delim == "rect"){
          data = read_table(fl,col_names = F, skip = b-1, na = info$missing_value,col_types = cols())
        }else{
          if(is.na(info$missing_value)){data = as.data.frame(fread(fl,stringsAsFactors = F, skip = b-1,strip.white = T, header = F, keepLeadingZeros = T))
          if(any(grepl("POSIXct",sapply(data,class)))){
            data = as.data.frame(fread(fl,stringsAsFactors = F, skip = b-1,strip.white = T, header = F, keepLeadingZeros = T, colClasses = list(character = grep("POSIXct",sapply(data,class)))))
          }
          }else{
            data = as.data.frame(fread(fl,stringsAsFactors = F, skip = b-1, na.strings = as.character(info$missing_value),strip.white = T, header = F, keepLeadingZeros = T))
            if(any(grepl("POSIXct",sapply(data,class)))){
              data = as.data.frame(fread(fl,stringsAsFactors = F, skip = b-1, na.strings = as.character(info$missing_value),strip.white = T, header = F, keepLeadingZeros = T, colClasses = list(character = grep("POSIXct",sapply(data,class)))))
              }
          }
        }

        # get the header line
        headers = as.character(fread(fl,stringsAsFactors = F, skip = n-1, nrows = 1, header = F))
        headers = gsub("/fields=","",headers) # in seabass format
        headers = headers[headers != "NA"]
        colnames(data)[1:length(headers)] = headers

        # remove empty columns
        if(length(which(!is.na(data[,ncol(data)]))) == 0){data = data[,-ncol(data)]}
        #data = data[rowSums(!is.na(data) | (data != "")) > 1,]
        # # remove empty rows dont do this with underway data - data frames are too big!
        # data = data[rowSums(matrix(unlist(lapply(as.matrix(data),is.empty)), ncol = ncol(data))) != ncol(data),]

        # reassign missing value
        for(cl in 1:ncol(data)){
          data[which(data[,cl] == info$missing_value),cl] <- NA}

        # get the variables listed in the metadata file and create a new reformatted data frame
        grabbed_vars = info[,..data_vars]
        idx = which(grabbed_vars %in% colnames(data))
        grabbed_vars = grabbed_vars[,..idx]
        colnames(grabbed_vars) = colnames(info[,..data_vars])[info[,..data_vars] %in% colnames(data)]
        dataf = data.frame(data[,as.character(grabbed_vars)],stringsAsFactors = F)
        colnames(dataf) = colnames(grabbed_vars)


        # 2 data variables to one header variable
         if(any(grepl("-",ID_info))){
          for(gb in which(grepl("-",ID_info))){
            ID_vars = unlist(strsplit(as.character(ID_info[,..gb]), "-"))
            dataf[,colnames(ID_info[,..gb])] = apply(data[, ID_vars],1, paste, collapse = "")
          }
        }
        ### date must be entered in as year-month-day
        if(grepl("-", info$DATE) & !grepl("header-", info$DATE)){
          date_vars = unlist(strsplit(as.character(info$DATE), "-"))
          dataf$DATE = sprintf("%04d-%02d-%02d", as.numeric(data[,date_vars[1]]), as.numeric(data[,date_vars[2]]), as.numeric(data[,date_vars[3]]))
          info$DATE_format = "%Y-%m-%d"
        }

        ### time must be hour-min or hour-min-sec
        if(grepl("-", info$TIME) & !grepl("header-", info$TIME)){
          time_vars = unlist(strsplit(as.character(info$TIME), "-"))
          if(length(time_vars) == 3){
            dataf$TIME = sprintf("%02d:%02d:%02d",as.numeric(data[,time_vars[1]]), as.numeric(data[,time_vars[2]]), as.numeric(data[,time_vars[3]]))
            info$TIME_format = "%H:%M:%S"}
          if(length(time_vars) == 2){
            dataf$TIME = sprintf("%02d:%02d",as.numeric(data[,time_vars[1]]), as.numeric(data[,time_vars[2]]))
            info$TIME_format = "%H:%M"}
        }


        f <- file( fl, open = "r" )

        # read header lines and get stored data to add to the new reformatted data table
        n2 = 0
        while( TRUE ){ # stops when at the end of the header lines
          n2 = n2+1
          if(n2 > n){break}
          line <- readLines( f, 1L )
          if(length(line)== 0){next}

          # get header data
          if(exists("header_vars")){

            # for each header variable store the associated data
            for(hd in header_vars){
              head_pattern = sub("header-","",info[,..hd])
              if(grepl(head_pattern, line )){
                dataf[,hd] = trimws(sub(info$header_sep,"",sub(head_pattern,"",line)))}
            }

          }

        }
        close( f )


        # append data frame
        if(fl == files[1]){
          data_appended = dataf
        }else{
          data_appended = rbind(data_appended,dataf)
          }

      }
      rm(data)
    }

    # create a full data frame, with consistant organisation
    data2 = data.frame(matrix(NA,nrow = nrow(data_appended),ncol = length(all_headers)))
    colnames(data2) = all_headers
    data2[,colnames(data_appended)] = data_appended

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
    #


    ### reformat date and time
    if(!is.empty(info$DATE) & !all(is.na(data2$DATE))){
      if(nchar(data2$DATE[which(!is.empty(data2$DATE))[1]])>12 & grepl("AADC",info$source)){data2$DATE = substr(data2$DATE,1,12)}} #problems with old AADC data
    if(info$DATE_format == "sec_since_year_start"){
      Y = substr(ex,5,8)
      if(!is.empty(info$DATE)){data2$DATE = as.POSIXct(paste(Y,"01","01", sep = "-"), tz = info$TZ) + as.numeric(data2$DATE)
      }

      if(!is.empty(info$TIME))
      {data2$TIME = as.POSIXct(paste(Y,"01","01", sep = "-"), tz = info$TZ) + as.numeric(data2$TIME)
      if(attributes(data2$TIME)$tzone != "UTC"){attributes(data2$TIME)$tzone = "UTC"}}


    }else{
      if(!is.empty(info$DATE) & !is.POSIXct(data2$DATE[1])){data2$DATE = as.Date(as.character(data2$DATE),format = info$DATE_format)
      }

      # times
      if(info$DATE_format != info$TIME_format){
        if(!is.empty(info$TIME))
        {dt = paste(data2$DATE, data2$TIME)
        data2$TIME = as.POSIXct(as.character(dt),format = paste("%Y-%m-%d",info$TIME_format),tz = info$TZ)
        if(info$TIME_format == "%H"){data2$TIME = data2$TIME + as.numeric(data2$TIME)%%1}
        if(attributes(data2$TIME)$tzone != "UTC"){attributes(data2$TIME)$tzone = "UTC"}
        }

      }else{
        if(!is.empty(info$TIME)  )
        {data2$TIME = as.POSIXct(as.character(data2$TIME),format = info$TIME_format,tz = info$TZ)
        if(attributes(data2$TIME)$tzone != "UTC"){attributes(data2$TIME)$tzone = "UTC"}}
      }}

    if(!is.empty(info$TIME)){

      for(dr in 1:nrow(data2)){
        data2$TIME2[dr] = sub(as.character(as.Date(data2$TIME[dr])),"",print.POSIXct2(data2$TIME[dr]))
        data2$TIME2[dr] = trimws(sub("UTC","",data2$TIME2[dr]))
      }
      data2$TIME = data2$TIME2
      data2 = data2[,-which(colnames(data2) == "TIME2")]
    }
    rm(dataf,data_appended)

    ### reformat position ###
    # positions
    if(info$LATITUDE == info$LONGITUDE & !is.empty(info$LATITUDE)){
      lat_lon = strsplit(data2$LATITUDE, split = " ")
      data2$LATITUDE = lat_lon[[1]][1]
      data2$LONGITUDE = lat_lon[[1]][2]
    }
    # if not in %deg
    # %deg %min %sec %pos format or %deg %min %pos format
    if(grepl("%pos",info$POSITION_format) & grepl("%min",info$POSITION_format)){
      # grab the delimiter
      delim = substr(sub("%pos","",sub("%deg","",sub("%min","",sub("%sec","",info$POSITION_format)))),1,1)

      ### start pos
      # split based on delimiter
      lat = gsub("[[:alpha:]]","",data2$LATITUDE)
      lat = strsplit(lat,split = delim)
      lon = gsub("[[:alpha:]]","",data2$LONGITUDE)
      lon = strsplit(lon,split = delim)
      # turn to numeric
      if(length(lat[[1]]) == 2){
        lat = as.numeric(lat[[1]][1]) + as.numeric(lat[[1]][2])/60
        lon = as.numeric(lon[[1]][1]) + as.numeric(lon[[1]][2])/60
      }
      if(length(lat[[1]]) == 3){
        lat = as.numeric(lat[[1]][1]) + as.numeric(lat[[1]][2])/60 + as.numeric(lat[[1]][3])/(60*60)
        lon = as.numeric(lon[[1]][1]) + as.numeric(lon[[1]][2])/60 + as.numeric(lon[[1]][3])/(60*60)
      }
      # turn to negative based on %pos
      data2$LATITUDE = as.numeric(lat)
      data2$LONGITUDE = as.numeric(lon)
      if(gsub("[^[:alpha:]]+","", data2$LATITUDE) == "S"){data2$LATITDUE = -as.numeric(lat)}
      if(gsub("[^[:alpha:]]+","", data2$LONGITUDE) == "W"){data2$LONGITUDE = -as.numeric(lon)}
    }

    if(info$POSITION_format != "%deg" & grepl("%deg",info$POSITION_format) & length(strsplit(info$POSITION_format,split = "%")[[1]]) == 2 ){
      for(dr in 1:nrow(data2)){
        data2$LATITUDE[dr] = sub(paste("\\",sub("%deg","",info$POSITION_format),sep = ""),"" ,data2$LATITUDE[dr])
        data2$LONGITUDE[dr] = sub(paste("\\",sub("%deg","",info$POSITION_format),sep = ""),"" ,data2$LONGITUDE[dr])}
    }
    data2$LATITUDE = as.numeric(data2$LATITUDE)
    data2$LONGITUDE = as.numeric(data2$LONGITUDE)

    #### Write the new file
    # new file name and path
    new_file = paste(paste(ex,"UWY",sep = "_"),".csv",sep = "")

    new_file_path = file.path(out_dir,new_file)
    # if the file exists, delete it
    if(file.exists(new_file_path)){unlink(new_file_path)}

    # if(exists("sub_meta")){
    #   files = list()
    #     for(sb in 1:nrow(sub_meta)){
    #     files = c(files, list.files(path = sub_meta[sb,]$path, pattern = sub_meta[sb,]$extention,full.names = T,recursive = F))}
    #   files = unlist(files)
    #   }else{
    # files = list.files(path = info$path, pattern = info$extention,full.names = T,recursive = F)}
    subsource = source_info %>% filter(source == info$source)
    cite_tags = c(unlist(strsplit(info$citation,";" )), unlist(strsplit(subsource$citations,";" )))
    cite_tags = cite_tags[is.character(cite_tags)]

    # reformat and write new file
    fd <- file(new_file_path, open = "wt")
    writeLines(paste(info$Data_type,gsub("-","",Sys.Date()),userID, sep = ""), fd)
    writeLines("#semi-manual exchange", fd)
    writeLines(paste("#ORIGINAL_UWYFILE/S:", paste(files, collapse="; ")), fd)
    writeLines(paste("#UWYFILE_MOD_DATE:",Sys.time(),"AEST"), fd)
    writeLines(paste("#SOURCED_FROM: ", info$source, "(",subsource$url,")", sep = ""), fd)

    if(is.empty(info$contact)){writeLines(paste("#DATASET_CONTACT:", info$PI),fd)}else{
      writeLines(paste("#DATASET_CONTACT: ", info$PI,"(",info$contact,")", sep = ""),fd)}

    writeLines(paste("#DOI/s:", paste(bib[trimws(unlist(strsplit(info$citation,";" )))]$doi, collapse = ",")), fd)
    writeLines(paste("#BIOMATE_CITE_TAGS:", paste(cite_tags, collapse = ",")), fd)

    dcite = format(bib[cite_tags], style = "text", .bibstyle = "BIOMATE")
    dcite = gsub(pattern = "\n", replacement = " ", dcite)
    writeLines(paste("#DATA_CITATION/S:", paste(dcite , collapse = "\n# and ")), fd)

    if(!is.empty(info$Notes)){writeLines(paste("#NOTE:", info$Notes), fd)}

     final_data = data2
     ######### NEED TO fix units
    ordered_names = c("DATE","TIME","LATITUDE","LONGITUDE",pig_names)
    ordered_units = c("YYYY-mm-dd","HH:MM:SS","DEGREES NORTH","DEGREES EAST")

    final_data = final_data[,ordered_names]
    rm(data2)

    writeLines(paste("NUMBER_HEADERS = 6"), fd)
    writeLines(paste("EXPOCODE =", ex), fd)
    writeLines(paste("SHIP =",platforms$`Platform Name`[match(substr(ex,1,4), platforms$`NODC code`)]),fd)
    writeLines(paste("TIMEZONE = UTC"), fd)
    writeLines("missing_value = -999", fd)
    writeLines("not_detected = -888", fd)


    # assign missing value number
    final_data$DATE = as.character(final_data$DATE)

    # remove rows with no pigment data
    final_data = final_data[rowSums(matrix(unlist(lapply(as.matrix(final_data[,pig_names]),is.empty)), ncol = length(pig_names))) != length(pig_names),]


    # remove columns with entire missing values
    for(cl in ncol(final_data):19){
      if(all(is.na(final_data[,cl]))){colnames(final_data)[cl] = NA
      ordered_units = ordered_units[-cl]
      }
    }
    final_data = final_data[,which(!is.na(colnames(final_data)))]
    writeLines(toString(colnames(final_data)), fd)
    heads = paste(colnames(final_data),"_u",sep = "")
    writeLines(toString(info[,..heads]), fd)
    units = lapply(colnames(final_data),function(x){if(paste(x,"_u",sep = "") %in% colnames(info)){info[,paste(x,"_u",sep = "")]}else{NA}})
    writeLines(toString(c(ordered_units,units)), fd)

    final_data[is.na(final_data)] <- -999
    final_data[final_data == "NA"] <- -999
    # reassign below detection limits value
    if(!is.empty(info$not_detected_value)){
      for(col in pig_names){
        cdx = which(colnames(final_data) == col)
        if(any(final_data[,cdx] == info$not_detected_value)){final_data[which(final_data[,cdx] == info$not_detected_value),cdx] = -888}
      }}
    for(row in 1:nrow(final_data))
    {writeLines(toString(final_data[row,]),fd)}
    close(fd)

    print(paste("row",rw,"EXPOCODE",ex,"successfully run"))
    ### checks
    # final_data$DEPTH = as.numeric(final_data$DEPTH)
    # final_data$DEPTH[final_data$CTD_IDs == "U"] = "U"
     #if(any(final_data$LATITUDE == -999 & final_data$STNNBR != "U" & final_data$DEPTH > 10 & !is.empty(final_data$DEPTH))){print(paste("Warning:There is some ctd data not matched, or missing",rw,ex,"for this many measurements", length(which(final_data$LATITUDE == -999 & final_data$STNNBR != "U" & final_data$DEPTH > 10 & !is.empty(final_data$DEPTH)))))}
    rm(final_data)
    if(exists("header_vars")){rm(header_vars)}
    if(exists("CTD_IDs")){rm(CTD_IDs)}
    if(exists("sub_meta")){rm(sub_meta)}
    rm(data_vars,header_in_line,header_line,headers,line)
  }

}

