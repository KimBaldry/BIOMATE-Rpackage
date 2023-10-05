#' @title Reformat particulate organic carbon (poc) data to WHP-Exchange format with MAREDAT headers
#'
#' @author Kimberlee Baldry
#' @description BIOMATE reformats particulate organic carbon (poc) files and saves them as a file that is tagged by EXPOCODE, METHOD and SOURCE. The resulting file follows WHP-Exchange format guidelines with the exception that MAREDAT pigment headers are used. Information is required from the user, stored within a user-entered metadata file particularly for batch runs, or provided via the BIOMATE submission app. The metadata file must be named "POC_meta.csv". All poc files must be split by EXPOCODE (can be multiple files in the same folder with identical extensions). Use the split_delim_file function to convert files if there are multiple voyages in a single file.
#' @note This function only supports well structured text delimited and NetCDF file formats. It is almost identical yo PIG_to_WHPE, despite alterations to pigment/poc header names, file paths and error messages.
#'
#' @return returns a successful run message and saves reformatted data to a "poc" subdirectory in path_out.
#'
#' @param file_path The path to the processing metadata file, titled "POC_meta.csv".
#' @param path_out The path where reformatted files are to be sent. A subdirectory "poc" will be created in here, if one does not already exist.
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


POC_to_WHPE = function(file_path, path_out,userID = "IMASUTASKB",row_start = 1,row_end = NA, t_thresh = 6, d_thresh = 1000, version = "1.0"){

  # this small function prevents the drop of midnight 00:00:00
  print.POSIXct2 <- function(x){format(x,"%Y-%m-%d %H:%M:%S %Z")}
  calc.dist.min = function(x){
    mins = which.min(distGeo(c(unmatched_df_p$LON[x], unmatched_df_p$LAT[x]), as.matrix(CTD_df[,c("LON_s","LAT_s")])))[1]
    dists = distGeo(c(unmatched_df_p$LON[x], unmatched_df_p$LAT[x]), as.matrix(CTD_df[mins,c("LON_s","LAT_s")]))
    if(any(!is.na(CTD_df$LON_b))){
      minb = which.min(distGeo(c(unmatched_df_p$LON[x], unmatched_df_p$LAT[x]), as.matrix(CTD_df[,c("LON_b","LAT_b")])))[1]
      mins = c(mins, minb)
      distb = distGeo(c(unmatched_df_p$LON[x], unmatched_df_p$LAT[x]), as.matrix(CTD_df[minb,c("LON_b","LAT_b")]))
      dists = c(dists, distb)}
    if(any(!is.na(CTD_df$LON_e))){
      mine = which.min(distGeo(c(unmatched_df_p$LON[x], unmatched_df_p$LAT[x]), as.matrix(CTD_df[,c("LON_e","LAT_e")])))[1]
      mins = c(mins, mine)
      diste = distGeo(c(unmatched_df_p$LON[x], unmatched_df_p$LAT[x]), as.matrix(CTD_df[mine,c("LON_e","LAT_e")]))
      dists = c(dists, diste)}
    data.frame(min = mins[which.min(dists)],diff = dists[which.min(dists)])
  }
  #calc.dist.diff = function(x){distGeo(c(unmatched_df_p$LON[x], unmatched_df_p$LAT[x]), as.matrix(CTD_df[unmatched_df_p$dist_min[x],c("LON","LAT")]))}
  calc.time.min = function(x){which.min(abs(unlist(lapply(col_t, FUN = function(y){difftime(CTD_df[,y], unmatched_df_t$time_m[x], units = "secs")}))))[1]}
  calc.time.diff =  function(x){abs(unlist(lapply(col_t, FUN = function(y){difftime(CTD_df[,y], unmatched_df_t$time_m[x], units = "secs")})))[unmatched_df_t$closest_t[x]]}

  ### directories ###
  # if the output directory does not exist for this layer, create it
  out_dir = file.path(path_out,"poc")
  if(!file.exists(out_dir))
  {dir.create(out_dir)}

  # check there is a link to a profiling_sensor layer
  ctd_path = file.path(path_out,"profiling_sensors")
  if(!file.exists(ctd_path))
  {stop(paste(ctd_path, "does not exist. Please link to the profiling sensors layer"))}

  #to link to UWY - to be implemented ***
  #uwy_path = file.path(path_out,"underway_sensors")
  #if(!file.exists(uwy_path))
  #{stop(paste(ctd_path, "does not exist. Please link to the underway sensors layer"))}

  ### metadata file ###
  # check that the file path exists
  meta_path = file.path(file_path,"POC_meta.csv")
  if(!file.exists(meta_path)){stop(paste("Error: There is no file named POC_meta.csv located in the dirctory",file_path))}
  # read it in
  meta = fread(file = meta_path,header = T,strip.white = T,stringsAsFactors = F)
  if(length(which(!is.na(meta[,ncol(meta)]))) == 0){meta = meta[,-ncol(meta)]}
  if(is.na(row_end)){row_end = nrow(meta)}

  # set up poc names and the final data frame
  all_headers = c("CTD_match","DATE","TIME_s","TIME_b","TIME_e","LATITUDE","LONGITUDE","STNNBR",	"CASTNO",	"DATE_analyser",	"TIME_analyser",	"LAT_analyser",	"LON_analyser","Sample_ID",	"BOTTLE",	"DEPTH")
  pig_names = colnames(meta)[-c(1:which(colnames(meta) == "PON_u"))]
  pig_names = pig_names[pig_names != "Notes"] # dont require Notes
  all_headers = c(all_headers, pig_names)
  pig_id_data = c("STNNBR",	"CASTNO","DATE_analyser","TIME_analyser",	"LAT_analyser",	"LON_analyser")

  ### reformat files ###
  # loop through the entries in the metadata file
  for(rw in row_start:row_end)
  {
    # relevent metadata for the row rw
    info = meta[rw,]
    if(nrow(meta %>% filter(source == info$source, Method == info$Method, EXPOCODE == info$EXPOCODE))>1){
      sub_meta = meta %>% filter(source == info$source, Method == info$Method, EXPOCODE == info$EXPOCODE)
      duplicate = T
      dup_id = which(unlist(lapply(1:nrow(sub_meta), function(x){r = all_equal(info[,1:11], sub_meta[x,1:11])
      if(r != T){r = F}
      return(r)})))
    }else{duplicate = F}

    ex = info$EXPOCODE
    if(!is.empty(info$header_sep)){ # convert words to punctuation
      if(info$header_sep == "colon"){info$header_sep = ":"}
      if(info$header_sep == "comma"){info$header_sep = ","}
      if(info$header_sep == "dash"){info$header_sep = "-"}
      if(info$header_sep == "equals"){info$header_sep = "="}
      if(info$header_sep == "space"){info$header_sep = " "}}
    if(is.empty(info$header_sep)){info$header_sep = ""}

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
    if(length(files) == 0){stop(paste("file_path in POC_meta.csv is wrong:",info$path,"It does not exist or the file is not here"))}

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
        #if(info$DATE_analyser =="year-month-day"){n_headers = n_headers - 1
        #if(info$TIME_analyser == "hour-minute-second" | info$TIME_analyser == "hour-minute"){n_headers = n_headers - 1}} # these wont be matched


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
          missing = missing[missing %in% pig_names[pig_names != "Notes"]]
          msg = paste("Header assignments might not be right. Missing pigment headers in",fl,"for",ex, "are:", paste(missing, collapse = ","))
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
        headers = gsub("/fields=","",headers)
        headers = headers[headers != "NA"]
        colnames(data)[1:length(headers)] = headers

        # remove empty columns
        if(length(which(!is.na(data[,ncol(data)]))) == 0){data = data[,-ncol(data)]}
        #data = data[rowSums(!is.na(data) | (data != "")) > 1,]
        # remove empty rows
        data = data[rowSums(matrix(unlist(lapply(as.matrix(data),is.empty)), ncol = ncol(data))) != ncol(data),]

        # reassign missing value
        for(cl in 1:ncol(data)){
          data[which(data[,cl] == info$missing_value),cl] <- NA}
        # identify underway data
        if(!is.empty(info$Underway_ID) & info$Underway_ID != "all"){
          uwy_var = strsplit(info$Underway_ID,split = "-")[[1]][1]
          uwy_val = sub(paste(uwy_var,"-",sep = ""),"",info$Underway_ID)
          uwy_ref = data[,uwy_var]}

        # get the variables listed in the metadata file and create a new reformatted data frame
        grabbed_vars = info[,..data_vars]
        idx = which(grabbed_vars %in% colnames(data))
        grabbed_vars = grabbed_vars[,..idx]
        colnames(grabbed_vars) = colnames(info[,..data_vars])[info[,..data_vars] %in% colnames(data)]
        dataf = data.frame(data[,as.character(grabbed_vars)],stringsAsFactors = F)
        colnames(dataf) = colnames(grabbed_vars)


        # 2 data variables to one header variable
        ID_info = ID_info[,-c("DATE_analyser", "TIME_analyser")]
        if(any(grepl("-",ID_info))){
          for(gb in which(grepl("-",ID_info))){
            ID_vars = unlist(strsplit(as.character(ID_info[,..gb]), "-"))
            dataf[,colnames(ID_info[,..gb])] = apply(data[, ID_vars],1, paste, collapse = "")
          }
        }
        ### date must be entered in as year-month-day
        if(grepl("-", info$DATE_analyser) & !grepl("header-", info$DATE_analyser)){
          date_vars = unlist(strsplit(as.character(info$DATE_analyser), "-"))
          dataf$DATE_analyser = sprintf("%04d-%02d-%02d", as.numeric(data[,date_vars[1]]), as.numeric(data[,date_vars[2]]), as.numeric(data[,date_vars[3]]))
          info$DATE_analyser_format = "%Y-%m-%d"
        }

        ### time must be hour-min or hour-min-sec
        if(grepl("-", info$TIME_analyser) & !grepl("header-", info$TIME_analyser)){
          time_vars = unlist(strsplit(as.character(info$TIME_analyser), "-"))
          if(length(time_vars) == 3){
            dataf$TIME_analyser = sprintf("%02d:%02d:%02d",as.numeric(data[,time_vars[1]]), as.numeric(data[,time_vars[2]]), as.numeric(data[,time_vars[3]]))
            info$TIME_analyser_format = "%H:%M:%S"}
          if(length(time_vars) == 2){
            dataf$TIME_analyser = sprintf("%02d:%02d",as.numeric(data[,time_vars[1]]), as.numeric(data[,time_vars[2]]))
            info$TIME_analyser_format = "%H:%M"}
        }

        # unit conversions
        conversions = colnames(info)[grepl("_c",colnames(info))]
        if(info$POC_u == "NG/L"){
          idx = colnames(dataf) %in% pig_names
          dataf[,idx] = dataf[,idx]*0.001
          # old code that specifies conversion in meta file
          #for(cn in conversions){
          #  if(!is.empty(info[,cn])){if(unlist(strsplit(cn,"_c")) %in% colnames(dataf)){dataf[,unlist(strsplit(cn,"_c"))] = dataf[,unlist(strsplit(cn,"_c"))]*info[,cn]}}
        }
        if(info$POC_u == "MG/L"){
          idx = colnames(dataf) %in% pig_names
          dataf[,idx] = dataf[,idx]*1000
        }
        if(info$POC_u == "UMOL/L"){
          if(any(colnames(dataf) == "POC")){
          dataf$POC =  dataf$POC*12.0107}
          if(any(colnames(dataf) == "DOC")){
          dataf$DOC =  dataf$DOC*12.0107}
          if(any(colnames(dataf) == "PON")){
          dataf$PON =  dataf$PON*28.02}
          if(any(colnames(dataf) == "DON")){
          dataf$DON =  dataf$DON*28.02}
          if(any(colnames(dataf) == "TPP")){
          dataf$TPP =  dataf$TPP*30.97}
        }
        if(!is.empty(info$PON_u) & info$POC_u %in% c("MG/M3",c("UG/L"))){
          dataf$PON = dataf$PON*28.02
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
          if(exists("uwy_ref")){uwy_data = uwy_ref}
        }else{
          data_appended = rbind(data_appended,dataf)
          if(exists("uwy_ref")){uwy_data = c(uwy_data,uwy_ref)}}

        #print(paste(fl,"read"))
        if(exists("uwy_ref")){rm(uwy_ref)}
      }
      rm(data)
    }

    # create a full data frame, with consistant organisation
    data2 = data.frame(matrix(NA,nrow = nrow(data_appended),ncol = length(all_headers)))
    colnames(data2) = all_headers
    data2[,colnames(data_appended)] = data_appended

    ### CASTNO and STNNBR refoematting ###
    if(!is.na(info$CASTNO) & !is.na(info$STNNBR)){
      if(info$CASTNO == info$STNNBR & !is.empty(info$STNNBR)){ # if stnnbr and castno are in the same line split the information. At the moment only works with castno and stnnbrs separated with a [.] SEABASS format. change to split by punctuation?
        data2$CASTNO = strsplit(as.character(data2$CASTNO),split = "[.]")[[1]][2]
        data2$STNNBR = strsplit(as.character(data2$STNNBR),split = "[.]")[[1]][1]}}
    # remove all special symbols and convert to an alphanumeric code
    data2$STNNBR = gsub("[^[:alnum:]]","",data2$STNNBR)
    data2$STNNBR = stri_replace_all_regex(data2$STNNBR, "\\b0*(\\d+)\\b", "$1")
    data2$CASTNO = gsub("[^[:alnum:]]","",data2$CASTNO)
    data2$CASTNO = stri_replace_all_regex(data2$CASTNO, "\\b0*(\\d+)\\b", "$1")

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

    ### Underway data assignment ### This process seems round-about. Can Mike improve it?
    if(!is.empty(info$Underway_ID)){
      if(info$Underway_ID == "all"){data2$STNNBR = "U"}else{
        if(length(grep(uwy_val,uwy_data))>0){
          data2$STNNBR[grep(uwy_val,uwy_data)] = "U"}}}
    if(exists("uwy_data")){rm(uwy_data)}

    # if station numbers are not indicated in the file assign unique numbers based on lat/lon. if underway indicated dont replace.
    if(is.empty(info$STNNBR)){
      idx = is.na(data2$STNNBR)
      temp = data2[idx,] %>% mutate(STNNBR = group_indices(., LAT_analyser, LON_analyser))
      data2$STNNBR[idx] = paste("chl",temp$STNNBR, sep = "_")
      rm(temp)
    }

    # replace missing cast numbers with 1
    if(length(which(data2$STNNBR != "U" & !is.na(data2$STNNBR) & is.na(data2$CASTNO)))>0){
      data2$CASTNO[which(data2$STNNBR != "U" & !is.na(data2$STNNBR) & is.na(data2$CASTNO))] = 1}

    ### reformat date and time
    if(!is.empty(info$DATE_analyser) & !all(is.na(data2$DATE_analyser))){
      if(nchar(data2$DATE_analyser[which(!is.empty(data2$DATE_analyser))[1]])>12 & grepl("AADC",info$source)){data2$DATE_analyser = substr(data2$DATE_analyser,1,12)}} #problems with old AADC data
    if(info$DATE_analyser_format == "sec_since_year_start"){
      Y = substr(ex,5,8)
      if(!is.empty(info$DATE_analyser)){data2$DATE_analyser = as.POSIXct(paste(Y,"01","01", sep = "-"), tz = info$TZ) + as.numeric(data2$DATE_analyser)
      }

      if(!is.empty(info$TIME_analyser))
      {data2$TIME_analyser = as.POSIXct(paste(Y,"01","01", sep = "-"), tz = info$TZ) + as.numeric(data2$TIME_analyser)
      if(attributes(data2$TIME_analyser)$tzone != "UTC"){attributes(data2$TIME_analyser)$tzone = "UTC"}}


    }else{
      if(!is.empty(info$DATE_analyser) & !is.POSIXct(data2$DATE_analyser[1])){data2$DATE_analyser = as.Date(as.character(data2$DATE_analyser),format = info$DATE_analyser_format)
      }

      # times
      if(info$DATE_analyser_format != info$TIME_analyser_format){
        if(!is.empty(info$TIME_analyser))
        {dt = paste(data2$DATE_analyser, data2$TIME_analyser)
        data2$TIME_analyser = as.POSIXct(as.character(dt),format = paste("%Y-%m-%d",info$TIME_analyser_format),tz = info$TZ)
        if(info$TIME_analyser_format == "%H"){data2$TIME_analyser = data2$TIME_analyser + as.numeric(data2$TIME_analyser)%%1}
        if(attributes(data2$TIME_analyser)$tzone != "UTC"){attributes(data2$TIME_analyser)$tzone = "UTC"}
        }

      }else{
        if(!is.empty(info$TIME_analyser)  )
        {data2$TIME_analyser = as.POSIXct(as.character(data2$TIME_analyser),format = info$TIME_analyser_format,tz = info$TZ)
        if(attributes(data2$TIME_analyser)$tzone != "UTC"){attributes(data2$TIME_analyser)$tzone = "UTC"}}
      }}

    # if(!is.empty(info$DATE_analyser) & info$DATE_analyser_format != "year-month-day"){
    # data2$DATE_analyser = as.Date(data2$DATE_analyser,format = info$DATE_analyser_format,tz = info$TZ)
    #   if(attributes(data2$DATE_analyser)$tzone != "UTC"){attributes(data2$DATE_analyser)$tzone = "UTC"}
    # }
    if(!is.empty(info$TIME_analyser)){

      for(dr in 1:nrow(data2)){
        data2$TIME_analyser2[dr] = sub(as.character(as.Date(data2$TIME_analyser[dr])),"",print.POSIXct2(data2$TIME_analyser[dr]))
        data2$TIME_analyser2[dr] = trimws(sub("UTC","",data2$TIME_analyser2[dr]))
      }
      data2$TIME_analyser = data2$TIME_analyser2
      data2 = data2[,-which(colnames(data2) == "TIME_analyser2")]
    }
    rm(dataf,data_appended)

    ### reformat position ###
    # positions
    if(info$LAT_analyser == info$LON_analyser & !is.empty(info$LAT_analyser)){
      lat_lon = strsplit(data2$LAT_analyser, split = " ")
      data2$LAT_analyser = lat_lon[[1]][1]
      data2$LON_analyser = lat_lon[[1]][2]
    }
    # if not in %deg
    # %deg %min %sec %pos format or %deg %min %pos format
    if(grepl("%pos",info$POSITION_format) & grepl("%min",info$POSITION_format)){
      # grab the delimiter
      delim = substr(sub("%pos","",sub("%deg","",sub("%min","",sub("%sec","",info$POSITION_format)))),1,1)

      ### start pos
      # split based on delimiter
      lat = gsub("[[:alpha:]]","",data2$LAT_analyser)
      lat = strsplit(lat,split = delim)
      lon = gsub("[[:alpha:]]","",data2$LON_analyser)
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
      data2$LAT_analyser = as.numeric(lat)
      data2$LON_analyser = as.numeric(lon)
      if(gsub("[^[:alpha:]]+","", data2$LAT_analyser) == "S"){data2$LAT_analyser = -as.numeric(lat)}
      if(gsub("[^[:alpha:]]+","", data2$LON_analyser) == "W"){data2$LON_analyser = -as.numeric(lon)}
    }

    if(info$POSITION_format != "%deg" & grepl("%deg",info$POSITION_format) & length(strsplit(info$POSITION_format,split = "%")[[1]]) == 2 ){
      for(dr in 1:nrow(data2)){
        data2$LAT_analyser[dr] = sub(paste("\\",sub("%deg","",info$POSITION_format),sep = ""),"" ,data2$LAT_analyser[dr])
        data2$LON_analyser[dr] = sub(paste("\\",sub("%deg","",info$POSITION_format),sep = ""),"" ,data2$LON_analyser[dr])}
    }
    data2$LAT_analyser = as.numeric(data2$LAT_analyser)
    data2$LON_analyser = as.numeric(data2$LON_analyser)
    # transfer to analyser info and wipe info to be replaced by matching CTD data
    data2$STNNBR_analyser = data2$STNNBR
    data2$CASTNO_analyser = data2$CASTNO
    data2$STNNBR[data2$STNNBR != "U"] = NA
    data2$CASTNO[data2$CASTNO != "U"] = NA

    #### Write the new file
    # new file name and path
    if(duplicate){new_file = paste(paste(ex,"POC",info$source,info$Method,dup_id,sep = "_"),".csv",sep = "")}else{
      new_file = paste(paste(ex,"POC",info$source,info$Method,sep = "_"),".csv",sep = "")}

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
    subsource$citations = trimws(subsource$citations)
    submethod = method_info %>% filter(Method == info$Method, analysis_type == info$analysis_type)
    submethod$citation = trimws(submethod$citation)
    cite_tags = c(unlist(strsplit(info$citation,";" )), unlist(strsplit(subsource$citations,";" )))
    cite_tags = trimws(cite_tags[is.character(cite_tags)])

    # reformat and write new file
    fd <- file(new_file_path, open = "wt")
    writeLines(paste(info$Data_type,gsub("-","",Sys.Date()),userID, sep = ""), fd)
    writeLines(paste("#version=",version), fd)
    writeLines("#semi-manual exchange", fd)
    writeLines(paste("#ORIGINAL_CHLFILE/S:", paste(files, collapse="; ")), fd)
    writeLines(paste("#CHLFILE_MOD_DATE:",Sys.time(),"AEST"), fd)
    writeLines(paste("#SOURCED_FROM: ", info$source, "(",subsource$url,")", sep = ""), fd)
    writeLines(paste("#ANALYSIS_METHOD:", paste(submethod$long_AT, "according to the", submethod$Method, "method. See BIOMATE supplementary information and method citations for more information.")), fd)

    if(is.empty(info$contact)){writeLines(paste("#DATASET_CONTACT:", info$PI),fd)}else{
      writeLines(paste("#DATASET_CONTACT: ", info$PI,"(",info$contact,")", sep = ""),fd)}

    writeLines(paste("#DOI/s:", paste(bib[trimws(unlist(strsplit(info$citation,";" )))]$doi, collapse = ",")), fd)
    writeLines(paste("#URL:", paste(bib[trimws(unlist(strsplit(info$citation,";" )))]$url, collapse = ",")), fd)
    writeLines(paste("#BIOMATE_CITE_TAGS:", paste(cite_tags, collapse = ",")), fd)

    dcite = format(bib[cite_tags], style = "text", .bibstyle = "BIOMATE")
    dcite = gsub(pattern = "\n", replacement = " ", dcite)
    writeLines(paste("#DATA_CITATION/S:", paste(dcite , collapse = "\n# and ")), fd)
    if(!is.empty(submethod$citation)){
      mcite = format(bib[submethod$citation], style = "text", .bibstyle = "BIOMATE")
      mcite = gsub(pattern = "\n", replacement = " ", mcite)
      writeLines(paste("#METHOD_CITATION/S:", paste(mcite, collapse = "\n# and ")), fd)}

    if(!is.empty(info$Notes)){writeLines(paste("#NOTE:", info$Notes), fd)}


    #### CTD_data ###
    # get matches and append data
    data2$CTD_IDs = NA
    # list available CTD_IDs
    if(any(is.na(data2$STNNBR))){
      ctd_files = list.files(ctd_path, pattern = ".csv", full.names = F)

      data2$CTD_IDs[data2$STNNBR_analyser == "U"] = "U"
      # ctd_ids from CTD files
      CTD_IDs = unlist(strsplit(ctd_files,"_ctd1.csv"))
      # create a data frame with split info
      ctd_split = strsplit(CTD_IDs, split = "_")
      CTD_info = data.frame("CTD_ID" = CTD_IDs, "EXPOCODE" = sapply(ctd_split, "[[", 1), "STNNBR" = sapply(ctd_split, "[[", 3), "CASTNO"= sapply(ctd_split, "[[", 4), stringsAsFactors = F)

      ### match based on STNNBR and CASTNO ###
      # CTD stnnbr can be longer but pigment cant be (this accounts for when a cast number is specified in the station number - easy fix)
      # need to remove exact matches so they dont pop up in inexact matches
      CTD_info = CTD_info %>% mutate(STNCAST = paste(STNNBR, CASTNO))
      data2_stncast = paste(data2$STNNBR_analyser,data2$CASTNO_analyser)
      CTD_info = CTD_info %>% filter(EXPOCODE == ex)
      # exact matches
      CTD_info_exact = CTD_info %>% filter(EXPOCODE == ex, STNCAST %in% data2_stncast)
      sdx = data2_stncast %in% CTD_info_exact$STNCAST


      ### check that the CTD matches match in date and position - if they don't different station numbers have been used for POC and PROF
      if(nrow(CTD_info_exact) > 0){

        for(idx in 1:nrow(CTD_info_exact)){
          ctd_file = file.path(ctd_path,paste(CTD_info_exact$CTD_ID[idx],"_ctd1.csv",sep = ""))
          sub_data2 = data2 %>% filter(STNNBR_analyser == CTD_info_exact$STNNBR[idx],CASTNO_analyser == CTD_info_exact$CASTNO[idx], !is.na(LAT_analyser))
          if(nrow(sub_data2) == 0){next}

          # open file and read relevent lines
          f <- file( ctd_file, open = "r" )
          time_s = NA
          time_b = NA
          time_e = NA
          while( TRUE ){
            line <- readLines( f, 1L ,skipNul = T)
            if( grepl( "DATE =", line ) ){
              CTD_info_exact$DATE[idx] <- trimws(sub("DATE =", "", line ))
              date = CTD_info_exact$DATE[idx]
            }
            if( grepl( "CTD_START_TIME =", line ) ){
              t_line = line
              time_s<- trimws(sub("UTC","",sub("CTD_START_TIME =", "", line )))

            }
            if( grepl( "CTD_BOTTOM_TIME =", line ) ){
              time_b <- trimws(sub("UTC","",sub( "CTD_BOTTOM_TIME =", "", line )))

            }
            if( grepl( "CTD_END_TIME =", line ) ){
              time_e <- trimws(sub("UTC","",sub( "CTD_END_TIME =", "", line )))

            }
            if( grepl( "CTD_START_LATITUDE =", line ) ){
              CTD_info_exact$LATITUDE_s[idx] <- as.numeric(sub("CTD_START_LATITUDE =", "", line ))
            }
            if( grepl( "CTD_START_LONGITUDE =", line ) ){
              CTD_info_exact$LONGITUDE_s[idx] <- as.numeric(sub("CTD_START_LONGITUDE =", "", line ))
            }
            if( grepl( "CTD_BOTTOM_LATITUDE =", line ) ){
              CTD_info_exact$LATITUDE_b[idx] <- as.numeric(sub("CTD_BOTTOM_LATITUDE =", "", line ))
            }
            if( grepl( "CTD_BOTTOM_LONGITUDE =", line ) ){
              CTD_info_exact$LONGITUDE_b[idx] <- as.numeric(sub("CTD_BOTTOM_LONGITUDE =", "", line ))
            }
            if( grepl( "CTD_END_LATITUDE =", line ) ){
              CTD_info_exact$LATITUDE_e[idx] <- as.numeric(sub("CTD_END_LATITUDE =", "", line ))
            }
            if( grepl( "CTD_END_LONGITUDE =", line ) ){
              CTD_info_exact$LONGITUDE_e[idx] <- as.numeric(sub("CTD_END_LONGITUDE =", "", line ))
            }
            if(grepl("CTDPRS", line)){break}
          }
          close(f)

          # get closest time difference
          times = as.POSIXct(c(paste(date,time_s), paste(date,time_b), paste(date,time_e)), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
          # check forward in time, if not add 24 hrs as likely crossed days.
          if(!is.na(times[2]) & !is.na(times[1])){
            t_diff_check = difftime(times[2], times[1], units = "secs")
            if(t_diff_check < 0){times[2] = times[2] + (3600*24)}
          }
          if(!is.na(times[3]) & !is.na(times[1])){
            t_diff_check = difftime(times[3],  times[1], units = "secs")
            if(t_diff_check < 0){times[3] = times[3] + (3600*24)}
          }

          CTD_info_exact$TIME_s[idx] = times[1]
          CTD_info_exact$TIME_b[idx] = times[2]
          CTD_info_exact$TIME_e[idx] = times[3]
          times = times[!is.na(times)]


          if(length(times) > 0 & !is.na(sub_data2$DATE_analyser[1])){
            t_diffs = abs(difftime(times, as.POSIXct(paste(sub_data2$DATE_analyser[1],sub_data2$TIME_analyser[1]),tz = "UTC"), units = "secs"))
            CTD_info_exact$t_diff[idx] = t_diffs[which.min(t_diffs)]}else{
              CTD_info_exact$t_diff[idx] = NA}
          # if the same date is recorded come up as match.
          if(length(times) == 0 & CTD_info_exact$DATE[idx] == sub_data2$DATE_analyser[1]){CTD_info_exact$t_diff[idx] = 1*60*60*24}

          # get position difference
          pos = matrix(c(CTD_info_exact$LONGITUDE_s[idx], CTD_info_exact$LATITUDE_s[idx],CTD_info_exact$LONGITUDE_b[idx], CTD_info_exact$LATITUDE_b[idx],CTD_info_exact$LONGITUDE_e[idx], CTD_info_exact$LATITUDE_e[idx]) , ncol = 2, byrow =T)
          pos_diff = distGeo(pos, c(sub_data2$LON_analyser[1], sub_data2$LAT_analyser[1]))
          CTD_info_exact$p_diff[idx] = pos_diff[which.max(pos_diff)]

          # identify if missmatch
          if(!is.empty(CTD_info_exact$p_diff[idx]) & !is.na(CTD_info_exact$t_diff[idx])){
            #note here t_diff is in days
            # use Johnson 2017 suggestion 8km and 1 day
            if(CTD_info_exact$p_diff[idx] > 8000 | CTD_info_exact$t_diff[idx] > 1*60*60*24){CTD_info_exact$nomatch[idx] = T}else{CTD_info_exact$nomatch[idx] = F}
          }else{CTD_info_exact$nomatch[idx] = T}

        }

        if(any(CTD_info_exact$nomatch)){print("PROF and POC station mismatch detected. Moving to closest PROF detection method.")}else{

          # asign CTD_IDs
          data2$CTD_IDs[sdx] = paste(ex,"CTD",data2$STNNBR_analyser[sdx],data2$CASTNO_analyser[sdx],sep = "_")
        }
        rm(sub_data2, times)
      }

      #was trying to find similar matches.... waste of time
      # CTD_info_notexact = CTD_info %>% filter(EXPOCODE == ex, !(STNCAST %in% data2_stncast), grepl(paste(unique(data2$STNNBR), collapse="|"),STNNBR), )
      # CTD_info_exact$STNCAST_2 = CTD_info_exact$STNCAST
      # if(nrow(CTD_info_notexact) > 0){
      # CTD_info_notexact$STNCAST_2 = NA
      # remaining_stnnbrs = unique(data2$STNNBR)[-grep(paste(CTD_info_exact$STNNBR, collapse="|"), unique(data2$STNNBR))]
      # CTD_info_notexact$STNCAST_2 = unlist(lapply(CTD_info_notexact$STNNBR,FUN = function(x){a =paste(remaining_stnnbrs[stri_detect_fixed(x,remaining_stnnbrs)])
      #                                                                                         if(length(a)== 0){return(NA)}else{return(a)}
      #                                                                                       }))
      # CTD_info_f = rbind(CTD_info_exact, CTD_info_notexact)}
      # else{CTD_info_f = CTD_info_exact}


      ### unmatched CTD casts ###
      un = is.na(data2$CTD_IDs)
      if(any(un)){
        # unmatched data
        unmatched_df = data.frame("STNCAST" = data2_stncast[un],"Time" = as.character(data2$TIME_analyser[un]), "Date" = as.character(data2$DATE_analyser[un]), "LAT" = data2$LAT_analyser[un],  "LON" = data2$LON_analyser[un], stringsAsFactors = F)
        unmatched_df = unmatched_df[!duplicated(unmatched_df$STNCAST),]

        # can only match when time AND position info available
        unmatched_df_t = unmatched_df %>% filter(!is.na(Time), !is.empty(Time), !is.empty(Date)) %>%
          mutate(time_m = as.POSIXct(paste(Date,Time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))
        notime = ifelse(nrow(unmatched_df_t) == 0,T,F)
        unmatched_df_p = unmatched_df %>% filter(!is.na(LAT))
        nopos = ifelse(nrow(unmatched_df_p) == 0,T,F)

        if(!nopos){

          # get ctd basic information
          CTD_df = data.frame("CTD_ID" = character(),"DATE" = character(), "TIME_s" = as.POSIXct(character()), "TIME_b"= as.POSIXct(character()), "TIME_e"= as.POSIXct(character()), "LAT_s"=numeric(), "LON_s" = numeric(), "LAT_b"=numeric(), "LON_b" = numeric(), "LAT_e"=numeric(), "LON_e" = numeric())
          if(nrow(unmatched_df) > 0){

            if(nrow(CTD_info) == 0){print(paste("There is no CTD data for", ex))}else{
              for(ctd in CTD_info$CTD_ID){
                ctd_file = file.path(ctd_path,paste(ctd,"_ctd1.csv",sep = ""))

                time_s = NA
                time_b = NA
                time_e = NA
                # open file and read relevent lines
                f <- file( ctd_file, open = "r" )

                while( TRUE ){
                  line <- readLines( f, 1L ,skipNul = T)
                  if( grepl( "DATE =", line ) ){
                    date <- trimws(sub("DATE =", "", line ))
                  }
                  if( grepl( "CTD_START_TIME =", line ) ){
                    t_line = line
                    time_s<- trimws(sub("UTC","",sub("CTD_START_TIME =", "", line )))

                  }
                  if( grepl( "CTD_BOTTOM_TIME =", line ) ){
                    time_b <- trimws(sub("UTC","",sub( "CTD_BOTTOM_TIME =", "", line )))

                  }
                  if( grepl( "CTD_END_TIME =", line ) ){
                    time_e <- trimws(sub("UTC","",sub( "CTD_END_TIME =", "", line )))

                  }
                  if( grepl( "CTD_START_LATITUDE =", line ) ){
                    lat_s <- as.numeric(sub("CTD_START_LATITUDE =", "", line ))
                  }
                  if( grepl( "CTD_START_LONGITUDE =", line ) ){
                    lon_s <- as.numeric(sub("CTD_START_LONGITUDE =", "", line ))
                  }
                  if( grepl( "CTD_BOTTOM_LATITUDE =", line ) ){
                    lat_b <- as.numeric(sub("CTD_BOTTOM_LATITUDE =", "", line ))
                  }
                  if( grepl( "CTD_BOTTOM_LONGITUDE =", line ) ){
                    lon_b <- as.numeric(sub("CTD_BOTTOM_LONGITUDE =", "", line ))
                  }
                  if( grepl( "CTD_END_LATITUDE =", line ) ){
                    lat_e <- as.numeric(sub("CTD_END_LATITUDE =", "", line ))
                  }
                  if( grepl( "CTD_END_LONGITUDE =", line ) ){
                    lon_e <- as.numeric(sub("CTD_END_LONGITUDE =", "", line ))
                  }
                  if(grepl("CTDPRS", line)){break}
                }
                close(f)

                times = as.POSIXct(c(paste(date,time_s), paste(date,time_b), paste(date,time_e)), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
                # check forward in time, if not add 24 hrs as likely crossed days.
                if(!is.na(times[2]) & !is.na(times[1])){
                  t_diff_check = difftime(times[2], times[1], units = "secs")
                  if(t_diff_check < 0){times[2] = times[2] + (3600*24)}
                }
                if(!is.na(times[3]) & !is.na(times[1])){
                  t_diff_check = difftime(times[3], times[1], units = "secs")
                  if(t_diff_check < 0){times[3] = times[3] + (3600*24)}
                }

                # gather ctd information
                CTD_df = CTD_df %>% add_row(CTD_ID = ctd, DATE = date, TIME_s = times[1], TIME_b= times[2], TIME_e= times[3], LAT_s = lat_s, LON_s = lon_s, LAT_b = lat_b, LON_b = lon_b, LAT_e = lat_e, LON_e = lon_e)


                #

                # intvls1 = c(270,270,270)*60
                # intvls2 = c(270,270,270)*60
                # tdx = which(is.finite(times))
                # if(all(is.na(times))){print(paste("There is no time information in ctd file",ctd))}
                # if(any(unmatched_df$Date == date)){
                #   unmatched_df_sub = unmatched_df %>% filter(Date == date)
                #   unmatched_df_sub_t = unmatched_df_sub %>% mutate(time_m = as.POSIXct(paste(Date,Time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))
                #   unmatched_df_sub_p
                #
                #   # match closest times
                #     max_t = max(times[tdx] + intvls1[tdx])
                #     min_t = min(times[tdx] - intvls2[tdx])
                #     res = unmatched_df_sub_t$time_m > min_t & unmatched_df_sub_t$time_m < max_t
                #     if(length(which(res)) > 1){res = which.min(abs(unmatched_df_sub_t$time_m - min_t))} # get minimum time difference
                #     if(any(res)){
                #       data2$CTD_IDs[data2_stncast == unmatched_df_sub_t$STNCAST[res]] = ctd
                #       data2$STNNBR[data2_stncast == unmatched_df_sub_t$STNCAST[res]] = stn
                #       data2$CASTNO[data2_stncast == unmatched_df_sub_t$STNCAST[res]] = cast
                #     }
                #     }
              }

              ### match closest times
              col_t = c(3:5)[!colSums(is.na(CTD_df[3:5])) == nrow(CTD_df)]
              if(!notime){
                #calculate time differences
                unmatched_df_t$closest_t = unlist(lapply(1:nrow(unmatched_df_t), calc.time.min))
                unmatched_df_t$t_diff = unlist(lapply(1:nrow(unmatched_df_t), calc.time.diff))

                # old code
                # if(length(col_t) > 1){
                #
                # unmatched_df_t = unmatched_df_t %>% mutate(closest_t = which.min(abs(unlist(lapply(col_t, FUN = function(x){CTD_df[,x] - time_m}))))[1]) %>%
                #                                   mutate(t_diff = abs(gather(CTD_df[,col_t])$value[closest_t] - time_m))}else{
                # unmatched_df_t = unmatched_df_t %>% mutate(closest_t = which.min(abs(unlist(lapply(col_t, FUN = function(x){CTD_df[,x] - time_m}))))[1]) %>%
                #                                       mutate(t_diff = abs(CTD_df[,col_t][closest_t] - time_m))
                #                                   }

                # check it is within 3 hours and assign CTD_ID
                unmatched_df_t = unmatched_df_t %>% filter(!is.na(closest_t), t_diff < t_thresh*60*60)

                # Any matches?
                if(nrow(unmatched_df_t) > 0){
                  # get CTD_IDs
                  unmatched_df_t$CTD_ID = unlist(lapply(1:nrow(unmatched_df_t),function(x){r = CTD_df$CTD_ID[unmatched_df_t$closest_t[x]%%nrow(CTD_df)]
                  if(is.empty(r)){r = CTD_df$CTD_ID[nrow(CTD_df)]}
                  return(r)}))
                  # join with data 2 IDs
                  joined_t = left_join(data.frame("STNCAST" = data2_stncast[un]),unmatched_df_t , by = "STNCAST")
                }}


              ### match closest positions

              # match with closest lat/lon on same date
              d_calc = lapply(1:nrow(unmatched_df_p), calc.dist.min)
              unmatched_df_p$dist_min  = unlist(lapply(d_calc,"[","min"))
              unmatched_df_p$d_diff = unlist(lapply(d_calc,"[","diff"))
              # check within distance threshold
              unmatched_df_p  = unmatched_df_p %>% filter(d_diff < d_thresh)


              # Any matches?

              if(nrow(unmatched_df_p) > 0 & notime){
                # if no time matches consider date of position matches.
                unmatched_df_p = unmatched_df_p %>% mutate(Date_match = ifelse(Date == CTD_df$DATE[dist_min],T,F))
                unmatched_df_p = unmatched_df_p %>% filter(Date_match)
              }
              if(nrow(unmatched_df_p) > 0){
                # get CTD_ID
                unmatched_df_p$CTD_ID = unlist(lapply(1:nrow(unmatched_df_p),function(x){CTD_df$CTD_ID[unmatched_df_p$dist_min[x]]}))
                # join with data 2 IDs
                joined_p = left_join(data.frame("STNCAST" = data2_stncast[un]),unmatched_df_p , by = "STNCAST")
                if(exists("joined_t")){
                  # join with data 2 IDs
                  joined_t_p = left_join(joined_t,unmatched_df_p[,c("STNCAST","dist_min","d_diff")] , by = "STNCAST")
                  rm(joined_t)
                }else{joined_t_p = joined_p
                rm(joined_p)}
                data2[un,c("CTD_IDs")] = joined_t_p$CTD_ID
                rm(joined_t_p)}
            }

            rm(unmatched_df, CTD_df, unmatched_df_p, unmatched_df_t)
          }
        }

      }
      rm(un)

      # append CTD_ID, STNNBR, CASTNO, LAT, LON, DATE, TIME
      # Now get the data from files that have identified matches
      CTD_list = unique(data2$CTD_IDs)
      CTD_list = CTD_list[!is.na(CTD_list)]
      CTD_list = CTD_list[CTD_list != "U"]
      # set up CTD data frame
      if(length(CTD_list) > 0){
        for(ctd in CTD_list){
          ctd_file = file.path(ctd_path,paste(ctd,"_ctd1.csv",sep = ""))

          idx = which(data2$CTD_IDs == ctd)

          # open file and read relevent lines
          f <- file( ctd_file, open = "r" )

          while( TRUE ){
            line <- readLines( f, 1L ,skipNul = T)
            if( grepl( "DATE =", line ) ){
              data2$DATE[idx] <- trimws(sub("DATE =", "", line ))
            }
            if( grepl( "CTD_START_TIME =", line ) ){
              t_line = line
              data2$TIME_s[idx] <- trimws(sub("UTC","",sub("CTD_START_TIME =", "", line )))
              if(is.empty(data2$TIME_s[idx])){data2$TIME_s[idx] = NA}
            }
            if( grepl( "CTD_BOTTOM_TIME =", line ) ){
              data2$TIME_b[idx] <- trimws(sub("UTC","",sub( "CTD_BOTTOM_TIME =", "", line )))
              if(is.empty(data2$TIME_b[idx])){data2$TIME_b[idx] = NA}
            }
            if( grepl( "CTD_END_TIME =", line ) ){
              data2$TIME_e[idx] <- trimws(sub("UTC","",sub( "CTD_END_TIME =", "", line )))
              if(is.empty(data2$TIME_e[idx])){data2$TIME_e[idx] = NA}
            }
            if( grepl( "CTD_START_LATITUDE =", line ) ){
              data2$LATITUDE[idx] <- sub("CTD_START_LATITUDE =", "", line )
            }
            if( grepl( "CTD_START_LONGITUDE =", line ) ){
              data2$LONGITUDE[idx] <- sub("CTD_START_LONGITUDE =", "", line )
            }
            if( grepl( "CTD_BOTTOM_LATITUDE =", line ) ){
              if(is.na(data2$LATITUDE[idx[1]])){data2$LATITUDE[idx] <- sub("CTD_BOTTOM_LATITUDE =", "", line )}
            }
            if( grepl( "CTD_BOTTOM_LONGITUDE =", line ) ){
              if(is.na(data2$LONGITUDE[idx[1]])){data2$LONGITUDE[idx] <- sub("CTD_BOTTOM_LONGITUDE =", "", line )}
            }
            if( grepl( "CTD_END_LATITUDE =", line ) ){
              if(is.na(data2$LATITUDE[idx[1]])){data2$LATITUDE[idx] <- sub("CTD_END_LATITUDE =", "", line )}
            }
            if( grepl( "CTD_END_LONGITUDE =", line ) ){
              if(is.na(data2$LONGITUDE[idx[1]])){data2$LONGITUDE[idx] <- sub("CTD_END_LONGITUDE =", "", line )}
            }
            if( grepl( "STNNBR or EVENTNBR =", line ) ){
              data2$STNNBR[idx] = trimws(sub("STNNBR or EVENTNBR =", "", line ))
            }
            if( grepl( "CASTNO =", line ) ){
              data2$CASTNO[idx] = trimws(sub("CASTNO =", "", line ))
            }
            # if( grepl( "STNNBR or EVENTNBR =", line ) ){
            #   if(as.character(data2$STNNBR[idx[1]]) != trimws(sub("STNNBR or EVENTNBR =", "", line ))){stop("station numbers dont match")}
            # }
            # if( grepl( "CASTNO =", line ) ){
            #   if(as.character(data2$CASTNO[idx[1]]) != trimws(sub("CASTNO =", "", line ))){stop("station numbers dont match")}
            # }
            if(grepl("CTDPRS", line)){break}
          }
          close( f )
          data2$CTD_match[idx] = 1

        }}
      rm(CTD_info_exact,CTD_info, CTD_list)
    }
    #data2 = data2[,which(colnames(data2) != "CTD_IDs")]
    ## something here for pulling from underway data
    final_data = data2
    ordered_names = c("CTD_IDs","DATE","TIME_s","TIME_b","TIME_e","LATITUDE","LONGITUDE","STNNBR",	"CASTNO",	"DATE_analyser",	"TIME_analyser",	"LAT_analyser",	"LON_analyser","STNNBR_analyser","CASTNO_analyser",	"Sample_ID",	"BOTTLE",	"DEPTH",pig_names)
    ordered_units = c("NA","YYYY-mm-dd","HH:MM:SS","HH:MM:SS","HH:MM:SS","DEGREES NORTH","DEGREES EAST","NA",	"NA",	"YYYY-mm-dd","HH:MM:SS","DEGREES NORTH","DEGREES EAST","NA","NA",	"NA",	"NA",	"METERS",rep("MG/M3",length(pig_names)))

    final_data = final_data[,ordered_names]
    if(any(is.na(final_data$CTD_match))){
      final_data$CTD_match[is.na(final_data$CTD_match)] = 0
    }
    rm(data2)

    writeLines(paste("NUMBER_HEADERS = 6"), fd)
    writeLines(paste("EXPOCODE =", ex), fd)
    writeLines(paste("SHIP =",platforms$`Platform Name`[match(substr(ex,1,4), platforms$`NODC code`)]),fd)
    writeLines(paste("TIMEZONE = UTC"), fd)
    writeLines("missing_value = -999", fd)
    writeLines("not_detected = -888", fd)


    # assign missing value number
    final_data$DATE = as.character(final_data$DATE)
    final_data$DATE_analyser = as.character(final_data$DATE_analyser)

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
    # units = lapply(colnames(final_data),function(x){if(paste(x,"_u",sep = "") %in% colnames(info)){info[,paste(x,"_u",sep = "")]}else{NA}})
    writeLines(toString(ordered_units), fd)

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

    print(paste("row",rw,"EXPOCODE",ex,"SOURCE",info$source,"METHOD",info$Method,"successfully run"))
    ### checks
    # final_data$DEPTH = as.numeric(final_data$DEPTH)
    # final_data$DEPTH[final_data$CTD_IDs == "U"] = "U"
    if(any(final_data$STNNBR == -999)){print(paste("Warning:There is no station number assigned for a data point in",rw,ex," Is this point underway?"))}
    #if(any(final_data$LATITUDE == -999 & final_data$STNNBR != "U" & final_data$DEPTH > 10 & !is.empty(final_data$DEPTH))){print(paste("Warning:There is some ctd data not matched, or missing",rw,ex,"for this many measurements", length(which(final_data$LATITUDE == -999 & final_data$STNNBR != "U" & final_data$DEPTH > 10 & !is.empty(final_data$DEPTH)))))}
    rm(final_data)
    if(exists("header_vars")){rm(header_vars)}
    if(exists("CTD_IDs")){rm(CTD_IDs)}
    if(exists("sub_meta")){rm(sub_meta)}
    rm(data_vars,header_in_line,header_line,headers,line)
  }

}
