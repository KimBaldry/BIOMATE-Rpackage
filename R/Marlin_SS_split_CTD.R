#' @title Split old Southern Surveyor files
#'
#' @author Kimberlee Baldry
#' @description This script splits old southern surveyor text files from ctds, and stors tab separated files (.tsv) in the same directory as the original file. It references a "Marlin_SS_split_CTD_meta.csv" file. This file must be stored in the path directory. It will output a data amendment file for records.
#' @note This function is only useful to reformat very old Southern Surveyor files.
#'
#' @return returns a "successful run" message and saves split files in the same directory as the original data file.
#'
#' @param path The path where the Marlin data is stored. This path must include sub-folders for expocodes
#'
#' @export

SS_CTDsplit = function(path){
  # read in the metadata file
  meta = read.csv(file = file.path(path, "Marlin_SS_split_CTD_meta.csv"), header = T, strip.white = T,stringsAsFactors = F)

  ### loop through and create the files for each EXPOCODE
  for(ex in meta$EXPOCODE)
  {
    fname = file.path(path,ex,meta$fname[which(meta$EXPOCODE == ex)])

    # open file and read relevent lines
    f <- file(fname, open = "r" )
    # assign a file ID for each file. Note that this isnt necesisarily the CTD station ID, just a file ID
    fileID = 1
    while( TRUE ){
      line <- readLines( f, 1L )
      # "EEEE" signals the end of teh compiled file. "SSSS" indicates the begining of a new CTD profile
      if(grepl( "EEEEEEE", line )){break}
      if( grepl( "SSSSS", line ) )
      {
        # close an open connection - if there is one
        if(exists("fd")){close(fd)}
        # delete the file if it already exists
        if(file.exists(file.path(path,ex,paste(ex,fileID,"ctd1.tsv",sep = "_")))){unlink(file.path(path,ex,paste(ex,fileID,"ctd1.tsv",sep = "_")))}
        # open a new connection for the new profile data to be written to
        fd <- file(file.path(path,ex,paste(ex,fileID,"ctd1.tsv",sep = "_")), open = "wt")
        fileID = fileID + 1
      }else{
        writeLines(line, fd)
      }
    }
    close(f)
    close(fd)
    rm(f, fileID, fd, line)
    created_files = list.files(file.path(path, ex),pattern = ".tsv")

    tf = file(file.path(path,ex,paste("Data_amendments",Sys.Date(),".txt", sep = "")),open = "wt")
    writeLines(paste("The following files were created from",fname,"using the function SS_CTDsplit() from Marlin_SS_split_CTD.R:"),tf)
    for(sp in created_files)
    {
      writeLines(paste("-",sp),tf)
    }
    close(tf)
  }

  return("successful run")
}
