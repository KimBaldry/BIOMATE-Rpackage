#' @title Reference BIO-MATE right
#'
#' @author Kimberlee Baldry
#' @description Export the BIO-MATE citation, along with individual data citations in on BibTeX file to make referencing easy. Upload this file into your reference manager or text editing software. A handy reference table is also compiles. Remember all individual data sets that are downloaded with BIO-MATE need to be cited appropriately.
#'
#'  @return A BibTeX file, saved to the directed location, and a data.table containing reference codes for datasets.
#'
#'  @param EXPOCODE A list of the EXPOCODES for the data you have used. Deefault is "all".
#'  @param streams The data streams to be cited. Can be "all" or a combination of "PIG", "PROF","POC' and "UWY".
#'  @param reports Are the cruise report citations needed: Y/N.
#'  @param path_out The path where the BibTeX file will be saved.
#'
#' @export



export_ref <- function(EXPOCODES = "all", streams = "all",  reports = "Y", path_out = "./"){

  # direct to citation files
  meta_path = file.path(.libPaths()[1],"BIOMATE","citations")
  file_names = list.files(meta_path)
  # select files besed on user-defined EXPOCODES
  if(any(EXPOCODES != "all")){
    EX = unlist(strsplit(file_names,".bib"))
    check = lapply(EXPOCODES, function(x){match(x, EX)})
    missing = is.na(check)
    selected_files = file_names[unlist(check[!is.na(check)])]
  }else{selected_files = file_names}

  # compile citations into one BibTEX file for user
  outName = file.path(path_out,"BIO-MATE_references.bib")
  if(file.exists(outName)){file.remove(outName)}
  outFile <- file(outName, "w")
  for (i in selected_files){
    lines <- readLines(file.path(meta_path,i))
    nlines = length(lines)

    ### remove unwanted BIO-MATE data streams and reports
    # create a reference index for the start of each reference using the @ symbol
    ref_idx = grep("@", lines)
    # get index of reports citations if required
    if(reports == "Y"){
      report_idx = grep("report", lines[ref_idx])
      if(length(report_idx) == 0){rm(report_idx)}
    }
    # get index of perticular data stream citations if required
    if(any(streams != "all")){
      ## change prof to ctd (mismatch in references)
      #if(any(grepl("PROF", streams, ignore.case = T))){stream[grep("PROF", stream, ignore.case = T)] = "CTD"}
      stream_idx = unlist(lapply(streams, function(x){grep(x,lines[ref_idx], ignore.case = T)}))
      # for now we dont care if a citation is missing for a stream. When the metadata table is complete we can check data existance against citation existance
      # cobine report idx and select stream idx if needed
      if(exists("report_idx")){final_idx = c(report_idx, stream_idx)}else{final_idx = stream_idx}
      # if all streams are selected, but no report citations are wanted, the final index is the complement of the report index
    }else{
      if(reports == "N"){
        report_idx = grep("report", lines[ref_idx])
        final_idx = setdiff(1:length(ref_idx), report_idx)

        }
      }
    if(reports == "Y" & !exists("final_idx")){final_idx = report_idx}
    # if we have made a new index (ie reports != "Y" or streams != "all") then remove the unwanted lines
    if(exists("final_idx")){
      if(length(final_idx) > 0){
        final_idx = final_idx[order(final_idx, decreasing = T)]
        keep_lines = NULL
        for(r in final_idx){
          if(r == length(ref_idx)){
            keep_lines = c(lines[ref_idx[r]:nlines], keep_lines)
          }else{
            keep_lines = c(lines[ref_idx[r]:(ref_idx[r+1]-1)], keep_lines)

          }
        }
        lines = keep_lines
      }else{print(paste("No references selected for", unlist(strsplit(i,".bib"))))
        next}
      rm(keep_lines, final_idx)
      if(exists("report_idx")){rm(report_idx)}
      if(exists("stream_idx")){rm(stream_idx)}
    }
    ### write the file
    writeLines(lines, outFile)
    if(i != selected_files[length(selected_files)]){writeLines("", outFile)}
    rm(lines, nlines, ref_idx)
  }
  close(outFile)
  ### print a message if some EXPOCODES are not matched to files.
  if(any(missing)){print(paste("There are no BibTEX files in BIOMATE for the following voyages:", paste(EXPOCODES[missing], collapse = ", "), sep = " "))}

}

