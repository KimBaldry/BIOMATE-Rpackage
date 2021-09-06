#' @title Create a new bibstyle for use in BIOMATE.
#'
#' @author Kimberlee Baldry
#' @description Executing this function will add a new "BIOMATE" style to a local library. This style is called within PIG_to_WHPE() and PROF_to_WHPE().
#'
#' @return returns message letting the user know the task was completed.
#'
#' @import tools
#'
#'
#' @export

add_bibstyle <- function(){
  
  BM = bibstyle(style = "text", .init = T)
  BM$formatMisc <- function (paper) 
  {
    result = paste(authorList(paper), " ", fmtYear(paper$year), " ", fmtTitle(paper$title), " ", paper$publisher, ". ", sep = "")
    if(!is.null(paper$organization)){result = paste(result,  paper$organization,". ", sep = "")}
    if(!is.null(paper$doi)){result = paste(result,  paper$doi,". ", sep = "")}
    if(!is.null(paper$url)){result = paste(result,  paper$url,". ", sep = "")}
    if(!is.null(paper$urldate)){result = paste(result, "Accessed on ",  paper$urldate,". ", sep = "")}
    result
  }
  environment(BM$formatMisc) = BM
  bibstyle("BIOMATE",BM)
  print("BIOMATE style has been added.")
}


