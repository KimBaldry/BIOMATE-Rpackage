#'  @title Check if an object is empty
#'
#'  @author Kimberlee Baldry
#'  @description This is a function that checks if an object (charactor, vector, integer) has no information contained with it (including null, empty, length of 0).
#'
#'  @return Logical result or vector of logical
#'
#'  @param x an object
#'  @param first.only Are we only investigating the fist entry of a character vector, or all entries of a character vector (in this case a vector is returned)?
#'  @param all.na.empty is a vector with all NA values counted as empty?
#'
#'  @import purrr
#'
#'  @export



is.empty <- function(x, first.only = TRUE, all.na.empty = TRUE) {
  # do we have a valid vector?
  if (!is.null(x)) {
    # if it's a character, check if we have only one element in that vector
    if (is.character(x)) {
      # characters may also be of length 0
      x = trimws(x)
      if (length(x) == 0) return(TRUE)
      # else, check all elements of x
      zero_len <- nchar(x) == 0
      # return result for multiple elements of character vector
      if (first.only) {
        zero_len <- is.logical(zero_len[1]) && length(zero_len[1]) == 1L && !is.na(zero_len[1]) && zero_len[1]
        if (length(x) > 0) x <- x[1]
      } else {
        return(unname(zero_len))
      }
      # we have a non-character vector here. check for length
    } else if (is.list(x)) {
      x <- purrr::compact(x)
      zero_len <- length(x) == 0
    } else {
      zero_len <- length(x) == 0
    }
  }

  return(any(is.null(x) || zero_len || (all.na.empty && all(is.na(x)))))
}

