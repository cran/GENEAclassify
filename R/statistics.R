
#' Calculates the ratio of the principle step frequencies.
#' 
#' Note this function is intended to be called from within 
#' an apply call in \code{\link{segmentation}}, hence the 
#' unusual syntax of defining the sampling frequency, freq
#' as an object in the parent frame
#'
#' @title Principal Frequency Ratio
#' @param principals principal frequency
#' @param freq single character naming single numeric variable 
#'     in parent environment stating sampling frequency
#' @param nfr single integer number of frames to search
#' @param \dots other arguments to be swallowed
#' @return single numeric
#' @export
#' @keywords internal
#' @examples
#'     Freq <- 100
#'     GENEAratio(principals = c(0.1, 3.2, 0.3, 0.0, 0.2))

GENEAratio <- function(principals, freq = "Freq", nfr = 3L, ...) { 
    
    freq <- get(freq, envir = parent.frame(n = nfr))
    
    return(sum(principals < (freq / 4), na.rm = TRUE) / 
        sum(principals > (freq / 4), na.rm = TRUE)) 
}


#' Step count is the number of step lengths divided by 2.
#' 
#' @title Count Steps
#' @param x vector
#' @param \dots other arguments to be swallowed
#' @return single numeric
#' @export
#' @keywords internal
#' @examples
#'    x1 <- c(20, 15, 10)
#'    GENEAcount(x1)
#'    x2 <- c(300, 255, 111)
#'    GENEAcount(x2)

GENEAcount <- function(x, ...) { return(round(length(x) / 2)) }

#' @title Skewness, a measure of centredness
#' 
#' @description Many datasets do not have the highest density 
#' of measurements in the middle of the distribution of the data.
#' Skewness is a measure of the asymmetry of a probability distribution.
#' 
#' @param x numeric vector
#' @param na.rm single logical should missing values be removed
#' @param \dots other arguments to be swallowed
#' @return single numeric
#' @export
#' @keywords internal
#' @examples
#' GENEAskew(1:10)
#' GENEAskew((1:10)^2)
#' GENEAskew((1:10)^0.5)

GENEAskew <- function(x, na.rm = TRUE, ...) {
    
    if (na.rm) { x <- x[!is.na(x)] }
    
    nn <- length(x)
    
    xbar <- mean(x)
    
    top <- sum((x - xbar)^3) / nn
    
    bottom <- sum((x - xbar)^2) / nn
    
    sk <- top / bottom^(3/2)
    
    return(sk)
}


#' Adding in a feedback call sumdiff that finds the difference of a function and then sums this difference. 
#' Called by \code{segmentation}.
#' @title Finding the sum of the differences
#' @param x vector of numeric values
#' @param na.rm single logical should missing values be removed
#' @return A single value data.
#' @export
#' @keywords internal
#' @examples
#'    tmp1 <- c(1,3,2,6,4,5,3,9,10)
#'    sumdiff(x = tmp1)

sumdiff <- function(x, na.rm=TRUE){
  tmp1 <- diff(na.omit(x))
  return(sum(tmp1))
}

#' Finding the mean of the differences
#' Called by \code{segmentation}.
#' @title Finding the mean of the differences
#' @param x vector of numeric values
#' @param na.rm Remove the na values in the vector
#' @return A single value data.
#' @export
#' @keywords internal
#' @examples
#'    tmp1 <- c(1,3,2,6,4,5,3,9,10)
#'    meandiff(x = tmp1)

meandiff <- function(x, na.rm=TRUE){
  tmp1 <- diff(na.omit(x))
  return(mean(tmp1))
}

#' Finding the absolute sum of the differences
#' called by \code{segmentation}.
#' @title Finding the sum of the differences
#' @param x vector of numeric values
#' @param na.rm single logical should missing values be removed
#' @return A single value data.
#' @export
#' @keywords internal
#' @examples
#'    tmp1 <- c(1,3,2,6,4,5,3,9,10)
#'    abssumdiff(x = tmp1)

abssumdiff <- function(x, na.rm=TRUE){
  tmp1 <- diff(na.omit(x))
  return(sum(abs(tmp1)))
}

#' Called by \code{segmentation}.
#' @title Finding the standard deviation of the differences
#' @param x vector of numeric values
#' @param na.rm single logical should missing values be removed
#' @return A single value data.
#' @export
#' @keywords internal
#' @examples
#'    tmp1 <- c(1,3,2,6,4,5,3,9,10)
#'    sddiff(x = tmp1)

sddiff <- function(x, na.rm=TRUE){
  tmp1 <- diff(na.omit(x))
  return(sd(tmp1))
}



