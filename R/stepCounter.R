
#' Function to calculate the number and variance of the steps in the data.
#'
#' @title Step Counter
#' @param data The data to use for calculating the steps. This should be a matrix
#' with columns time, x, y, z.
#' @param samplefreq The sampling frequency of the data, in hertz,
#' when calculating the step number (default 100).
#' @param smlen single integer number of data points used as a window
#' when counting zero crossing events
#' @param filterorder single integer, order of the Butterworth bandpass filter,
#' passed to argument n of \code{\link[signal]{butter}}.
#' @param threshold A threshold that acts as a proxy for magnitude on the filtered xz signal for calculating steps. 
#' @param boundaries length 2 numeric vector specifying lower and upper bounds
#' of Butterworth or Chebychev filter (default \code{c(0.04, 1)}),
#' passed to argument W of \code{\link[signal]{butter}} or \code{\link[signal]{cheby1}}.
#' @param stepmethod select method to use: \enumerate{
#'     \item 'Chebyfilter' (default)
#'     \item 'Butterfilter'
#'     \item 'longrun'
#'     \item 'none'
#' }
#' @param Rp the decibel level that the cheby filter takes, see \code{\link[signal]{cheby1}}.
#' @param plot.it single logical create plot of data and zero crossing points (default \code{FALSE}).
#' @param Centre If Centre set to true (default) then the step counter zeros the xz series before filtering.
#' @param STFT If STFT is TRUE then the Step Counter uses the STFT function to find the length of the window for each segment.
#' @param fun character vector naming functions by which to summarize steps.
#' "count" is an internally implemented summarizing function that returns step count.
#' @return Returns a vector with length fun.
#' @export
#' @importFrom signal butter cheby1
#' @examples
#' d1 <- matrix(c(100, 101, -0.79, -0.86,
#'         -0.17, -0.14, 0.53, 0.46),
#'     nrow = 2, ncol = 4)
#' stepCounter(data = d1, stepmethod = "longrun")

stepCounter <- function(data,
                        samplefreq = 100,
                        smlen = 20L,
                        filterorder = 4L,
                        threshold = 0.01,
                        boundaries = c(0.15, 1.0),
                        stepmethod = c("Chebyfilter","Butterfilter","longrun","none"),
                        Rp = 0.5,
                        plot.it = FALSE,
                        Centre = TRUE,
                        STFT = FALSE,
                        fun = c("count","mean", "sd", "mad")) {
    
    if (missing(data)) {stop("data is missing") }
    if (missing(stepmethod)){stepmethod = "Chebyfilter"}
    if (!is.character(fun)) { stop("fun must be character vector of function names") }
    if (length(fun) < 1L) { stop("fun must name at least one function") }

    method <- match.arg(arg = stepmethod)
    #Create the data for calculating the steps by adding x and z columns
    xzSeries <- data[, 2] + data[, 4]
    
    # Centre the xz series about 0.
    if (Centre == TRUE){xzSeries = xzSeries - mean(xzSeries)}
    
    # Calculate the mean principal frequency and use this for the window length.
    if (STFT == TRUE){
      # Calculate the stft.
      ft <- (stft(as.matrix(data), quiet = TRUE,
                  reassign = TRUE, date.col = FALSE,
                  freq = samplefreq))
      
      samplefreq <- ft$sampling.frequency
      # If any nas created 
      if (is(ft, class2 = "try-error")) { ft <- NA }
       
      if (is.na(ft) == TRUE){smlen = 25}
      
      # Create the moving window.
      smlen = round(samplefreq/mean(ft$values))
      
      if (smlen > 25){smlen = 25}
      
      print("")
    }
    
    centreData <- switch(stepmethod,

        "Butterfilter" = {
            ##calculate the normalised frequency
            normFreq <- boundaries / (0.5 * samplefreq)
            ##calculate butterworth filter parameter
            butterFilter <- butter(n = filterorder,
                                    W = normFreq,
                                    type = "pass",
                                    plane = "z")
            ##apply the bandpass filter
            signal::filter(butterFilter, xzSeries / (0.5 * samplefreq)) },

        "Chebyfilter" = {
            normFreq <- boundaries / (0.5 * samplefreq)
            ##calculate Chebyfilter type1.
            ChebyFilter <- cheby1( n = filterorder,
                                   Rp = Rp,
                                   W = normFreq,
                                   type = "pass",
                                   plane = "z")
            ##apply the bandpass filter
            signal::filter(ChebyFilter, xzSeries / (0.5 * samplefreq)) },

        "longrun" = {
            debias(x = xzSeries, fun = runmean, k = samplefreq)},

        "none" = {
            xzSeries})

    #Find the zero crossings
    zeroCrossing <- getZeros(x = centreData - threshold, len = smlen)

    #Find the number of zero crossings
    sumZeroCrossings <- sum(zeroCrossing)

    # find the time of the zero crossings - add one because of diff
    zeroCrossingPoints <- which(zeroCrossing) + 1

    if (! 1 %in% zeroCrossingPoints) {
        zeroCrossingPoints <- c(1, zeroCrossingPoints)
    }

    if (!length(xzSeries) %in% zeroCrossingPoints) {
        zeroCrossingPoints <- c(zeroCrossingPoints, length(xzSeries))
    }

    if (plot.it) {
        plot(x = data[, 1], y = centreData,
            type = "l",
            main = paste("Counted Crossings:", sumZeroCrossings),
            xlab = "", ylab = "Centered X-Z Signal")
        abline(v = data[zeroCrossing, 1], col = 2)
    }

    #find duration of steps
    stepDuration <- diff(data[zeroCrossingPoints, 1, drop = TRUE] * samplefreq)
    #find Segment duration
    SegmentDuration <- data[length(data[,1]),1]-data[1,1]
    # Calculate in minutes
    SegmentDuration=SegmentDuration/60
    res <- numeric(length(fun))
    names(res) <- fun

    if ("count" %in% fun) {
        res["count"] <- sumZeroCrossings
        fun <- fun[fun != "count"]
    }
    if ("mean" %in% fun){
      res["mean"] <- (sumZeroCrossings)/SegmentDuration
      fun <- fun[fun != "mean"]
    }
    for (i in fun) {
        val <- try(get(x = i, mode = "function")(stepDuration))
        if (is(val, class2 = "try-error")) { val <- NA }
        res[i] <- val
    }
    return(res)
}

#' calculate the running mean, bumping the window against the ends of the segment
#' borrowed in part from \pkg{caTools}
#' 
#' @title calculate the running mean
#' @param x numeric vector
#' @param k single integer window size
#' @return numeric vector, debiased
#' @export
#' @keywords internal
#' @examples
#'     x1 <- c(-1, 1, 2, 2, 3, 3, 3, 4, 5)
#'     runmean(x1, k = 3)

runmean <- function(x, k = 2) {
  
  n <- length(x)
  
  if (k <= 1) {
    
    y <- x
    
  } else {
    
    if (k > n) { k <- n }
    
    k2 <- k %/% 2
    
    k1 <- k - k2 - 1
    
    y <- c(sum(x[seq_len(k)]), diff(x, k))
    
    y <- cumsum(y) / k
    
    y <- c(rep(y[1], k1), y, rep(y[length(y)], k2))
    
    y[seq.int(from = k2 + 1, to = n)] <- y[seq_len(n - k2)]
    
    for (i in seq_len(k - 1)) { y[i] <- mean(x[seq_len(i)], na.rm = TRUE) }
  }
  return(y)
}

#' Centres a vector at zero with the same distribution
#'
#' @title debias a vector
#' @param x numeric vector
#' @param fun function to calculate centre (default \code{mean})
#' @param \dots additional arguments to pass to fun
#' @return numeric vector with centre zero
#' @keywords internal
#' @export
#' @examples
#'     x1 <- c(-1, 1, 2, 2, 3, 3, 3, 4, 5)
#'     debias(x = x1)
#'     x2 <- c(-10, -10, -9, -8, -8, -7, -5, -4)
#'     debias(x = x2)
#'     debias(x = x2, fun = runmean, k = 3)
#'     x3 <- c(4, 2, 10, 2, 3 ,6, 7, 6, 1)
#'     debias(x = x3, fun = runmean, k = 4)

debias <- function(x, fun = mean, ...) {
    middle <- fun(x, ...)
    sn <- sign(middle)
    sn[sn == 1] <- -1
    out <- x + sn * middle
    return(out)
}


#' Find positions in a numeric vector where sign has changed on average for more than n elements.
#'
#' @title Find Zero Crossing Points
#' @param x numeric vector
#' @param len natural number window length for sign comparison
#' @return logical vector of length zero
#' @export
#' @keywords internal
#' @examples
#'    x1 <- -7:8
#'    getZeros(x = x1)
#'    x2 <- rep(-1:1, times = 6)
#'    getZeros(x = x2)

getZeros <- function(x, len = 3) {

    if (len < 1) { stop("len must be a natural number greater than 0") }

    matched <- zeros <- rep(FALSE, times = length(x))

    if (length(x) > len) {

        dx <- abs(c(rep(0, times = len), diff(sign(x), lag = len)))

        zeros <- dx == max(dx, na.rm = TRUE)

        if (all(zeros)) {

            zeros <- !zeros

        } else {
            # check for min length len change

            signature <- c(rep(TRUE, times = len - 1), FALSE)

            for (r in len + seq_len(length(zeros) - len)) {

                if (all(zeros[seq.int(from = r - len + 1, to = r)] == signature)) {

                    matched[r - len + 1] <- TRUE
                }
            }
            
            zeros <- matched
        }

    } else { warning("x is shorter than len") }

    return(zeros)
}
