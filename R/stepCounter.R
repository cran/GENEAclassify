
#' Function to calculate the number and variance of the steps in the data.
#'
#' @title Step Counter
#' @param AccData The data to use for calculating the steps. This should either an AccData object or a vector.
#' @param samplefreq The sampling frequency of the data, in hertz,
#' when calculating the step number (default 100).
#' @param filterorder single integer, order of the Chebyshev bandpass filter,
#' passed to argument n of \code{\link[signal]{cheby1}}.
#' @param boundaries length 2 numeric vector specifying lower and upper bounds
#' of Chebychev filter (default \code{c(0.5, 5)} Hz),
#' passed to argument W of \code{\link[signal]{butter}} or \code{\link[signal]{cheby1}}.
#' @param Rp the decibel level that the cheby filter takes, see \code{\link[signal]{cheby1}}.
#' @param plot.it single logical create plot of data and zero crossing points (default \code{FALSE}).
#' @param hysteresis The hysteresis applied after zero crossing. (default 100mg)
#' @param fun character vector naming functions by which to summarize steps.
#' "count" is an internally implemented summarizing function that returns step count.
#' @param verbose single logical should additional progress reporting be printed at the console? (default FALSE).
#' @return Returns a vector with length fun.
#' @export
#' @importFrom signal cheby1
#' @examples
#' d1 <- sin(seq(0.1, 100, 0.1))/2 + rnorm(1000)/10 + 1
#' Steps4 = stepCounter(d1)
#' length(Steps4)
#' mean(Steps4)
#' sd(Steps4)
#' plot(Steps4)

stepCounter <- function(AccData, 
                        samplefreq = 100,
                        filterorder = 2,
                        boundaries = c(0.5, 5), # 
                        Rp = 3,
                        plot.it = FALSE, 
                        hysteresis = 0.1, 
                        verbose = verbose,
                        fun = c("GENEAcount","mean", "sd", "mad")) {
  
  if (missing(AccData)) { stop("data is missing") }
  if (!is.character(fun)) { stop("fun must be character vector of function names") }
  if (length(fun) < 1L) { stop("fun must name at least one function") }
  
  # Going to remove the amplitude variable from this step counter. Can come back to this later. 
  res <- numeric(length(fun))
  names(res) <- fun
  
  #### Check whether an AccData obejct or a vector ####
  
  # If Not an AccData object is it a numerical vector (Can be timestamps and a vector)
  if (class(AccData) == "AccData"){
    StepData = AccData$data.out[,3]
    samplefreq = AccData$freq
  } else if (class(AccData) == "numeric"){
    StepData = AccData
    if (missing(samplefreq)){
      warning("No samplefreq is given. samplefreq set to default, samplefreq = 100")
      samplefreq = 100
    }
  } else if (class(AccData) == "matrix" & dim(AccData)[2] == 2 | 
             class(AccData) == "data.frame" & dim(AccData)[2] == 2){
    StepData = AccData[,2]
    if (missing(samplefreq)){
      warning("No samplefreq is given. samplefreq set to default, samplefreq = 100")
      samplefreq = 100
    }
  } else {
    stop("Step Counter must use either an AccData object, Numerical Vector or a 2D Matrix of time and StepData")
  }
  
  Filter <- cheby1(n = filterorder,                               # order of filter
                   Rp = Rp,                                       # ripple of passband
                   W = boundaries/samplefreq,                     # lower then upper frequencies of bandpass
                   type = "pass",
                   plane = "z")
  
  #### Apply the bandpass filter ####
  filteredData = signal::filter(Filter, StepData) 
  
  state = -1                                                       # initialise step state
  interval = 0                                                     # initialise the interval counter
  cadence = numeric(0)                                             # initialise first element of array for intervals
  samples = length(filteredData)                                   # loop through all samples
  
  for (a in 1:samples) {
    if ((filteredData[a] > hysteresis) && (state < 0)){            # new step started
      state = 1                                                    # set the state
      cadence[length(cadence)] = interval + 1                      # write the step interval
      cadence[length(cadence)+1] = 0                               # initialise to record the next step    
      interval = 0                                                 # reset the step counter
    } else if ((-1*filteredData[a] > hysteresis) && (state > 0)) { # hysteresis reset condition met
      state = -1                                                   # reset the state
      interval = interval + 1                                      # increment the interval
    } else {
      interval = interval + 1                                      # increment the interval
    }
    cadence[length(cadence)] = interval                            # capture last part step
  }
  
  cadence = cadence/samplefreq                                     # divide by the sample frequency to get seconds
  
  if ("GENEAcount" %in% fun) {
    res["GENEAcount"] <- length(cadence)
    fun <- fun[fun != "GENEAcount"]
  }
  
  if ("mean" %in% fun) {
    res["mean"] <- 60 /mean(cadence)
    fun <- fun[fun != "mean"]
  }
  
  if ("median" %in% fun) {
    res["median"] <- 60 /median(cadence)
    fun <- fun[fun != "median"]
  }
  
  for (i in fun) {
    val <- try(get(x = i, mode = "function")(cadence))
    if (is(val, class2 = "try-error")) { val <- NA }
    res[i] <- val
  }
  
  return(res)
  
}

#' Peak detect function 
#' 
#' @title Find_Peaks
#' @param x timeseries signal 
#' @param m The number of points either side of the peak to required to be a peak. 
#' @details Finds the peaks and valleys within the signal passed to the function.  
#' @keywords Internal
#' 

find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  return(pks)
}



