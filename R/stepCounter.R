
#' Function to calculate the number and variance of the steps in the data.
#'
#' @title Step Counter
#' @param data The data to use for calculating the steps. This should be a matrix
#' with columns time, x, y, z.
#' @param samplefreq The sampling frequency of the data, in hertz,
#' when calculating the step number (default 100).
#' @param smlen single integer number of data points used as a window
#' when counting zero crossing events
#' @param  AxesMethod Select which axes to count the steps. \enumerate{
#'     \item 'X'
#'     \item 'Y' (default)
#'     \item 'Z'
#'     \item 'XY'
#'     \item 'XZ'
#'     \item 'YZ'
#'     \item 'XYZ'
#' }
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
#' @param verbose single logical should additional progress reporting be printed at the console? (default FALSE).
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
                        AxesMethod = c("X","Y","Z","XZ","XY","YZ","XYZ"), 
                        filterorder = 4L,
                        threshold = 0.01,
                        boundaries = c(0.15, 1.0),
                        stepmethod = c("Chebyfilter","Butterfilter","longrun","none"),
                        Rp = 0.5,
                        plot.it = FALSE,
                        Centre = TRUE,
                        STFT = FALSE,
                        verbose = FALSE,
                        fun = c("GENEAcount","mean", "GENEAamplitude", "sd", "mad", "GENEAwavelength")) {
    
    if (missing(data)) {stop("data is missing") }
    if (missing(stepmethod)){stepmethod = "Chebyfilter"}
    if (missing(AxesMethod)){AxesMethod = "XZ"}
    if (!is.character(fun)) { stop("fun must be character vector of function names") }
    if (length(fun) < 1L) { stop("fun must name at least one function") }

    # Going to remove the amplitude variable from this step counter. Can come back to this later. 
    res <- numeric(length(fun))
    names(res) <- fun
  
    if ("GENEAamplitude" %in% fun){
      res["GENEAamplitude"] <- 0
      fun <- fun[fun != "GENEAamplitude"]
    }
  
    if ("GENEAwavelength" %in% fun){
      res["GENEAwavelength"] <- 0
      fun <- fun[fun != "GENEAwavelength"]
    }
  
    if ("GENEAdistance" %in% fun){
      res["GENEAdistance"] <- 0
      fun <- fun[fun != "GENEAdistance"]
    }
  
    AxesMethod <- match.arg(arg = AxesMethod)
    method <- match.arg(arg = stepmethod)
    #Create the data for calculating the steps by adding x and z columns
    xzSeries <- switch(AxesMethod,
                       "X" = {data[, 2]},
                       "Y" = {data[, 3]},
                       "Z" = {data[, 4]},
                       "XZ" = {data[, 2] + data[, 4]},
                       "XY" = {data[, 2] + data[, 3]},
                       "YZ" = {data[, 3] + data[, 4]},
                       "XYZ" = {data[, 2] + data[, 3] + data[, 4]}
    )
    
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
      
      if (verbose == TRUE){
        print("smlen value is:");print(smlen)
        }
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
        p <- plot(x = data[, 1], y = centreData,
            type = "l",
            main = paste("Counted Crossings:", sumZeroCrossings),
            xlab = "", ylab = "Centered X-Z Signal")
        p <- p + abline(v = data[zeroCrossing, 1], col = 2)
        print(p)
    }

    #find duration of steps
    stepDuration <- diff(data[zeroCrossingPoints, 1, drop = TRUE] * samplefreq)
    #find Segment duration
    SegmentDuration <- data[length(data[,1]),1]-data[1,1]
    # Calculate in minutes
    SegmentDuration <- SegmentDuration/60
    
    if ("GENEAcount" %in% fun) {
        res["GENEAcount"] <- sumZeroCrossings
        fun <- fun[fun != "GENEAcount"]
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

#' Second Step counting method that uses the Peak and Valley detector after a moving average has been applied to the data. 
#' 
#' @title Step Counter 2

#' @param data The data to use for calculating the steps. This should be a matrix
#' with columns time, x, y, z.
#' @param samplefreq The sampling frequency of the data, in hertz,
#' when calculating the step number (default 100).
#' @param smlen single integer number of data points used as a window
#' when counting zero crossing events
#' @param  AxesMethod Select which axes to count the steps. \enumerate{
#'     \item 'X'
#'     \item 'Y' (default)
#'     \item 'Z'
#'     \item 'XY'
#'     \item 'XZ'
#'     \item 'YZ'
#'     \item 'XYZ'
#' }
#' @param ma.smooth Should a moving average filter be applied to the data. 
#' @param filterorder Order of the moving average filter to apply to the data. This variable will be passed to \code{\link[forecast]{ma}}
#' @param Peak_Threshold Number of values either side of the peak/valley that are higher/lower for the value to qualify as a peak/valley 
#' @param Step_Threshold The difference between a peak, valley then peak or valley, peak then valley to constitute a step.
#' @param sd_Threshold A Threshold used to determine when to calculate steps based on the standard deviation between potential steps. If the standard deviation is above this threshold then the algorithm does not calculate steps.
#' @param magsa_Threshold A Threshold used to determine when to calculate steps based on the mean magnitude of acceleration of a segmentd. If the mean magnitude of the segment is below this threshold then the algorithm does not calculate steps.
#' @param plot.it single logical create plot of data and peak/valley detection points (default \code{FALSE}).
#' @param Centre If Centre set to true (default) then the step counter zeros the xz series before filtering.
#' @param fun character vector naming functions by which to summarize steps.
#' "count" is an internally implemented summarizing function that returns step count.
#' @param verbose single logical should additional progress reporting be printed at the console? (default FALSE).

#' @return Returns a vector with length fun.
#' @export
#' @importFrom signal butter cheby1
#' @examples
#' d1 <- matrix(c(100, 101, -0.79, -0.86,
#'         -0.17, -0.14, 0.53, 0.46),
#'     nrow = 2, ncol = 4)
#' stepCounter(data = d1, stepmethod = "longrun")

stepCounter2 = function(data,
                        samplefreq = 100,
                        smlen = 20L,
                        AxesMethod = c("X","Y","Z","XZ","XY","YZ","XYZ"), 
                        ma.smooth = TRUE,
                        filterorder = 4L,
                        Peak_Threshold = 5, 
                        Step_Threshold = 0.5,
                        sd_Threshold = 150,
                        magsa_Threshold = 0.15,
                        plot.it = FALSE,
                        Centre = TRUE,
                        verbose = FALSE,
                        fun = c("GENEAcount", "GENEAamplitude", "GENEAwavelength",
                                "mean", "sd", "mad")){
  
  if (verbose){print("Stepcounter2 initiated")}
  if (missing(data)) {stop("data is missing") }
  if (missing(AxesMethod)){AxesMethod = "XZ"}
  if (!is.character(fun)) { stop("fun must be character vector of function names") }
  if (length(fun) < 1L) { stop("fun must name at least one function") }
  
  AxesMethod <- match.arg(arg = AxesMethod)
  # Decide on an axes to use. 
  xzSeries <- switch(AxesMethod,
                     "X" = {data[, 2]},
                     "Y" = {data[, 3]},
                     "Z" = {data[, 4]},
                     "XZ" = {data[, 2] + data[, 4]},
                     "XY" = {data[, 2] + data[, 3]},
                     "YZ" = {data[, 3] + data[, 4]},
                     "XYZ" = {data[, 2] + data[, 3] + data[, 4]}
  )
  
  # Centre the xz series about 0.
  if (Centre == TRUE){xzSeries = xzSeries - mean(xzSeries)}
  
  if (ma.smooth == TRUE){
    # xzSeries = (ma(xzSeries, filterorder, centre = TRUE))
    
    # Using the ma function from forecast and filter from stats
    if (abs(filterorder - round(filterorder)) > 1e-08) 
      stop("order must be an integer")
    if (filterorder%%2 == 0) 
      w <- c(0.5, rep(1, filterorder - 1), 0.5)/filterorder
    else w <- rep(1, filterorder)/filterorder
    xzSeries = stats::filter(xzSeries, w)
    
    }
  
  # Now I need to find the peaks and valleys. Plotting these as well
  mavals    = na.omit(xzSeries) # Will this cause issues? - Need to replace with a 0 ideally
  mapeaks   = find_peaks( mavals, m = Peak_Threshold)
  mavalleys = find_peaks(-mavals, m = Peak_Threshold)
 
  # This is to prevent the step counting attempting to count steps when none exsist. 
  if (length(mapeaks) == 0 | length(mavalleys) == 0 ){
    if (verbose == TRUE){print("No peaks or valleys")}
    res <- numeric(length(fun))
    names(res) <- fun
    
    if ("GENEAcount" %in% fun) {
      res["GENEAcount"] <- 0
      fun <- fun[fun != "GENEAcount"]
    }
    if ("mean" %in% fun){
      res["mean"] <- 0
      fun <- fun[fun != "mean"]
    }
    if ("GENEAamplitude" %in% fun){
      res["GENEAamplitude"] <- 0
      fun <- fun[fun != "GENEAamplitude"]
    }
    if ("GENEAwavelength" %in% fun){
      res["GENEAwavelength"] <- 0
      fun <- fun[fun != "GENEAwavelength"]
    }
    
    for (i in fun) {
      val <- try(get(x = i, mode = "function")(0))
      if (is(val, class2 = "try-error")) { val <- NA }
      res[i] <- val
    }
    return(res)
  }
  
  if (plot.it == TRUE){
    p <- plot(mavals, type = "l")
    p <- p + points(mapeaks   + Peak_Threshold , mavals[(mapeaks)], col = "blue", pch = 1, lwd = 5)
    p <- p + points(mavalleys + Peak_Threshold , mavals[(mavalleys)], col = "red" , pch = 1, lwd = 5)
    cat(p)
  }
  
  # Firstly deciding where to start 
  Steps = lastpeak = 0; n = 1;  StepIn = StepAmp = StepInDiff = PeakUsed = ValleysUsed = c()
  
  # If FALSE Valleys first, TRUE Peaks first. Checking the indices of both the peaks found and valleys. 
  # If the other way around then swap mapeaks with mavalleys 
  if (mapeaks[1] > mavalleys[1]){
    Tmp4 = mapeaks
    Tmp5 = mavalleys 
    mapeaks = Tmp5
    mavalleys = Tmp4
  }
  
  # Running through the peaks that have now been found as potential markers for step detection.
  for (i in 1:length(mapeaks)){
    
    # Routine to skip over mapeaks where a peak has already been found 
    if (!is.null(PeakUsed) & !is.null(mapeaks)){
      if (any(PeakUsed == mapeaks[i])){
        # Need to account for a breakpoint here. 
        # Once the peak has been used I need to move to the next one. 
        next
      }
    }
    
    # Next if the last peak was greater than the current peak
    if (lastpeak > mapeaks[i]){
      next
    }
    
    posvalleys = mavalleys[mavalleys < (mapeaks[i] + smlen) & mavalleys > mapeaks[i]]
    
    if (verbose == TRUE){print(posvalleys)}
    
    if (length(posvalleys) == 0){
      if (verbose == TRUE){
        print("posvalleys has length 0")
      }
      next
    }
    
    # Checking the corresponding valley/peak to check whether the step has a large enough amplitude. 
    for (j in 1:length(posvalleys)){
      # Need this to move back down a for loop 
      if (!is.null(PeakUsed) & !is.null(mapeaks)){
        if (any(PeakUsed == mapeaks[i])){
          # Need to account for a breakpoint here. 
          # Once the peak has been used I need to move to the next one. 
          next
        }
      }
      
      # Has this valley been used already
      if (!is.null(ValleysUsed)){
        if (any(ValleysUsed == posvalleys[j])){
          next
        }
      }
      
      # Check to see that the difference between the peak and valley is above the threshold
      if (abs(mavals[posvalleys[j]] - mavals[mapeaks[i]]) > Step_Threshold){
        
        # Now find all the points that are the window away from the mapeaks to give pospeaks
        pospeaks = mapeaks[mapeaks < (posvalleys[j] + smlen) & mapeaks > (posvalleys[j] + 10) & mapeaks > lastpeak]
        
        # Need to identify which valley I had got to here.
        if (length(pospeaks) == 0){
          next
        } 
        
        # Run through the posvalleys 
        for (k in 1:length(pospeaks)){
          
          if (!is.null(PeakUsed) & !is.null(mapeaks)){
            if (any(PeakUsed == mapeaks[i])){
              # Need to account for a breakpoint here. 
              # Once the peak has been used I need to move to the next one. 
              next
            }
          }
          
          # Checking that the next peak has not already been used. 
          if (!is.null(PeakUsed)){
            if (any(PeakUsed == pospeaks[k])){
              # print("Missing these numbers");print(i)
              next
            }
          }
          
          # Checking that the next valley has not been used. 
          if (!is.null(ValleysUsed)){
            if (any(ValleysUsed == posvalleys[j])){
              # print("Missing these numbers");print(i)
              next
            }
          }
          
          # Check to see wether the difference is above the threshold
          if (abs(mavals[pospeaks[k]] - mavals[posvalleys[j]]) > Step_Threshold){
            Steps  = Steps + 1
            StepIn[n] = posvalleys[j]
            # Needs to be something here to move the valley position along accordingly - So double steps aren't counted
            
            PeakUsed[n] = mapeaks[i] # Add the starting peak. 
            
            # Now need to set the next peak as the end peak used. In this case pospeak is 41 where k = 1
            LastPeak = pospeaks[k]
            
            ValleysUsed[n] = posvalleys[j]
            
            # This gives the difference between recordings in terms of number of recordings. 
            StepInDiff[n] = LastPeak - PeakUsed[n] # Finds the length of the step recorded.
            
            # Now calculating the amplitude of this step between these points. 
            StepAmp[n] = (abs(mavals[(PeakUsed[n])] - mavals[(ValleysUsed[n])]) + abs(mavals[LastPeak] - mavals[(ValleysUsed[n])]))/4
            
            n = n + 1
            next
          }
        }
        next
      }
    }
  }
  
  # Dividing the steps be two
  if (Steps > 1){
    Steps <- floor(Steps/2)
    len = length(StepIn)
    seq_rm = seq(1, len, by = 2)
    # Only remove if there is a vector longer than 1.
    if (length(seq_rm) > 0){
      StepIn = StepIn[-seq_rm]
      StepInDiff = StepInDiff[-seq_rm]
    }
  }
  
  stepDuration <- diff(data[StepIn, 1, drop = TRUE] * samplefreq, na.rm = T)
  # Find Segment duration
  SegmentDuration <- data[length(data[,1]),1] - data[1,1]
  # Calculate in minutes
  SegmentDurationMin <- SegmentDuration / 60
  res <- numeric(length(fun))
  names(res) <- fun
  
  # Finding step.sd and magsa from data.
  if (length(stepDuration) > 0){
    # Looking at step.sd to check whether over 75. 
    
    sd_lim   <- try(sd(stepDuration, na.rm = T))
    magsa_lim <- try(mean((data[,2]^2 + data[,3]^2 + data[,4]^2)^0.5, na.rm = T))
    
    # Invidiaully pull these errors out. 
    if (is(sd_lim, class2 = "try-error") | is.na(sd_lim)){
      sd_lim <- 0
    }
    
    if (is(magsa_lim, class2 ="try-error") | is.na(magsa_lim)){
      magsa_lim <- 0
    }
    
    if (sd_lim > sd_Threshold | magsa_lim < magsa_Threshold){
      if ("GENEAcount" %in% fun) {
        res["GENEAcount"] <- 0
        fun <- fun[fun != "GENEAcount"]
      }
      if ("mean" %in% fun){
        res["mean"] <- 0
        fun <- fun[fun != "mean"]
      }
      if ("GENEAamplitude" %in% fun){
        res["GENEAamplitude"] <- 0
        fun <- fun[fun != "GENEAamplitude"]
      }
      if ("GENEAwavelength" %in% fun){
        res["GENEAwavelength"] <- 0
        fun <- fun[fun != "GENEAwavelength"]
      }
      
      for (i in fun) {
        val <- try(get(x = i, mode = "function")(0))
        if (is(val, class2 = "try-error")) { val <- NA }
        res[i] <- val
      }
    } else{
      # Run over the steps here. 
      if ("GENEAcount" %in% fun) {
        res["GENEAcount"] <- Steps
        fun <- fun[fun != "GENEAcount"]
      }
      if ("mean" %in% fun){
        res["mean"] <- (Steps)/SegmentDurationMin
        fun <- fun[fun != "mean"]
      }
      
      if ("GENEAamplitude" %in% fun){
        res["GENEAamplitude"] <- mean(StepAmp, na.rm = T)
        fun <- fun[fun != "GENEAamplitude"]
      }
      
      if ("GENEAwavelength" %in% fun){
        res["GENEAwavelength"] <- mean(StepInDiff, na.rm = T)
        fun <- fun[fun != "GENEAwavelength"]
      }
      
      if ("GENEAdistance" %in% fun){
        res["GENEAdistance"] <- mean(diff(StepIn, na.rm = T), na.rm =T)
        fun <- fun[fun != "GENEAdistance"]
      }
      
      for (i in fun) {
        val <- try(get(x = i, mode = "function")(stepDuration))
        if (is(val, class2 = "try-error")) { val <- NA }
        res[i] <- val
      }
    } 
  }
  
  return(res)
}

#' Peak detector function 
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
  pks
}



