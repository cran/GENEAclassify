

#' @title Select Times as Segment Changepoints
#' 
#' @description Trims down the number of changepoints in time series 
#' data to give segments that can reasonably be classified into 
#' discrete periods of activity. These will not normally be less than
#' 30 seconds in duration.
#' 
#' @details The \pkg{changepoint} package provides tools for optimally 
#' segmenting time series data.
#' @param time numeric vector
#' @param intervalseconds An integer number of seconds (usually greater than 5) 
#' during which at most one changepoint may occur (typically between 5 and 30).
#' If intervalseconds is NA, all times will be returned.
#' If intervalseconds <= mininterval it will be ignored 
#' (with a warning if verbose is \code{TRUE})
#' @param changeupdown cpt object
#' @param changedegrees cpt object
#' @param mininterval single numeric
#' @param verbose single logical should all warnings be reported? (default TRUE)
#' @return numeric vector of times
#' @export
#' @keywords internal
#' @import changepoint
#' @importFrom methods is
#' @importFrom stats as.formula filter mad median na.omit predict sd setNames
#' @examples
#' 
#' library(changepoint)
#' d1 <- c(54, 30, 27, 53, 100, 204, 197)
#' d2 <- c(67, 64, 70, 79, 69, 60, 54)
#' c1 <- cpt.var(d1, penalty = "SIC", pen.value = 1e-3, method = "PELT")
#' c2 <- cpt.var(d2, penalty = "SIC", pen.value = 1e-3, method = "PELT")
#' changeTimes(time = 0:6, intervalseconds = 30, 
#'     changeupdown = c1, changedegrees = c2)

changeTimes <- function(time, intervalseconds = 30, 
    changeupdown, changedegrees, mininterval = 5, verbose = TRUE) {
    
    if (missing(time)) { stop("time is missing") }
    if (!is.numeric(time)) { 
        stop("time must be numeric, but was ", class(time)) }
    if (!is(object = changeupdown, class2 = "cpt")) {
        stop("changeupdown must be a cpt object, but was ", 
            class(changeupdown))
    }
    if (!is(object = changedegrees, class2 = "cpt")) {
        stop("changedegrees must be a cpt object, but was ", 
            class(changedegrees))
    }
    ## extract change points for both updown and rotation
    cptUpDown <- cpts(changeupdown)
    cptDegrees <- cpts(changedegrees)
    
    allChanges <- unique(sort(c(cptUpDown, cptDegrees)))
    
    ## get time of changes
    timeUpDown <- c(time[1], time[sort(cptUpDown)], time[length(time)])
    
    timeDegrees <- c(time[1], time[sort(cptDegrees)], time[length(time)])

    allTimes <- c(time[1], time[allChanges], time[length(time)])

    durationUpDown <- diff(timeUpDown)
    durationDegrees <- diff(timeDegrees)
    
    if (is.null(intervalseconds)) { intervalseconds <- NA }
    
    if (length(intervalseconds) != 1L) {
        warning("intervalseconds must be length 1, but was length ", 
            length(intervalseconds))
        intervalseconds <- intervalseconds[1]
    }
    if (!is.na(intervalseconds)) {
        
        if (intervalseconds > mininterval) {
            
            shortdurationUpDown <- which(as.numeric(durationUpDown) < intervalseconds)
            
            shortdurationDegrees <- which(as.numeric(durationDegrees) < intervalseconds)
            
            upDownSegments <- removeShortSegments(
                shortduration = shortdurationUpDown, 
                changes = cptUpDown, 
                variance = param.est(changeupdown)$variance, 
                time = time)
                
            degreesSegments <- removeShortSegments(
                shortduration = shortdurationDegrees, 
                changes = cptDegrees, 
                variance = param.est(changedegrees)$variance, 
                time = time)
            
            allChanges <- sort(c(
                    setNames(upDownSegments$cpts, rep("U", length = length(upDownSegments$cpts))), 
                    setNames(degreesSegments$cpts, rep("D", length = length(degreesSegments$cpts)))))
            
            allChanges <- setNames(unique(allChanges), names(allChanges)[!duplicated(allChanges)])
            allTimes <-  c(time[1], time[allChanges], time[length(time)])
            allDurations <- as.numeric(allTimes[-1] - allTimes[-length(allTimes)] )
            
            shortduration <- which(as.numeric(allDurations) < intervalseconds)
            shortduration <- unique(ifelse(test = names(allChanges)[shortduration] == "D", 
                yes = shortduration, 
                no = shortduration + 1))
            shortduration <- na.omit(shortduration)
            
            if (length(shortduration) > 0) { 
                allChanges <- allChanges[-shortduration] }
            
            allTimes <-  c(time[1], time[allChanges], time[length(time)])
            
            # If allChanges is not null
            
            if (length(allChanges) > 0 ){
            
                allTimeDiff <- diff(allTimes)
                
                shortalltimes <- which(as.numeric(allTimeDiff) < intervalseconds)
                
                FinalVariance = length(allChanges) - 1
                
                for (i in 1:length(allChanges)){
                  if (names(allChanges[i]) == "D"){
                    FinalVariance[i] <- param.est(changedegrees)$variance[i]
                  }
                  if (names(allChanges[i]) == "U"){
                    FinalVariance[i] <- param.est(changeupdown)$variance[i]
                  }
                }
                
                allSegments <- removeShortSegments(
                    shortduration = shortalltimes,
                    changes = unname(allChanges), 
                    variance = param.est(changeupdown)$variance, 
                    time = time)
                
                allTimes = allSegments$times
            }
        } else {
            if (verbose) { 
                warning("ignoring intervalseconds <= ", mininterval, " and returning all segments") }
        }
    }
    
    return(allTimes)
}

