

#' Loads the data into R and creates format required for segmentation.
#' 
#' @title Data import function
#' @param bindata File path to binary data to be segmented.
#' @param downsample Rate to downsample the data, defaults to every 100th observation. For no downsampling set NULL.
#' @param start Where to start reading observations.
#' @param end Where to end reading observations.
#' @param blocksize data chunk size for \code{read.bin}.
#' @param ... additional arguments passed through.
#' @details Reads in the binary data file and extracts the information required for the segmentation procedure.
#' @return Returns a list containing a matrix of the data including the x, y and z axis data, vectors of the up down (elevation) 
#' and degrees (rotation), a vector of time stamps, a vector of vector magnitudes and the serial number of the device. 
#' @export
#' @import GENEAread
#' @examples
#' ##    segData <- dataImport(bindata = list.files("RunWalk.bin", full = TRUE)[1])
#' ##     names(segData)


dataImport <- function(bindata, downsample = 100, start = NULL, end = NULL, blocksize = 500, ...) {

    binaryData <- read.bin(binfile = bindata, calibrate = TRUE, downsample = downsample, 
        start = start, end = end, blocksize = blocksize)

    if (is.null(downsample)) {
        
        binaryDataOut <- binaryData$data.out
        
    } else {
        
        binaryDataFULL <- read.bin(bindata, calibrate = TRUE, start = start, end = end, 
            blocksize = blocksize)
        
        binaryDataOut <- binaryDataFULL$data.out
        
        rm(binaryDataFULL)
    }
    
    serial <- binaryData$header["Device_Unique_Serial_Code", ][[1]]

    rightWrist <- grepl("right wrist", binaryData$header["Device_Location_Code", ][[1]])
    leftWrist <- grepl("left wrist", 
            binaryData$header["Device_Location_Code", ][[1]]) | grepl("[[:blank:]]", 
            gsub("", " ", binaryData$header["Device_Location_Code", ][[1]]))

    if(!leftWrist & !rightWrist){
        warning("Note: data assumed to be left wrist")
    }

    Intervals <- get.intervals(binaryData, length = NULL, 
        incl.date = TRUE, size = 1e6)

    ## Extract time, light and temp data
    Time <- Intervals[, "timestamp"]
    Light <- binaryData$data.out[, "light"]
    Temp <- binaryData$data.out[, "temperature"]

    rm(binaryData)

    ## extract the up/down and rotation data
    dataUpDown <- updown(Intervals)
    dataDegrees <- degrees(Intervals)

    if (rightWrist) {
        dataUpDown <- -1 * dataUpDown
    }

    
    vecMagnitude <- abs(sqrt(rowSums((Intervals[, c("x", "y", "z")])^2)) - 1)
    
    geneaBin <- list(Data = Intervals, Freq = downsample, UpDown = dataUpDown, Degrees = dataDegrees, 
        Time = Time, Light = Light, Temp = Temp, Magnitude = vecMagnitude, 
        RawData = binaryDataOut, Serial = serial)
    
    class(geneaBin) <- c(class(geneaBin), "GENEAbin")
    
    return(geneaBin)
}


#' Extract up/down time series
#' 
#' Extract code from \code{positionals} to perform data conversion to up/down time series.
#' Input is expected to be result of \code{get.intervals}.
#' @title Extract data relating to the up/down component
#' @param x data output from \code{get.intervals}
#' @return The up/down vertical elevation data (y-axis)
#' @export
#' @keywords internal
#' @examples
#'     d1 <- matrix(c(100, 101, -0.79, -0.86, -0.17, -0.14, 0.53, 0.46), 
#'         nrow = 2, ncol = 4)
#'     colnames(d1) <- c("timestamp", "x", "y", "z")
#'     updown(x = d1)

updown <- function(x) {
    
    numerator <- x[, "y"]
    
    # magnitude
    denominator <- sqrt(rowSums(x[, c("x", "y", "z")]^2)) 
    
    ud <- (-acos(numerator / denominator) * 180 / pi + 90)
    
    return(ud)
}


#' Extract data relating to the rotation component.
#' 
#' Called by \code{dataImport}.
#' Note: the "+ 1" has been removed from the original implementation.
#' @title Extract rotation time series
#' @param x data output from get.intervals
#' @return The degrees (rotation) data.
#' @export
#' @keywords internal
#' @examples
#'    d1 <- matrix(c(100, 101, -0.79, -0.86, -0.17, -0.14, 0.53, 0.46), 
#'         nrow = 2, ncol = 4)
#'    colnames(d1) <- c("timestamp", "x", "y", "z")
#'    degrees(x = d1)

degrees <- function(x) {
    
    magnitude <- sqrt(rowSums(x[, c("x", "z")]^2))
    
    deg <- 361 * (sign(-x[, "x"]) * 180 * acos(-x[, "z"] / magnitude) / pi + 180) / 360
    
    deg <- floor(deg) - 45 
    
    return(deg)
}
