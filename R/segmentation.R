
#' @name GENEAclassify-package
#' @description This package provides tools to perform the segmentation
#' and classification of GENEActiv accelorometer data. The high frequency
#' time series data is split into segments based on activity changepoints.
#' An \pkg{rpart} fit is trained against known activities at each segment.
#' This fit can then then be used to guess behaviours from test data when
#' activity at each time point has not been reported. This allows detailed
#' behaviour profiles to be created for the wearer.
#' @docType package
#' @title Classification of accelorometer data
#' @keywords package
#' @import MASS 

NULL

#' @title Perform Segmentation on GENEActiv Accelerometer Data
#'
#' @description Perform segmentation of Activinsights accelerometer data.
#' Data are smoothed by the second, or by 10 data points,
#' whichever number of records is greater.
#'
#' Filtering is performed by tools from \pkg{waveslim}.
#' Options are passed to \code{\link[waveslim]{wavelet.filter}}.
#'
#' @param data the GENEActiv bin object to be segmented which should be the output
#' of the \code{\link{dataImport}} function.
#' @param outputfile single character, file name for saving the segmentation output as CSV
#' (and if plot.it is TRUE, corresponding plot PNG). If NULL, create no files.
#' @param outputdir single character, the absolute or relative path to directory in which
#' plot and changes files) should be created, or NULL
#' (default "GENEAclassification"). Ignored if outputfile is NULL.
#' @param datacols character vector constructed 'column.summary'.
#' This object specifies the data and summary to output for the classification.
#' The first part of each element must name column in the GENEAbin datasets specified by filePath.
#' Derived columns may also be selected:\itemize{
#'     \item Step (zero-crossing step counter method),
#'     \item Principal.Frequency.
#' }
#' The second should be the name of a function that evaluates to lenth one.
#' The functions must contain only alphabetical characters
#' (no numbers, underscores or punctuation).
#' The default matrix, specified using the length 1 character vector
#' \code{'default'} is: \itemize{
#'     \item UpDown.mean
#'     \item UpDown.sd
#'     \item UpDown.mad
#'     \item Degrees.mean
#'     \item Degrees.var
#'     \item Degrees.sd
#'     \item Light.mean
#'     \item Light.max
#'     \item Temp.mean
#'     \item Temp.sumdiff
#'     \item Temp.meandiff
#'     \item Temp.abssumdiff
#'     \item Temp.sddiff
#'     \item Magnitude.mean
#'     \item Step.GENEAcount
#'     \item Step.sd
#'     \item Step.mean
#'     \item Step.GENEAamplitude
#'     \item Step.GENEAwavelength
#'     \item Principal.Frequency.median
#'     \item Principal.Frequency.mad
#'     \item Principal.Frequency.ratio
#'     \item Principal.Frequency.sumdiff
#'     \item Principal.Frequency.meandiff
#'     \item Principal.Frequency.abssumdiff
#'     \item Principal.Frequency.sddiff
#' }
#' @param decimalplaces named numeric vector of decimal places with which to
#' round summary columns. \code{NULL} returns unrounded values.
#' The length 1 character vector 'default' applies default roundings: \itemize{
#'     \item Start.Time = 0,
#'     \item Degrees.mean = 3,
#'     \item Degrees.median = 3,
#'     \item Degrees.var = 3,
#'     \item Degrees.sd = 3,
#'     \item Degrees.mad = 3,
#'     \item Magnitude.mean = 3,
#'     \item UpDown.mean = 3,
#'     \item UpDown.median = 3,
#'     \item UpDown.var = 3
#'     \item UpDown.sd = 3,
#'     \item UpDown.mad = 3,
#'     \item Principal.Frequency.median = 3,
#'     \item Principal.Frequency.mad = 3,
#'     \item Principal.Frequency.ratio = 3,
#'     \item Principal.Frequency.sumdiff = 3,
#'     \item Principal.Frequency.meandiff = 3,
#'     \item Principal.Frequency.abssumdiff = 3,
#'     \item Principal.Frequency.sddiff = 3,
#'     \item Light.mean = 0,
#'     \item Light.max = 0,
#'     \item Temp.mean = 1,
#'     \item Temp.sumdiff = 3
#'     \item Temp.meandiff = 3
#'     \item Temp.abssumdiff = 3
#'     \item Temp.sddiff = 3
#'     \item Step.GENEAcount = 0
#'     \item Step.sd = 1
#'     \item Step.mean = 0
#'     \item Step.GENEAamplitude = 3
#'     \item Step.GENEAwavelength = 3
#'     \item Step.GENEAdistance = 3
#' }
#' @param filterWave single logical, should a smoothing filter from \code{\link[waveslim]{wave.filter}} be applied? (default FALSE).
#' @param filtername single character, the name of the wavelet to use for smoothing
#' when filter is TRUE. (default "haar") Passed to \code{link[waveslim]{wave.filter}}.
#' @param j single numeric, the level to which to smooth. Passed to \code{link[waveslim]{wave.filter}} (default 8).
# Changepoint varaibles 
#' @param penalty single characgter, the penalty to use for changepoint detection. default ("SIC").
#' @param pen.value1 Value of the type 1 error required when penalty is "Asymptotic".
#' @param pen.value2 Default set as NULL and so equals pen.value1 if no input. 
#' @param intervalseconds An integer number of seconds between 5 and 30 during which at most one changepoint may occur.
#' @param plot.it single logical, Creates a plot showing the zero crossings counted by the step counting algorithm#' @param Centre Centres the xz signal about 0 when set to True.
#' @param mininterval single numeric that defines the smallest changepoint initially found. Passed to \code{\link[changepoint]{cpt.var}} as the variable minseglen
#' @param changepoint defines the change point analysis to use. UpDownDegrees performs the change point analysis on the variance of arm elevation and wrist rotation. 
#' TempFreq performs a change point on the variance in the temeprature and frequency (Typically better for sleep behaviours).

# Step Counter variables 
#' @param  AxesMethod Select which axes to count the steps. \enumerate{
#'     \item 'X'
#'     \item 'Y' (default)
#'     \item 'Z'
#'     \item 'XY'
#'     \item 'XZ'
#'     \item 'YZ'
#'     \item 'XYZ'
#' }
#' @param Centre single logical, Centres the xz signal about 0 (default TRUE) when counting the zero crossings within the step counting algorithm.
#' @param STFT If STFT is TRUE then the Step Counter uses the STFT function to find the length of the window for each segment.
#' @param win The window length at which to compute the STFT for the changepoint analysis. See \code{\link[GENEAread]{stft}}.
#' @param smlen defines the window length used within the step counting alogirthm.
#' @param threshold Threshold for the step counter to register a step. 
#' @param stepmethod defines the method used by the step counting algoirthm, see \code{\link[GENEAclassify]{stepCounter}} for details.
#' @param samplefreq The sampling frequency of the data, in hertz,
#' when calculating the step number. (default 100).
#' @param boundaries to pass to the filter in the step counting algorithm.
#' @param Rp the decibel level that the cheby filter takes. See \code{\link[signal]{cheby1}}.
#' @param filterorder The order of the filter applied with respect to the butter or cheby options if stepCounter is used. The order of the moving average filter if step counter 2 is used.
#' See \code{\link[signal]{cheby1}} or \code{\link[signal]{butter}}.
#' @param peaks single logical to indicate which step counter to use. If TRUE \code{\link[GENEAclassify]{stepCounter2}} will be used,
#' if FALSE \code{\link[GENEAclassify]{stepCounter}} will be used. (default TRUE).
#' @param ma.smooth Should a moving average filter be applied to the data. 
#' @param Peak_Threshold Number of values either side of the peak/valley that are higher/lower for the value to qualify as a peak/valley 
#' @param Step_Threshold The difference between a peak, valley then peak or valley, peak then valley to constitute a step.
#' @param sd_Threshold A Threshold used to determine when to calculate steps based on the standard deviation between potential steps. If the standard deviation is above this threshold then the algorithm does not calculate steps.
#' @param magsa_Threshold A Threshold used to determine when to calculate steps based on the mean magnitude of acceleration of a segmentd. If the mean magnitude of the segment is below this threshold then the algorithm does not calculate steps.
#' @param verbose single logical should additional progress reporting be printed at the console? (default TRUE).
#' @param ... other arguments to be passed to \code{\link{dataImport}},
#' \code{\link{segmentation}} and other functions with these functions.
#' @details Performs the segmentation procedure on the provided elevation data.
#' Optionally a wavelet filter is first applied to smooth the data.
#' The number of changes occuring in a given number of seconds may be controlled using the
#' intervalseconds argument. Changes will be removed based on which segments are the closest match in terms of variance.
#' A series of features for each of the segments will then be calculated and returned as a csv file.
#' @return The segment data are returned invisibly. This data frame has columns:\itemize{
#'     \item Serial.Number
#'     \item Start.Time
#'     \item Segment.Start.Time
#'     \item Segment.Duration
#'     \item UpDown.median
#'     \item UpDown.var
#'     \item Degrees.median
#'     \item Degrees.mad
#' }
#' In addition, the requested columns are included.
#' Optionally, as a side effect a csv file is returned listing all of the segments
#' found in the data along with a variety of features for that segment.
#' Optionally a png file plotting the data and the
#' detected changes can also be produced.
#' @export
#' @importFrom grDevices dev.off png
#' @examples
#' ### Load the data to segment keeping only the first quarter of the data
#' ## library(GENEAread)
#' ## testfile = file.path(system.file(package = "GENEAread"),
#' ##                                  "binfile",
#' ##                                  "TESTfile.bin")
#' ## segData <- dataImport(binfile = testfile,
#' ##     downsample = 100, start = 0, end = 0.25)
#' ## head(segData)
#' ### Run loaded data through segmentation function
#' ## segment <- segmentation(data = segData, outputfile = NULL)
#' ## head(segment)
#' ## segment2 <- segmentation(data = segData, outputfile = NULL,
#' ##     datacols = "Degrees.skew")
#' ## head(segment2)

segmentation <- function(data,
                         outputfile = "detectedChanges",
                         outputdir = "GENEAclassification",
                         datacols = "default",
                         decimalplaces = "default",
                         filterWave = FALSE,
                         filtername = "haar",
                         j = 8,
                         # Changepoint variables 
                         changepoint = c("UpDownDegrees", "TempFreq", "UpDownFreq", 
                                         "UpDownMean", "UpDownVar", "UpDownMeanVar",
                                         "DegreesMean", "DegreesVar", "DegreesMeanVar", 
                                         "UpDownMeanVarDegreesMeanVar", 
                                         "UpDownMeanVarMagMeanVar"), 
                         penalty = "Manual",
                         pen.value1 = 40,
                         pen.value2 = 400,
                         intervalseconds = 30,
                         mininterval = 5,
                         # Step Coutner 2 variables 
                         peaks = TRUE,
                         AxesMethod = c("X","Y","Z","XZ","XY","YZ","XYZ"), 
                         ma.smooth = TRUE,
                         Peak_Threshold = 3, 
                         Step_Threshold = 0.5,
                         sd_Threshold = 150,
                         magsa_Threshold = 0.15,
                         # Step Counter 1 Variables
                         stepmethod = c("Chebyfilter","Butterfilter","longrun","none"),
                         boundaries = c(0.15, 1.0), 
                         samplefreq = 100,
                         smlen = 100L,
                         threshold = 0.001,
                         filterorder = 4L,  
                         Rp = 0.5, 
                         plot.it = FALSE,
                         Centre = TRUE,
                         STFT = FALSE,
                         win = 10,
                         verbose = FALSE,
                         ...) {
  
    if (missing(data)) { stop("data is missing") }
    if (missing(stepmethod)) {stepmethod = "Chebyfilter"} # Set Chebyfilter as the default.
    if (missing(AxesMethod)) {AxesMethod = "XZ"} # Set the Y-axis as the default
    if (missing(changepoint)) {changepoint = "UpDownMeanVarDegreesMeanVar"} 
    if (is.null(pen.value2)) {pen.value2 = pen.value1}
  
    changepoint <- match.arg(arg = changepoint)
    
  
    if (!is(data, class2 = "list")) {
        stop("data should be a GENEAbin list") }

    dataNames <- c("Data", "UpDown", "Degrees", "Time", "Light", "Temp",
        "Magnitude", "Serial", "RawData", "Freq")

    missNames <- !(dataNames %in% names(data))

    if (any(missNames) && verbose) {
        warning("expected columns not present: ",
            paste(dataNames[missNames], collapse = ", ")) }

    ## Create the outputs directory
    if (!is.null(outputfile)) {
        if (is.null(outputdir)) {

            outputdir <- tempdir()

            if (verbose) { cat("file output location: ", outputdir, "\n") }
        }

        if (!file.exists(outputdir)) {

            dir.create(outputdir)
        }
    }
    requiredCols <- c("UpDown.median", "UpDown.mad",
        "Degrees.median", "Degrees.mad")

    if (identical(datacols, "default")) {

      dataCols <- c("UpDown.mean",
                    "UpDown.var",
                    "UpDown.sd",
                    "Degrees.mean",
                    "Degrees.var",
                    "Degrees.sd",
                    "Magnitude.mean",
                    # Frequency Variables
                    "Principal.Frequency.median",
                    "Principal.Frequency.mad",
                    "Principal.Frequency.GENEAratio",
                    "Principal.Frequency.sumdiff",
                    "Principal.Frequency.meandiff",
                    "Principal.Frequency.abssumdiff",
                    "Principal.Frequency.sddiff",
                    # Light Variables
                    "Light.mean", 
                    "Light.max",
                    # Temperature Variables
                    "Temp.mean",
                    "Temp.sumdiff",
                    "Temp.meandiff",
                    "Temp.abssumdiff",
                    "Temp.sddiff",
                    # Step Variables
                    "Step.GENEAcount", 
                    "Step.sd",
                    "Step.mean", 
                    "Step.GENEAamplitude", 
                    "Step.GENEAwavelength",
                    "Step.GENEAdistance")
  
    } else {

        if (!is.character(datacols)) {

            stop("datacols must be a character vector")
        }
        dataCols <- datacols
        
    }
    
    dataCols <- c(requiredCols, dataCols)

    dataFixed <- c("Serial.Number", "Start.Time", "Segment.Start.Time",
        "Segment.Duration")
    
    dataCols <- dataCols[!dataCols %in% dataFixed]
    
    dataCols <- dataCols[!duplicated(dataCols)]
    
    whereFuns <- regexpr("\\.[A-Za-z]+$", dataCols)

    #### Check Functions #######################################################################

    funs <- substr(dataCols, start = whereFuns + 1, stop = nchar(dataCols))
    
    uFuns <- unique(funs)
    
    missFuns <- !sapply(uFuns, exists)
    
    if (any(missFuns)) {
        stop("functions requested by datacols do not exist: ",
            paste(uFuns[missFuns], collapse = ", ")) }

    #### Check Columns #######################################################################

    cols <- substr(dataCols, start = 1, stop = whereFuns - 1)

    uCols <- unique(cols)

    expectCols <- c(names(data), "Principal.Frequency", "Step")

    missCols <- !uCols %in% expectCols

    if (any(missCols) && verbose) {
        warning("columns requested by datacols not present in data: ",
            paste(uCols[missCols], collapse = ", ")) }

    #### Create Summary matrix #######################################################################

    # these matrices control the data objects that are used for analysis
    # and the functions that should be applied to them

    # split into three matices based on whether the columns exist
    # or are derived from raw data

    dataColsMat <- cbind(cols, funs)

    der <- cols %in% c("Principal.Frequency", "Step")

    dataColsMatstep <- dataColsMat[cols == "Step", , drop = FALSE]

    dataColsMatFreq <- dataColsMat[cols == "Principal.Frequency", , drop = FALSE]

    dataColsMatExist <- dataColsMat[!der, , drop = FALSE]

    ##  split up data

    UpDown <- data$UpDown

    Degrees <- data$Degrees

    Time <- data$Time

    Light <- data$Light

    Temp <- data$Temp

    Magnitude <- data$Magnitude

    Serial <- data$Serial

    # Raw Data for Step Counts
    xyzdata <- as.data.frame(data$RawData)

    Freq <- data$Freq

    rm(data)

    #### QC time #######################################################################

    # QC time

    if (!all(Time > 0)) { stop("not all Time values are positive") }

    diffTime <- diff(Time)

    typical <- median(diffTime, na.rm = TRUE)

    isFarFromTypical <- diffTime > typical * 1.01 | diffTime < typical / 1.01

    if (any(isFarFromTypical)  && verbose) {
        warning("time increments vary by more than 1% (",
            sum(isFarFromTypical, na.rm = TRUE), " records)") }

    #### Wavelet filtering #######################################################################

    ## apply optional wavelet filtering

    if (filterWave) { 

        if (requireNamespace(package = "waveslim")) {

            # decompose up-down time series
            modwtUpDown <- try(
                    modwt(x = UpDown,
                        wf = filtername,
                        n.levels = j),
                silent = !verbose)

            if (!is(modwtUpDown, class2 = "try-error")) {

                d <- modwtUpDown
                d[["d1"]] <- d[["d2"]] <- d[["d3"]] <- d[["d4"]] <- d[["d5"]] <- numeric(
                    length = length(UpDown))
                modwtUpDown <- d

                # reconstruct time series
                UpDown <- imodwt(y = modwtUpDown)

            } else {
                warning(
                    paste("filtering will be skipped: error during filtering\n",
                        modwtUpDown))
            }

        } else {
            warning("filtering will be skipped: waveslim not loaded")
        }
    }

    #### Change point method #######################################################################

    ## find changepoints based on updown data and rotation data

    allTimes <- switch(changepoint,
                       
                       "UpDownDegrees" = {
                         
                         changeUpDown <- cpt.var(data = UpDown,
                                                 penalty = penalty,
                                                 pen.value = pen.value1,
                                                 minseglen = mininterval,
                                                 method = "PELT")
                         
                         changeDegrees <- cpt.var(data = Degrees,
                                                  penalty = penalty,
                                                  pen.value = pen.value2,
                                                  minseglen = mininterval,
                                                  method = "PELT")
                         
                         ## find times at segment boundaries
                         
                         changeTimes(time = Time,
                                     intervalseconds = intervalseconds,
                                     changeupdown = changeUpDown,
                                     changedegrees = changeDegrees,
                                     verbose = verbose)},
                       
                       
                       "TempFreq" = {
                         
                         spdata = xyzdata[, c("x", "y", "z")]
                         
                         PFrequency <- sapply(X = spdata, FUN = function(x) {
                           
                           ft <- try(stft(as.matrix(x), quiet = TRUE, win = win,
                                          reassign = TRUE, date.col = FALSE,
                                          freq = Freq)$principals, silent = TRUE)
                           
                           if (is(ft, class2 = "try-error")) { ft <- NA }
                           
                           return(ft)
                         })
                         
                         FreqData= c()
                         
                         for (i in 1:length(PFrequency[,1])){
                           FreqData[i] = mean(PFrequency[i,])
                         }
                         
                         ## Change point analysis using STFT and temperature
                         
                         changeTemp <- cpt.var(data = Temp,
                                               penalty = penalty,
                                               pen.value = pen.value1,
                                               minseglen = mininterval,
                                               method = "PELT")
                         
                         changeSTFT <- cpt.var(data =  FreqData,
                                               penalty = penalty,
                                               pen.value = pen.value2,
                                               minseglen = mininterval,
                                               method = "PELT")
                         
                         changeTimes(time = Time,
                                     intervalseconds = intervalseconds,
                                     changeupdown = changeTemp,
                                     changedegrees = changeSTFT,
                                     verbose = verbose)},
                       
                       "UpDownFreq" = {
                         
                         print("UpDownFreq applied")
                         
                         spdata = xyzdata[, c("x", "y", "z")]
                         
                         PFrequency <- sapply(X = spdata, FUN = function(x) {
                           
                           ft <- try(stft(as.matrix(x), quiet = TRUE, win = win,
                                          reassign = TRUE, date.col = FALSE,
                                          freq = Freq)$principals, silent = TRUE)
                           
                           if (is(ft, class2 = "try-error")) { ft <- NA }
                           
                           return(ft)
                         })
                         
                         FreqData= c()
                         
                         for (i in 1:length(PFrequency[,1])){
                           FreqData[i] = mean(PFrequency[i,])
                         }
                         
                         changeUpDown <- cpt.var(data = UpDown,
                                                 penalty = penalty,
                                                 pen.value = pen.value1,
                                                 minseglen = mininterval,
                                                 method = "PELT")
                         
                         changeSTFT <- cpt.var(data =  FreqData,
                                               penalty = penalty,
                                               pen.value = pen.value2,
                                               minseglen = mininterval,
                                               method = "PELT")
                         
                         changeTimes(time = Time,
                                     intervalseconds = intervalseconds,
                                     changeupdown = changeUpDown,
                                     changedegrees = changeSTFT,
                                     verbose = verbose)
                       },
                       
                       # Changepoint only on the Arm Elevation
                       
                       "UpDownMean" = {
                         UpDownMean = cpt.mean(data = UpDown,
                                          penalty = penalty,
                                          pen.value = pen.value1,
                                          minseglen = mininterval,
                                          method = "PELT")
                         
                         Time[cpts(UpDownMean)]
                         
                       }, 
                       "UpDownVar" = {
                         UpDownVar = cpt.var(data = UpDown,
                                        penalty = penalty,
                                        pen.value = pen.value1,
                                        minseglen = mininterval,
                                        method = "PELT")
                         
                         Time[cpts(UpDownVar)]
                         
                       },
                       
                       "UpDownMeanVar" = {
                         
                         UpDownMeanVar = cpt.meanvar(data = UpDown,
                                                penalty = penalty,
                                                pen.value = pen.value1,
                                                minseglen = mininterval,
                                                method = "PELT")
                         
                         Time[cpts(UpDownMeanVar)]
                         
                       },
                       
                       # Changepoint only on the Wrist Rotation
                       
                       "DegreesMean" = {
                         DegreesMean = cpt.mean(data = Degrees,
                                               penalty = penalty,
                                               pen.value = pen.value1,
                                               minseglen = mininterval,
                                               method = "PELT")
                         
                         Time[cpts(DegreesMean)]
                         
                       }, 
                       "DegreesVar" = {
                         DegreesVar = cpt.var(data = Degrees,
                                             penalty = penalty,
                                             pen.value = pen.value1,
                                             minseglen = mininterval,
                                             method = "PELT")
                         
                         Time[cpts(DegreesVar)]
                         
                       },
                       
                       "DegreesMeanVar" = {
                         
                         DegreesMeanVar = cpt.meanvar(data = Degrees,
                                                     penalty = penalty,
                                                     pen.value = pen.value1,
                                                     minseglen = mininterval,
                                                     method = "PELT")
                         
                         Time[cpts(DegreesMeanVar)]
                         
                       },
                       
                       # Changepoint I reqiure - Might need modification later on. Seems okay right now. 
                       
                       "UpDownMeanVarDegreesMeanVar" = {
                         
                         UpDownMeanVar = cpt.meanvar(data = UpDown,
                                                     penalty = penalty,
                                                     pen.value = pen.value1,
                                                     minseglen = mininterval,
                                                     method = "PELT")
                         
                         DegreesMeanVar = cpt.meanvar(data = Degrees,
                                                      penalty = penalty,
                                                      pen.value = pen.value2,
                                                      minseglen = mininterval,
                                                      method = "PELT")
                         
                         changeTimes(time = Time,
                                     intervalseconds = intervalseconds,
                                     changeupdown = UpDownMeanVar,
                                     changedegrees = DegreesMeanVar,
                                     verbose = verbose)
                         
                       }, 
                       
                       "UpDownMeanVarMagMeanVar" = {
                         
                         UpDownMeanVar = cpt.meanvar(data = UpDown,
                                                     penalty = penalty,
                                                     pen.value = pen.value1,
                                                     minseglen = mininterval,
                                                     method = "PELT")
                         
                         MagMeanVar = cpt.meanvar(data = Degrees,
                                                      penalty = penalty,
                                                      pen.value = pen.value2,
                                                      minseglen = mininterval,
                                                      method = "PELT")
                         
                         changeTimes(time = Time,
                                     intervalseconds = intervalseconds,
                                     changeupdown = UpDownMeanVar,
                                     changedegrees = MagMeanVar,
                                     verbose = verbose)
                       }
    )
    
    
    #### If AllTimes has length 0 then ####
    if (length(allTimes) == 0){
      allTimes = c(Time[1], Time[length(Time)])
      allDurations <- Time[length(Time)] - Time[1]
    } else{
      allDurations <- as.numeric(allTimes[-1] - allTimes[-length(allTimes)])
    }
    
    #### segment durations ####
    output <- data.frame(Serial.Number = Serial,
                         Start.Time = as.integer(allTimes[-length(allTimes)]),
                         Segment.Start.Time = format(x = convert.time(allTimes[-length(allTimes)]),
                                                                      format = "%H:%M:%S"),
                         Segment.Duration = allDurations,
                         stringsAsFactors = FALSE)

    # update allTimes to catch all Time if few allTimes are returned
    if (length(allTimes) < 3) { allTimes <- range(Time) }

    allTimes[c(which.min(allTimes), which.max(allTimes))] <- allTimes[
        c(which.min(allTimes), which.max(allTimes))] + c(-1, 1)

    if (plot.it) {

        ## format Time
        Timeplot <- convert.time(Time)

        TIME <- format(x = Timeplot, format = "%H:%M:%S")

        TIME <- as.POSIXct(x = TIME, format = "%H:%M:%S")

        timeUpDownPlot <- convert.time(sort(c(Time[1], allTimes, Time[length(Time)])))

        TIME2 <- format(x = timeUpDownPlot, format = "%H:%M:%S")

        TIME2 <- as.POSIXct(x = TIME2, format = "%H:%M:%S")

        if (!is.null(outputfile)) {
            outputf <- gsub(pattern = "\\.[CcPSsNVvG]$", replacement = ".png", x = outputfile)
            if (!grepl(pattern = "\\.png$", x = outputf)) { outputf <- paste0(outputf, ".png") }
            png(filename = file.path(outputdir, outputf))
        }

        plot(x = TIME, y = UpDown,
            type = "l",
            xlab = "Time", ylab = "")
        abline(v = as.numeric(TIME2), col = 2)

        if (!is.null(outputfile)) { dev.off() }
    }

    #### summarize existing cols ##############################################

    cutPointEx <- cut(Time, allTimes)

    # split variables

    UpDown <- split(x = UpDown, f = cutPointEx)

    Degrees <- split(x = Degrees, f = cutPointEx)

    Time <- split(x = Time, f = cutPointEx)

    Magnitude <- split(x = Magnitude, f = cutPointEx)

    Light <- split(x = Light, f = cutPointEx)

    Temp <- split(x = Temp, f = cutPointEx)

    existSummary <- summariseCols(colfun = dataColsMatExist)

    output <- cbind(output, existSummary)

    output$Cuts <- as.numeric(rownames(output))

    #### summarize derived cols ###############################################

    ## redefine cutpoint here

    cutPointDe <- cut(xyzdata$timestamp, allTimes)

    #### StepCounter ####

    if (nrow(dataColsMatstep) > 0) {

        spdata <- split(xyzdata[, c("timestamp", "x", "y", "z")], f = cutPointDe)
        # The max smlen can be is 30! Could change in future. 
        if (!"smlen" %in% names(match.call())) { smlen <- 30L }
        
        # Accounting for the smlen change in frequency. Optimised for 100Hz 
        smlen = (smlen * (100 / Freq))

        if (peaks == TRUE){
          stepNumber <- lapply(spdata, function(x, samplefreq, smlen) {
               stepCounter2(x[, c("timestamp", "x", "y", "z")], 
                            samplefreq = samplefreq,
                            smlen = smlen,
                            AxesMethod = AxesMethod, 
                            ma.smooth = ma.smooth,
                            filterorder = filterorder,
                            Peak_Threshold = Peak_Threshold, 
                            Step_Threshold = Step_Threshold,
                            sd_Threshold = sd_Threshold,
                            magsa_Threshold = magsa_Threshold,
                            plot.it = plot.it,
                            Centre = Centre,
                            verbose = verbose,
                            fun = dataColsMatstep[, "funs", drop = TRUE])},
                            samplefreq = max(10, Freq), smlen = smlen)
        }
        else{
          stepNumber <- lapply(spdata, function(x, samplefreq, smlen) {
                stepCounter(x[, c("timestamp", "x", "y", "z")],
                            samplefreq = samplefreq,
                            smlen = smlen,
                            AxesMethod = AxesMethod, 
                            stepmethod = stepmethod,
                            filterorder = filterorder,
                            threshold = threshold,
                            Rp = Rp,
                            boundaries = boundaries,
                            plot.it = plot.it,
                            Centre = Centre,
                            STFT = STFT,
                            fun = dataColsMatstep[, "funs", drop = TRUE])},
                            samplefreq = max(10, Freq), smlen = smlen)
        }

        stepNumber <- as.data.frame(do.call("rbind", stepNumber),
           stringsAsFactors = FALSE)
        colnames(stepNumber) <- paste0("Step.", colnames(stepNumber))

        # TODO check this is always okay
        if (nrow(stepNumber) != nrow(output)) { stop("missmatch in rows when creating step") }
        stepNumber$Cuts <- output$Cuts

        output <- merge(x = output, y = stepNumber)
    }

    #### Principal.Frequency ####

    if (nrow(dataColsMatFreq) > 0) {

        spdata <- split(x = xyzdata[, c("x", "y", "z")], f = cutPointDe)

        Principal.Frequency <- sapply(X = spdata, FUN = function(x) {

                ft <- try(stft(as.matrix(x), quiet = TRUE, type = "mv",
                            reassign = TRUE, date.col = FALSE,
                            freq = Freq)$principals, silent = TRUE)

                if (is(ft, class2 = "try-error")) { ft <- NA }

                return(ft)
            })

        principalSummary <- summariseCols(colfun = dataColsMatFreq)

        principalSummary$Cuts <- as.numeric(rownames(principalSummary))

        output <- merge(x = output, y = principalSummary)
    }

    output$Cuts <- NULL

    if (!is.null(decimalplaces)) {

        if (identical(decimalplaces, "default")) {

            decimalplaces <- c(Start.Time = 0,
                               Degrees.mean = 3, Degrees.median = 3,
                               Degrees.var = 3, Degrees.sd = 3,
                               Degrees.mad = 3, Magnitude.mean = 3,
                               UpDown.mean = 3, UpDown.median = 3,
                               UpDown.var = 3, UpDown.sd = 3, UpDown.mad = 3,
                               Principal.Frequency.median = 3,
                               Principal.Frequency.mad = 3,
                               Principal.Frequency.ratio = 3,
                               Principal.Frequency.sumdiff = 3,
                               Principal.Frequency.meandiff = 3,
                               Principal.Frequency.abssumdiff = 3,
                               Principal.Frequency.sddiff = 3,
                               Light.mean = 0, Light.max = 0,
                               Temp.mean = 2, 
                               Temp.sumdiff = 3,
                               Temp.meandiff = 3,
                               Temp.abssumdiff = 3,
                               Temp.sddiff = 3,
                               Step.count = 0, 
                               Step.sd = 1,
                               Step.mean = 3, 
                               Step.GENEAamplitude = 3, 
                               Step.GENEAwavelength = 3, 
                               Step.GENEAdistance = 3)
        }

        if (!(is(object = decimalplaces, class2 = "numeric") &&
                !is.null(names(decimalplaces)))) {

            stop("decimalplaces must be a named numeric vector")
        }

        decimalplaces <- decimalplaces[decimalplaces %in% names(output)]

        for (nn in seq_along(decimalplaces)) {

            column <- names(decimalplaces)[nn]

            output[, column] <- round(x = output[, column, drop = TRUE],
                digits = decimalplaces[nn])
        }
    }

    if (!is.null(outputfile)) {
        # write out
        if (!grepl(pattern = ".\\.[CcSsVv]$", x = outputfile)) {
            outputfile <- paste0(outputfile, ".csv")
        }

        write.csv(output, file = file.path(outputdir, outputfile), row.names = FALSE)
    }
    if (verbose) { cat("Analysis Complete!\n") }

    return(invisible(output))
}

# appease R CMD check for suggested functions
globalVariables(c("modwt", "imodwt"))

