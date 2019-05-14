
#' Import and summarize GENEActiv bin data for manual classification.
#'
#' @title import and segment one or more bin files.
#' @param testfile character vector stating path to a GENEActiv bin file, or a folder containing GENEActiv bin files.
#' @param start Where to start reading observations.
#' @param end Where to end reading observations.
#' @param Use.Timestamps To use timestamps as the start and end time values this has to be set to TRUE. (Default FALSE)
#' @param mmap.load Default is (.Machine$sizeof.pointer >= 8). see \code{\link[mmap]{mmap}} for more details
#' @param outputtoken single character string to be appended to the file name
#' for saving the segmenation output (default '_segmentated').
#' @param outputdir The absolute or relative path to directory in which artifacts (plot and changes files) should be created, or NULL
#' (default "GENEAclassification").
#' @param datacols a vector constructed 'column.summary' or 'default'. See \code{\link{segmentation}} for details.
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
#'     \item Step.count = 0
#'     \item Step.sd = 1
#'     \item Step.mean = 0
#'     \item Step.GENEAamplitude = 3
#'     \item Step.GENEAwavelength = 3
#'     \item Step.GENEAdistance = 3
#' }
#' This can be changed by using a named list. e.g decimalplaces = c(Start.Time = 2, Degrees.mean = 4).
#' @param filterWave single logical, should a smoothing filter from \code{\link[waveslim]{wave.filter}} be applied? (default FALSE).
#' @param filtername single character, the name of the wavelet to use for smoothing
#' when filter is TRUE. (default "haar") Passed to \code{link[waveslim]{wave.filter}}.
#' @param j single numeric, the level to which to smooth. Passed to \code{link[waveslim]{wave.filter}} (default 8).
#' @param penalty single characgter, the penalty to use for changepoint detection. default ("SIC").
#' @param pen.value1 Value of the type 1 error required when penalty is "Asymptotic".
#' @param pen.value2 Default set as NULL and so equals pen.value1 if no input. 
#' @param intervalseconds An integer number of seconds between 5 and 30 during which at most one changepoint may occur.
#' @param plot.it single logical, Creates a plot showing the zero crossings counted by the step counting algorithm#' @param Centre Centres the xz signal about 0 when set to True.
#' @param mininterval single numeric that defines the smallest changepoint initially found. Passed to \code{\link[changepoint]{cpt.var}} as the variable minseglen
#' @param plot.seg single logical, Creates a plot displaying the changepoint locations.
#' @param plot.seg.outputfile The name of the png file created that shows the change points on a positionals plots.
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
#' @param changepoint defines the change point analysis to use. UpDownDegrees performs the change point analysis on the variance of arm elevation and wrist rotation. 
#' TempFreq performs a change point on the variance in the temeprature and frequency (Typically better for sleep behaviours).
#' @param samplefreq The sampling frequency of the data, in hertz,
#' when calculating the step number. (default 100).
#' @param boundaries to pass to the filter in the step counting algorithm.
#' @param Rp the decibel level that the cheby filter takes. See \code{\link[signal]{cheby1}}.
#' @param filterorder The order of the filter applied with respect to the butter or cheby options. 
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
#' @return segmented data are returned
#' @export
#' @importFrom grDevices hcl
#' @importFrom graphics axis points
#' @importFrom stats quantile
#' @seealso The returned object can be interrogated with \code{\link[=head.GENEAbin]{head}}.
#' @examples
#' ## segData <- getGENEAsegments(
#' ##     testfile = "subj02_left_wrist_2013-01-03_17-21-06_Being_Driven3.bin",
#' ##     outputdir = file.path(tempdir(), "GENEAclassification"),
#' ##    filter = TRUE, j = 5)
#' ## head(segData)
#' ## list.files(file.path(tempdir(), "GENEAclassification"))

getGENEAsegments <- function(testfile, 
                             start = NULL, 
                             end = NULL, 
                             Use.Timestamps = FALSE,
                             mmap.load = (.Machine$sizeof.pointer >= 8),
                             outputtoken = "_segmented",
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
                             plot.seg = FALSE,
                             plot.seg.outputfile = "Changepoint",
                             Centre = TRUE,
                             STFT = FALSE,
                             win = 10,
                             verbose = FALSE,
                             ...) {

  #### Error Catching ####
    if (!(length(verbose) == 1 && is.logical(verbose))) { stop("verbose should be a single logical") }

    # Ensure variables are being passed correctly
    if (missing(stepmethod)) {stepmethod = "Chebyfilter"} # Set Chebyfilter as the default.
    if (missing(AxesMethod)) {AxesMethod = "XZ"}
    if (missing(changepoint)) {changepoint = "UpDownMeanVarDegreesMeanVar"}
    if (is.null(pen.value2)) {pen.value2 = pen.value1}
  
    # files should exist

    info <- file.info(testfile)

    noInfo <- is.na(info[, "isdir"])

    if (any(noInfo)) { stop("testfile not found: ", paste(testfile[noInfo], collapse = ", ")) }

    aDir <- info[, "isdir"]

    # get files in named folder(s)
    # only import *.bin files

    binPattern <- "*\\.[bB][iI][nN]$"

    if (any(aDir)) {

        tst <- testfile[!aDir]

        for (dd in testfile[aDir]) {

            tst <- c(tst, list.files(path = dd, pattern = binPattern, full.names = TRUE))
        }

        testfile <- tst
    }

    isBin <- grepl(pattern = binPattern, x = testfile)

    if (any(!isBin)) {

        notBin <- "1 file"

        if (sum(!isBin) > 1) { notBin <- paste(sum(!isBin), "files") }

        if (verbose) {
            warning("only GENEA *.bin files can be processed, ignoring ", notBin) }

        if (verbose) {
            apply(matrix(testfile[!isBin]), margin = 1, FUN = cat, "\n") }

        testfile <- testfile[isBin]
    }

    nfile <- length(testfile)

    if (nfile < 1) { stop("testfile refers to ", nfile, " files") }

    if (verbose) {

        ntxt <- paste(nfile, "files...")

        if (nfile == 1) { ntxt <- "1 file..." }

        cat("reading", ntxt, "\n")
    }

    # check output names

    if (!(length(outputtoken) == 1 && is.character(outputtoken))) {
        stop("outputtoken should be a single character") }

    if (!(length(outputdir) == 1 && is.character(outputdir))) {
        stop("outputdir should be a single character") }


    #### check datacols ####

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

    # collect data and perform segmentation

    output <- vector(mode = "list", length = length(testfile))

    names(output) <- testfile

    for (ff in testfile) {

        inDat <- try(dataImport(bindata = ff, 
                                start = start, 
                                end = end, 
                                Use.Timestamps = Use.Timestamps,
                                mmap.load = mmap.load, ...))

        if (!is(inDat, "try-error")) {

            # give segmentation outputs the root file name in the output folder

            whereDirChars <- unlist(gregexpr("(\\\\|/)", ff))

            whereStart <- 1

            if (length(whereDirChars) > 1) {
                whereDirChars <- whereDirChars[length(whereDirChars)] }

            if (whereDirChars > 0) {
                whereStart <- whereDirChars +  1 }

            shortName <- substr(ff, start = whereStart, stop = nchar(ff) - 4)
           
            outName <- paste0(shortName, outputtoken)

            
            #### perform segmentation ####
            segData <- try(segmentation(data = inDat,
                                        outputfile = outName,
                                        outputdir = outputdir,
                                        datacols = dataCols,
                                        decimalplaces = decimalplaces,
                                        filterWave = filterWave,
                                        filtername = filtername,
                                        j = j, 
                                        changepoint = changepoint,
                                        penalty = penalty,
                                        pen.value1 = pen.value1,
                                        pen.value2 = pen.value2,
                                        intervalseconds = intervalseconds,
                                        mininterval = mininterval,
                                        peaks = peaks,
                                        AxesMethod = AxesMethod,
                                        ma.smooth = ma.smooth,
                                        Peak_Threshold = Peak_Threshold, 
                                        Step_Threshold = Step_Threshold,
                                        sd_Threshold = sd_Threshold,
                                        magsa_Threshold = magsa_Threshold,
                                        stepmethod = stepmethod,  
                                        boundaries= boundaries,
                                        samplefreq = samplefreq,
                                        smlen = smlen,
                                        threshold = threshold,
                                        filterorder = filterorder,  
                                        Rp = Rp, 
                                        plot.it = plot.it,
                                        Centre = Centre,
                                        STFT = STFT,
                                        win = win,
                                        verbose = verbose,...))

            
            if (is(segData, "try-error")) {

                warning("error during segmentation of ", ff)

                segData <- NULL
            }
            
        } else {

            warning("error during import of ", ff)

            segData <- NULL
        }

        output[[ff]] <- segData
        
        #### Plot Changepoints ####
        if (plot.seg == TRUE){
          AccData = read.bin(ff, start = start, end = end, ...)
          tmp2 = get.intervals(AccData, start = 0, end = 1, incl.date = T)
          ind = rep(T, nrow(tmp2))
          col = hcl(0:360)
          max.points = 1e6
          
          if(!is.null(plot.seg.outputfile)){
            png(paste0(plot.seg.outputfile,".png"))}
          
          else{
            png(paste0(shortName,"_changepoints",".png"))
            }
  
          plot(convert.time(tmp2[ind,1]) -> x,
               -acos(tmp2[ind,3] / sqrt(rowSums(tmp2[ind,-1]^2)) ) *180/pi +90 ->y, 
               col = col[ floor( length(col)* (sign(-tmp2[ind,2]) * 180 *acos(-tmp2[ind,4] / sqrt(rowSums(tmp2[ind,-c(1,3)]^2)))/pi +180)/360   ) + 1 ],
               ylim = c(-90, 100),xlim=c(min(convert.time(tmp2[,1])),max(convert.time(tmp2[,1]))), 
               xlab = "Time", ylab="Up/Down", pch=".", cex= 2, yaxt = "n"); abline(h = c(-2:2) * 45, lty= 2); axis(2, at = 45 * -2:2)
          
          for (j in 1:length(segData$Start.Time)){
            abline(v = convert.time(segData$Start.Time[j]),col="red")}
          
          points(convert.time(seq(tmp2[1,1] , quantile(tmp2[,1], 0.2), len = 361) )-> tmp, rep(95, 361), col =col[ floor( length(col)* seq(0.999, 0 , len = 361)) +1 ] , pch = "|")
          text(tmp[1], 95, "CCW")
          text(tmp[361], 95, "CW")
          points(tmp[c( 90, 180, 270) +1], rep(95, 3), pch = "|")
          dev.off()
          rm(list=c("AccData","col","max.points","tmp2","ind"))
        }
    }
    
    out <- do.call(what = "rbind", args = output)

    out$Source <- rownames(out)

    rownames(out) <- seq_len(nrow(out))

    out <- out[, c("Source", colnames(out)[-ncol(out)])]

    return(out)
}
