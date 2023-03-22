
#' @title Classify Data into Categories defined in an rpart GENEA fit
#' 
#' @description Perform classification on segmented GENEActiv bin data using an 
#' rpart GENEA training fit.
#' 
#' @param testfile character string stating path to a GENEActiv bin file, or a folder 
#' containing GENEActiv bin files.
#' @param start Where to start reading observations.
#' @param end Where to end reading observations.
#' @param Use.Timestamps To use timestamps as the start and end time values this has to be set to TRUE. (Default FALSE)
#' @param radians calculate degrees rotation in radians.
#' @param mmap.load Default is (.Machine$sizeof.pointer >= 8). see \code{\link[mmap]{mmap}} for more details
#' @param trainingfit a GENEA rpart object created by \code{\link{createGENEAmodel}} 
#' that gives the decision tree that was fitted from the training data. 
#' These are the parameters used to predict the new data.
#' @param newdata a new data frame that is to be classified (provide instead of testfile). 
#' The data must contain the \code{\link{features}} named in the trainingfit.
#' @param outputname file name root (excluding extension) for saving the 
#' classification output (default "classified").
#' @param outputdir absolute or relative path to directory in which artifacts 
#' (plot and changes files) should be created or \code{NULL} 
#' (default "GENEAclassification"). 
#' @param verbose single logical should additional progress reporting be 
#' printed at the console (default \code{TRUE}).
#' @param allprobs single logical should all estimated probabilities be 
#' reported rather than probability of selected class (default \code{FALSE}).
#' @param setinf single numeric an arbitrary value to replace Inf in calculated 
#' columns or NA to ignore Inf values. (default 100).
#' -setinf is used to replace -Inf. Alternatively, use setinf NULL to leave Inf as is.
#' @param ... other arguments to be passed to \code{\link{dataImport}}, 
#' \code{\link{segmentation}} and other functions
#' @return The function will return the data frame that was provided as newdata with 
#' additional columns. \enumerate{
#'     \item \code{Class}, a factor indicating that the predicted category of the segment
#'     \item \code{p.Class}, estimated probability that the prediction is correct
#' }
#' Alternatively, by setting argument allprobs to TRUE, a column constructed 
#' 'p.level' containing the estimated probability of each possible class 
#' will be returned instead.
#' @param datacols a vector constructed 'column.summary' or 'default'. See \code{\link{segmentation}} for details.
#' @param penalty single characgter, the penalty to use for changepoint detection. default ("SIC")
#' @param pen.value1 Value of the type 1 error required when penalty is "Asymptotic".
#' @param pen.value2 Default set as NULL and so equals pen.value1 if no input. 
#' @param mininterval single numeric that defines the smallest changepoint initially found. Passed to \code{\link[changepoint]{cpt.var}} as the variable minseglen
#' @param intervalseconds An integer number of seconds between 5 and 30 during which at most one changepoint may occur.
#' @param plot.it (logical) Creates a plot showing the zero crossings counted by the step counting algorithm#' @param Centre Centres the xz signal about 0 when set to True.
#' @param plot.seg (logical) Creates a plot displaying the changepoint locations
#' @param plot.seg.outputfile The name of the png file created that shows the change points on a positionals plots.
#' @param changepoint defines the change point analysis to use. UpDownDegrees performs the change point analysis on the variance of arm elevation and wrist rotation. TempFreq performs a change point on the variance in the temeprature and Frequency (Typically better for sleep behaviours) 
#' @param samplefreq The sampling frequency of the data, in hertz,
#' when calculating the step number. (default 100)
#' @param boundaries to passed to the filter in the step counting algorithm.
#' @param Rp the decibel level that the cheby filter takes. see \code{\link[signal]{cheby1}}
#' @param filterorder The order of the filter applied with respect to the cheby options. 
#' @param hysteresis The hysteresis applied after zero crossing. (default 100mg)
#' @param stft_win numeric for the window to calculate the frequency of an event using the \code{\link[GENEAread]{stft}} function.
#' @details This function will apply the rules determined by the rpart GENEA 
#' decision tree passed to argument trainingfit to the columns 
#' of newdata to classify into classes 
#' (view using \code{"\link[=levels.GENEA]{levels}"}).
#' @export
#' @examples 
#' ## segData <- read.csv(system.file(package = "GENEAclassify", 
#' ##       "testdata", "trainingData9.csv"))
#' ## The training fit can be created by provided the file path to the training data
#' ## in the function getTrainingData - see the help file for more details
#' ## Uses the fitted decision tree to predict the segmented data
#' ## class9 <- classifyGENEA(testfile = "trainingData9.csv",
#' ##                         newdata = segData, 
#' ##                         trainingfit = trainingFit)
#' ## head(class9)
#' ## table(class9$Class)

classifyGENEA <- function(testfile, 
                          start = NULL,
                          end = NULL,
                          Use.Timestamps = FALSE,
                          radians = FALSE,
                          mmap.load = (.Machine$sizeof.pointer >= 8),
                          trainingfit = trainingFit, 
                          newdata, 
                          allprobs = FALSE, 
                          setinf = 100,
                          outputname = "_classified", 
                          outputdir = "GENEAclassification", 
                          datacols = "default",
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
                          # Step Counter Variables 
                          samplefreq = 100,
                          filterorder = 2,
                          boundaries = c(0.5, 5), 
                          Rp = 3,
                          plot.it = FALSE,
                          hysteresis = 0.1, 
                          stft_win = 10, 
                          # Plots
                          plot.seg = FALSE, 
                          plot.seg.outputfile = "Changepoint",
                          verbose = TRUE,
                          ...) {
  
  #### 1. Exceptions ####
  
  if (!is(trainingfit, "rpart")) { stop("trainingfit should be class rpart") } 
  
  if (!is.null(outputname)) {
    if (length(outputname) != 1L && !is.character(outputname)) { 
      stop("outputname should be single character")
    }
    if (length(outputdir) != 1L && !is.character(outputdir)) { 
      stop("outputdir should be single character")
    }
  }
  if (length(verbose) != 1L && !is.logical(verbose)) { 
    stop("verbose should be single logical") }
  
  if (length(allprobs) != 1L && !is.logical(allprobs)) { 
    stop("allprobs should be single logical") }
  
  if (!is.null(setinf)) {
    if (length(setinf) != 1L && !(is.numeric(setinf) || is.logical(setinf))) { 
      stop("setinf should be single numeric") }
  }
  
  if (is.null(pen.value2)) {pen.value2 = pen.value1}
  
  fnames <- features(trainingfit)
  
  #### 2. getGENEAsegments ####
  
  if (missing(newdata)) {
    
    if (!missing(testfile)) {
      
      if (!all(is.character(testfile))) { 
        stop("testfile should be class character") }
      
      # Ensure variables are being passed correctly
      if (missing(changepoint)) {changepoint = "UpDownMeanVarDegreesMeanVar"}
      
      newData <- getGENEAsegments(testfile = testfile,
                                  start = start, 
                                  end = end, 
                                  Use.Timestamps = Use.Timestamps,
                                  radians = radians, 
                                  mmap.load = mmap.load,
                                  outputtoken = "_segmented", 
                                  outputdir = outputdir, 
                                  datacols = datacols,
                                  # Changepoint variables 
                                  changepoint = changepoint,
                                  penalty = penalty,
                                  pen.value1 = pen.value1,
                                  pen.value2 = pen.value2,
                                  intervalseconds = intervalseconds,
                                  mininterval = mininterval,
                                  # StepCounter Variables
                                  samplefreq = samplefreq,
                                  filterorder = filterorder,
                                  boundaries = boundaries, 
                                  Rp = Rp,
                                  plot.it = plot.it,
                                  hysteresis = hysteresis,
                                  stft_win = stft_win, 
                                  # Plots 
                                  plot.seg = plot.seg,
                                  plot.seg.outputfile = "Changepoint",
                                  verbose = verbose)
      
    } else { 
      stop("testfile and newdata are missing but one should be provided") }
    
  } else { 
    
    if (!is(newdata, "data.frame")) { 
      stop("newdata must be a data.frame with columns: ", 
           paste(fnames, collapse = ", ")) }
    
    newData <- newdata 
  }
  
  dnames <- colnames(newData)
  
  missingcols <- !fnames %in% dnames
  
  if (any(missingcols)) { 
    stop("trainingfit features missing from newdata: ", 
         paste(fnames[missingcols], collapse = ", "))
  }
  
  # note that columns required by the fit must be provided, 
  # even if they are not used in the final model
  
  provided <- names(trainingfit$variable.importance)
  
  unused <- provided[!provided %in% fnames]
  
  missingprov <- !unused %in% dnames
  
  if (any(missingprov)) {
    
    dummyRows <- matrix(as.numeric(NA), 
                        nrow = nrow(newData), 
                        ncol = length(missingprov))
    
    colnames(dummyRows) <- unused
    
    newData <- cbind(newData, as.data.frame(dummyRows))
  }
  
  # Create the outputs directory 
  if (is.null(outputdir)) {
    
    outputdir <- tempdir()
    
    if (verbose) { cat("file output location: ", outputdir, "\n") }
  }
  if (!file.exists(outputdir)) {
    
    dir.create(outputdir)
  }
  
  print(outputname)
  
  outputFile <- NULL
  if (outputname != "_classified"){
    # Append the file path to the output file names
    outputname = paste0(outputname ,".csv")
    outputFile <- file.path(outputdir, outputname)
  } else{
    outputname = "_classifed.csv"
  }
  print(outputFile)
  
  classes <- levels(trainingfit)
  
  nLevs <- length(classes)
  
  if (verbose) { 
    cat(paste0("Classes identified in trainingfit\n", 
               paste(paste("    *", classes), collapse = "\n"), "\n"))
  }
  
  # recode Inf
  if (!is.null(setinf)) {
    
    newData[is.na(newData)] <- -99999999
    
    newData[newData == Inf] <- setinf
    
    newData[newData == -Inf] <- -setinf
    
    newData[newData == -99999999] <- NA
  }
  
  #### 3. predict new data from training data tree fit ####
  pred <- predict(object = trainingfit, newdata = newData, type = "class")
  
  ## find probabilities for prediction based on tree fit
  predProb <- predict(object = trainingfit, newdata = newData, type = "prob")
  
  if (allprobs) {
    
    colnames(predProb) <- paste0("p.", colnames(predProb))
    
  } else {
    
    predProb <- cbind(predProb, as.numeric(pred))
    
    predProb <- apply(X = predProb, 
                      MARGIN = 1, 
                      FUN = function(x) { x[x[length(x)]] })
    
    predProb <- data.frame('p.Class' = predProb)
  }
  
  ## create return by adding new columns to input data
  classifiedData <- cbind(newData, Class = pred, predProb)
  
  # remove dummy columns
  if (any(missingprov)) {
    classifiedData[unused] <- NULL
  }
  
  ## Taking the file name and adding the token to the end
  ff <- testfile 
  whereDirChars <- unlist(gregexpr("(\\\\|/)", ff))
  
  whereStart <- 1
  
  if (length(whereDirChars) > 1) {
    whereDirChars <- whereDirChars[length(whereDirChars)] }
  
  if (whereDirChars > 0) {
    whereStart <- whereDirChars +  1 }
  
  shortName <- substr(ff, start = whereStart, stop = nchar(ff) - 4)
  
  outName <- paste0(shortName, outputname)
  
  if (!is.null(outputFile)){
    write.csv(x = classifiedData, file = file.path(outputdir, outputname))
  } else {
    write.csv(x = classifiedData, file = file.path(outputdir, outName))
  }
  
  ## return classified data
  return(classifiedData)
}

