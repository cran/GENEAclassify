## ----global_options, warning = FALSE, eval = FALSE, echo = FALSE--------------
#  knitr::opts_chunk$get("root.dir")

## ----installing the dependencies, eval = FALSE--------------------------------
#  
#  install.packages("GENEAread", repos = "http://cran.us.r-project.org")
#  install.packages("changepoint", repos = "http://cran.us.r-project.org")
#  install.packages("signal", repos = "http://cran.us.r-project.org")
#  install.packages("mmap", repos = "http://cran.us.r-project.org")
#  
#  # Load in the libraries
#  library(GENEAread)
#  library(changepoint)
#  library(signal)
#  library(mmap)

## ----Installing from Source, eval = FALSE-------------------------------------
#  # You will need to change the folder location inside setwd("") to the directory where you saved the tar.gz file
#  # Note that R only uses / not \ when refering to a file/directory location
#  setwd("/Users/owner/Documents/GENEActiv")
#  install.packages("GENEAclassify_1.5.1.tar.gz", repos=NULL, type="source")

## ----loading in the GENEAclassify library, eval = FALSE-----------------------
#  library(GENEAclassify)

## ----installing from GitHub, eval = FALSE-------------------------------------
#  install.packages("devtools",repos = "http://cran.us.r-project.org")
#  library(devtools)
#  
#  install_github("https://github.com/Langford/GENEAclassify_1.41.git",
#                 auth_token = "7f0051aaca453eaabf0e60d49bcf752c0fea0668")
#  

## ----Run library function again GENEAclassify library,eval=FALSE--------------
#  
#  library(GENEAclassify)
#  

## ----run the vignette,eval = FALSE--------------------------------------------
#  
#  vignette("GENEAclassifyDemo", package = NULL, lib.loc = NULL, all = TRUE)
#  

## ----Loading Data then Segmenting, eval = FALSE-------------------------------
#   # Name of the file to analyse
#  DataFile = "DataDirectory/jl_left wrist_010094_2012-01-30 20-39-54.bin"
#  ImportedData = dataImport(DataFile, downsample = 100, start = 0, end = 0.1)
#  head(ImportData)

## ---- eval = FALSE------------------------------------------------------------
#  # These are some of the output variables from segmentation and getGENEAsegments
#   dataCols <- c("UpDown.mean",
#                  "UpDown.var",
#                  "UpDown.sd",
#                  "Degrees.mean",
#                  "Degrees.var",
#                  "Degrees.sd",
#                  "Magnitude.mean",
#                  # Frequency Variables
#                  "Principal.Frequency.median",
#                  "Principal.Frequency.mad",
#                  "Principal.Frequency.GENEAratio",
#                  "Principal.Frequency.sumdiff",
#                  "Principal.Frequency.meandiff",
#                  "Principal.Frequency.abssumdiff",
#                  "Principal.Frequency.sddiff",
#                  # Light Variables
#                  "Light.mean",
#                  "Light.max",
#                  # Temperature Variables
#                  "Temp.mean",
#                  "Temp.sumdiff",
#                  "Temp.meandiff",
#                  "Temp.abssumdiff",
#                  "Temp.sddiff",
#                  # Step Variables
#                  "Step.GENEAcount",
#                  "Step.sd",
#                  "Step.mean")
#  
#  # Performing the segmentation now given the dataCols we want to find.
#  
#  SegDataFile = segmentation(ImportedData, dataCols)
#  # View the data from the segmentation
#  head(SegDataFile)

## ----segment a datafile, eval = FALSE-----------------------------------------
#   # Name of the file to analyse
#  DataFile = "DataDirectory/jl_left wrist_010094_2012-01-30 20-39-54.bin"
#  SegDataFile = getGENEAsegments(DataFile, dataCols, start = 0, end = 0.1)

## ----Displaying varying step counting alogrithms, eval = FALSE----------------
#  
#  WalkingData = "TrainingData/Walking/walking_jl_right wrist_024603_2015-12-12 15-36-47.bin"
#  
#  # Starting with default filter
#  W1 = getGENEAsegments(WalkingData, plot.it = TRUE)
#  
#  # plot.it Shows the crossing points. Turn this on for all plots to see how each filter works
#  # List the step outputs here.
#  W1$Step.GENEAcount; W1$Step.sd; W1$Step.mean
#  
#  W2 = getGENEAsegments(WalkingData, filteroder = 4)
#  # Changing the filterorder changes the order of the chebyshev filter applied.
#  W2$Step.GENEAcount; W2$Step.sd; W2$Step.mean
#  
#  W3 = getGENEAsegments(WalkingData, boundaries = c(0.15, 1))
#  # List the step outputs here.
#  W3$Step.GENEAcount; W3$Step.sd; W3$Step.mean
#  
#  # Changing the deicbel paramter
#  W4 = getGENEAsegments(WalkingData, Rp = 3)
#  W4$Step.GENEAcount; W4$Step.sd; W4$Step.mean
#  
#  # Increasing the hystersis
#  W5 = getGENEAsegments(WalkingData, hysteresis = 0.1)
#  W5$Step.GENEAcount; W5$Step.sd; W5$Step.mean

## ----loading TrainingData.csv, eval = FALSE-----------------------------------
#  # Change the file path to the location of GENEAclassify.
#  setwd("/Users/owner/Documents/GENEActiv/GENEAclassify_1.41/Data")
#  TrainingData = read.table("TrainingData.csv", sep = ",")
#  
#  # The data can also be called through from the package.
#  data(TrainingData)
#  TrainingData

## ---- eval = FALSE------------------------------------------------------------
#  ClassificationModel = createGENEAmodel(TrainingData,
#                                         features = c("Segment.Duration",
#                                                      "UpDown.mean", "UpDown.sd",
#                                                      "Degrees.mean", "Degrees.sd",
#                                                      "Magnitude.mean",
#                                                      "Light.mean",
#                                                      "Temp.mean",
#                                                      "Step.sd", "Step.count", "Step.mean",
#                                                      "Principal.Frequency.median", "Principal.Frequency.mad")
#                     )

## ---- eval = FALSE------------------------------------------------------------
#  ClassificationModel = createGENEAmodel(TrainingData,
#                                         features = c("UpDown.mean", "UpDown.sd",
#                                                      "Degrees.mean", "Degrees.sd",
#                                                      "Magnitude.mean",
#                                                      "Step.sd", "Step.mean",
#                                                      "Principal.Frequency.median",
#                                                      "Principal.Frequency.mad"))

## ----classifying a File, eval = FALSE-----------------------------------------
#  DataFile = "jl_left wrist_010094_2012-01-30 20-39-54.bin" # Change to the file to classify
#  ClassifiedFile = classifyGENEA(DataFile,
#                                 trainingfit = ClassificationModel,
#                                 start = "3:00",
#                                 end = "1 3:00")

## ----classifying a Directory, eval = FALSE------------------------------------
#  ClassifiedDirectory = classifyGENEA(DataDirectory,
#                                      trainingfit = ClassificationModel,
#                                      start = "3:00",
#                                      end = "1 3:00")

## ----Segmentation RunWalk file, echo = FALSE, eval = FALSE--------------------
#  SegData = getGENEAsegments("RunWalk.bin", end = "9:23")
#  head(SegData)

## ----List creation,eval = FALSE-----------------------------------------------
#  Activity = c("Running",
#               "Running",
#               "Walking")

## ----Attaching Activities, eval = FALSE---------------------------------------
#  SegData = cbind(SegData, ActivitiesListed)

## ---- eval = FALSE------------------------------------------------------------
#  SegData$Activity[1:2] = "Running"
#  SegData$Activity[3] = "Walking"

## ---- eval = FALSE------------------------------------------------------------
#  Cycling = getGENEAsegments("TrainingData/Cycling")
#  Cycling$Activity = "Cycling"
#  
#  NonWear = getGENEAsegments("TrainingData/NonWear")
#  NonWear$Activity = "NonWear"
#  
#  onthego = getGENEAsegments("TrainingData/onthego")
#  onthego$Activity = "onthego"
#  
#  Running = getGENEAsegments("TrainingData/Running")
#  Running$Activity = "Running"
#  
#  Sitting = getGENEAsegments("TrainingData/Sitting")
#  Sitting$Activity = "Sitting"
#  
#  Sleep = getGENEAsegments("TrainingData/Sleep")
#  Sleep$Activity = "Sleep"
#  
#  Standing = getGENEAsegments("TrainingData/Standing")
#  Standing$Activity = "Standing"
#  
#  Swimming = getGENEAsegments("TrainingData/Swimming")
#  Swimming$Activity = "Swimming"
#  
#  Transport = getGENEAsegments("TrainingData/Transport")
#  Transport$Activity = "Transport"
#  
#  Walking = getGENEAsegments("TrainingData/Walking")
#  Walking$Activity = "Walking"
#  
#  Workingout = getGENEAsegments("TrainingData/Workingout")
#  Workingout$Activity = "Workingout"

## ---- Combining Segments, eval = FALSE----------------------------------------
#  TrainingData = rbind(Cycling,
#                       NonWear,
#                       onthego,
#                       Running,
#                       Sitting,
#                       Sleep,
#                       Standing,
#                       Swimming,
#                       Transport,
#                       Walking,
#                       Workingout)

## ---- eval = FALSE------------------------------------------------------------
#  ClassificationModel = createGENEAmodel(TrainingData,
#                     features = c("UpDown.mean",
#                                  "UpDown.sd","Degrees.mean",
#                                  "Degrees.sd","Magnitude.mean",
#                                  "Step.sd","Step.mean",
#                                  "Principal.Frequency.median",
#                                  "Principal.Frequency.mad"))

