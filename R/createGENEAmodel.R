#' @name trainingFit
#' @title Example classification tree
#' @description rpart object declaring the decision tree for the classification of GENEActiv bin data into typical groups.
#' @docType data
#' @format An \code{rpart} object with class \code{GENEA} containing the decision 
#' tree information required for prediction.
#' @source Output of \code{createGENEAmodel} on experimental training data.
#' @keywords datasets
#' @seealso \code{"\link[=levels.GENEA]{levels}"} \code{\link{features}}
#' @examples
#' data(trainingFit)
#' class(trainingFit)
#' levels(trainingFit)
#' features(trainingFit)
#' plot(trainingFit)
#' text(trainingFit, cex = 0.5)

globalVariables("trainingFit")

#' @name TrainingData
#' @title Example Training Data set
#' @description A manually classified training data set provided with the package. 
#' @docType data
#' @format Manually classified segmented data that can be used to build a classification model.
#' @source Manually classified by Actvinsights.
#' @keywords datasets
#' @seealso \code{"\link[=levels.GENEA]{levels}"} \code{\link{features}}
#' @examples
#' data(TrainingData)
#' class(TrainingData)

globalVariables("TrainingData")

#' @title Create training data decision tree model
#' 
#' @description From data frame create a decision tree that can be used for 
#' classifying data into specified categories. The data frame may optionally 
#' contain a reserved column Source, specifying the provenance of the record.
#' The data frame must contain a column, by default named Activity, specifing 
#' the classes into which the model fit should be classified.
#' 
#' @param data data frame containing segmented GENEActiv bin data.
#' @param outputtree name of the png file that shows the classification tree plot.
#' @param plot a logical value indicating whether a plot of the classification tree
#' should be plotted. The default is TRUE.
#' @param features character vector naming independent variables to use in classification.
#'     Alternatively, a numeric vector specifying the variables to pass to the classification 
#'     or NULL, in which case all variables are used in the order of the supplied training dataset.
#'     Note that including large numbers of variables (>7) may result in long run times.
#' @param category single character naming the dependent variable to use (default 'Activity').
#' @param verbose single logical should additional progress reporting be printed 
#' at the console? (default \code{TRUE})
#' @param \dots other arguments for \code{rpart}
#' @details The function will create an rpart classification tree for the training data based 
#' upon the parameters passed to features. The model created, an GENEA rpart object can be used 
#' within the function \code{"\link[=classifyGENEA]{classifyGENEA}"} to classify GENEA bin files.
#' 
#' @return A GENEA rpart fit.
#' @seealso The returned object can be interrogated with 
#'     \code{features}, the variables used in defining the model,
#'     and \code{"\link[=levels.GENEA]{levels}"}, the response categories predicted by the model.
#' @export
#' @import rpart
#' @importFrom graphics abline par text
#' @examples 
#' ## dataPath <- file.path(system.file(package = "GENEAclassify"),
#' ##                                   "testdata",
#' ##                                   "trainingData9.csv")
#' ##
#' ## t1 <- read.csv(file = dataPath)
#' ## 
#' ## f1 <- createGENEAmodel(data = t1,
#' ##                        features = c("Degrees.var",
#' ##                                     "UpDown.mad",
#' ##                                     "Magnitude.mean"),
#' ##                        category = "Activity")
#' ## 
#' ## class(f1)
#' ## levels(f1)
#' ## features(f1)
#' ## plot(f1)
#' ## text(f1)

createGENEAmodel  <- function(data, 
                              outputtree = NULL,
                              features = c("Segment.Duration", 
                                           "Principal.Frequency.mad", 
                                           "UpDown.sd", 
                                           "Degrees.sd"),
                              category = "Activity", 
                              plot = TRUE, 
                              verbose = TRUE, 
                              ...) {
    
    xVars <- colnames(data)
    
    if (!(length(category) == 1 && (is.character(category) || is.numeric(category)))) { 
        stop("one column in data must be specified to category which contains the training classes") }
    
    if (is.character(category)) { if (!category %in% xVars) { stop("category not found in data") } }
    
    if (is.character(numeric)) { if (!category %in% seq_along(xVars)) { stop("category not found in data") } }
    
    if (!length(unique(data[, category])) > 1) { 
        stop("more than one class in column ", category, " required for meaningful classification") }
    
    if (is.null(features)) {
        
        features <- seq_along(xVars)
        
        features <- features[!xVars %in% c(category, "Source")]
        
        xVars <- xVars[features]
        
    } else {
        
        if (length(features) == 0) { 
            stop("features is length 0, at least one variable is needed by classifier!") }
        
        if (is.factor(features)) { features <- levels(features)[features] }
        
        if (is.numeric(features)) {
            
            if (any(features < 1)) { 
                stop("negative indexing of features is not permitted") }
            
            if (any(features > ncol(data))) { 
                stop("features must be selected from columns provided") }
            
            xVars <- xVars[features]
            
            if (any(xVars %in% c(category, "Source"))) { 
                stop("features must not include 'Source' or '", category, "'") }
            
        } else {
            
            if (is.character(features)) {
                
                notFound <- !features %in% xVars
                
                if (any(notFound)) { 
                    stop("features not found in data: ", paste(features[notFound], collapse = ", ")) }
                
                xVars <- features
                
            } else { 
                stop("features object class not recognized, please provide a character or numeric vector") }
        }
    }
    
    if (verbose) { 
        cat(paste0("columns provided for training fit:\n", 
                paste(paste("    *", xVars), collapse = "\n"), "\n")) }
    
    formula <- as.formula(paste(category, "~", paste(xVars, collapse = " + ")))
    
	fit <- rpart(formula, method = "class", data = data, ...)
    
    selected <- features(fit)
    
    if (verbose) { 
        cat(paste0("columns used by training fit:\n", 
                paste(paste("    *", selected), collapse = "\n"), "\n")) }  
    
    if (plot) {
        
    	if (nrow(fit$cptable) > 1) {
            
            px <- par("xpd")
            
            par(xpd = TRUE)
            
            plot(fit, uniform = TRUE, main = "Classification Tree")
            
            text(fit, cex = 0.6)
            
            par(xpd = px)
            
            if (is.null(outputtree) == FALSE){
              # Check name of outputtree 
              
              if (grepl(outputtree, ".png") == FALSE){
                stop("outputtree must have the suffix .png to save the classification tree plot")
              }
              
              png(paste0(outputtree))
              
              px <- par("xpd")
              
              par(xpd = TRUE)
              
              plot(fit, uniform = TRUE, main = "Classification Tree")
              
              text(fit, cex = 0.6)
              
              par(xpd = px)
              
              dev.off()
            }
        
        } else { 
            warning("no classification to plot")
        }
    }
    
    fit$call$formula <- formula
    
    class(fit) <- c(class(fit), "GENEA")
    
    attr(x = fit, which = "features") <- selected
    
	return(fit)
}

