
#' @title Get features of a GENEA rpart fit
#' 
#' @description Shows the response features specified by a GENEA (rpart) object.
#' 
#' @param fit GENEA rpart object
#' @return Character vector naming features selected by GENEA rpart fit.
#' @export
#' @keywords internal
#' @examples
#' data(trainingFit)
#' features(fit = trainingFit)

features <- function(fit) {
    
    featr <- attr(fit, "features")
    
    if (is.null(featr)) {
        
        if (is(fit, "rpart")) {
            
            featr <- rownames(fit$splits)[inGroup(fit$splits[, "count"])]
            
            featr <- featr[!duplicated(featr)]
            
        } else { 
            warning("fit is not an rpart object; features could not be identified")
        }
    }
    
    return(featr)
}

#' @title Which values are in a new group?
#' 
#' @description Choose elements at each step.
#' Used by \code{\link{features}} to identify unique components of splits table.
#' 
#' @param x numeric vector
#' @return logical vector with length x
#' @keywords internal
#' @examples 
#' GENEAclassify:::inGroup(x = rep(-1, 5))
#' x1 <- c(1, 1, 1, 2, 2, 2, 2, 3, 3)
#' x1t <- GENEAclassify:::inGroup(x = x1)
#' x1[x1t]

inGroup <- function(x) {
    
    x <- diff(x = x) != 0
    
    x <- c(TRUE, x)
    
    return(x)
}


#' @title Get levels of a GENEA rpart response
#' 
#' @description Shows the levels used in a GENEA (rpart) fit object.
#' 
#' @param x GENEA rpart object
#' @return character vector naming levels classified by GENEA rpart fit
#' @method levels GENEA
#' @export
#' @keywords internal
#' @examples
#' data(trainingFit)
#' levels(x = trainingFit)

levels.GENEA <- function(x) {
    
    return(attr(x = x, which = "ylevels"))
}


#' @title Get features of a GENEActiv bin data object.
#' 
#' @description Shows the first n elements of a GENEActiv bin (list) data object.
#' 
#' @param x GENEActiv bin list object
#' @param \dots additional arguments
#' @return list showing the first part of each list element.
#' @method head GENEAbin
#' @export
#' @keywords internal
#' @import utils
#' @examples
#' temp <- list(Data = matrix(rnorm(40), ncol = 4), 
#'     UpDown = rnorm(10), 
#'     Degrees = sample(-90:90, 10),
#'     Time = 1:10, 
#'     Light = runif(10),
#'     Temp = rep(20, 10),
#'     Magnitude = runif(10), 
#'     Serial = "012345",
#'     RawData = matrix(rnorm(700), ncol = 7),
#'     Freq = 100)
#' class(temp) <- c("list", "GENEAbin")
#' head(x = temp)

head.GENEAbin <- function(x, ...) {

    lapply(X = x, FUN = head, ...)
}
