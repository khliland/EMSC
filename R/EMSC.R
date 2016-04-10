#' Extended multiplicative signal correction (EMSC)
#'
#' Performs model-based background correction and normalisation of spectra. EMSC handles
#' variations in scaling, polynomial baselines and interferents. Parameters for corrections
#' are stored for further analysis, and spectra are corrected accordingly.
#'
#' This is the main EMSC function performing all calculations. It can be run with
#' no parameters (defaults are used), with a predefined EMSC model object or with
#' parameters that are passed on to the EMSC model building function \code{\link{EMSC_model}}.
#'
#' @param X \code{matrix} containing spectra as rows.
#' @param model an EMSC model to use instead of the other parameters.
#' @param ... named model parameters for EMSC_model.
#'
#' @return An object of class EMSC is returned. This contains:
#' \itemize{
#'  \item{\code{corrected}:}{ \code{matrix} of corrected spectra.}
#'  \item{\code{parameters}:}{ \code{matrix} of fitted parameter values.}
#'  \item{\code{model}:}{ object containing input all input parameters.}
#' }
#'
#' @seealso \code{\link{EMSC_model}} \code{\link{predict.EMSC}}
#' @references H. Martens, E. Stark, Extended multiplicative signal correction and spectral
#'  interference subtraction: new preprocessing methods for near infrared spectroscopy.
#'  J Pharm Biomed Anal. 1991; 9(8):625-35.
#'
#' @examples
#' data(milk)
#' Raman      <- milk$Raman[, 850:3300]
#' EMSC.basic <- EMSC(Raman)
#' EMSC.poly6 <- EMSC(Raman, degree = 6)
#' EMSC.ref   <- EMSC(Raman, degree = 6, reference = Raman[30, ])
#'
#' \dontrun{
#' old.par  <- par(mfrow = c(2,2), mar = c(4,4,1,1))
#' xlim     <- rev(as.numeric(range(colnames(Raman))))
#' matplot(colnames(Raman), t(Raman), type = 'l', xlim = xlim,
#'         ylab = 'Relative intensity', xlab = 'Raw spectra')
#' matplot(colnames(Raman), t(EMSC.basic$corrected), type = 'l', xlim = xlim,
#'         ylab = 'Relative intensity', xlab = 'Corrected (basic)')
#' matplot(colnames(Raman), t(EMSC.poly6$corrected), type = 'l', xlim = xlim,
#'         ylab = 'Relative intensity', xlab = 'Corrected (6th degree polynomial)')
#' matplot(colnames(Raman), t(EMSC.ref$corrected),   type = 'l', xlim = xlim,
#'         ylab = 'Relative intensity', xlab = 'Corrected (reference = spec. #30)')
#' par(old.par)
#' }
#'
#' @importFrom pracma mldivide
#' @export
EMSC <- function(X, model = NULL, ...){
  n <- dim(X)[1]

  # Make EMSC model if not supplied
  if(is.null(model)){
    mf <- match.call(expand.dots = TRUE)
    if(is.null(x <- mf$x)){
      x <- X
    }
    model <- EMSC_model(x = X, ...)
  }

  # Extract from model object
  terms <- model$sizes[2] + 2 + model$sizes[3] # Polynomial + interferents

  # Apply weights (if included)
  if(!is.null(model$weights)){
    Xw     <- X * rep(model$weights, each = n)
    modelw <- model$model * rep(model$weights, each = n)
  } else {
    Xw     <- X
    modelw <- model$model
  }

  # Handling of 0 objects
  PP <- apply(Xw==0, 1, all)

  # Perform parameter estimation
  if(any(PP)){
    Parameters0  <- mldivide(t(modelw), t(Xw[!PP,]))

    k <- 2:terms
    Corrected0 <- Xw[!PP,] - crossprod(Parameters0[k,, drop = FALSE], modelw[k,, drop = FALSE])

    corrected        <- Xw
    corrected[!PP,]  <- Corrected0
    parameters       <- matrix(0, nrow(Parameters0), nrow(Xw))
    parameters[,!PP] <- Parameters0

  } else {
    Parameters0  <- mldivide(t(modelw), t(Xw))

    k <- 2:terms
    Corrected0 <- Xw - crossprod(Parameters0[k,, drop = FALSE], modelw[k,, drop = FALSE])
    Corrected0 <- Corrected0 * 1./Parameters0[1,] # correct multipl. eff.

    corrected  <- Corrected0
    parameters <- Parameters0
  }

  # Return
  object <- list(corrected = corrected, parameters = parameters, model = model)
  class(object) <- c("EMSC", "list")
  object
}


#' Model object for extended multiplicative signal correction (EMSC)
#'
#' Sets up an EMSC model to be applied to one or more set of spectra.
#'
#' @param x \code{numeric} vector containing abcissas of spectra to be corrected, reference spectrum
#' with/without names or matrix to be corrected.
#' @param reference \code{numeric} vector containing the reference spectrum.
#' @param degree \code{integer} giving the polynomial degree of the baseline; 0 or higher, default = 2.
#' @param interferent \code{numeric} vector containing a spectral component to remove.
#' @param constituent \code{numeric} vector containing a spectral component to include.
#' @param weights \code{numeric} vector of abcissas weights.
#' @param rep_corr proportion of variance or number of subspace components in replicate space (not implemented yet).
#'
#' @return An EMSC model is returned containing all parameters.
#'
#' @seealso \code{\link{EMSC}} \code{\link{predict.EMSC}}
#' @export
EMSC_model <- function(x, reference = NA, degree = 2,
                       interferent = NULL, constituent = NULL, weights = NULL, rep_corr = NULL){

  # Handle abcissas
  if(is.data.frame(x))
    x <- as.matrix(x)
  if(is.matrix(x)){
    if(!is.ordered(abcissas <- colnames(x))){
      abcissas <- 1:dim(x)[2]
    }
    if(length(reference) == 1 && is.na(reference))
      reference <- colMeans(x)
  } else {
    if(!is.ordered(abcissas <- names(x))){
      if(!is.ordered(abcissas <- x)){
        abcissas <- 1:length(x)
      }
    }
  }

  # Sizes
  n.i <- ifelse(is.matrix(interferent), nrow(interferent), ifelse(is.null(interferent),0, 1))
  n.c <- ifelse(is.matrix(constituent), nrow(constituent), ifelse(is.null(constituent),0, 1))
  p <- length(abcissas)

  Start   <- abcissas[1]
  End     <- abcissas[p]
  C  <- 0.5*(Start+End)
  M  <- 2.0/(Start-End)

  # Construct polynomials
  if(is.null(degree)){
    degree <- -1
  }
  model <- matrix(0, degree+2, p)
  model[1,] <- reference
  mod_names <- "Reference"
  if(!is.null(degree)){ # Contains at least baseline
    model[2,] <- 1
    mod_names <- c(mod_names, "Baseline")
    if(degree > 0){ # Contains polynomial
      for(i in 1:degree){
        model[i+2,] <- (M*(abcissas-C))^i
        mod_names <- c(mod_names, paste("Degree", i))
      }
      model[3,] <- M*(Start-abcissas)-1
    }
  }

  # Add interferents and constituents
  model <- rbind(model, interferent, constituent)
  if(!is.null(interferent)){
    mod_names <- c(mod_names, paste("Interferent", 1:n.i))
  }
  if(!is.null(constituent)){
    mod_names <- c(mod_names, paste("Constituent", 1:n.c))
  }

  # Handle replicate correction
  # (not implemented yet)

  # Add dimnames
  dimnames(model) <- list(term = mod_names, abcissas = abcissas)

  # Return
  list(model = model, abcissas = abcissas, sizes = c(p, degree, n.i, n.c), weights = weights)
}


#' Predict Method for EMSC
#' 
#' Prediction for \code{EMSC} ojects. Corrections are calculated for the new
#' \code{matrix} based on the EMSC model used in the input object.
#' 
#' @param object An object fitted by the \code{EMSC} function.
#' @param newdata A \code{matrix} or object convertable to a matrix containing observations as rows.
#' 
#' @seealso \code{\link{EMSC}} \code{\link{EMSC_model}}
#' 
#' @examples
#' data(milk)
#' Raman.cal <- milk$Raman[  1:90,  850:3300]
#' Raman.val <- milk$Raman[-(1:90), 850:3300]
#' EMSC.cal  <- EMSC(Raman.cal)
#' EMSC.val  <- predict(EMSC.cal, Raman.val)
#' identical(EMSC.cal$model, EMSC.val$model) # Same model, reference spectrum, etc.
#' 
#' matplot(t(EMSC.cal$corrected), type = 'l', col = 'black', lty = 1, ylab = 'Intensity')
#' matplot(t(EMSC.val$corrected), type = 'l', col = 'red', lty = 2, add = TRUE)
#' legend('topleft', legend = c('Calibration','Validation'), lty = 1:2, col = 1:2)
#' 
#' @export
predict.EMSC <- function(object, newdata = NULL){
  if(is.null(newdata))
    return(object)
  
  newdata <- unclass(as.matrix(newdata))
  EMSC(newdata, object$model)
}