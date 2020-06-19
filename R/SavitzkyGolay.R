#' Savitzky-Golay filtering and derivatives
#'
#' @param X \code{matrix} containing spectra as rows.
#' @param poly Polynomial degree of smoother.
#' @param width Window width of smoother, default = 11, must be an odd number.
#' @param deriv Derivative degree, can be 0, default = 2. 
#' @param ends Handling of spectrum ends, i.e. first and last (width-1)/2 points. Default is "cut", i.e. remove ends, "extrapolate" copies the first/last estimable point, while "zeros" is included for backward compatilibity (fill with 0).
#'
#' @return A matrix of filtered spectra (possibly with derivatives)
#' @export
#'
#' @examples
#' data(fishoil)
#' Raman    <- fishoil$Raman[, 850:3300]
#' SavGol   <- SavitzkyGolay(Raman)
#' old.par  <- par(mfrow = c(2,1), mar = c(4,4,1,1))
#' matplot(colnames(Raman), t(Raman), type = 'l',
#'         ylab = 'Relative intensity', xlab = 'Raw spectra')
#' matplot(colnames(SavGol), t(SavGol), type = 'l',
#'         ylab = 'Relative intensity', xlab = 'Smoothed 2nd derivative')
#' par(old.par)
SavitzkyGolay <- function(X, poly=3, width=11, deriv=2, ends=c("cut","extrapolate","zeros")){
  # Strip X
  X <- unclass(as.matrix(X))
  
  if(width%%2 == 0)
    stop("'width' must be an odd number!")
  
  S <- outer((-(width-1)/2):((width-1)/2), 0:poly,"^") # Flipped Vandermonde matrix
  R <- qr.R(qr(S))
  g <- S %*% tcrossprod(solve(R))
  
  M  <- X*0 # matrix(0, nrow(X), ncol(X))
  F1 <- (width+1)/2
  F2 <- (-F1+1):(F1-1)
  if(deriv == 0){
    for(j in F1:(ncol(X)-(width-1)/2)){ # Calculate the n-th derivative of the i-th spectrum
      M[,j] <- X[,j + F2] %*% g[,1]
    }
  } else {
    for(j in F1:(ncol(X)-(width-1)/2)){ # Calculate the n-th derivative of the i-th spectrum
      M[,j] <- X[,j + F2] %*% (deriv * g[,deriv+1])
    }
  }
  if(ends[1] == "cut")
    M <- M[,((width-1)/2+1) : (ncol(M)-(width-1)/2+1)]
  if(ends[1] == "extrapolate")
    M[,1:(width-1)/2] <- M[,(width-1)/2+1]; M[,(ncol(M)-(width-1)/2+1):ncol(M)] <- M[,ncol(M)-(width-1)/2]
  return(M)
}
