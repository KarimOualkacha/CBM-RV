##' Second Main function of the CBM-RV package
##'
##' The CBM.RV.null function is the second main function of the CBM-RV
##' used package.
##'
##' @title CBM.RV.null main function
##' @param y matrix of K phenotypes data (one row-entry per individual), of
##'     row dimension \eqn{n}, with \eqn{n} the number of subjects)
##' @param covariates matrix of covariates NOT-including intercept (dimension:
##'     \eqn{n \times p}, with \eqn{p} the number of covariates)
##' @param kin represents the kinship matrix (dimension:
##'     \eqn{n \times n})
##' @return A list containing varaince components estimates under the null model. The
##'     output list contains the following arguments:
##'     \itemize{
##'     \item \code{inv.sqrt.Omega0}: the inverse of the phenotypes varainces-covariances matrix under the null model
##'     \item \code{PolygenicVC}: the estimate of the polygenic varainces-covariances matrix
##'     \item \code{Env.VC}: the estimate of the Environmental varainces-covariances matrix
##'     \item \code{U}: the orthogonal matrix obtained form the spectral decomposition of the knship matrix: kinship = U S U^t
##'     \item \code{S}: the diagonal matrix obtained form the spectral decomposition of the knship matrix: kinship = U S U^t
##'     }
##' @author Karim Oualkacha
##' @export
CBM.RV.Null <- function(Y, 
                        covariates = covariates, 
                        kin = kin, 
                        nb.running.guess = 3){
  #---- need to estimate parameters/VCs under the null model----# 
  ss = eigen(kin)
  U = ss$vectors
  S = ss$values
  vc.res = REML.VCs.under.null.BFGS(Y, covariates = covariates, U=U, S=S, nb.running.guess = 3)
  inv.sqrt.Omega0<-eigen.Omega0(vc.res$Polygenic.VC, vc.res$Env.VC, kin)
  res.null <- list(inv.sqrt.Omega0 = inv.sqrt.Omega0, PolygenicVC = vc.res$Polygenic.VC, Env.VC = vc.res$Env.VC, U=U, S=S)
  return(res.null)
}
