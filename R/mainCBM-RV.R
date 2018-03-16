##' Main function of the CBM-RV package
##'
##' The CBM.RV function is the main function of the CBM-RV
##' used package.
##'
##' @title CBM.RV main function
##' @param y matrix of K phenotypes data (one row-entry per individual), of
##'     row dimension \eqn{n}.
##' @param G matrix of the SNPs within the genomic regions in
##'     which the association test specified by the \code{method}
##'     parameter should be run. The matrix has a dimension:
##'     \eqn{n \times r}, with \eqn{r} the number of SNPs/variants.
##' @param covariates matrix of covariates NOT-including intercept (dimension:
##'     \eqn{n \times p}, with \eqn{p} the number of covariates)
##' @param null.obj a list returned by the function CBM.RV.Null(). The list contains the inverse of the variances-covariances matrix under the null model.
##' The list contains also the spectrale decomposition of the of kinship matrix of the data under study. See function CBM.RV.Null() for more details about the ouptut of this function.
##' @param weights optional numeric vector of genotype weights. If
##'     this opting is not specified, the beta distribution is used
##'     for weighting the variants, with each weight given by
##'     \eqn{w_i = dbeta(f_i, 1, 25)^2}, with \eqn{f_i} the minor
##'     allele frequency (MAF) of variant \eqn{i}. This default is the
##'     same as used by the
##'     \href{https://cran.r-project.org/web/packages/SKAT/}{\code{SKAT}
##'     package}. This vector is used as the diagonal of the
##'     \eqn{r \times r} matrix \eqn{W}, with \eqn{r} the number of
##'     variants.
##' @param method character, selects the method to use for the
##'     association testing. Can be one of the following:
##' \itemize{
##' \item \code{"optimal.p"} (default)
##' \item \code{"minimum.p"}, minimum p-value
##' \item \code{"Qscore"}, p-value based on the Q score statistic 
##' \item \code{"QFisher"}, p-value based on the Fisher combination of p-values derived from Q1 and Q2 
##' }
##' @param copfit character, selects the copula to use for modelling the (Q1, Q2) dependence. 
##' Can be one of the following:
##' \itemize{
##' \item \code{"Gaussian"} (default)
##' \item \code{"Student"}, Student copula with degree of freedom equals 10
##' \item \code{"Clayton"}, Clayton copula 
##' \item \code{"Frank"}, Frank copula
##' }
##' @return A data frame containing results of the association test
##'     specified by the \code{method} parameter using the copula model for modelling (Q1,Q2) dependency 
##'     by the \code{copfit} parameter. The
##'     output data frame contains the following columns:
##'     \itemize{
##'     \item \code{Score.Test}: the score of the given association test
##'     \item \code{P.value}: the p-value of the association test
##'     \item \code{alpha}: the estimated dependence parameter
##'     }
##' @author Karim Oualkacha
##' @export
CBM.RV <- function(y=NULL,
                   G=NULL,
                   covariates=NULL,
                   null.obj=NULL,
                   weights=NULL,
                   copfit="Gaussian",
                   method="optimal.p"
)
{
  ## Parameter checks
  y <- check_pheno(y)
  check_covariates(covariates, y)
  check_copula(copfit)
  check_method(method)
  check_weights(weights)
  
  inv.sqrt.Omega0 = null.obj$inv.sqrt.Omega0
  U = null.obj$U
  S = null.obj$S
  
  #----- calculation of p-values
  message("Starting association analysis of the genomic region...")
  message("Step 1: eigenvalue/eigenvector kernel matrices calculations...")
  res12 <- K1K2(Y, G, inv.sqrt.Omega0, covariates = covariates)
  K1.tilde = res12$K1.tilde
  K2.tilde = res12$K2.tilde
  Y.tilde = res12$Y.tilde
  param.K12 = eig.val.K12(Y.tilde,K1.tilde,K2.tilde)
  a = param.K12$a; b = param.K12$b; c = param.K12$c; V.Q12 = param.K12$V.Q12; cor.12 = param.K12$cor.12;
  values1 = param.K12$values1; values2 = param.K12$values2;
  exp.product = param.K12$exp.product
  Q1 = param.K12$Q1; Q2 = param.K12$Q2; K1 = param.K12$K1; K2 = param.K12$K2; 
  #----------- Q.1 and Q2 ---------#
  p.Q1 = davies(Q1,values1)$Qq
  p.Q2 = davies(Q2,values2)$Qq
  #-----------
  if (copfit == "Gaussian") { family=1; borne.inf=-0.9; borne.sup=0.9}
  if (copfit == "Student") {family=2; borne.inf=-0.9; borne.sup=0.9} 
  if (copfit == "Clayton") {family=3; borne.inf=.00001; borne.sup=5}
  if (copfit == "Frank") {family=5; borne.inf=.00001; borne.sup=10}
  
  message("Step 2: Starting estimation of the depenence copula parameter...")
  paramCop = param.cop(values1,values2,cop=copfit,borne.inf,borne.sup,exp.product,par2=10) # compute the copula parameter using formulas from MASKAT.pdf
  
  message("Step 3: Starting p.value calculation...")
  switch(method,
         optimal.p={
           R = try(Compute.Q.al(Q1,Q2,values1,values2,K1,K2,family=family,par=paramCop,par2=10),TRUE)
           if(class(R)=="try-error"){
             p.value=NA
           } else {
             p.value = try(R$p.value,TRUE)
           }
           p.value
         },
         minimum.p={
           Q <- Compute.Q.max(Q1,Q2)
           p.value = Compute.pvalue.Q.max(Q,values1,values2,family=family,par=paramCop,par2=10)
           p.value
         },
         Qscore={
           Q <- compute.Q.opt(Q1,Q2,V.Q12)
           p.value = try(Compute.pvalue.Q.opt(Q,a,b,values1,values2,family=family,par=paramCop,par2=10),TRUE) # compute analytical p-value using formulas from MASKAT.pdf
           if(class(p.value)=="try-error"){
             p.value=NA
           }
           p.value
         },
         QFisher={
           Q <- Compute.Q.Fisher(Q1,Q2,values1,values2)
           p.value = Compute.pvalue.Q.Fisher.Copula(Q,values1,values2,family=family,par= paramCop,par2=10)
           p.value
         },
         QMFKM={
           Q <- Q1 + Q2
           K.sum = K1 + K2
           values.sum = eigs(K.sum,matrix.rank(K.sum,method="chol"))$values
           p.value = davies(Q,values.sum)$Qq
           p.value
         }
  )
  
  return(list(Score.Test=Q, p.value=p.value, alpha=paramCop))
}




