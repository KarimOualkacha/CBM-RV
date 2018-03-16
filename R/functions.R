# Prog_LLC_KO file: new version 28 april 2017
# Problem 1 reported by JS: davies.SURV.2 approximates some p-values near to zero by a small negative value 
# Problem 1 solved see davies.SURV.2() function  
##########################
## Procedures generales ##
##########################

## Calcul de la trace d'un produit de matrices

mult.i=function(i,A,B) {(A[i,]%*%B[,i])[1,1]}
mult = function(A,B,n) {sum(sapply(1:n,mult.i,A=A,B=B))}

## Calcul de la fonction de survie des formes quadratiques par l'approximation de davies en forme matricielle pour les integrales 

davies.SURV.2 = function(i,q,values) {
                                     qqq = davies(q[i],lambda=values)$Qq
  if (qqq > 1) qqq = 1
  if (qqq < 0) qqq =0 
    return(qqq)
                                     }
davies.SURV = function(q,values) {sapply(1:(length(q)),davies.SURV.2,q=q,values=values)}

## Calcul des quantiles des formes quadratiques

davies.SURV.3 = function(q,values,res) {davies(q,lambda=values)$Qq-res}
davies.QUANTILE.2 = function(i,q,values) {uniroot(davies.SURV.3,lower=0,upper=100,values=values,res=q[i])$root}
davies.QUANTILE = function(q,values) {sapply(1:(length(q)),davies.QUANTILE.2,q=q,values=values)}


#################################################################
## Calcul du parametre de la copule a partir de la correlation ##
#################################################################

integrant.expectation.product.quadratic.forms=function(x,y,values1,values2,family,par,h=1e-10)
{
S1=davies.SURV(x,values1)
S2=davies.SURV(y,values2)
S1.h=davies.SURV(x+h,values1)
S2.h=davies.SURV(y+h,values2)
f1=(S1-S1.h)/h
f2=(S2-S2.h)/h
cop=BiCopPDF(u1=S1,u2=S2,family=family,par=par)
return(x*y*cop*f1*f2)
}

compute.expectation.product.Q1.Q2 = function(par,values1,values2,family,res=0,Q1max,Q2max)
{
I=integral2(integrant.expectation.product.quadratic.forms,xmin=0,xmax=Q1max,ymin=0,ymax=Q2max,values1=values1,values2=values2,family=family,par=par)
return(I$Q-res)
}

compute.copula.parameter = function(values1,values2,family,borne.inf,borne.sup,exp.product)
{
Q1max=davies.QUANTILE(1e-20,values1) # earch for q1 satisfies P(Q1>q1) approximatly 0
Q2max=davies.QUANTILE(1e-20,values2)
U=uniroot(compute.expectation.product.Q1.Q2,lower=borne.inf,upper=borne.sup,values1=values1,values2=values2,family=family,res=exp.product,Q1max=Q1max,Q2max=Q2max)
return(U$root)
}

#------------------------------------------------------------------------------------#
#
#------------------------------------------------------------------------------------#
param.cop = function(values1,values2,cop,borne.inf,borne.sup,exp.product,par2=10){ 
  
  if(cop=="Gaussian"){
    res<-compute.copula.parameter(values1,values2,family=1,borne.inf=borne.inf,borne.sup=borne.sup,exp.product=exp.product)
  }
else if(cop=="Student") {
  res<-compute.copula.parameter2(values1,values2,family=2,borne.inf=borne.inf,borne.sup=borne.sup,exp.product=exp.product,par2=10)
}
else if(cop=="Clayton") {
  res<-compute.copula.parameter(values1,values2,family=3,borne.inf=borne.inf,borne.sup=borne.sup,exp.product=exp.product)
}
#else if(cop=="Gumbel") {
#  res<-compute.copula.parameter(values1,values2,family=4,borne.inf=??,borne.sup=??,exp.product=exp.product)
#}
else if(cop=="Frank") {
  res<-compute.copula.parameter(values1,values2,family=5,borne.inf=borne.inf,borne.sup=borne.sup,exp.product=exp.product)
}
#else if(cop=="Joe"){
#  res<-compute.copula.parameter(values1,values2,family=6,borne.inf=??,borne.sup=??,exp.product=exp.product)
#}
  return(res)
  }
#------------------------------------------------------------------------------------#
#####################################################################################
## Calcul du parametre de la copule a partir de la correlation (pour une copule t) ##
#####################################################################################

integrant.expectation.product.quadratic.forms2=function(x,y,values1,values2,family,par,par2,h=1e-10)
{
S1=davies.SURV(x,values1)
S2=davies.SURV(y,values2)
S1.h=davies.SURV(x+h,values1)
S2.h=davies.SURV(y+h,values2)
f1=(S1-S1.h)/h
f2=(S2-S2.h)/h
cop=BiCopPDF(u1=S1,u2=S2,family=family,par=par,par2=par2)
return(x*y*cop*f1*f2)
}

compute.expectation.product.Q1.Q2.2 = function(par,par2,values1,values2,family,res=0,Q1max,Q2max)
{
I=integral2(integrant.expectation.product.quadratic.forms2,xmin=0,xmax=Q1max,ymin=0,ymax=Q2max,values1=values1,values2=values2,family=family,par=par,par2=par2)
return(I$Q-res)
}

compute.copula.parameter2 = function(values1,values2,family,borne.inf,borne.sup,exp.product,par2)
{
Q1max=davies.QUANTILE(1e-20,values1)
Q2max=davies.QUANTILE(1e-20,values2)
U=uniroot(compute.expectation.product.Q1.Q2.2,lower=borne.inf,upper=borne.sup,values1=values1,values2=values2,family=family,res=exp.product,Q1max=Q1max,Q2max=Q2max,par2=par2)
return(U$root)
}

#####################################
## Calcul des statistiques de test ##
#####################################

### Q.opt

compute.Q.opt = function(Q1,Q2,V)
{
Q=c(Q1,Q2)
stat = (t(Q)%*%V%*%Q)[1,1]
return(stat)
}

to.integrate=function(q2,q,a,b,values1,values2,family,par,par2=10,h=1e-10)
{
Q1.sup = (sqrt(q2^2*(b^2-a^2)+a*q)-b*q2)/a
U1 = davies.SURV(Q1.sup,values1)
U2 = davies.SURV(q2,values2)
U2.h = davies.SURV(q2+h,values2)
f2 = (U2-U2.h)/h
if (family == 2){  
  Cop = BiCopHfunc2(U1,U2,family=family,par=par,par2=par2)
  } else if (family != 2){ Cop = BiCopHfunc2(U1,U2,family=family,par=par) }
res = (1-Cop)*f2
return(res)
}


Compute.pvalue.Q.opt=function(Q.opt,a,b,values1,values2,family,par, par2 = 10)
{
lim.sup=sqrt(a*Q.opt/(a^2-b^2))
if (family == 2){  
I=integrate(to.integrate,lower=0,upper=lim.sup,q=Q.opt,a=a,b=b,values1=values1,values2=values2,family=family,par=par,par2=par2)
} else if (family != 2){ 
  I=integrate(to.integrate,lower=0,upper=lim.sup,q=Q.opt,a=a,b=b,values1=values1,values2=values2,family=family,par=par)  
}  
return(1-I$value)
}

############################################

### Q.star

Compute.Q.star=function(Q1,Q2,values1,values2,V.star)
{
Q1.star = qnorm(davies.SURV(Q1,values1))
Q2.star = qnorm(davies.SURV(Q2,values2))
Q.star = c(Q1.star,Q2.star)
return((t(Q.star)%*%V.star%*%Q.star)[1,1])
}


############################################

### Q.max

Compute.Q.max = function(Q1,Q2) {pmax(Q1,Q2)}

Compute.pvalue.Q.max = function(Q.max,values1,values2,family,par,par2=10)
{
S1=davies.SURV(Q.max,values1)
S2=davies.SURV(Q.max,values2)
if (family == 2){  
  Cop = BiCopCDF(u1=S1,u2=S2,family=family,par=par,par2=par2)
} else if (family != 2){ Cop = BiCopCDF(u1=S1,u2=S2,family=family,par=par) }
pvalue = S1 + S2 - Cop
return(pvalue)
}


############################################

###### Fisher's combination method

Compute.Q.Fisher = function(Q1,Q2,values1,values2)
{
p1 = davies.SURV(Q1,values1)
p2 = davies.SURV(Q2,values2)
res = -2*log(p1) - 2*log(p2)
return(res)
}


to.integrate.F = function(v2,q,values1,values2,family,par,par2=10)
{
U1 = 1 - pchisq(q-v2,df=2)
U2 = 1 - pchisq(v2,df=2)
f2 = dchisq(v2,df=2)
if (family == 2){  
  Cop = BiCopHfunc2(U1,U2,family=family,par=par, par2 = par2)
} else if (family != 2){ Cop = BiCopHfunc2(U1,U2,family=family,par=par) }
res = (1-Cop)*f2
return(res)
}

Compute.pvalue.Q.Fisher.Copula = function(Q.Fisher,values1,values2,family,par,par2=10)
{
  if (family == 2){  
    I = integrate(to.integrate.F,lower=0,upper=Q.Fisher,q=Q.Fisher,values1=values1,values2=values2,family=family,par=par,par2=par2)
  } else if (family != 2){ I = integrate(to.integrate.F,lower=0,upper=Q.Fisher,q=Q.Fisher,values1=values1,values2=values2,family=family,par=par) }
    return((1-I$value))
}

Compute.pvalue.Q.Fisher.Approx = function(Q.Fisher,f,c) {1-pchisq(Q.Fisher/c,df=f)}

############################################

###### Q.al = al*Q1 + (1-al)*Q2
to.integrate3 = function(q2,P1,P2,P3,P4,P5,al1,al2,al3,al4,al5,values1,values2,family,par,par2=10,h=1e-10)
{
R1 = (P1-(1-al1)*q2)/al1
R2 = (P2-(1-al2)*q2)/al2
R3 = (P3-(1-al3)*q2)/al3
R4 = (P4-(1-al4)*q2)/al4
R5 = (P5-(1-al5)*q2)/al5
R=apply(cbind(R1,R2,R3,R4,R5),1,min)

U1=davies.SURV(R,values1)
U2=davies.SURV(q2,values2)
U2.h = davies.SURV(q2+h,values2)
f2 = (U2-U2.h)/h
if (family == 2){  
  Cop = BiCopHfunc2(U1,U2,family=family,par=par, par2 = par2)
} else if (family != 2){ Cop = BiCopHfunc2(U1,U2,family=family,par=par) }
res = (1-Cop)*f2
return(res)
}

#-------------------------------------------------------------------------------------------#
#
#-------------------------------------------------------------------------------------------#
Compute.Q.al = function(Q1,Q2,values1,values2,K1,K2,family,par,par2=10)
{
  
p0 = davies.SURV(Q1,values1)

al1 = 0.2
Q.al1 = al1*Q1 + (1-al1)*Q2
K.al1 = al1*K1 + (1-al1)*K2
values.al1 = eigs(K.al1,matrix.rank(K.al1,method="chol"))$values
p1 = davies.SURV(Q.al1,values.al1)

al2 = 0.4
Q.al2 = al2*Q1 + (1-al2)*Q2
K.al2 = al2*K1 + (1-al2)*K2
values.al2 = eigs(K.al2,matrix.rank(K.al2,method="chol"))$values
p2 = davies.SURV(Q.al2,values.al2)

al3 = 0.6
Q.al3 = al3*Q1 + (1-al3)*Q2
p3 = davies.SURV(Q.al3,values.al2)

al4 = 0.8
Q.al4 = al4*Q1 + (1-al4)*Q2
p4 = davies.SURV(Q.al4,values.al1)

al5=1
p5 = davies.SURV(Q2,values2)

p = min(c(p0,p1,p2,p3,p4,p5))

P0 = davies.QUANTILE(0,values1)
P1 = davies.QUANTILE(p,values.al1)
P2 = davies.QUANTILE(p,values.al2)
P3 = P2
P4 = P1
P5 = davies.QUANTILE(p,values2)

al=c(0.2,0.4,0.6,0.8)
P=c(P1,P2,P3,P4)

lim.sup.Q2 = min(c(P0,P/(1-al)))

if (family == 2){  
  I=integrate(to.integrate3,lower=0,upper=lim.sup.Q2,P1=P1,P2=P2,P3=P3,P4=P4,P5=P5,
              al1=al1,al2=al2,al3=al3,al4=al4,al5=al5,values1=values1,values2=values2,family,par,par2=par2)
} else if (family != 2) {   I=integrate(to.integrate3,lower=0,upper=lim.sup.Q2,P1=P1,P2=P2,P3=P3,P4=P4,P5=P5,
                                        al1=al1,al2=al2,al3=al3,al4=al4,al5=al5,values1=values1,values2=values2,family,par)
}
result = data.frame(stat=p,p.value=(1-I$value))
return(result)
}

#------------------------------------------------------------------#
#
#------------------------------------------------------------------#
perm.Q.opt <- function(Q.opt, K1, K2, V.Q12, nb.perm=1000){
  
  MCMC.Q.opt = NULL
  for (s.MCMC in 1:nb.perm){
    MCMC.Z = rnorm(dim(K1)[1],0,1)
    Q1.perm <- t(MCMC.Z) %*% K1 %*% MCMC.Z
    Q2.perm <- t(MCMC.Z) %*% K2 %*% MCMC.Z
    Q.opt.perm = compute.Q.opt(Q1.perm,Q2.perm,V.Q12)
    MCMC.Q.opt = c(MCMC.Q.opt, Q.opt.perm)
  }
  p.Q.opt <- mean( MCMC.Q.opt > as.numeric(Q.opt) )

  m1=mean(MCMC.Q.opt)
  m2=mean((MCMC.Q.opt-m1)^2)
  m4=mean((MCMC.Q.opt-m1)^4)
  g=(m4/(m2^2))-3
  df=12/g
  stat=df+(Q.opt-m1)*sqrt(2*df/m2)
  corrected.p.value = 1-pchisq(stat,df=df)
  results<-list(p.Q.opt=p.Q.opt, p.Q.opt.perm = corrected.p.value)
  results
}

#------------------------------------------------------------------#
#
#------------------------------------------------------------------#
eig.val.K12<- function(Y.tilde, K1.tilde,K2.tilde){
  
  n=dim(K1.tilde)[1]
  values1.tilde = eigs(K1.tilde,matrix.rank(K1.tilde,method="chol"))$values
  values2.tilde = eigs(K2.tilde,matrix.rank(K2.tilde,method="chol"))$values
  
  K1 = K1.tilde / sum(values1.tilde)
  K2 = K2.tilde / sum(values2.tilde)
  
  values1 = values1.tilde / sum(values1.tilde)
  values2 = values2.tilde / sum(values2.tilde)
  
  v.Q = 2*mult(K1,K1,n)
  cov.12=2*mult(K1,K2,n)
  cor.12=cov.12/sqrt(v.Q*v.Q)
  
  V=solve(matrix(c(v.Q,cov.12,cov.12,v.Q),ncol=2))
  
  a=V[1,1]
  b=V[1,2]
  c=V[2,2]
  
  Q1 = t(Y.tilde) %*% K1 %*% Y.tilde
  Q2 = t(Y.tilde) %*% K2 %*% Y.tilde
  
  exp.product = cov.12+sum(values1)*sum(values2)
  res<-list(Q1=Q1, Q2=Q2, K1=K1, K2=K2, a=a, b=b,c=c,values1=values1,values2=values2,V.Q12=V, exp.product = exp.product, cor.12 = cor.12)
  
  res
}

#-------------------------------------------------------------------------------------------#
#                                                                                           #
#-------------------------------------------------------------------------------------------#
Kernel.Geno <- function(G,freq.MAF){
  if( length(freq.MAF) == 1){
    w = (dbeta(freq.MAF, 1, 25))^2
    K = w * G %*% t(G)
  } else
  {
    w = vector(length = length(freq.MAF))
    for (i in 1:length(freq.MAF)){
      w[i] = (dbeta(freq.MAF[i], 1, 25))^2
    }
    w = diag(w)
    K = G %*% w %*% t(G)
  }
  return(K)
}

#-------------------------------------------------------------------------------------------#
#                                                                                           #
#-------------------------------------------------------------------------------------------#
eigen.Omega0 = function(Sigma.RG, Sigma.e, kin){
  
  Omega.0 = Sigma.RG %x% kin + Sigma.e %x% diag(1, dim(kin)[1], dim(kin)[1])
  eigen.Omega0 = eigen(Omega.0)
  V.0 = eigen.Omega0$vectors 
  D.0 = 1 / sqrt(eigen.Omega0$values)
  inv.sqrtD.0 = diag(D.0)
  inv.sqrt.Omega0 = V.0 %*% inv.sqrtD.0 %*% t(V.0)
  return(inv.sqrt.Omega0)
}

#-------------------------------------------------------------------------------------------#
#                                                                                           #
#-------------------------------------------------------------------------------------------#
K1K2 = function(Y.multi, G, inv.sqrt.Omega0, covariates = NULL){
  # Y.multi is a matrix N  x 2 with rows designate subjects and columns are pheno1 and pheno 2
  # G: N x r matrix of genotypes  
  # inv.sqrt.Omega0: the sqrt of the Variance-covariance matrix of (Y.multi[,1]^t, Y.multi[,2]^t)^t under the null  
  
  N = dim(Y.multi)[1]
  d = dim(Y.multi)[2]
  #---------------------------------------
  freq.MAF = apply(G, 2, mean)/2
  K <- Kernel.Geno(G,freq.MAF)
  #---------------------------------------
  Y.Cov = c(Y.multi) 
  un.n = c(rep(1,N))
  if (is.null(covariates)){ X = kronecker(diag(1,d,d),un.n) } 
  else if(!is.null(covariates)){
    x = cbind(un.n,covariates)
    X = kronecker(diag(1,d,d),x) 
  } 
  #---------------------------------------
  X.tilde = inv.sqrt.Omega0 %*% X
  Y.tilde = inv.sqrt.Omega0 %*% Y.Cov
  
  inv.tXX = solve( t(X.tilde) %*% X.tilde )
  M.0 = diag(1, (d*N), (d*N)) - ( X.tilde %*% inv.tXX %*% t(X.tilde) )
  
  M.sqrtOmega0 = M.0 %*% inv.sqrt.Omega0 ##t
  
  K.1 = matrix(c(1,0,0,0),nrow=2) %x% K ##t
  K1.tilde = M.sqrtOmega0 %*% K.1  
  K1.tilde = K1.tilde %*% t(M.sqrtOmega0)/2
  
  K.2 = matrix(c(0,0,0,1),nrow=2) %x% K 
  K2.tilde = M.sqrtOmega0 %*% K.2 
  K2.tilde = K2.tilde %*% t(M.sqrtOmega0)/2
  
  res <- list(Y.tilde = Y.tilde, K1.tilde = K1.tilde, K2.tilde = K2.tilde)
  res
}

## Mar 27, 2015 %% modified by KO on 09 Jan 2017:
# 1) I added p-value calculations from different methods as separate functions 
## subfunction code for MASKAT
## in simulation10

#---------------------------------------------#
#
#---------------------------------------------#
Neg.LogLikelihood.ASKAT <- function(delta, S, Ut.y, Ut.x, n) {
  W    <- diag(1/(delta+S))
  beta <- solve(t(Ut.x) %*% W %*% Ut.x) %*% t(Ut.x) %*% W %*% Ut.y
  s.g  <- mean((Ut.y-Ut.x %*% beta)^2 / (delta+S))
  LL   <- n * log(s.g) + sum(log(S+delta))
  return(LL)
}
#---------------------------------------------#
#
#---------------------------------------------#
Estim.H0.ASKAT <- function(y, X, S, U)
{
  Ut.x <- t(U) %*% X
  Ut.y <- t(U) %*% y
  delta <- optimize(Neg.LogLikelihood.ASKAT,
                    interval=c(1e-4, 100),
                    S=S,
                    Ut.y=Ut.y,
                    Ut.x=Ut.x,
                    n=length(y))$min
  
  W <- diag(1/(delta + S))
  beta <- solve(t(Ut.x) %*% W %*% Ut.x) %*% t(Ut.x) %*% W %*% Ut.y
  s.g <- mean((Ut.y - Ut.x %*% beta)^2 / (delta + S))
  s.e <- delta * s.g
  return(c(s.e, s.g, beta))
}
#---------------------------------------------#
#
#---------------------------------------------#
guess.values <- function(Y,S,U,covariates=NULL){
  
  N = dim(Y)[1]
  y1 = Y[,1]  
  y2 = Y[,2]
  un.n = c(rep(1,N))
  if (is.null(covariates)){
    X = as.matrix(un.n) 
  } else if(!is.null(covariates)){
    X = as.matrix(cbind(un.n,covariates))
  }
  
  intival.1 <- Estim.H0.ASKAT(y1, X, S, U)
  init.s.e.1 <- intival.1[1]
  init.s.g.1 <- intival.1[2]
  beta.cov.1 <- intival.1[3:length(intival.1)]
  
  intival.2 <- Estim.H0.ASKAT(y2, X, S, U)
  init.s.e.2 <- intival.2[1]
  init.s.g.2 <- intival.2[2]
  beta.cov.2 <- intival.2[3:length(intival.2)] 
  
  res1 <- y1 - X %*% beta.cov.1  
  res2 <- y2 - X %*% beta.cov.2  
  
  init.s.g.12 <- init.s.e.12 <- cor(res1,res2) 
  
 all.intival <- c(init.s.g.1,init.s.g.2,init.s.g.12,init.s.e.1,init.s.e.2,init.s.e.12) 
  #list(init.s.g.11 = init.s.g.1,
#                      init.s.g.22 = init.s.g.2,
#                      init.s.e.11 = init.s.e.1,
#                      init.s.e.22 = init.s.e.2,
#                      init.s.g.12 = init.s.g.12,
#                      init.s.e.12 = init.s.e.12)
  return(all.intival)
}
#---------------------------------------------#
#
#---------------------------------------------#
inv.V.22 <- function(d, S.g, S.e){
  a<-d*S.g+S.e
  a <- c(a)
  aa <- (1/(a[1]*a[4] - a[2]*a[3]) + 1e-4) * matrix(c(a[4],-a[2],-a[3],a[1]),2,2)
  return(aa)
}
#---------------------------------------------#
#
#---------------------------------------------#
det.V.22 <- function(d, S.g, S.e){
  a<-d*S.g+S.e
  det.a = determinant(a, logarithm=T)[1:2]
  log_det_a <- det.a$modulus[1] * det.a$sign
  return(log_det_a)
}
#---------------------------------------------#
#
#---------------------------------------------#
reml.lik.Nelder<-function(theta,y.tilde,X.tilde,U,S,I){
  #theta[1]: polygenic variance for pheno1
  #theta[2]: polygenic variance for pheno2
  #theta[3]: polygenic covariance between pheno1 and pehno2
  #theta[4]: residual variance for pheno1
  #theta[5]: residual variance for pheno2
  #theta[6]: residual covariance between pheno1 and pehno2
  log.s.g.11 <- theta[1]
  log.s.g.22 <- theta[2]
  s.g.12 <- theta[3]
  log.s.e.11 <- theta[4]
  log.s.e.22 <- theta[5]
  s.e.12 <- theta[6]
  ## added 18th april 
  k11 = exp(log.s.g.11) + abs(s.g.12)
  k22 = exp(log.s.g.22) + abs(s.g.12)
  K_est = matrix(c(k11,s.g.12,s.g.12,k22),nrow=2)
  e11 = exp(log.s.e.11) + abs(s.e.12)
  e22 = exp(log.s.e.22) + abs(s.e.12)
  e_est = matrix(c(e11,s.e.12,s.e.12,e22),nrow=2)
  V_est_inv <- bdiag(lapply(S, inv.V.22, S.g=K_est, S.e=e_est))
  log_det_V_est <- sum(sapply(S, det.V.22, S.g=K_est, S.e=e_est))
  XV_estX=t(X.tilde)%*%V_est_inv%*%X.tilde + diag(1e-4,dim(X.tilde)[2],dim(X.tilde)[2])
  determinant2 = determinant(XV_estX, logarithm=T)[1:2]
  log_det_XV_estX <- determinant2$modulus[1] * determinant2$sign
  beta=solve(XV_estX)%*%t(X.tilde)%*%V_est_inv%*%y.tilde
  mu_hat=X.tilde%*%beta
  logl<- log_det_V_est + log_det_XV_estX + t(y.tilde-mu_hat)%*%V_est_inv%*%(y.tilde-mu_hat)
  return(as.numeric(logl))
}

reml.lik<-function(theta,y.tilde,X.tilde,U,S,I){
  #theta[1]: polygenic variance for pheno1
  #theta[2]: polygenic variance for pheno2
  #theta[3]: polygenic covariance between pheno1 and pehno2
  #theta[4]: residual variance for pheno1
  #theta[5]: residual variance for pheno2
  #theta[6]: residual covariance between pheno1 and pehno2
  log.s.g.11 <- theta[1]
  log.s.g.22 <- theta[2]
  s.g.12 <- theta[3] * sqrt(exp(log.s.g.11) *exp(log.s.g.22))  
  log.s.e.11 <- theta[4]
  log.s.e.22 <- theta[5]
  s.e.12 <- theta[6] * sqrt(exp(log.s.e.11) * exp(log.s.e.22))
  ## added 18th april 
  #k11 = exp(log.s.g.11) + abs(s.g.12)
  #k22 = exp(log.s.g.22) + abs(s.g.12)
  K_est = matrix(c(exp(log.s.g.11),s.g.12,s.g.12,exp(log.s.g.22)),nrow=2)
  #e11 = exp(log.s.e.11) + abs(s.e.12)
  #e22 = exp(log.s.e.22) + abs(s.e.12)
  e_est = matrix(c(exp(log.s.e.11),s.e.12,s.e.12,exp(log.s.e.22)),nrow=2)
  V_est_inv <- bdiag(lapply(S, inv.V.22, S.g=K_est, S.e=e_est))
  log_det_V_est <- sum(sapply(S, det.V.22, S.g=K_est, S.e=e_est))
  XV_estX=t(X.tilde)%*%V_est_inv%*%X.tilde + diag(1e-4,dim(X.tilde)[2],dim(X.tilde)[2])
  determinant2 = determinant(XV_estX, logarithm=T)[1:2]
  log_det_XV_estX <- determinant2$modulus[1] * determinant2$sign
  beta=solve(XV_estX)%*%t(X.tilde)%*%V_est_inv%*%y.tilde
  mu_hat=X.tilde%*%beta
  logl<- log_det_V_est + log_det_XV_estX + t(y.tilde-mu_hat)%*%V_est_inv%*%(y.tilde-mu_hat)
  return(as.numeric(logl))
}

reml.lik.BFGS<-function(theta,y.tilde,X.tilde,U,S,I){
  #theta[1]: polygenic variance for pheno1
  #theta[2]: polygenic variance for pheno2
  #theta[3]: polygenic covariance between pheno1 and pehno2
  #theta[4]: residual variance for pheno1
  #theta[5]: residual variance for pheno2
  #theta[6]: residual covariance between pheno1 and pehno2
  s.g.11 <- theta[1]
  s.g.22 <- theta[2]
  s.g.12 <- theta[3] * sqrt(s.g.11 * s.g.22)  
  s.e.11 <- theta[4]
  s.e.22 <- theta[5]
  s.e.12 <- theta[6] * sqrt(s.e.11 * s.e.22)
  K_est = matrix(c(s.g.11,s.g.12,s.g.12,s.g.22),nrow=2)
  e_est = matrix(c(s.e.11,s.e.12,s.e.12,s.e.22),nrow=2)
  V_est_inv <- bdiag(lapply(S, inv.V.22, S.g=K_est, S.e=e_est))
  log_det_V_est <- sum(sapply(S, det.V.22, S.g=K_est, S.e=e_est))
  XV_estX=t(X.tilde)%*%V_est_inv%*%X.tilde + diag(1e-4,dim(X.tilde)[2],dim(X.tilde)[2])
  determinant2 = determinant(XV_estX, logarithm=T)[1:2]
  log_det_XV_estX <- determinant2$modulus[1] * determinant2$sign
  beta=solve(XV_estX)%*%t(X.tilde)%*%V_est_inv%*%y.tilde
  mu_hat=X.tilde%*%beta
  logl<- log_det_V_est + log_det_XV_estX + t(y.tilde-mu_hat)%*%V_est_inv%*%(y.tilde-mu_hat)
  return(as.numeric(logl))
}
#-------------------------------------------------------------------------------------------#
#
#-------------------------------------------------------------------------------------------#
REML.VCs.under.null <- function(Y, 
                                covariates = NULL, 
                                U,
                                S,
                                nb.running.guess = 1,
                                Initial.values = NULL){ 
  #----- some check ----------------------#
  if (!is.null(Initial.values)){
    if (length(Initial.values) != 6){
      warning("VCs have 6 degrees of freedom. If it is not NULL, Initial.values should be a vector of length 6 ")
      stop()
    }
  }
  #---------------------------------------#
  Y.Cov = c(t(Y)) 
  N = dim(Y)[1]
  d = dim(Y)[2]
  if (is.null(covariates)){
    un.n = c(rep(1,N))
    X = kronecker(un.n,diag(1,d,d)) 
  }
  else if(!is.null(covariates)){
    un.n = c(rep(1,N))
    x = cbind(un.n,covariates)
    X = kronecker(x,diag(1,d,d)) 
  } 
  
  if (is.null(Initial.values)){
    init <- guess.values(Y,S,U,covariates)
    init.para <- c(log(init$init.s.g.11),log(init$init.s.g.22),init$init.s.g.12,
                   log(init$init.s.e.11),log(init$init.s.e.22),init$init.s.e.12)  
  }
  else   if (!is.null(Initial.values)){
    init.para <- c(log(Initial.values[1]),log(Initial.values[2]),Initial.values[3],log(Initial.values[4]),log(Initial.values[5]),Initial.values[6])
  }
  
  VC.all = NULL
  LL = NULL
  
  I = diag(1,N,N)  
  y.tilde = kronecker(t(U),diag(1,2,2)) %*% Y.Cov
  X.tilde = kronecker(t(U),diag(1,2,2)) %*% X
  
  if (nb.running.guess >= 1){ 
    print("waiting for null model VCs estimations: first guess run of optim function for the null model is running")
    all = optim(init.para, reml.lik.1, y=y.tilde, X=X.tilde, U=U, S=S, I=I, method = "L-BFGS-B")
    VC.all = cbind(VC.all, c(exp(all$par[1]),exp(all$par[2]),all$par[3],exp(all$par[4]),exp(all$par[5]),all$par[6]))
    LL = c(LL,0.5*all$value)
    
    if (nb.running.guess == 1){
    message("Estimation of VCs under H0 is run only one time with one vector of initial values")
    VCs = VC.all
    
    sg11 = exp(VCs[1]) + abs(VCs[3])
    sg22 = exp(VCs[2]) + abs(VCs[3])
    sg12 = VCs[3]
    S.g = matrix(c(sg11,sg12,sg12,sg22),ncol = 2)
    
    se11 = exp(VCs[4]) + abs(VCs[4])
    se22 = exp(VCs[5]) + abs(VCs[4])
    se12 = VCs[4]
    S.e = matrix(c(se11,se12,se12,se22),ncol = 2)
    #S.g12 = VCs[3] * sqrt(VCs[1] * VCs[2])
    #S.g = matrix(c(VCs[1],S.g12,S.g12,VCs[2]),ncol = 2)
    #S.e12 = VCs[6] * sqrt(VCs[4] * VCs[5])
    #S.e = matrix(c(VCs[4],S.e12,S.e12,VCs[5]),ncol = 2)
    
    }
    } else if (nb.running.guess > 1){
  for (i in 2:nb.running.guess){
    print(paste("running optim with the ", i, "-th guess values of variance components", sep=""))
    print(runif(1,0.1,1)*init.para)
    all = optim(runif(1,0.1,1)*init.para, reml.lik, y=y.tilde, X=X.tilde, U=U, S=S, I=I, method = "L-BFGS-B")
    VC.all = cbind(VC.all,c(exp(all$par[1]),exp(all$par[2]),all$par[3],exp(all$par[4]),exp(all$par[5]),all$par[6]))
    LL = c(LL,all$value) 
    }
    ind = which(LL>0)
    LL = LL[ind]
    VC.all = VC.all[,ind]
    VCs = VC.all[,sort.int(LL, index.return=TRUE)$ix[1]]
    S.g12 = VCs[3] * sqrt(VCs[1] * VCs[2])
    S.g = matrix(c(VCs[1],S.g12,S.g12,VCs[2]),ncol = 2)
    S.e12 = VCs[6] * sqrt(VCs[4] * VCs[5])
    S.e = matrix(c(VCs[4],S.e12,S.e12,VCs[5]),ncol = 2)
    }
  
  results <- list(VCs = VCs, Polygenic.VC = S.g, Env.VC = S.e, VC.all=VC.all, LL=LL)
}
#-------------------------------------------------------------------------------------------#
#
#-------------------------------------------------------------------------------------------#
REML.VCs.under.null.Nelder <- function(Y, 
                                covariates = NULL, 
                                U,
                                S,
                                nb.running.guess = 1,
                                Initial.values = NULL){ 
  #----- some check ----------------------#
  if (!is.null(Initial.values)){
    if (length(Initial.values) != 6){
      warning("VCs have 6 degrees of freedom. If it is not NULL, Initial.values should be a vector of length 6 ")
      stop()
    }
  }
  #---------------------------------------#
  Y.Cov = c(t(Y)) 
  N = dim(Y)[1]
  d = dim(Y)[2]
  if (is.null(covariates)){
    un.n = c(rep(1,N))
    X = kronecker(un.n,diag(1,d,d)) 
  }
  else if(!is.null(covariates)){
    un.n = c(rep(1,N))
    x = cbind(un.n,covariates)
    X = kronecker(x,diag(1,d,d)) 
  } 
  
  if (is.null(Initial.values)){
    init <- guess.values(Y,S,U,covariates)
    init.para <- c(log(1-abs(init[3])),log(1-abs(init[3])),init[3],
                   log(1-abs(init[3])),log(1-abs(init[6])),init[6])  
  }
  else   if (!is.null(Initial.values)){ # needs to be cheked
    init.para <- c(log(Initial.values[1]),log(Initial.values[2]),Initial.values[3],log(Initial.values[4]),log(Initial.values[5]),Initial.values[6])
  }
  
  VC.all = NULL
  LL = NULL
  
  I = diag(1,N,N)  
  y.tilde = kronecker(t(U),diag(1,2,2)) %*% Y.Cov
  X.tilde = kronecker(t(U),diag(1,2,2)) %*% X
  
  if (nb.running.guess >= 1){ 
    print("waiting for null model VCs estimations: first guess run of optim function for the null model is running")
    all = optim(init.para, reml.lik.Nelder, y=y.tilde, X=X.tilde, U=U, S=S, I=I, method = "Nelder-Mead")
    VC.all = cbind(VC.all, all$par)
    LL = c(LL,0.5*all$value)
    
    if (nb.running.guess == 1){
      message("Estimation of VCs under H0 is run only one time with one vector of initial values")
      VCs = VC.all
    }
  } else if (nb.running.guess > 1){
    for (i in 2:nb.running.guess){
      print(paste("running optim with the ", i, "-th guess values of variance components", sep=""))
      print(runif(1,0.1,1)*init.para)
      all = optim(runif(1,0.1,1)*init.para, reml.lik.Nelder, y=y.tilde, X=X.tilde, U=U, S=S, I=I, method = "Nelder-Mead")
      VC.all = cbind(VC.all, all$par)
      LL = c(LL,all$value) 
    }
    ind = which(LL>0)
    LL = LL[ind]
    VC.all = VC.all[,ind]
    VCs = VC.all[,sort.int(LL, index.return=TRUE)$ix[1]]
    
    }
  sg11 = exp(VCs[1]) + abs(VCs[3])
  sg22 = exp(VCs[2]) + abs(VCs[3])
  sg12 = VCs[3]
  S.g = matrix(c(sg11,sg12,sg12,sg22),ncol = 2)
  
  se11 = exp(VCs[4]) + abs(VCs[6])
  se22 = exp(VCs[5]) + abs(VCs[6])
  se12 = VCs[6]
  S.e = matrix(c(se11,se12,se12,se22),ncol = 2)
  
  results <- list(VCs = VCs, Polygenic.VC = S.g, Env.VC = S.e, VC.all=VC.all, LL=LL)
}
#-------------------------------------------------------------------------------------------#
#
#-------------------------------------------------------------------------------------------#
REML.VCs.under.null.BFGS <- function(Y, 
                                covariates = NULL, 
                                U,
                                S,
                                nb.running.guess = 1,
                                Initial.values = NULL){ 
  #----- some check ----------------------#
  if (!is.null(Initial.values)){
    if (length(Initial.values) != 6){
      warning("VCs have 6 degrees of freedom. If it is not NULL, Initial.values should be a vector of length 6 ")
      stop()
    }
  }
  #---------------------------------------#
  Y.Cov = c(t(Y)) 
  N = dim(Y)[1]
  d = dim(Y)[2]
  if (is.null(covariates)){
    un.n = c(rep(1,N))
    X = kronecker(un.n,diag(1,d,d)) 
  }
  else if(!is.null(covariates)){
    un.n = c(rep(1,N))
    x = cbind(un.n,covariates)
    X = kronecker(x,diag(1,d,d)) 
  } 
  
  if (is.null(Initial.values)){
    init <- guess.values(Y,S,U,covariates)
    init.para <- init
    }
  else   if (!is.null(Initial.values)){
    init.para <- Initial.values
    }
  
  VC.all = NULL
  LL = NULL
  
  I = diag(1,N,N)  
  y.tilde = kronecker(t(U),diag(1,2,2)) %*% Y.Cov
  X.tilde = kronecker(t(U),diag(1,2,2)) %*% X

  if (nb.running.guess >= 1){ 
    print("waiting for null model VCs estimations: first guess run of optim function for the null model is running")
    all = optim(init.para, reml.lik.BFGS, y=y.tilde, X=X.tilde, U=U, S=S, I=I, method = "L-BFGS-B", lower=c(0.001,0.001,-0.98,0.001,0.001,-0.98), upper=c(Inf,Inf,0.98,Inf,Inf,0.98))
    VC.all = cbind(VC.all, all$par)
    LL = c(LL,0.5*all$value)
    
  if (nb.running.guess == 1){
      message("Estimation of VCs under H0 is run only one time with one vector of initial values")
      VCs = VC.all
    }
  if (nb.running.guess > 1){
    for (i in 2:nb.running.guess){
      print(paste("running optim with the ", i, "-th guess values of variance components", sep=""))
      init.par2 = runif(1,0.1,1)*init.para
      all = optim(init.par2, reml.lik.BFGS, y=y.tilde, X=X.tilde, U=U, S=S, I=I, method = "L-BFGS-B", lower=c(0.001,0.001,-0.98,0.001,0.001,-0.98), upper=c(Inf,Inf,0.98,Inf,Inf,0.98))
      VC.all = cbind(VC.all,all$par)
      LL = c(LL,all$value) 
    }
    ind = which(LL>0)
    LL = LL[ind]
    VC.all = VC.all[,ind]
    VCs = VC.all[,sort.int(LL, index.return=TRUE)$ix[1]]
  }
  }
  S.g12 = VCs[3] * sqrt(VCs[1] * VCs[2])
  S.g = matrix(c(VCs[1],S.g12,S.g12,VCs[2]),ncol = 2)
  S.e12 = VCs[6] * sqrt(VCs[4] * VCs[5])
  S.e = matrix(c(VCs[4],S.e12,S.e12,VCs[5]),ncol = 2)
  
  results <- list(VCs = VCs, Polygenic.VC = S.g, Env.VC = S.e, VC.all=VC.all, LL=LL)
}
#######################################################################################################################
#######################################################################################################################
####  copPar2KT: given dependence parameter of some copula, this function evaluate the Kendall's tau in term of alpha 
####  using functions of Vinecopula package
########################################################################################################################
copPar2KT<-function(alpha,cop,par2=0){
  
  if(cop=="Gaussian")     
    res<-BiCopPar2Tau(1,alpha)
  else if(cop=="Student")  
    res<-BiCopPar2Tau(2,alpha,par2=par2)
  else if(cop=="Clayton")  
    res<-copClayton@tau(alpha)
  else if(cop=="Gumbel") 
    res<-copGumbel@tau(alpha+1)
  else if(cop=="Frank")  
    res<-copFrank@tau(alpha)
  else if(cop=="Joe") 
    res<-copJoe@tau(alpha+1)
  res
}

#copPar2KT(alpha=.5,cop="Student",par2=3)


#######################################################################################################################
####  given kendall's tau of some copula, this function return the correponding dependence parameter 
####  using functions of Vinecopula package
########################################################################################################################
ftau<-function(tau,alpha,cop,par2){tau-copPar2KT(alpha=alpha,cop=cop,par2=par2)}

copKT2Par<-function(tau,cop,par2=0){
  
  if(cop=="Gaussian")     
    res<-BiCopTau2Par(1,tau)
  else if(cop=="Student")  
    res<-uniroot(ftau,lower=-.9,upper=.9,extendInt="yes",tau=tau,cop="Student",par2=par2)$root
  else if(cop=="Clayton"){
    #tau<-max(1e-7,tau)
    res<-copClayton@iTau(tau) }
  else if(cop=="Gumbel"){ 
    #tau<-max(1,tau)
    res<-copGumbel@iTau(2-1/tau) }
  else if(cop=="Frank") {
    #tau<-min(.99999,max(1e-7,tau)) 
    res<-copFrank@iTau(tau) }
  else if(cop=="Joe"){ 
    tau<-max(1e-10,tau)
    res<-uniroot(ftau,lower=1e-9,upper=700,extendInt="yes",tau=tau,cop="Joe",par2=par2)$root
  }
  res
}
