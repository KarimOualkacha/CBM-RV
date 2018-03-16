#----------------------------------------------#
library(CompQuadForm)
library(MASS)
library(Matrix)
library(mvtnorm)
library(mnormt)
library(numDeriv)
library(kinship2)
library(VineCopula)
library(copula)
#----------------------------------------------#
library(pracma)
library(matrixcalc)
library(rARPACK)
#----------------------------------------------#
#----------------------------------------------#
#setwd("/Users/KOualkachaUQAM/Dropbox/Etudiants/Jianping/SSCcode/Codes4JS/CBMRV/R")
source("functions.R")
source("inputChecksCBMRV.R")
source("mainCBM-RV.R")
source("CBM-RV-Null.R")
#setwd("/Users/KOualkachaUQAM/Dropbox/Etudiants/Jianping/SSCcode/Codes4JS/CBMRV/data")
load("kin1.Rdata")
kin <- 2*kin1
Geno = as.matrix(read.table("genotypes.dat",sep = " "))
Y = as.matrix(read.table("phenotypes.dat",sep = " "))
covariates = as.matrix(read.table("covariates.dat",sep = " "))

#----------------------------------------------#

null.obj = CBM.RV.Null(Y, covariates = covariates, kin, nb.running.guess = 3)
CBM.RV(y=Y, G=Geno, covariates=covariates, null.obj=null.obj, copfit="Gaussian", method="optimal.p")

