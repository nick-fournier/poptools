#### Installing packages ####
if("mipfp" %in% rownames(installed.packages()) == FALSE) {install.packages("mipfp")}
if("dplyr" %in% rownames(installed.packages()) == FALSE) {install.packages("dplyr")}
if("plyr" %in% rownames(installed.packages()) == FALSE) {install.packages("plyr")}
if("reshape2" %in% rownames(installed.packages()) == FALSE) {install.packages("reshape2")}
if("devtools" %in% rownames(installed.packages()) == FALSE) {install.packages("devtools")}
if("RCurl" %in% rownames(installed.packages()) == FALSE) {install.packages("RCurl")}
library(devtools)
if("aws.s3" %in% rownames(installed.packages()) == FALSE) {install.packages("aws.s3")}
if("glmnet" %in% rownames(installed.packages()) == FALSE) {install.packages("glmnet")}
if("doParallel" %in% rownames(installed.packages()) == FALSE) {install.packages("doParallel")}
if("parallel" %in% rownames(installed.packages()) == FALSE) {install.packages("parallel")}
if("pryr" %in% rownames(installed.packages()) == FALSE) {install.packages("pryr")}
if("filehash" %in% rownames(installed.packages()) == FALSE) {install.packages("filehash")}
if("data.table" %in% rownames(installed.packages()) == FALSE) {install.packages("data.table")}
if("extrafont" %in% rownames(installed.packages()) == FALSE) {install.packages("extrafont")}
if("ggplot2" %in% rownames(installed.packages()) == FALSE) {install.packages("ggplot2")}
if("abind" %in% rownames(installed.packages()) == FALSE) {install.packages("abind")}
if("mlogit" %in% rownames(installed.packages()) == FALSE) {install.packages("mlogit")}
if("mnlogit" %in% rownames(installed.packages()) == FALSE) {install.packages("mnlogit")}
if("pbapply" %in% rownames(installed.packages()) == FALSE) {install.packages("pbapply")}
if("Rcpp" %in% rownames(installed.packages()) == FALSE) {install.packages("Rcpp")}
if("RcppArmadillo" %in% rownames(installed.packages()) == FALSE) {install.packages("RcppArmadillo")}
if("RcppProgress" %in% rownames(installed.packages()) == FALSE) {install.packages("RcppProgress")}
if("Matrix" %in% rownames(installed.packages()) == FALSE) {install.packages("Matrix")}
if("nnls" %in% rownames(installed.packages()) == FALSE) {install.packages("nnls")}
if("limSolve" %in% rownames(installed.packages()) == FALSE) {install.packages("limSolve")}
#if("clpAPI" %in% rownames(installed.packages()) == FALSE) {install.packages("clpAPI")}
#install.packages("clpAPI", configure.args = c("--with-clp-lib=/home/nick/coin-Clp/lib", "--with-clp-include=/home/nick/coin-Clp/include/coin"))
#install.packages("clpAPI", configure.args = "--prefix=/home/nick/coin-Clp")
if("gridExtra" %in% rownames(installed.packages()) == FALSE) {install.packages("gridExtra")}
if("lme4" %in% rownames(installed.packages()) == FALSE) {install.packages("lme4")}
if("ROCR" %in% rownames(installed.packages()) == FALSE) {install.packages("ROCR")}
if("RcppGSL" %in% rownames(installed.packages()) == FALSE) {install.packages("RcppGSL")}
if("survey" %in% rownames(installed.packages()) == FALSE) {install.packages("survey")}
if("bnlearn" %in% rownames(installed.packages()) == FALSE) {install.packages("bnlearn")}
if("visNetwork" %in% rownames(installed.packages()) == FALSE) {install.packages("visNetwork")}
if("HydeNet" %in% rownames(installed.packages()) == FALSE) {install.packages("HydeNet")}

library(glmnet)
library(HydeNet)
library(visNetwork)
library(bnlearn)
library(survey)
library(RcppGSL)
library(RcppProgress)
library(mnlogit)
library(mipfp)
library(dplyr)
library(plyr)
library(reshape2)
library(aws.s3)
library(foreach)
library(doParallel)
library(parallel)
library(pryr)
library(filehash)
library(data.table)
library(RCurl)
library(extrafont)
library(ggplot2)
library(mlogit)
library(pbapply)
library(abind)
library(Rcpp)
library(RcppArmadillo)
library(Matrix)
library(limSolve)
library(nnls)
#library(clpAPI)
library(gridExtra)
library(lme4)
library(ROCR)
#library(CppHHassignment) #custom package
