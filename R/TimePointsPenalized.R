library(expm)
library(glmnet)
library(Rcpp)
library(pROC)
library(stringr)
library(tidyverse)
library(ggpubr)
library(ComplexHeatmap)
require(doParallel)
library(RcppArmadillo)
library(Rcpp)

system("export OPENBLAS_NUM_THREADS=1")
system("export GOTO_NUM_THREADS=1")
system("export OMP_NUM_THREADS=1")
compileAttributes()
# Rcpp::sourceCpp("/home/misha/Documents/Development/TimePointsPenalized/src/FittingFunctions.cpp")
sourceCpp("src/FittingFunctions.cpp")

fitTimePointsPenalized <- function(y0, x0, FollowUp, lam1V, gamma, tV, standardize=TRUE, Clinilal0=data.frame(case_control0=y0), cores=1)
{     
  registerDoParallel(cores = cores)
  if (standardize) for (i in 1:ncol(x0)) x0[,i] <- (x0[,i] - mean(x0[,i]))/sd(x0[,i])
  Intercept <- 0
  beta <- rep(0,ncol(x0)*length(tV))
  y <- c()
  Clinilal <- data.frame()
  samplesT <- 1:(nrow(x0)*length(tV))
  GenesT <- rep("",ncol(x0)*length(tV))
  Clinilal0$sample <- rownames(x0)
  Clinilal0$time <- FollowUp
  
  for (it in 1:length(tV))
  {
    t <- tV[it]
    case_controlT <- ifelse(FollowUp>t,0,ifelse(y0==1,1,-1))
    y <- c(y,case_controlT)
    ClinilalT <- Clinilal0
    ClinilalT$StatusT <- case_controlT
    ClinilalT$FollowUp <- FollowUp
    Clinilal <- rbind(Clinilal,ClinilalT)
    GenesT[(1+(it-1)*ncol(x0)):(it*ncol(x0))] <- paste0(colnames(x0),"_t_",it)
    samplesT[(1+(it-1)*nrow(x0)):(it*nrow(x0))] <- paste0(rownames(x0),"_t_",it)
  }
  names(beta) <- GenesT
  names(y) <- samplesT
  Ind <- which(y %in% c(0,1))
  Clinilal <- Clinilal[Ind,]
  y <- y[Ind]
  
  IndFor0 <- c()
  IndTFor0 <- c()
  w <- y*0
  for (it in 1:length(tV))
  {
    IndT <- which(Clinilal$time==tV[it])
    w[IndT][which(y[IndT]==0)] <- 1/sum(y[IndT]==0)
    w[IndT][which(y[IndT]==1)] <- 1/sum(y[IndT]==1)
    IndFor0 <- c(IndFor0,which(rownames(x0) %in% samples[IndT]))
    IndTFor0 <- c(IndTFor0,which(rownames(x0) %in% samples[IndT])*0+it)
  }
  fits <- list()
  for (ilam1 in 1:length(lam1V))
  {
    lam1 <- lam1V[ilam1]
    lam2 <- gamma*lam1
    fits <- list(fits,.Call(`_TimePointsPenalized_FitRound`, x0, y, tV, lam1, lam2, beta, Intercept, w, IndFor0, IndTFor0));
    
    beta <- fits[[ilam1]]$beta
    Intercept <- fits[[ilam1]]$Intercept
  }
}       


