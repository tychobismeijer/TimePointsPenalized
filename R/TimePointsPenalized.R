#' fits lasso with penalized differences between adjacent time pounts coefficients 
#'
#' @param y0 case/control vector (no time iformation - "naive" approach)
#' @param x0 gene expression matrix (rows-samples by columns-genes)
#' @param FollowUp follow-up times (recurrence time for recurrences and follow-up for patients with no recurrences)
#' @param lam1V array of lasso penalty prefactor
#' @param gamma prefactor of the second penalty term - differences between adjacent time points coefficients
#' @param tV array of time points
#' @param standardize TRUE/FALSE standardization of the x0 columns (zero mean, unit variance)
#' @param Clinilal0 dataframe with clinical information (same order as rows of x0)
#' @param cores number of cores for parallelization (using foreach)
#' @export
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


