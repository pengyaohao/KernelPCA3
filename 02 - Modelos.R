setwd("C:/Users/b05652877465/Desktop/KERNEL PCA/ROBUST")

rm(list = ls());gc()

eff.frontier <- function (covariance,vec_mean, short="no", max.allocation=NULL,
                          risk.premium.up=.5, risk.increment=.005){
  # return argument should be a m x n matrix with one column per security
  # short argument is whether short-selling is allowed; default is no (short
  # selling prohibited)max.allocation is the maximum % allowed for any one
  # security (reduces concentration) risk.premium.up is the upper limit of the
  # risk premium modeled (see for loop below) and risk.increment is the
  # increment (by) value used in the for loop
  
  
  n <- ncol(covariance)
  
  # Create initial Amat and bvec assuming only equality constraint
  # (short-selling is allowed, no allocation constraints)
  Amat <- matrix (1, nrow=n)
  bvec <- 1
  meq <- 1
  
  # Then modify the Amat and bvec if short-selling is prohibited
  if(short=="no"){
    Amat <- cbind(1, diag(n))
    bvec <- c(bvec, rep(0, n))
  }
  
  # And modify Amat and bvec if a max allocation (concentration) is specified
  if(!is.null(max.allocation)){
    if(max.allocation > 1 | max.allocation <0){
      stop("max.allocation must be greater than 0 and less than 1")
    }
    if(max.allocation * n < 1){
      stop("Need to set max.allocation higher; not enough assets to add to 1")
    }
    Amat <- cbind(Amat, -diag(n))
    bvec <- c(bvec, rep(-max.allocation, n))
  }
  
  # Calculate the number of loops
  loops <- risk.premium.up / risk.increment + 1
  loop <- 1
  
  # Initialize a matrix to contain allocation and statistics
  # This is not necessary, but speeds up processing and uses less memory
  eff <- matrix(nrow=loops, ncol=n+3)
  # Now I need to give the matrix column names
  colnames(eff) <- c(names(vec_mean), "Std.Dev", "Exp.Return", "sharpe")
  i <- 0 
  # Loop through the quadratic program solver
  for (i in seq(from=0, to=risk.premium.up, by=risk.increment)){
    dvec <- vec_mean * i # This moves the solution along the EF
    sol <- solve.QP(covariance, dvec=dvec, Amat=Amat, bvec=bvec, meq=meq)
    eff[loop,"Std.Dev"] <- sqrt(sum(sol$solution*colSums((covariance*sol$solution))))
    eff[loop,"Exp.Return"] <- as.numeric(sol$solution %*% vec_mean)
    eff[loop,"sharpe"] <- eff[loop,"Exp.Return"] / eff[loop,"Std.Dev"]
    eff[loop,1:n] <- sol$solution
    loop <- loop+1
  }
  
  return(as.data.frame(eff))
}

library(data.table)
library(dplyr)
library(bit64)
library(stringr)
library(readxl)
library(reshape2)
library(knitr)
library(ggplot2)
library(plyr)
library(gplots)
library(lubridate)
library(tcltk)
library(googlesheets)
library(gsheet)
library(corrplot)
library(xtable)
library(tidyverse)

library(magrittr)
library(RnavGraphImageData)
library(dplyr)
library(e1071)
library(kernlab)
library(quadprog)
library(robust)
library(PerformanceAnalytics)


vec_parms <- NULL
arquivos <- list.files(path = "data/global",full.names = T)
nomes_arquivos <- gsub(".rds","",list.files(path = "data/global"),fixed = T)
arquivos <- list.files(path = "C:/Users/b05652877465/Desktop/KERNEL PCA/ROBUST/dados",full.names = T)
nomes_arquivos <- gsub(".rds","",list.files(path = "C:/Users/b05652877465/Desktop/KERNEL PCA/ROBUST/dados"),fixed = T)
i <- 4

n_in_sample <- 100/15 # 1/10 in sample
# n_in_sample <- 10 # 1/10 in sample

# for (i in 1:7) {
for (i in 1:length(nomes_arquivos)) {
# for (i in 1:1) {
  base <- readRDS(arquivos[i])
  data_in_sample <- rev(base$DT)[trunc(nrow(base)/n_in_sample)] 
  var_modelo <- setdiff(names(base),"DT")
  dados_in_sample <- base %>%
    filter(DT<=data_in_sample)
  
  dados_out_of_sample <- base %>%
    filter(DT>data_in_sample)
  
  RF_mean <- mean(dados_out_of_sample$RF,na.rm=T)
  
  ### PEARSON ----
    
  mpearson  <- cor(dados_in_sample %>% select(var_modelo))
  covariance <- var(dados_in_sample %>% select(var_modelo))
  vec_mean <- colMeans(dados_in_sample %>% select(var_modelo))
  

  
  # Run the eff.frontier function based on no short and 50% alloc. restrictions
  eff <- eff.frontier(covariance,vec_mean, short="no", max.allocation=NULL,
                      risk.premium.up=1, risk.increment=.1)
  
  # Find the optimal portfolio
  eff.optimal.point <- eff %>% filter(sharpe==max(eff$sharpe)) 
  
  # # Conferindo valores
  wm <- eff.optimal.point %>% select(var_modelo)
  # sum(round(wm,10))
  

  A <-  as.matrix(data.matrix(dados_out_of_sample %>% select(var_modelo)),
                  ncol = ncol(dados_out_of_sample)) %*%  t(as.matrix(round(wm,10)))
  
  
  dados_out_of_sample <- dados_out_of_sample %>%
    mutate(rs = cumsum(r))
  
  
  
  ## Não vamos comparar shape ratio de diferentes frequências
  q <- 1
  
  
  ### Risk Free
  RF_mean <- mean(dados_out_of_sample$RF,na.rm=T)
  
  
  # Sp500
  auto_cov_sp <- acf(dados_out_of_sample$r,lag=1)$acf[2]
  n_q_sp <- q^(1/2)*(1 + (2*auto_cov_sp)/(1- auto_cov_sp)*(1 - (1-auto_cov_sp^q)/(q*(1-auto_cov_sp))))^(-1/2)
  
  sr_sp <- (mean(dados_out_of_sample$rs)-RF_mean)/sd(dados_out_of_sample$rs)
  sr_sp_month <-  n_q_sp*sr_sp 
  
  # Carteira
  auto_cov_out <- acf(A,lag=1)$acf[2]
  n_q <- q^(1/2)*(1 + (2*auto_cov_out)/(1- auto_cov_out)*(1 - (1-auto_cov_out^q)/(q*(1-auto_cov_out))))^(-1/2)
  sr <- (mean(A)-RF_mean)/sd(A)
  
  sr_month <- n_q*sr 

  sr_month_0 <- sr_month 
  VIID_q <- (n_q^2)*(1+1/2*(sr^2))
  VGMM_add <- ifelse(q==1,0,q*(sr)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q <- VIID_q   + VGMM_add 
  
  sd_sr_q <- (VGMM_q/nrow(dados_out_of_sample))^(1/2)
  
  
  library(PerformanceAnalytics)
  
  sor <- n_q*(mean(A)-RF_mean)/DownsideDeviation(A,MAR = RF_mean)
  sor_sp <- (mean(dados_out_of_sample$rs)-RF_mean)/DownsideDeviation(dados_out_of_sample$rs,MAR = RF_mean)
  
    
  VIID_q2 <- (n_q^2)*(1+1/2*(sor^2))
  VGMM_add2 <- ifelse(q==1,0,q*(sor)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q2 <- VIID_q2   + VGMM_add2 
  
  sd_sr_q2 <- (VGMM_q2/nrow(dados_out_of_sample))^(1/2)
  
  vec_parms_temp <- data.frame(data = nomes_arquivos[i],
                               cov_method = "Pearson",
                               initial_in = min(dados_in_sample$DT),
                               end_in = max(dados_in_sample$DT),
                               initial_out = min(dados_out_of_sample$DT),
                               end_out = max(dados_out_of_sample$DT),
                               N_in = nrow(dados_in_sample),
                               N_out = nrow(dados_out_of_sample),
                               X = length(var_modelo),
                               X_pos_n = length(which(wm>0.00)),
                               X_pos_valid_n = length(which(wm>0.01)),
                               X_pos_valid_sum = paste0(round(sum(wm[(which(wm>0.01))])*100,2),"%"),
                               w_max = paste0(round(max(wm[(which(wm>0.01))])*100,2),"%"),
                               w_min = paste0(round(min(wm[(which(wm>0.01))])*100,2),"%"),
                               r_final = paste0(round(cumsum(A)[nrow(dados_out_of_sample)]*100,4),"%"),
                               Q = NA, 
                               lamda_max = NA,
                               eigen_up = NA,
                               eigen_v_up = NA,
                               eigen_max = NA,
                               eigen_v_max = NA,
                               st_dev=sd(A),
                               down_dev=DownsideDeviation(A,MAR = RF_mean),
                               sharpe_ratio = sr_month ,
                               p_sharp_ratio_sp = 1-pnorm(sr_month-sr_sp_month,0,sd_sr_q),
                               p_sharp_ratio_0 = 1-pnorm(sr_month,0,sd_sr_q),
                               rmt_improvement = NA ,
                               sortino_ratio = sor ,
                               p_sort_ratio_sp = 1-pnorm(sor-sor_sp,0,sd_sr_q2),
                               # rmt_improvement_so = NA ,
                               p_sort_ratio_0 = 1-pnorm(sor,0,sd_sr_q2),
                               auto_cov_ar1 = auto_cov_out,
                               p_sharp_ratio_rmt = NA,
                               p_sort_ratio_rmt = NA)
  

  
  vec_parms <- rbind(vec_parms,vec_parms_temp)
  
  ### RMT ----
  ## ref: J. Bouchaud, M. Potters, Theory of Financial Risks—From Statistical Physics to Risk Management, Cambridge University Press, UK, 2000.
  Q <- nrow(dados_in_sample)/(length(var_modelo))
  ## falta entender sigma!!!
  eigen0 <- eigen(mpearson)
  # eigen0 <- eigen(covariance)
  lamda_max <- (1 + 1/Q + (1/Q)^.5)
  
  
  g1 <- eigen0$values[eigen0$values>=lamda_max]
  g2 <- eigen0$values[eigen0$values<lamda_max]
  eigen1<- diag(c(g1,rep(mean(eigen0$values),length(g2))))
  cor1 <- eigen0$vectors %*% eigen1 %*% t(eigen0$vectors)
  covariance_rmt <-  cor1
  for(ii in 1:ncol(cor1)){
    for(jj in 1:nrow(cor1)){
      covariance_rmt[ii,jj] <- cor1[ii,jj] * (covariance[ii,ii]*covariance[jj,jj])^.5 
    }
  }
  
  # Run the eff.frontier function based on no short and 50% alloc. restrictions
  eff <- eff.frontier(covariance_rmt,vec_mean, short="no", max.allocation=NULL,
                      risk.premium.up=1, risk.increment=.1)
  
  # Find the optimal portfolio
  eff.optimal.point <- eff %>% filter(sharpe==max(eff$sharpe)) 
  
  # Conferindo valores
  wm <- eff.optimal.point %>% select(var_modelo)
  
  B <-  as.matrix(data.matrix(dados_out_of_sample %>% select(var_modelo)),
                  ncol = ncol(dados_out_of_sample)) %*%
    t(as.matrix(round(wm,10)))
  
  
  
  # Carteira
  auto_cov_out <- acf(B,lag=1)$acf[2]
  n_q <- q^(1/2)*(1 + (2*auto_cov_out)/(1- auto_cov_out)*(1 - (1-auto_cov_out^q)/(q*(1-auto_cov_out))))^(-1/2)
  sr <- (mean(B)-RF_mean)/sd(B)
  
  sr_month <- n_q*sr 
  VIID_q <- (n_q^2)*(1+1/2*(sr^2))
  VGMM_add <- ifelse(q==1,0,q*(sr)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q <- VIID_q   + VGMM_add 
  
  sd_sr_q <- (VGMM_q/nrow(dados_out_of_sample))^(1/2)
  
  
  sor <- n_q*(mean(B)-RF_mean)/DownsideDeviation(B,MAR = RF_mean)
  sor_sp <- (mean(dados_out_of_sample$rs)-RF_mean)/DownsideDeviation(dados_out_of_sample$rs,MAR = RF_mean)
  
  
  VIID_q2 <- (n_q^2)*(1+1/2*(sor^2))
  VGMM_add2 <- ifelse(q==1,0,q*(sor)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q2 <- VIID_q2   + VGMM_add2 
  
  sd_sr_q2 <- (VGMM_q2/nrow(dados_out_of_sample))^(1/2)

  
  vec_parms_temp <- data.frame(data = nomes_arquivos[i],
                               cov_method = "Pearson RMT",
                               initial_in = min(dados_in_sample$DT),
                               end_in = max(dados_in_sample$DT),
                               initial_out = min(dados_out_of_sample$DT),
                               end_out = max(dados_out_of_sample$DT),
                               N_in = nrow(dados_in_sample),
                               N_out = nrow(dados_out_of_sample),
                               X = length(var_modelo),
                               X_pos_n = length(which(wm>0.00)),
                               X_pos_valid_n = length(which(wm>0.01)),
                               X_pos_valid_sum = paste0(round(sum(wm[(which(wm>0.01))])*100,2),"%"),
                               w_max = paste0(round(max(wm[(which(wm>0.01))])*100,2),"%"),
                               w_min = paste0(round(min(wm[(which(wm>0.01))])*100,2),"%"),
                               r_final = paste0(round(cumsum(B)[nrow(dados_out_of_sample)]*100,4),"%"),
                               Q = Q, 
                               lamda_max = lamda_max,
                               eigen_up = paste0(length(g1)," (",round(length(g1)/(ncol(dados_in_sample) - 1)*100,2),"%)"),
                               eigen_v_up = paste0(round(sum(g1)/(sum(g1)+sum(g2))*100,2),"%"),
                               eigen_max = g1[1]/lamda_max,
                               eigen_v_max = paste0(round(sum(g1[1])/(sum(g1)+sum(g2))*100,2),"%"),
                               st_dev=sd(B),
                               down_dev=DownsideDeviation(B,MAR = RF_mean),
                               sharpe_ratio = sr_month ,
                               p_sharp_ratio_sp = 1-pnorm(sr_month-sr_sp_month,0,sd_sr_q),
                               p_sharp_ratio_0 = 1-pnorm(sr_month,0,sd_sr_q),
                               rmt_improvement =((mean(B-A)+1)^(252)-1)*100,
                               sortino_ratio = sor ,
                               p_sort_ratio_sp = 1-pnorm(sor-sor_sp,0,sd_sr_q2),
                               p_sort_ratio_0 = 1-pnorm(sor,0,sd_sr_q2),
                               auto_cov_ar1 = auto_cov_out,
                               p_sharp_ratio_rmt =  1-pnorm(sr_month-sr_month_0,0,sd_sr_q),
                               p_sort_ratio_rmt =  1-pnorm(sor-sor_sp,0,sd_sr_q2))
  
  
  
  
  
  
  
  
  vec_parms <- rbind(vec_parms,vec_parms_temp)
  
  
  ## saving results
  dados_out_of_sample$pearson <- cumsum(A)
  dados_out_of_sample$RMT <- cumsum(B)
  
  ### ROBUST 1 ----
  
  library(robust)
  
  # covRob(dd,estim = 'mcd')
  # covRob(dd,estim = 'weighted')
  # covRob(dd,estim = 'donostah')
  # covRob(dd,estim = 'pairwiseGK')
  covariance <- covRob(dados_in_sample %>% select(var_modelo),estim = 'mcd')$cov
  mpearson  <- cov2cor(covRob(dados_in_sample %>% select(var_modelo),estim = 'mcd')$cov)
  vec_mean <- colMeans(dados_in_sample %>% select(var_modelo))
  
  
  
  # Run the eff.frontier function based on no short and 50% alloc. restrictions
  eff <- eff.frontier(covariance,vec_mean, short="no", max.allocation=NULL,
                      risk.premium.up=1, risk.increment=.1)
  
  # Find the optimal portfolio
  eff.optimal.point <- eff %>% filter(sharpe==max(eff$sharpe)) 
  
  # # Conferindo valores
  wm <- eff.optimal.point %>% select(var_modelo)
  # sum(round(wm,10))
  
  
  A <-  as.matrix(data.matrix(dados_out_of_sample %>% select(var_modelo)),
                  ncol = ncol(dados_out_of_sample)) %*%  t(as.matrix(round(wm,10)))
  
  
  dados_out_of_sample <- dados_out_of_sample %>%
    mutate(rs = cumsum(r))
  
  
  
  ## Não vamos comparar shape ratio de diferentes frequências
  q <- 1
  
  
  ### Risk Free
  RF_mean <- mean(dados_out_of_sample$RF,na.rm=T)
  
  
  # Sp500
  auto_cov_sp <- acf(dados_out_of_sample$r,lag=1)$acf[2]
  n_q_sp <- q^(1/2)*(1 + (2*auto_cov_sp)/(1- auto_cov_sp)*(1 - (1-auto_cov_sp^q)/(q*(1-auto_cov_sp))))^(-1/2)
  
  sr_sp <- (mean(dados_out_of_sample$rs)-RF_mean)/sd(dados_out_of_sample$rs)
  sr_sp_month <-  n_q_sp*sr_sp 
  
  # Carteira
  auto_cov_out <- acf(A,lag=1)$acf[2]
  n_q <- q^(1/2)*(1 + (2*auto_cov_out)/(1- auto_cov_out)*(1 - (1-auto_cov_out^q)/(q*(1-auto_cov_out))))^(-1/2)
  sr <- (mean(A)-RF_mean)/sd(A)
  
  sr_month <- n_q*sr 
  
  sr_month_0 <- sr_month 
  VIID_q <- (n_q^2)*(1+1/2*(sr^2))
  VGMM_add <- ifelse(q==1,0,q*(sr)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q <- VIID_q   + VGMM_add 
  
  sd_sr_q <- (VGMM_q/nrow(dados_out_of_sample))^(1/2)
  
  
  
  sor <- n_q*(mean(A)-RF_mean)/DownsideDeviation(A,MAR = RF_mean)
  sor_sp <- (mean(dados_out_of_sample$rs)-RF_mean)/DownsideDeviation(dados_out_of_sample$rs,MAR = RF_mean)
  
  
  VIID_q2 <- (n_q^2)*(1+1/2*(sor^2))
  VGMM_add2 <- ifelse(q==1,0,q*(sor)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q2 <- VIID_q2   + VGMM_add2 
  
  sd_sr_q2 <- (VGMM_q2/nrow(dados_out_of_sample))^(1/2)
  
  
  
  vec_parms_temp <- data.frame(data = nomes_arquivos[i],
                               cov_method = "MCD",
                               initial_in = min(dados_in_sample$DT),
                               end_in = max(dados_in_sample$DT),
                               initial_out = min(dados_out_of_sample$DT),
                               end_out = max(dados_out_of_sample$DT),
                               N_in = nrow(dados_in_sample),
                               N_out = nrow(dados_out_of_sample),
                               X = length(var_modelo),
                               X_pos_n = length(which(wm>0.00)),
                               X_pos_valid_n = length(which(wm>0.01)),
                               X_pos_valid_sum = paste0(round(sum(wm[(which(wm>0.01))])*100,2),"%"),
                               w_max = paste0(round(max(wm[(which(wm>0.01))])*100,2),"%"),
                               w_min = paste0(round(min(wm[(which(wm>0.01))])*100,2),"%"),
                               r_final = paste0(round(cumsum(A)[nrow(dados_out_of_sample)]*100,4),"%"),
                               Q = NA, 
                               lamda_max = NA,
                               eigen_up = NA,
                               eigen_v_up = NA,
                               eigen_max = NA,
                               eigen_v_max = NA,
                               st_dev=sd(A),
                               down_dev=DownsideDeviation(A,MAR = RF_mean),
                               sharpe_ratio = sr_month ,
                               p_sharp_ratio_sp = 1-pnorm(sr_month-sr_sp_month,0,sd_sr_q),
                               p_sharp_ratio_0 = 1-pnorm(sr_month,0,sd_sr_q),
                               rmt_improvement = NA ,
                               sortino_ratio = sor ,
                               p_sort_ratio_sp = 1-pnorm(sor-sor_sp,0,sd_sr_q2),
                               # rmt_improvement_so = NA ,
                               p_sort_ratio_0 = 1-pnorm(sor,0,sd_sr_q2),
                               auto_cov_ar1 = auto_cov_out,
                               p_sharp_ratio_rmt = NA,
                               p_sort_ratio_rmt = NA)
  
  
  
  vec_parms <- rbind(vec_parms,vec_parms_temp)
  
  ### RMT ROBUST 1 ----
  ## ref: J. Bouchaud, M. Potters, Theory of Financial Risks—From Statistical Physics to Risk Management, Cambridge University Press, UK, 2000.
  Q <- nrow(dados_in_sample)/(length(var_modelo))
  ## falta entender sigma!!!
  eigen0 <- eigen(mpearson)
  # eigen0 <- eigen(covariance)
  lamda_max <- (1 + 1/Q + (1/Q)^.5)
  
  
  g1 <- eigen0$values[eigen0$values>=lamda_max]
  g2 <- eigen0$values[eigen0$values<lamda_max]
  eigen1<- diag(c(g1,rep(mean(eigen0$values),length(g2))))
  cor1 <- eigen0$vectors %*% eigen1 %*% t(eigen0$vectors)
  covariance_rmt <-  cor1
  for(ii in 1:ncol(cor1)){
    for(jj in 1:nrow(cor1)){
      covariance_rmt[ii,jj] <- cor1[ii,jj] * (covariance[ii,ii]*covariance[jj,jj])^.5 
    }
  }
  
  # Run the eff.frontier function based on no short and 50% alloc. restrictions
  eff <- eff.frontier(covariance_rmt,vec_mean, short="no", max.allocation=NULL,
                      risk.premium.up=1, risk.increment=.1)
  
  # Find the optimal portfolio
  eff.optimal.point <- eff %>% filter(sharpe==max(eff$sharpe)) 
  
  # Conferindo valores
  wm <- eff.optimal.point %>% select(var_modelo)
  
  B <-  as.matrix(data.matrix(dados_out_of_sample %>% select(var_modelo)),
                  ncol = ncol(dados_out_of_sample)) %*%
    t(as.matrix(round(wm,10)))
  
  
  
  # Carteira
  auto_cov_out <- acf(B,lag=1)$acf[2]
  n_q <- q^(1/2)*(1 + (2*auto_cov_out)/(1- auto_cov_out)*(1 - (1-auto_cov_out^q)/(q*(1-auto_cov_out))))^(-1/2)
  sr <- (mean(B)-RF_mean)/sd(B)
  
  sr_month <- n_q*sr 
  VIID_q <- (n_q^2)*(1+1/2*(sr^2))
  VGMM_add <- ifelse(q==1,0,q*(sr)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q <- VIID_q   + VGMM_add 
  
  sd_sr_q <- (VGMM_q/nrow(dados_out_of_sample))^(1/2)
  
  
  
  sor <- n_q*(mean(B)-RF_mean)/DownsideDeviation(B,MAR = RF_mean)
  sor_sp <- (mean(dados_out_of_sample$rs)-RF_mean)/DownsideDeviation(dados_out_of_sample$rs,MAR = RF_mean)
  
  
  VIID_q2 <- (n_q^2)*(1+1/2*(sor^2))
  VGMM_add2 <- ifelse(q==1,0,q*(sor)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q2 <- VIID_q2   + VGMM_add2 
  
  sd_sr_q2 <- (VGMM_q2/nrow(dados_out_of_sample))^(1/2)
  
  
  
  vec_parms_temp <- data.frame(data = nomes_arquivos[i],
                               cov_method = "MCD RMT",
                               initial_in = min(dados_in_sample$DT),
                               end_in = max(dados_in_sample$DT),
                               initial_out = min(dados_out_of_sample$DT),
                               end_out = max(dados_out_of_sample$DT),
                               N_in = nrow(dados_in_sample),
                               N_out = nrow(dados_out_of_sample),
                               X = length(var_modelo),
                               X_pos_n = length(which(wm>0.00)),
                               X_pos_valid_n = length(which(wm>0.01)),
                               X_pos_valid_sum = paste0(round(sum(wm[(which(wm>0.01))])*100,2),"%"),
                               w_max = paste0(round(max(wm[(which(wm>0.01))])*100,2),"%"),
                               w_min = paste0(round(min(wm[(which(wm>0.01))])*100,2),"%"),
                               r_final = paste0(round(cumsum(B)[nrow(dados_out_of_sample)]*100,4),"%"),
                               Q = Q, 
                               lamda_max = lamda_max,
                               eigen_up = paste0(length(g1)," (",round(length(g1)/(ncol(dados_in_sample) - 1)*100,2),"%)"),
                               eigen_v_up = paste0(round(sum(g1)/(sum(g1)+sum(g2))*100,2),"%"),
                               eigen_max = g1[1]/lamda_max,
                               eigen_v_max = paste0(round(sum(g1[1])/(sum(g1)+sum(g2))*100,2),"%"),
                               st_dev=sd(B),
                               down_dev=DownsideDeviation(B,MAR = RF_mean),
                               sharpe_ratio = sr_month ,
                               p_sharp_ratio_sp = 1-pnorm(sr_month-sr_sp_month,0,sd_sr_q),
                               p_sharp_ratio_0 = 1-pnorm(sr_month,0,sd_sr_q),
                               rmt_improvement =((mean(B-A)+1)^(252)-1)*100,
                               sortino_ratio = sor ,
                               p_sort_ratio_sp = 1-pnorm(sor-sor_sp,0,sd_sr_q2),
                               p_sort_ratio_0 = 1-pnorm(sor,0,sd_sr_q2),
                               auto_cov_ar1 = auto_cov_out,
                               p_sharp_ratio_rmt =  1-pnorm(sr_month-sr_month_0,0,sd_sr_q),
                               p_sort_ratio_rmt =  1-pnorm(sor-sor_sp,0,sd_sr_q2))

  
  vec_parms <- rbind(vec_parms,vec_parms_temp)
  
  
  ## saving results
  dados_out_of_sample$MCD <- cumsum(A)
  dados_out_of_sample$MCD_RMT <- cumsum(B)
  
  
  
  ### ROBUST 2 ----
  
  library(robust)
  
  # covRob(dd,estim = 'mcd')
  # covRob(dd,estim = 'weighted')
  # covRob(dd,estim = 'donostah')
  # covRob(dd,estim = 'pairwiseGK')
  covariance <- covRob(dados_in_sample %>% select(var_modelo),estim = 'weighted')$cov
  mpearson  <- cov2cor(covRob(dados_in_sample %>% select(var_modelo),estim = 'weighted')$cov)
  vec_mean <- colMeans(dados_in_sample %>% select(var_modelo))
  
  
  
  # Run the eff.frontier function based on no short and 50% alloc. restrictions
  eff <- eff.frontier(covariance,vec_mean, short="no", max.allocation=NULL,
                      risk.premium.up=1, risk.increment=.1)
  
  # Find the optimal portfolio
  eff.optimal.point <- eff %>% filter(sharpe==max(eff$sharpe)) 
  
  # # Conferindo valores
  wm <- eff.optimal.point %>% select(var_modelo)
  # sum(round(wm,10))
  
  
  A <-  as.matrix(data.matrix(dados_out_of_sample %>% select(var_modelo)),
                  ncol = ncol(dados_out_of_sample)) %*%  t(as.matrix(round(wm,10)))
  
  
  dados_out_of_sample <- dados_out_of_sample %>%
    mutate(rs = cumsum(r))
  
  
  
  ## Não vamos comparar shape ratio de diferentes frequências
  q <- 1
  
  
  ### Risk Free
  RF_mean <- mean(dados_out_of_sample$RF,na.rm=T)
  
  
  # Sp500
  auto_cov_sp <- acf(dados_out_of_sample$r,lag=1)$acf[2]
  n_q_sp <- q^(1/2)*(1 + (2*auto_cov_sp)/(1- auto_cov_sp)*(1 - (1-auto_cov_sp^q)/(q*(1-auto_cov_sp))))^(-1/2)
  
  sr_sp <- (mean(dados_out_of_sample$rs)-RF_mean)/sd(dados_out_of_sample$rs)
  sr_sp_month <-  n_q_sp*sr_sp 
  
  # Carteira
  auto_cov_out <- acf(A,lag=1)$acf[2]
  n_q <- q^(1/2)*(1 + (2*auto_cov_out)/(1- auto_cov_out)*(1 - (1-auto_cov_out^q)/(q*(1-auto_cov_out))))^(-1/2)
  sr <- (mean(A)-RF_mean)/sd(A)
  
  sr_month <- n_q*sr 
  
  sr_month_0 <- sr_month 
  VIID_q <- (n_q^2)*(1+1/2*(sr^2))
  VGMM_add <- ifelse(q==1,0,q*(sr)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q <- VIID_q   + VGMM_add 
  
  sd_sr_q <- (VGMM_q/nrow(dados_out_of_sample))^(1/2)
  
  
  
  sor <- n_q*(mean(A)-RF_mean)/DownsideDeviation(A,MAR = RF_mean)
  sor_sp <- (mean(dados_out_of_sample$rs)-RF_mean)/DownsideDeviation(dados_out_of_sample$rs,MAR = RF_mean)
  
  
  VIID_q2 <- (n_q^2)*(1+1/2*(sor^2))
  VGMM_add2 <- ifelse(q==1,0,q*(sor)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q2 <- VIID_q2   + VGMM_add2 
  
  sd_sr_q2 <- (VGMM_q2/nrow(dados_out_of_sample))^(1/2)
  
  
  
  vec_parms_temp <- data.frame(data = nomes_arquivos[i],
                               cov_method = "WMCD",
                               initial_in = min(dados_in_sample$DT),
                               end_in = max(dados_in_sample$DT),
                               initial_out = min(dados_out_of_sample$DT),
                               end_out = max(dados_out_of_sample$DT),
                               N_in = nrow(dados_in_sample),
                               N_out = nrow(dados_out_of_sample),
                               X = length(var_modelo),
                               X_pos_n = length(which(wm>0.00)),
                               X_pos_valid_n = length(which(wm>0.01)),
                               X_pos_valid_sum = paste0(round(sum(wm[(which(wm>0.01))])*100,2),"%"),
                               w_max = paste0(round(max(wm[(which(wm>0.01))])*100,2),"%"),
                               w_min = paste0(round(min(wm[(which(wm>0.01))])*100,2),"%"),
                               r_final = paste0(round(cumsum(A)[nrow(dados_out_of_sample)]*100,4),"%"),
                               Q = NA, 
                               lamda_max = NA,
                               eigen_up = NA,
                               eigen_v_up = NA,
                               eigen_max = NA,
                               eigen_v_max = NA,
                               st_dev=sd(A),
                               down_dev=DownsideDeviation(A,MAR = RF_mean),
                               sharpe_ratio = sr_month ,
                               p_sharp_ratio_sp = 1-pnorm(sr_month-sr_sp_month,0,sd_sr_q),
                               p_sharp_ratio_0 = 1-pnorm(sr_month,0,sd_sr_q),
                               rmt_improvement = NA ,
                               sortino_ratio = sor ,
                               p_sort_ratio_sp = 1-pnorm(sor-sor_sp,0,sd_sr_q2),
                               # rmt_improvement_so = NA ,
                               p_sort_ratio_0 = 1-pnorm(sor,0,sd_sr_q2),
                               auto_cov_ar1 = auto_cov_out,
                               p_sharp_ratio_rmt = NA,
                               p_sort_ratio_rmt = NA)
  
  
  
  vec_parms <- rbind(vec_parms,vec_parms_temp)
  
  ### RMT ROBUST 2 ----
  ## ref: J. Bouchaud, M. Potters, Theory of Financial Risks—From Statistical Physics to Risk Management, Cambridge University Press, UK, 2000.
  Q <- nrow(dados_in_sample)/(length(var_modelo))
  ## falta entender sigma!!!
  eigen0 <- eigen(mpearson)
  # eigen0 <- eigen(covariance)
  lamda_max <- (1 + 1/Q + (1/Q)^.5)
  
  
  g1 <- eigen0$values[eigen0$values>=lamda_max]
  g2 <- eigen0$values[eigen0$values<lamda_max]
  eigen1<- diag(c(g1,rep(mean(eigen0$values),length(g2))))
  cor1 <- eigen0$vectors %*% eigen1 %*% t(eigen0$vectors)
  covariance_rmt <-  cor1
  for(ii in 1:ncol(cor1)){
    for(jj in 1:nrow(cor1)){
      covariance_rmt[ii,jj] <- cor1[ii,jj] * (covariance[ii,ii]*covariance[jj,jj])^.5 
    }
  }
  
  # Run the eff.frontier function based on no short and 50% alloc. restrictions
  eff <- eff.frontier(covariance_rmt,vec_mean, short="no", max.allocation=NULL,
                      risk.premium.up=1, risk.increment=.1)
  
  # Find the optimal portfolio
  eff.optimal.point <- eff %>% filter(sharpe==max(eff$sharpe)) 
  
  # Conferindo valores
  wm <- eff.optimal.point %>% select(var_modelo)
  
  B <-  as.matrix(data.matrix(dados_out_of_sample %>% select(var_modelo)),
                  ncol = ncol(dados_out_of_sample)) %*%
    t(as.matrix(round(wm,10)))
  
  
  
  # Carteira
  auto_cov_out <- acf(B,lag=1)$acf[2]
  n_q <- q^(1/2)*(1 + (2*auto_cov_out)/(1- auto_cov_out)*(1 - (1-auto_cov_out^q)/(q*(1-auto_cov_out))))^(-1/2)
  sr <- (mean(B)-RF_mean)/sd(B)
  
  sr_month <- n_q*sr 
  VIID_q <- (n_q^2)*(1+1/2*(sr^2))
  VGMM_add <- ifelse(q==1,0,q*(sr)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q <- VIID_q   + VGMM_add 
  
  sd_sr_q <- (VGMM_q/nrow(dados_out_of_sample))^(1/2)
  
  
  
  sor <- n_q*(mean(B)-RF_mean)/DownsideDeviation(B,MAR = RF_mean)
  sor_sp <- (mean(dados_out_of_sample$rs)-RF_mean)/DownsideDeviation(dados_out_of_sample$rs,MAR = RF_mean)
  
  
  VIID_q2 <- (n_q^2)*(1+1/2*(sor^2))
  VGMM_add2 <- ifelse(q==1,0,q*(sor)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q2 <- VIID_q2   + VGMM_add2 
  
  sd_sr_q2 <- (VGMM_q2/nrow(dados_out_of_sample))^(1/2)
  
  
  
  vec_parms_temp <- data.frame(data = nomes_arquivos[i],
                               cov_method = "WMCD RMT",
                               initial_in = min(dados_in_sample$DT),
                               end_in = max(dados_in_sample$DT),
                               initial_out = min(dados_out_of_sample$DT),
                               end_out = max(dados_out_of_sample$DT),
                               N_in = nrow(dados_in_sample),
                               N_out = nrow(dados_out_of_sample),
                               X = length(var_modelo),
                               X_pos_n = length(which(wm>0.00)),
                               X_pos_valid_n = length(which(wm>0.01)),
                               X_pos_valid_sum = paste0(round(sum(wm[(which(wm>0.01))])*100,2),"%"),
                               w_max = paste0(round(max(wm[(which(wm>0.01))])*100,2),"%"),
                               w_min = paste0(round(min(wm[(which(wm>0.01))])*100,2),"%"),
                               r_final = paste0(round(cumsum(B)[nrow(dados_out_of_sample)]*100,4),"%"),
                               Q = Q, 
                               lamda_max = lamda_max,
                               eigen_up = paste0(length(g1)," (",round(length(g1)/(ncol(dados_in_sample) - 1)*100,2),"%)"),
                               eigen_v_up = paste0(round(sum(g1)/(sum(g1)+sum(g2))*100,2),"%"),
                               eigen_max = g1[1]/lamda_max,
                               eigen_v_max = paste0(round(sum(g1[1])/(sum(g1)+sum(g2))*100,2),"%"),
                               st_dev=sd(B),
                               down_dev=DownsideDeviation(B,MAR = RF_mean),
                               sharpe_ratio = sr_month ,
                               p_sharp_ratio_sp = 1-pnorm(sr_month-sr_sp_month,0,sd_sr_q),
                               p_sharp_ratio_0 = 1-pnorm(sr_month,0,sd_sr_q),
                               rmt_improvement =((mean(B-A)+1)^(252)-1)*100,
                               sortino_ratio = sor ,
                               p_sort_ratio_sp = 1-pnorm(sor-sor_sp,0,sd_sr_q2),
                               p_sort_ratio_0 = 1-pnorm(sor,0,sd_sr_q2),
                               auto_cov_ar1 = auto_cov_out,
                               p_sharp_ratio_rmt =  1-pnorm(sr_month-sr_month_0,0,sd_sr_q),
                               p_sort_ratio_rmt =  1-pnorm(sor-sor_sp,0,sd_sr_q2))
  
  
  vec_parms <- rbind(vec_parms,vec_parms_temp)
  
  
  ## saving results
  dados_out_of_sample$WMCD <- cumsum(A)
  dados_out_of_sample$WMCD_RMT <- cumsum(B)
  
  
  
  
  
  
  ### ROBUST 4 ----
  
  library(robust)
  
  # covRob(dd,estim = 'mcd')
  # covRob(dd,estim = 'weighted')
  # covRob(dd,estim = 'donostah')
  # covRob(dd,estim = 'pairwiseGK')
  covariance <- covRob(dados_in_sample %>% select(var_modelo),estim = 'pairwiseGK')$cov
  mpearson  <- cov2cor(covRob(dados_in_sample %>% select(var_modelo),estim = 'pairwiseGK')$cov)
  vec_mean <- colMeans(dados_in_sample %>% select(var_modelo))
  
  
  
  # Run the eff.frontier function based on no short and 50% alloc. restrictions
  eff <- eff.frontier(covariance,vec_mean, short="no", max.allocation=NULL,
                      risk.premium.up=1, risk.increment=.1)
  
  # Find the optimal portfolio
  eff.optimal.point <- eff %>% filter(sharpe==max(eff$sharpe)) 
  
  # # Conferindo valores
  wm <- eff.optimal.point %>% select(var_modelo)
  # sum(round(wm,10))
  
  
  A <-  as.matrix(data.matrix(dados_out_of_sample %>% select(var_modelo)),
                  ncol = ncol(dados_out_of_sample)) %*%  t(as.matrix(round(wm,10)))
  
  
  dados_out_of_sample <- dados_out_of_sample %>%
    mutate(rs = cumsum(r))
  
  
  
  ## Não vamos comparar shape ratio de diferentes frequências
  q <- 1
  
  
  ### Risk Free
  RF_mean <- mean(dados_out_of_sample$RF,na.rm=T)
  
  
  # Sp500
  auto_cov_sp <- acf(dados_out_of_sample$r,lag=1)$acf[2]
  n_q_sp <- q^(1/2)*(1 + (2*auto_cov_sp)/(1- auto_cov_sp)*(1 - (1-auto_cov_sp^q)/(q*(1-auto_cov_sp))))^(-1/2)
  
  sr_sp <- (mean(dados_out_of_sample$rs)-RF_mean)/sd(dados_out_of_sample$rs)
  sr_sp_month <-  n_q_sp*sr_sp 
  
  # Carteira
  auto_cov_out <- acf(A,lag=1)$acf[2]
  n_q <- q^(1/2)*(1 + (2*auto_cov_out)/(1- auto_cov_out)*(1 - (1-auto_cov_out^q)/(q*(1-auto_cov_out))))^(-1/2)
  sr <- (mean(A)-RF_mean)/sd(A)
  
  sr_month <- n_q*sr 
  
  sr_month_0 <- sr_month 
  VIID_q <- (n_q^2)*(1+1/2*(sr^2))
  VGMM_add <- ifelse(q==1,0,q*(sr)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q <- VIID_q   + VGMM_add 
  
  sd_sr_q <- (VGMM_q/nrow(dados_out_of_sample))^(1/2)
  
  
  
  sor <- n_q*(mean(A)-RF_mean)/DownsideDeviation(A,MAR = RF_mean)
  sor_sp <- (mean(dados_out_of_sample$rs)-RF_mean)/DownsideDeviation(dados_out_of_sample$rs,MAR = RF_mean)
  
  
  VIID_q2 <- (n_q^2)*(1+1/2*(sor^2))
  VGMM_add2 <- ifelse(q==1,0,q*(sor)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q2 <- VIID_q2   + VGMM_add2 
  
  sd_sr_q2 <- (VGMM_q2/nrow(dados_out_of_sample))^(1/2)
  
  
  
  vec_parms_temp <- data.frame(data = nomes_arquivos[i],
                               cov_method = "OGK",
                               initial_in = min(dados_in_sample$DT),
                               end_in = max(dados_in_sample$DT),
                               initial_out = min(dados_out_of_sample$DT),
                               end_out = max(dados_out_of_sample$DT),
                               N_in = nrow(dados_in_sample),
                               N_out = nrow(dados_out_of_sample),
                               X = length(var_modelo),
                               X_pos_n = length(which(wm>0.00)),
                               X_pos_valid_n = length(which(wm>0.01)),
                               X_pos_valid_sum = paste0(round(sum(wm[(which(wm>0.01))])*100,2),"%"),
                               w_max = paste0(round(max(wm[(which(wm>0.01))])*100,2),"%"),
                               w_min = paste0(round(min(wm[(which(wm>0.01))])*100,2),"%"),
                               r_final = paste0(round(cumsum(A)[nrow(dados_out_of_sample)]*100,4),"%"),
                               Q = NA, 
                               lamda_max = NA,
                               eigen_up = NA,
                               eigen_v_up = NA,
                               eigen_max = NA,
                               eigen_v_max = NA,
                               st_dev=sd(A),
                               down_dev=DownsideDeviation(A,MAR = RF_mean),
                               sharpe_ratio = sr_month ,
                               p_sharp_ratio_sp = 1-pnorm(sr_month-sr_sp_month,0,sd_sr_q),
                               p_sharp_ratio_0 = 1-pnorm(sr_month,0,sd_sr_q),
                               rmt_improvement = NA ,
                               sortino_ratio = sor ,
                               p_sort_ratio_sp = 1-pnorm(sor-sor_sp,0,sd_sr_q2),
                               # rmt_improvement_so = NA ,
                               p_sort_ratio_0 = 1-pnorm(sor,0,sd_sr_q2),
                               auto_cov_ar1 = auto_cov_out,
                               p_sharp_ratio_rmt = NA,
                               p_sort_ratio_rmt = NA)
  
  
  
  vec_parms <- rbind(vec_parms,vec_parms_temp)
  
  ### RMT ROBUST 4 ----
  ## ref: J. Bouchaud, M. Potters, Theory of Financial Risks—From Statistical Physics to Risk Management, Cambridge University Press, UK, 2000.
  Q <- nrow(dados_in_sample)/(length(var_modelo))
  ## falta entender sigma!!!
  eigen0 <- eigen(mpearson)
  # eigen0 <- eigen(covariance)
  lamda_max <- (1 + 1/Q + (1/Q)^.5)
  
  
  g1 <- eigen0$values[eigen0$values>=lamda_max]
  g2 <- eigen0$values[eigen0$values<lamda_max]
  eigen1<- diag(c(g1,rep(mean(eigen0$values),length(g2))))
  cor1 <- eigen0$vectors %*% eigen1 %*% t(eigen0$vectors)
  covariance_rmt <-  cor1
  for(ii in 1:ncol(cor1)){
    for(jj in 1:nrow(cor1)){
      covariance_rmt[ii,jj] <- cor1[ii,jj] * (covariance[ii,ii]*covariance[jj,jj])^.5 
    }
  }
  
  # Run the eff.frontier function based on no short and 50% alloc. restrictions
  eff <- eff.frontier(covariance_rmt,vec_mean, short="no", max.allocation=NULL,
                      risk.premium.up=1, risk.increment=.1)
  
  # Find the optimal portfolio
  eff.optimal.point <- eff %>% filter(sharpe==max(eff$sharpe)) 
  
  # Conferindo valores
  wm <- eff.optimal.point %>% select(var_modelo)
  
  B <-  as.matrix(data.matrix(dados_out_of_sample %>% select(var_modelo)),
                  ncol = ncol(dados_out_of_sample)) %*%
    t(as.matrix(round(wm,10)))
  
  
  
  # Carteira
  auto_cov_out <- acf(B,lag=1)$acf[2]
  n_q <- q^(1/2)*(1 + (2*auto_cov_out)/(1- auto_cov_out)*(1 - (1-auto_cov_out^q)/(q*(1-auto_cov_out))))^(-1/2)
  sr <- (mean(B)-RF_mean)/sd(B)
  
  sr_month <- n_q*sr 
  VIID_q <- (n_q^2)*(1+1/2*(sr^2))
  VGMM_add <- ifelse(q==1,0,q*(sr)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q <- VIID_q   + VGMM_add 
  
  sd_sr_q <- (VGMM_q/nrow(dados_out_of_sample))^(1/2)
  
  
  
  sor <- n_q*(mean(B)-RF_mean)/DownsideDeviation(B,MAR = RF_mean)
  sor_sp <- (mean(dados_out_of_sample$rs)-RF_mean)/DownsideDeviation(dados_out_of_sample$rs,MAR = RF_mean)
  
  
  VIID_q2 <- (n_q^2)*(1+1/2*(sor^2))
  VGMM_add2 <- ifelse(q==1,0,q*(sor)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q2 <- VIID_q2   + VGMM_add2 
  
  sd_sr_q2 <- (VGMM_q2/nrow(dados_out_of_sample))^(1/2)
  
  
  
  vec_parms_temp <- data.frame(data = nomes_arquivos[i],
                               cov_method = "OGK RMT",
                               initial_in = min(dados_in_sample$DT),
                               end_in = max(dados_in_sample$DT),
                               initial_out = min(dados_out_of_sample$DT),
                               end_out = max(dados_out_of_sample$DT),
                               N_in = nrow(dados_in_sample),
                               N_out = nrow(dados_out_of_sample),
                               X = length(var_modelo),
                               X_pos_n = length(which(wm>0.00)),
                               X_pos_valid_n = length(which(wm>0.01)),
                               X_pos_valid_sum = paste0(round(sum(wm[(which(wm>0.01))])*100,2),"%"),
                               w_max = paste0(round(max(wm[(which(wm>0.01))])*100,2),"%"),
                               w_min = paste0(round(min(wm[(which(wm>0.01))])*100,2),"%"),
                               r_final = paste0(round(cumsum(B)[nrow(dados_out_of_sample)]*100,4),"%"),
                               Q = Q, 
                               lamda_max = lamda_max,
                               eigen_up = paste0(length(g1)," (",round(length(g1)/(ncol(dados_in_sample) - 1)*100,2),"%)"),
                               eigen_v_up = paste0(round(sum(g1)/(sum(g1)+sum(g2))*100,2),"%"),
                               eigen_max = g1[1]/lamda_max,
                               eigen_v_max = paste0(round(sum(g1[1])/(sum(g1)+sum(g2))*100,2),"%"),
                               st_dev=sd(B),
                               down_dev=DownsideDeviation(B,MAR = RF_mean),
                               sharpe_ratio = sr_month ,
                               p_sharp_ratio_sp = 1-pnorm(sr_month-sr_sp_month,0,sd_sr_q),
                               p_sharp_ratio_0 = 1-pnorm(sr_month,0,sd_sr_q),
                               rmt_improvement =((mean(B-A)+1)^(252)-1)*100,
                               sortino_ratio = sor ,
                               p_sort_ratio_sp = 1-pnorm(sor-sor_sp,0,sd_sr_q2),
                               p_sort_ratio_0 = 1-pnorm(sor,0,sd_sr_q2),
                               auto_cov_ar1 = auto_cov_out,
                               p_sharp_ratio_rmt =  1-pnorm(sr_month-sr_month_0,0,sd_sr_q),
                               p_sort_ratio_rmt =  1-pnorm(sor-sor_sp,0,sd_sr_q2))
  
  
  vec_parms <- rbind(vec_parms,vec_parms_temp)
  
  
  ## saving results
  dados_out_of_sample$OGK <- cumsum(A)
  dados_out_of_sample$OGK_RMT <- cumsum(B)
  
  
  
  
  

  ### KERNEL GAUSS  ----
  
  
  mpearson  <- cor(dados_in_sample %>% select(var_modelo))
  covariance <- var(dados_in_sample %>% select(var_modelo))
  vec_mean <- colMeans(dados_in_sample %>% select(var_modelo))

  rbf1 <- rbfdot(sigma = 1/(((sum(diag(covariance))/nrow(covariance))^2)*2))
  K11 <- kernelMatrix(rbf1,covariance)
  
  covariance_krmt <-  K11
  for(ii in 1:ncol(K11)){
    for(jj in 1:nrow(K11)){
      covariance_krmt[ii,jj] <- K11[ii,jj] * (covariance[ii,ii]*covariance[jj,jj])^.5 
    }
  }
  
  
  
  # Run the eff.frontier function based on no short and 50% alloc. restrictions
  eff <- eff.frontier(covariance_krmt,vec_mean, short="no", max.allocation=NULL,
                      risk.premium.up=1, risk.increment=.1)
  
  # Find the optimal portfolio
  eff.optimal.point <- eff %>% filter(sharpe==max(eff$sharpe)) 
  
  
  # Conferindo valores
  wm <- eff.optimal.point %>% select(var_modelo)
  
  C <-  as.matrix(data.matrix(dados_out_of_sample %>% select(var_modelo)),
                  ncol = ncol(dados_out_of_sample)) %*%
    t(as.matrix(round(wm,10)))
  
  # dados_out_of_sample$KERNEL <- cumsum(C)
  
  # Carteira
  auto_cov_out <- acf(C,lag=1)$acf[2]
  n_q <- q^(1/2)*(1 + (2*auto_cov_out)/(1- auto_cov_out)*(1 - (1-auto_cov_out^q)/(q*(1-auto_cov_out))))^(-1/2)
  sr <- (mean(C)-RF_mean)/sd(C)
  
  sr_month <- n_q*sr 
  sr_month_0 <- sr_month
  VIID_q <- (n_q^2)*(1+1/2*(sr^2))
  VGMM_add <- ifelse(q==1,0,q*(sr)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q <- VIID_q   + VGMM_add 
  
  sd_sr_q <- (VGMM_q/nrow(dados_out_of_sample))^(1/2)
  
  
  sor <- n_q*(mean(C)-RF_mean)/DownsideDeviation(C,MAR = RF_mean)
  sor_sp <- (mean(dados_out_of_sample$rs)-RF_mean)/DownsideDeviation(dados_out_of_sample$rs,MAR = RF_mean)
  
  
  VIID_q2 <- (n_q^2)*(1+1/2*(sor^2))
  VGMM_add2 <- ifelse(q==1,0,q*(sor)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q2 <- VIID_q2   + VGMM_add2 
  
  sd_sr_q2 <- (VGMM_q2/nrow(dados_out_of_sample))^(1/2)
  
  
  
  
  vec_parms_temp <- data.frame(data = nomes_arquivos[i],
                               cov_method = "GAUSS",
                               initial_in = min(dados_in_sample$DT),
                               end_in = max(dados_in_sample$DT),
                               initial_out = min(dados_out_of_sample$DT),
                               end_out = max(dados_out_of_sample$DT),
                               N_in = nrow(dados_in_sample),
                               N_out = nrow(dados_out_of_sample),
                               X = length(var_modelo),
                               X_pos_n = length(which(wm>0.00)),
                               X_pos_valid_n = length(which(wm>0.01)),
                               X_pos_valid_sum = paste0(round(sum(wm[(which(wm>0.01))])*100,2),"%"),
                               w_max = paste0(round(max(wm[(which(wm>0.01))])*100,2),"%"),
                               w_min = paste0(round(min(wm[(which(wm>0.01))])*100,2),"%"),
                               r_final = paste0(round(cumsum(C)[nrow(dados_out_of_sample)]*100,4),"%"),
                               Q = Q, 
                               lamda_max = lamda_max,
                               eigen_up = NA,
                               eigen_v_up = NA,
                               eigen_max = NA,
                               eigen_v_max = NA,
                               st_dev=sd(A),
                               down_dev=DownsideDeviation(C,MAR = RF_mean),
                               sharpe_ratio = sr_month ,
                               p_sharp_ratio_sp = 1-pnorm(sr_month-sr_sp_month,0,sd_sr_q),
                               p_sharp_ratio_0 = 1-pnorm(sr_month,0,sd_sr_q),
                               rmt_improvement = NA ,
                               sortino_ratio = sor ,
                               p_sort_ratio_sp = 1-pnorm(sor-sor_sp,0,sd_sr_q2),
                               # rmt_improvement_so = NA ,
                               p_sort_ratio_0 = 1-pnorm(sor,0,sd_sr_q2),
                               auto_cov_ar1 = auto_cov_out,
                               p_sharp_ratio_rmt = NA,
                               p_sort_ratio_rmt = NA)
  
  
  vec_parms <- rbind(vec_parms_temp,vec_parms)
  
  
  ### KRMT GAUSS ----
  
  eigen0 <- eigen(K11)
  g1 <- eigen0$values[eigen0$values>=lamda_max]
  g2 <- eigen0$values[eigen0$values<lamda_max]
  eigen1<- diag(c(g1,rep(mean(g2),length(g2))))
  cor1 <- eigen0$vectors %*% eigen1 %*% t(eigen0$vectors)
  covariance_krmt <-  cor1
  for(ii in 1:ncol(cor1)){
    for(jj in 1:nrow(cor1)){
      covariance_krmt[ii,jj] <- cor1[ii,jj] * (covariance[ii,ii]*covariance[jj,jj])^.5 
    }
  }
  
  
  # Run the eff.frontier function based on no short and 50% alloc. restrictions
  eff <- eff.frontier(covariance_krmt,vec_mean, short="no", max.allocation=NULL,
                      risk.premium.up=1, risk.increment=.1)
  
  # Find the optimal portfolio
  eff.optimal.point <- eff %>% filter(sharpe==max(eff$sharpe)) 
  
  
  # Conferindo valores
  wm <- eff.optimal.point %>% select(var_modelo)
  
  D <-  as.matrix(data.matrix(dados_out_of_sample %>% select(var_modelo)),
                  ncol = ncol(dados_out_of_sample)) %*%
    t(as.matrix(round(wm,10)))
  
  # dados_out_of_sample$KRMT <- cumsum(D)
  
  
  # Carteira
  auto_cov_out <- acf(D,lag=1)$acf[2]
  n_q <- q^(1/2)*(1 + (2*auto_cov_out)/(1- auto_cov_out)*(1 - (1-auto_cov_out^q)/(q*(1-auto_cov_out))))^(-1/2)
  sr <- (mean(D)-RF_mean)/sd(D)
  
  sr_month <- n_q*sr 
  VIID_q <- (n_q^2)*(1+1/2*(sr^2))
  VGMM_q <- VIID_q + (q*sr)^2*(sum((1- (1:(q-1))/q)^2))
  
  sd_sr_q <- (VGMM_q/nrow(dados_out_of_sample))^(1/2)
  
  
  
  sor <- n_q*(mean(D)-RF_mean)/DownsideDeviation(D,MAR = RF_mean)
  sor_sp <- (mean(dados_out_of_sample$rs)-RF_mean)/DownsideDeviation(dados_out_of_sample$rs,MAR = RF_mean)
  
  
  VIID_q2 <- (n_q^2)*(1+1/2*(sor^2))
  VGMM_add2 <- ifelse(q==1,0,q*(sor)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q2 <- VIID_q2   + VGMM_add2 
  
  sd_sr_q2 <- (VGMM_q2/nrow(dados_out_of_sample))^(1/2)
  
  
  
  vec_parms_temp <- data.frame(data = nomes_arquivos[i],
                               cov_method = "GAUSS RMT",
                               initial_in = min(dados_in_sample$DT),
                               end_in = max(dados_in_sample$DT),
                               initial_out = min(dados_out_of_sample$DT),
                               end_out = max(dados_out_of_sample$DT),
                               N_in = nrow(dados_in_sample),
                               N_out = nrow(dados_out_of_sample),
                               X = length(var_modelo),
                               X_pos_n = length(which(wm>0.00)),
                               X_pos_valid_n = length(which(wm>0.01)),
                               X_pos_valid_sum = paste0(round(sum(wm[(which(wm>0.01))])*100,2),"%"),
                               w_max = paste0(round(max(wm[(which(wm>0.01))])*100,2),"%"),
                               w_min = paste0(round(min(wm[(which(wm>0.01))])*100,2),"%"),
                               r_final = paste0(round(cumsum(D)[nrow(dados_out_of_sample)]*100,4),"%"),
                               Q = Q, 
                               lamda_max = lamda_max,
                               eigen_up = paste0(length(g1)," (",round(length(g1)/(ncol(dados_in_sample) - 1)*100,2),"%)"),
                               eigen_v_up = paste0(round(sum(g1)/(sum(g1)+sum(g2))*100,2),"%"),
                               eigen_max = g1[1]/lamda_max,
                               eigen_v_max = paste0(round(sum(g1[1])/(sum(g1)+sum(g2))*100,2),"%"),
                               st_dev=sd(B),
                               down_dev=DownsideDeviation(D,MAR = RF_mean),
                               sharpe_ratio = sr_month ,
                               p_sharp_ratio_sp = 1-pnorm(sr_month-sr_sp_month,0,sd_sr_q),
                               p_sharp_ratio_0 = 1-pnorm(sr_month,0,sd_sr_q),
                               rmt_improvement =((mean(D-C)+1)^(252)-1)*100,
                               sortino_ratio = sor ,
                               p_sort_ratio_sp = 1-pnorm(sor-sor_sp,0,sd_sr_q2),
                               p_sort_ratio_0 = 1-pnorm(sor,0,sd_sr_q2),
                               auto_cov_ar1 = auto_cov_out,
                               p_sharp_ratio_rmt =  1-pnorm(sr_month-sr_month_0,0,sd_sr_q),
                               p_sort_ratio_rmt =  1-pnorm(sor-sor_sp,0,sd_sr_q2))
  
  
  
  
  vec_parms <- rbind(vec_parms_temp,vec_parms)
  
  dados_out_of_sample$GAUSS <- cumsum(C)
  dados_out_of_sample$GAUSS_RMT <- cumsum(D)
  
  
  ### KERNEL POLY 2  ----
  
  
  # mpearson  <- cor(dados_in_sample %>% select(var_modelo))
  # covariance <- var(dados_in_sample %>% select(var_modelo))
  # vec_mean <- colMeans(dados_in_sample %>% select(var_modelo))
  # 
  # 
  # # poly <- polydot(degree=2,scale = sqrt(sum(diag(mpearson))*sum(diag(covariance))),offset = sqrt(mean(mpearson)*mean(covariance)))
  # poly <- polydot(degree=2,scale = 1/((sum(diag(covariance))/nrow(covariance))^2)*2)
  # K11 <- kernelMatrix(poly,covariance)
  # # K11 <- kernelMatrix(poly,mpearson)
  # # K11 <- K11-mean(K11)/((sd(K11)))
  # 
  # covariance_krmt <- K11
  # 
  # # vec_sd <- diag(K11)^.5
  # # for(ii in 1:ncol(K11)){
  # #   for(jj in 1:nrow(K11)){
  # #     covariance_krmt[ii,jj] <- K11[ii,jj] / (vec_sd[ii]*vec_sd[jj])
  # #     # covariance_krmt[ii,jj] <- K11[ii,jj] / (covariance[ii]*covariance[jj])
  # #   }
  # # }
  # 
  # 
  # # covariance_krmt_poly <-  K12_cor
  # for(ii in 1:ncol(K11)){
  #   for(jj in 1:nrow(K11)){
  #     covariance_krmt[ii,jj] <- covariance_krmt[ii,jj] * (covariance[ii,ii]*covariance[jj,jj])^.5
  #   }
  # }

  
  covariance_krmt <- matrix(NA, nrow = ncol(dados_in_sample %>% select(var_modelo)), 
                            ncol = ncol(dados_in_sample %>% select(var_modelo)))
  
  dd <- 1
  for (i in 1:ncol(dados_in_sample %>% select(var_modelo))){
    for (j in i:ncol(dados_in_sample %>% select(var_modelo))){
      covariance_krmt[i,j] <- covariance_krmt[j,i] <- 
        ((dados_in_sample %>% select(var_modelo))[,i]%*%(dados_in_sample %>% select(var_modelo))[,j]+dd)^2
    }
  }
  
  colnames(covariance_krmt) <- rownames(covariance_krmt) <- colnames(dados_in_sample %>% select(var_modelo))
  
  # covariance_krmt <- cov(covariance_krmt)
  
  # Run the eff.frontier function based on no short and 50% alloc. restrictions
  eff <- eff.frontier(covariance_krmt,vec_mean, short="no", max.allocation=NULL,
                      risk.premium.up=1, risk.increment=.1)
  
  # Find the optimal portfolio
  eff.optimal.point <- eff %>% filter(sharpe==max(eff$sharpe)) 
  
  
  # Conferindo valores
  wm <- eff.optimal.point %>% select(var_modelo)
  
  C <-  as.matrix(data.matrix(dados_out_of_sample %>% select(var_modelo)),
                  ncol = ncol(dados_out_of_sample)) %*%
    t(as.matrix(round(wm,10)))
  
  # dados_out_of_sample$KERNEL <- cumsum(C)
  
  # Carteira
  auto_cov_out <- acf(C,lag=1)$acf[2]
  n_q <- q^(1/2)*(1 + (2*auto_cov_out)/(1- auto_cov_out)*(1 - (1-auto_cov_out^q)/(q*(1-auto_cov_out))))^(-1/2)
  sr <- (mean(C)-RF_mean)/sd(C)
  
  sr_month <- n_q*sr 
  sr_month_0 <- sr_month
  VIID_q <- (n_q^2)*(1+1/2*(sr^2))
  VGMM_add <- ifelse(q==1,0,q*(sr)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q <- VIID_q   + VGMM_add 
  
  sd_sr_q <- (VGMM_q/nrow(dados_out_of_sample))^(1/2)
  
  
  sor <- n_q*(mean(C)-RF_mean)/DownsideDeviation(C,MAR = RF_mean)
  sor_sp <- (mean(dados_out_of_sample$rs)-RF_mean)/DownsideDeviation(dados_out_of_sample$rs,MAR = RF_mean)
  
  
  VIID_q2 <- (n_q^2)*(1+1/2*(sor^2))
  VGMM_add2 <- ifelse(q==1,0,q*(sor)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q2 <- VIID_q2   + VGMM_add2 
  
  sd_sr_q2 <- (VGMM_q2/nrow(dados_out_of_sample))^(1/2)
  
  
  
  
  vec_parms_temp <- data.frame(data = nomes_arquivos[i],
                               cov_method = "POLY 2",
                               initial_in = min(dados_in_sample$DT),
                               end_in = max(dados_in_sample$DT),
                               initial_out = min(dados_out_of_sample$DT),
                               end_out = max(dados_out_of_sample$DT),
                               N_in = nrow(dados_in_sample),
                               N_out = nrow(dados_out_of_sample),
                               X = length(var_modelo),
                               X_pos_n = length(which(wm>0.00)),
                               X_pos_valid_n = length(which(wm>0.01)),
                               X_pos_valid_sum = paste0(round(sum(wm[(which(wm>0.01))])*100,2),"%"),
                               w_max = paste0(round(max(wm[(which(wm>0.01))])*100,2),"%"),
                               w_min = paste0(round(min(wm[(which(wm>0.01))])*100,2),"%"),
                               r_final = paste0(round(cumsum(C)[nrow(dados_out_of_sample)]*100,4),"%"),
                               Q = Q, 
                               lamda_max = lamda_max,
                               eigen_up = NA,
                               eigen_v_up = NA,
                               eigen_max = NA,
                               eigen_v_max = NA,
                               st_dev=sd(A),
                               down_dev=DownsideDeviation(C,MAR = RF_mean),
                               sharpe_ratio = sr_month ,
                               p_sharp_ratio_sp = 1-pnorm(sr_month-sr_sp_month,0,sd_sr_q),
                               p_sharp_ratio_0 = 1-pnorm(sr_month,0,sd_sr_q),
                               rmt_improvement = NA ,
                               sortino_ratio = sor ,
                               p_sort_ratio_sp = 1-pnorm(sor-sor_sp,0,sd_sr_q2),
                               # rmt_improvement_so = NA ,
                               p_sort_ratio_0 = 1-pnorm(sor,0,sd_sr_q2),
                               auto_cov_ar1 = auto_cov_out,
                               p_sharp_ratio_rmt = NA,
                               p_sort_ratio_rmt = NA)
  
  
  vec_parms <- rbind(vec_parms_temp,vec_parms)
  
  
  ### KRMT POLY 2 ----
  
  eigen0 <- eigen(K11)
  g1 <- eigen0$values[eigen0$values>=lamda_max]
  g2 <- eigen0$values[eigen0$values<lamda_max]
  eigen1<- diag(c(g1,rep(mean(g2),length(g2))))
  cor1 <- eigen0$vectors %*% eigen1 %*% t(eigen0$vectors)
  covariance_krmt <-  cor1
  for(ii in 1:ncol(cor1)){
    for(jj in 1:nrow(cor1)){
      covariance_krmt[ii,jj] <- cor1[ii,jj] * (covariance[ii,ii]*covariance[jj,jj])^.5 
    }
  }
  
  
  # Run the eff.frontier function based on no short and 50% alloc. restrictions
  eff <- eff.frontier(covariance_krmt,vec_mean, short="no", max.allocation=NULL,
                      risk.premium.up=1, risk.increment=.1)
  
  # Find the optimal portfolio
  eff.optimal.point <- eff %>% filter(sharpe==max(eff$sharpe)) 
  
  
  # Conferindo valores
  wm <- eff.optimal.point %>% select(var_modelo)
  
  D <-  as.matrix(data.matrix(dados_out_of_sample %>% select(var_modelo)),
                  ncol = ncol(dados_out_of_sample)) %*%
    t(as.matrix(round(wm,10)))
  
  # dados_out_of_sample$KRMT <- cumsum(D)
  
  
  # Carteira
  auto_cov_out <- acf(D,lag=1)$acf[2]
  n_q <- q^(1/2)*(1 + (2*auto_cov_out)/(1- auto_cov_out)*(1 - (1-auto_cov_out^q)/(q*(1-auto_cov_out))))^(-1/2)
  sr <- (mean(D)-RF_mean)/sd(D)
  
  sr_month <- n_q*sr 
  VIID_q <- (n_q^2)*(1+1/2*(sr^2))
  VGMM_q <- VIID_q + (q*sr)^2*(sum((1- (1:(q-1))/q)^2))
  
  sd_sr_q <- (VGMM_q/nrow(dados_out_of_sample))^(1/2)
  
  
  
  sor <- n_q*(mean(D)-RF_mean)/DownsideDeviation(D[,1],MAR = RF_mean)
  sor_sp <- (mean(dados_out_of_sample$rs)-RF_mean)/DownsideDeviation(dados_out_of_sample$rs,MAR = RF_mean)
  
  
  VIID_q2 <- (n_q^2)*(1+1/2*(sor^2))
  VGMM_add2 <- ifelse(q==1,0,q*(sor)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q2 <- VIID_q2   + VGMM_add2 
  
  sd_sr_q2 <- (VGMM_q2/nrow(dados_out_of_sample))^(1/2)
  
  
  
  vec_parms_temp <- data.frame(data = nomes_arquivos[i],
                               cov_method = "POLY 2 RMT",
                               initial_in = min(dados_in_sample$DT),
                               end_in = max(dados_in_sample$DT),
                               initial_out = min(dados_out_of_sample$DT),
                               end_out = max(dados_out_of_sample$DT),
                               N_in = nrow(dados_in_sample),
                               N_out = nrow(dados_out_of_sample),
                               X = length(var_modelo),
                               X_pos_n = length(which(wm>0.00)),
                               X_pos_valid_n = length(which(wm>0.01)),
                               X_pos_valid_sum = paste0(round(sum(wm[(which(wm>0.01)),1])*100,2),"%"),
                               w_max = paste0(round(max(wm[(which(wm>0.01)),1])*100,2),"%"),
                               w_min = paste0(round(min(wm[(which(wm>0.01)),1])*100,2),"%"),
                               r_final = paste0(round(cumsum(D)[nrow(dados_out_of_sample)]*100,4),"%"),
                               Q = Q, 
                               lamda_max = lamda_max,
                               eigen_up = paste0(length(g1)," (",round(length(g1)/(ncol(dados_in_sample) - 1)*100,2),"%)"),
                               eigen_v_up = paste0(round(sum(g1)/(sum(g1)+sum(g2))*100,2),"%"),
                               eigen_max = g1[1]/lamda_max,
                               eigen_v_max = paste0(round(sum(g1[1])/(sum(g1)+sum(g2))*100,2),"%"),
                               st_dev=sd(B),
                               down_dev=DownsideDeviation(D[,1],MAR = RF_mean),
                               sharpe_ratio = sr_month ,
                               p_sharp_ratio_sp = 1-pnorm(sr_month-sr_sp_month,0,sd_sr_q),
                               p_sharp_ratio_0 = 1-pnorm(sr_month,0,sd_sr_q),
                               rmt_improvement =((mean(D[,1]-C)+1)^(252)-1)*100,
                               sortino_ratio = sor ,
                               p_sort_ratio_sp = 1-pnorm(sor-sor_sp,0,sd_sr_q2),
                               p_sort_ratio_0 = 1-pnorm(sor,0,sd_sr_q2),
                               auto_cov_ar1 = auto_cov_out,
                               p_sharp_ratio_rmt =  1-pnorm(sr_month-sr_month_0,0,sd_sr_q),
                               p_sort_ratio_rmt =  1-pnorm(sor-sor_sp,0,sd_sr_q2))
  
  
  
  
  vec_parms <- rbind(vec_parms_temp,vec_parms)
  
  
  dados_out_of_sample$POLY2 <- cumsum(C[,1])
  dados_out_of_sample$POLY2_RMT <- cumsum(D[,1])
  
  
  ### KERNEL POLY 3  ----
  
  
  # mpearson  <- cor(dados_in_sample %>% select(var_modelo))
  # covariance <- var(dados_in_sample %>% select(var_modelo))
  # vec_mean <- colMeans(dados_in_sample %>% select(var_modelo))
  # 
  # # poly <- polydot(degree=2,scale = sqrt(sum(diag(mpearson))*sum(diag(covariance))),offset = sqrt(mean(mpearson)*mean(covariance)))
  # poly <- polydot(degree=3,scale = 1/((sum(diag(covariance))/nrow(covariance))^2)*2)
  # K11 <- kernelMatrix(poly,covariance)
  # # K11 <- kernelMatrix(poly,mpearson)
  # # K11 <- K11-mean(K11)/((sd(K11)))
  # 
  # covariance_krmt <- K11
  # 
  # # vec_sd <- diag(K11)^.5
  # # for(ii in 1:ncol(K11)){
  # #   for(jj in 1:nrow(K11)){
  # #     covariance_krmt[ii,jj] <- K11[ii,jj] / (vec_sd[ii]*vec_sd[jj])
  # #     # covariance_krmt[ii,jj] <- K11[ii,jj] / (covariance[ii]*covariance[jj]) 
  # #   }
  # # }
  # 
  # 
  # # covariance_krmt_poly <-  K12_cor
  # for(ii in 1:ncol(K11)){
  #   for(jj in 1:nrow(K11)){
  #     covariance_krmt[ii,jj] <- covariance_krmt[ii,jj] * (covariance[ii,ii]*covariance[jj,jj])^.5
  #   }
  # }
  
  # covariance_krmt <- cov(covariance_krmt)
  
  
  covariance_krmt <- matrix(NA, nrow = ncol(dados_in_sample %>% select(var_modelo)), 
                            ncol = ncol(dados_in_sample %>% select(var_modelo)))
  
  dd <- 1
  for (i in 1:ncol(dados_in_sample %>% select(var_modelo))){
    for (j in i:ncol(dados_in_sample %>% select(var_modelo))){
      covariance_krmt[i,j] <- covariance_krmt[j,i] <- 
        ((dados_in_sample %>% select(var_modelo))[,i]%*%(dados_in_sample %>% select(var_modelo))[,j]+dd)^3
    }
  }
  
  colnames(covariance_krmt) <- rownames(covariance_krmt) <- colnames(dados_in_sample %>% select(var_modelo))
  
  
  
  # Run the eff.frontier function based on no short and 50% alloc. restrictions
  eff <- eff.frontier(covariance_krmt,vec_mean, short="no", max.allocation=NULL,
                      risk.premium.up=1, risk.increment=.1)
  
  # Find the optimal portfolio
  eff.optimal.point <- eff %>% filter(sharpe==max(eff$sharpe)) 
  
  
  # Conferindo valores
  wm <- eff.optimal.point %>% select(var_modelo)
  
  C <-  as.matrix(data.matrix(dados_out_of_sample %>% select(var_modelo)),
                  ncol = ncol(dados_out_of_sample)) %*%
    t(as.matrix(round(wm,10)))
  
  # dados_out_of_sample$KERNEL <- cumsum(C)
  
  # Carteira
  auto_cov_out <- acf(C,lag=1)$acf[2]
  n_q <- q^(1/2)*(1 + (2*auto_cov_out)/(1- auto_cov_out)*(1 - (1-auto_cov_out^q)/(q*(1-auto_cov_out))))^(-1/2)
  sr <- (mean(C)-RF_mean)/sd(C)
  
  sr_month <- n_q*sr 
  sr_month_0 <- sr_month
  VIID_q <- (n_q^2)*(1+1/2*(sr^2))
  VGMM_add <- ifelse(q==1,0,q*(sr)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q <- VIID_q   + VGMM_add 
  
  sd_sr_q <- (VGMM_q/nrow(dados_out_of_sample))^(1/2)
  
  
  sor <- n_q*(mean(C)-RF_mean)/DownsideDeviation(C,MAR = RF_mean)
  sor_sp <- (mean(dados_out_of_sample$rs)-RF_mean)/DownsideDeviation(dados_out_of_sample$rs,MAR = RF_mean)
  
  
  VIID_q2 <- (n_q^2)*(1+1/2*(sor^2))
  VGMM_add2 <- ifelse(q==1,0,q*(sor)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q2 <- VIID_q2   + VGMM_add2 
  
  sd_sr_q2 <- (VGMM_q2/nrow(dados_out_of_sample))^(1/2)
  
  
  
  
  vec_parms_temp <- data.frame(data = nomes_arquivos[i],
                               cov_method = "POLY 3",
                               initial_in = min(dados_in_sample$DT),
                               end_in = max(dados_in_sample$DT),
                               initial_out = min(dados_out_of_sample$DT),
                               end_out = max(dados_out_of_sample$DT),
                               N_in = nrow(dados_in_sample),
                               N_out = nrow(dados_out_of_sample),
                               X = length(var_modelo),
                               X_pos_n = length(which(wm>0.00)),
                               X_pos_valid_n = length(which(wm>0.01)),
                               X_pos_valid_sum = paste0(round(sum(wm[(which(wm>0.01))])*100,2),"%"),
                               w_max = paste0(round(max(wm[(which(wm>0.01))])*100,2),"%"),
                               w_min = paste0(round(min(wm[(which(wm>0.01))])*100,2),"%"),
                               r_final = paste0(round(cumsum(C)[nrow(dados_out_of_sample)]*100,4),"%"),
                               Q = Q, 
                               lamda_max = lamda_max,
                               eigen_up = NA,
                               eigen_v_up = NA,
                               eigen_max = NA,
                               eigen_v_max = NA,
                               st_dev=sd(A),
                               down_dev=DownsideDeviation(C,MAR = RF_mean),
                               sharpe_ratio = sr_month ,
                               p_sharp_ratio_sp = 1-pnorm(sr_month-sr_sp_month,0,sd_sr_q),
                               p_sharp_ratio_0 = 1-pnorm(sr_month,0,sd_sr_q),
                               rmt_improvement = NA ,
                               sortino_ratio = sor ,
                               p_sort_ratio_sp = 1-pnorm(sor-sor_sp,0,sd_sr_q2),
                               # rmt_improvement_so = NA ,
                               p_sort_ratio_0 = 1-pnorm(sor,0,sd_sr_q2),
                               auto_cov_ar1 = auto_cov_out,
                               p_sharp_ratio_rmt = NA,
                               p_sort_ratio_rmt = NA)
  
  
  vec_parms <- rbind(vec_parms_temp,vec_parms)
  
  
  ### KRMT POLY 3 ----
  
  eigen0 <- eigen(K11)
  g1 <- eigen0$values[eigen0$values>=lamda_max]
  g2 <- eigen0$values[eigen0$values<lamda_max]
  eigen1<- diag(c(g1,rep(mean(g2),length(g2))))
  cor1 <- eigen0$vectors %*% eigen1 %*% t(eigen0$vectors)
  covariance_krmt <-  cor1
  for(ii in 1:ncol(cor1)){
    for(jj in 1:nrow(cor1)){
      covariance_krmt[ii,jj] <- cor1[ii,jj] * (covariance[ii,ii]*covariance[jj,jj])^.5 
    }
  }
  
  
  # Run the eff.frontier function based on no short and 50% alloc. restrictions
  eff <- eff.frontier(covariance_krmt,vec_mean, short="no", max.allocation=NULL,
                      risk.premium.up=1, risk.increment=.1)
  
  # Find the optimal portfolio
  eff.optimal.point <- eff %>% filter(sharpe==max(eff$sharpe)) 
  
  
  # Conferindo valores
  wm <- eff.optimal.point %>% select(var_modelo)
  
  D <-  as.matrix(data.matrix(dados_out_of_sample %>% select(var_modelo)),
                  ncol = ncol(dados_out_of_sample)) %*%
    t(as.matrix(round(wm,10)))
  
  # dados_out_of_sample$KRMT <- cumsum(D)
  
  
  # Carteira
  auto_cov_out <- acf(D,lag=1)$acf[2]
  n_q <- q^(1/2)*(1 + (2*auto_cov_out)/(1- auto_cov_out)*(1 - (1-auto_cov_out^q)/(q*(1-auto_cov_out))))^(-1/2)
  sr <- (mean(D)-RF_mean)/sd(D)
  
  sr_month <- n_q*sr 
  VIID_q <- (n_q^2)*(1+1/2*(sr^2))
  VGMM_q <- VIID_q + (q*sr)^2*(sum((1- (1:(q-1))/q)^2))
  
  sd_sr_q <- (VGMM_q/nrow(dados_out_of_sample))^(1/2)
  
  
  
  sor <- n_q*(mean(D)-RF_mean)/DownsideDeviation(D[,1],MAR = RF_mean)
  sor_sp <- (mean(dados_out_of_sample$rs)-RF_mean)/DownsideDeviation(dados_out_of_sample$rs,MAR = RF_mean)
  
  
  VIID_q2 <- (n_q^2)*(1+1/2*(sor^2))
  VGMM_add2 <- ifelse(q==1,0,q*(sor)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q2 <- VIID_q2   + VGMM_add2 
  
  sd_sr_q2 <- (VGMM_q2/nrow(dados_out_of_sample))^(1/2)
  
  
  
  vec_parms_temp <- data.frame(data = nomes_arquivos[i],
                               cov_method = "POLY 3 RMT",
                               initial_in = min(dados_in_sample$DT),
                               end_in = max(dados_in_sample$DT),
                               initial_out = min(dados_out_of_sample$DT),
                               end_out = max(dados_out_of_sample$DT),
                               N_in = nrow(dados_in_sample),
                               N_out = nrow(dados_out_of_sample),
                               X = length(var_modelo),
                               X_pos_n = length(which(wm>0.00)),
                               X_pos_valid_n = length(which(wm>0.01)),
                               X_pos_valid_sum = paste0(round(sum(wm[(which(wm>0.01))])*100,2),"%"),
                               w_max = paste0(round(max(wm[(which(wm>0.01))])*100,2),"%"),
                               w_min = paste0(round(min(wm[(which(wm>0.01))])*100,2),"%"),
                               r_final = paste0(round(cumsum(D)[nrow(dados_out_of_sample)]*100,4),"%"),
                               Q = Q, 
                               lamda_max = lamda_max,
                               eigen_up = paste0(length(g1)," (",round(length(g1)/(ncol(dados_in_sample) - 1)*100,2),"%)"),
                               eigen_v_up = paste0(round(sum(g1)/(sum(g1)+sum(g2))*100,2),"%"),
                               eigen_max = g1[1]/lamda_max,
                               eigen_v_max = paste0(round(sum(g1[1])/(sum(g1)+sum(g2))*100,2),"%"),
                               st_dev=sd(B),
                               down_dev=DownsideDeviation(D[,1],MAR = RF_mean),
                               sharpe_ratio = sr_month ,
                               p_sharp_ratio_sp = 1-pnorm(sr_month-sr_sp_month,0,sd_sr_q),
                               p_sharp_ratio_0 = 1-pnorm(sr_month,0,sd_sr_q),
                               rmt_improvement =((mean(D[,1]-C)+1)^(252)-1)*100,
                               sortino_ratio = sor ,
                               p_sort_ratio_sp = 1-pnorm(sor-sor_sp,0,sd_sr_q2),
                               p_sort_ratio_0 = 1-pnorm(sor,0,sd_sr_q2),
                               auto_cov_ar1 = auto_cov_out,
                               p_sharp_ratio_rmt =  1-pnorm(sr_month-sr_month_0,0,sd_sr_q),
                               p_sort_ratio_rmt =  1-pnorm(sor-sor_sp,0,sd_sr_q2))
  
  
  
  
  vec_parms <- rbind(vec_parms_temp,vec_parms)
  
  
  dados_out_of_sample$POLY3 <- cumsum(C)
  dados_out_of_sample$POLY3_RMT <- cumsum(D[,1])
  
  
  ### KERNEL POLY 4  ----
  
  
  # mpearson  <- cor(dados_in_sample %>% select(var_modelo))
  # covariance <- var(dados_in_sample %>% select(var_modelo))
  # vec_mean <- colMeans(dados_in_sample %>% select(var_modelo))
  # 
  # # poly <- polydot(degree=2,scale = sqrt(sum(diag(mpearson))*sum(diag(covariance))),offset = sqrt(mean(mpearson)*mean(covariance)))
  # poly <- polydot(degree=4,scale = 1/((sum(diag(covariance))/nrow(covariance))^2)*2)
  # K11 <- kernelMatrix(poly,covariance)
  # # K11 <- kernelMatrix(poly,mpearson)
  # # K11 <- K11-mean(K11)/((sd(K11)))
  # 
  # covariance_krmt <- K11
  # 
  # # vec_sd <- diag(K11)^.5
  # # for(ii in 1:ncol(K11)){
  # #   for(jj in 1:nrow(K11)){
  # #     covariance_krmt[ii,jj] <- K11[ii,jj] / (vec_sd[ii]*vec_sd[jj])
  # #     # covariance_krmt[ii,jj] <- K11[ii,jj] / (covariance[ii]*covariance[jj]) 
  # #   }
  # # }
  # 
  # 
  # # covariance_krmt_poly <-  K12_cor
  # for(ii in 1:ncol(K11)){
  #   for(jj in 1:nrow(K11)){
  #     covariance_krmt[ii,jj] <- covariance_krmt[ii,jj] * (covariance[ii,ii]*covariance[jj,jj])^.5
  #   }
  # }
  
  # covariance_krmt <- cov(covariance_krmt)
  
  
  covariance_krmt <- matrix(NA, nrow = ncol(dados_in_sample %>% select(var_modelo)), 
                            ncol = ncol(dados_in_sample %>% select(var_modelo)))
  
  dd <- 1
  for (i in 1:ncol(dados_in_sample %>% select(var_modelo))){
    for (j in i:ncol(dados_in_sample %>% select(var_modelo))){
      covariance_krmt[i,j] <- covariance_krmt[j,i] <- 
        ((dados_in_sample %>% select(var_modelo))[,i]%*%(dados_in_sample %>% select(var_modelo))[,j]+dd)^4
    }
  }
  
  colnames(covariance_krmt) <- rownames(covariance_krmt) <- colnames(dados_in_sample %>% select(var_modelo))
  
  
  
  # Run the eff.frontier function based on no short and 50% alloc. restrictions
  eff <- eff.frontier(covariance_krmt,vec_mean, short="no", max.allocation=NULL,
                      risk.premium.up=1, risk.increment=.1)
  
  # Find the optimal portfolio
  eff.optimal.point <- eff %>% filter(sharpe==max(eff$sharpe)) 
  
  
  # Conferindo valores
  wm <- eff.optimal.point %>% select(var_modelo)
  
  C <-  as.matrix(data.matrix(dados_out_of_sample %>% select(var_modelo)),
                  ncol = ncol(dados_out_of_sample)) %*%
    t(as.matrix(round(wm,10)))
  
  # dados_out_of_sample$KERNEL <- cumsum(C)
  
  # Carteira
  auto_cov_out <- acf(C,lag=1)$acf[2]
  n_q <- q^(1/2)*(1 + (2*auto_cov_out)/(1- auto_cov_out)*(1 - (1-auto_cov_out^q)/(q*(1-auto_cov_out))))^(-1/2)
  sr <- (mean(C)-RF_mean)/sd(C)
  
  sr_month <- n_q*sr 
  sr_month_0 <- sr_month
  VIID_q <- (n_q^2)*(1+1/2*(sr^2))
  VGMM_add <- ifelse(q==1,0,q*(sr)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q <- VIID_q   + VGMM_add 
  
  sd_sr_q <- (VGMM_q/nrow(dados_out_of_sample))^(1/2)
  
  
  sor <- n_q*(mean(C)-RF_mean)/DownsideDeviation(C,MAR = RF_mean)
  sor_sp <- (mean(dados_out_of_sample$rs)-RF_mean)/DownsideDeviation(dados_out_of_sample$rs,MAR = RF_mean)
  
  
  VIID_q2 <- (n_q^2)*(1+1/2*(sor^2))
  VGMM_add2 <- ifelse(q==1,0,q*(sor)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q2 <- VIID_q2   + VGMM_add2 
  
  sd_sr_q2 <- (VGMM_q2/nrow(dados_out_of_sample))^(1/2)
  
  
  
  
  vec_parms_temp <- data.frame(data = nomes_arquivos[i],
                               cov_method = "POLY 4",
                               initial_in = min(dados_in_sample$DT),
                               end_in = max(dados_in_sample$DT),
                               initial_out = min(dados_out_of_sample$DT),
                               end_out = max(dados_out_of_sample$DT),
                               N_in = nrow(dados_in_sample),
                               N_out = nrow(dados_out_of_sample),
                               X = length(var_modelo),
                               X_pos_n = length(which(wm>0.00)),
                               X_pos_valid_n = length(which(wm>0.01)),
                               X_pos_valid_sum = paste0(round(sum(wm[(which(wm>0.01))])*100,2),"%"),
                               w_max = paste0(round(max(wm[(which(wm>0.01))])*100,2),"%"),
                               w_min = paste0(round(min(wm[(which(wm>0.01))])*100,2),"%"),
                               r_final = paste0(round(cumsum(C)[nrow(dados_out_of_sample)]*100,4),"%"),
                               Q = Q, 
                               lamda_max = lamda_max,
                               eigen_up = NA,
                               eigen_v_up = NA,
                               eigen_max = NA,
                               eigen_v_max = NA,
                               st_dev=sd(A),
                               down_dev=DownsideDeviation(C,MAR = RF_mean),
                               sharpe_ratio = sr_month ,
                               p_sharp_ratio_sp = 1-pnorm(sr_month-sr_sp_month,0,sd_sr_q),
                               p_sharp_ratio_0 = 1-pnorm(sr_month,0,sd_sr_q),
                               rmt_improvement = NA ,
                               sortino_ratio = sor ,
                               p_sort_ratio_sp = 1-pnorm(sor-sor_sp,0,sd_sr_q2),
                               # rmt_improvement_so = NA ,
                               p_sort_ratio_0 = 1-pnorm(sor,0,sd_sr_q2),
                               auto_cov_ar1 = auto_cov_out,
                               p_sharp_ratio_rmt = NA,
                               p_sort_ratio_rmt = NA)
  
  
  vec_parms <- rbind(vec_parms_temp,vec_parms)
  
  
  ### KRMT POLY 4 ----
  
  eigen0 <- eigen(K11)
  g1 <- eigen0$values[eigen0$values>=lamda_max]
  g2 <- eigen0$values[eigen0$values<lamda_max]
  eigen1<- diag(c(g1,rep(mean(g2),length(g2))))
  cor1 <- eigen0$vectors %*% eigen1 %*% t(eigen0$vectors)
  covariance_krmt <-  cor1
  for(ii in 1:ncol(cor1)){
    for(jj in 1:nrow(cor1)){
      covariance_krmt[ii,jj] <- cor1[ii,jj] * (covariance[ii,ii]*covariance[jj,jj])^.5 
    }
  }
  
  
  # Run the eff.frontier function based on no short and 50% alloc. restrictions
  eff <- eff.frontier(covariance_krmt,vec_mean, short="no", max.allocation=NULL,
                      risk.premium.up=1, risk.increment=.1)
  
  # Find the optimal portfolio
  eff.optimal.point <- eff %>% filter(sharpe==max(eff$sharpe)) 
  
  
  # Conferindo valores
  wm <- eff.optimal.point %>% select(var_modelo)
  
  D <-  as.matrix(data.matrix(dados_out_of_sample %>% select(var_modelo)),
                  ncol = ncol(dados_out_of_sample)) %*%
    t(as.matrix(round(wm,10)))
  
  # dados_out_of_sample$KRMT <- cumsum(D)
  
  
  # Carteira
  auto_cov_out <- acf(D,lag=1)$acf[2]
  n_q <- q^(1/2)*(1 + (2*auto_cov_out)/(1- auto_cov_out)*(1 - (1-auto_cov_out^q)/(q*(1-auto_cov_out))))^(-1/2)
  sr <- (mean(D)-RF_mean)/sd(D)
  
  sr_month <- n_q*sr 
  VIID_q <- (n_q^2)*(1+1/2*(sr^2))
  VGMM_q <- VIID_q + (q*sr)^2*(sum((1- (1:(q-1))/q)^2))
  
  sd_sr_q <- (VGMM_q/nrow(dados_out_of_sample))^(1/2)
  
  
  
  sor <- n_q*(mean(D)-RF_mean)/DownsideDeviation(D[,1],MAR = RF_mean)
  sor_sp <- (mean(dados_out_of_sample$rs)-RF_mean)/DownsideDeviation(dados_out_of_sample$rs,MAR = RF_mean)
  
  
  VIID_q2 <- (n_q^2)*(1+1/2*(sor^2))
  VGMM_add2 <- ifelse(q==1,0,q*(sor)^2*(sum((1- (1:(q-1))/q)^2)))
  ##eq26 lo2002
  VGMM_q2 <- VIID_q2   + VGMM_add2 
  
  sd_sr_q2 <- (VGMM_q2/nrow(dados_out_of_sample))^(1/2)
  
  
  
  vec_parms_temp <- data.frame(data = nomes_arquivos[i],
                               cov_method = "POLY 4 RMT",
                               initial_in = min(dados_in_sample$DT),
                               end_in = max(dados_in_sample$DT),
                               initial_out = min(dados_out_of_sample$DT),
                               end_out = max(dados_out_of_sample$DT),
                               N_in = nrow(dados_in_sample),
                               N_out = nrow(dados_out_of_sample),
                               X = length(var_modelo),
                               X_pos_n = length(which(wm>0.00)),
                               X_pos_valid_n = length(which(wm>0.01)),
                               X_pos_valid_sum = paste0(round(sum(wm[(which(wm>0.01)),1])*100,2),"%"),
                               w_max = paste0(round(max(wm[(which(wm>0.01)),1])*100,2),"%"),
                               w_min = paste0(round(min(wm[(which(wm>0.01)),1])*100,2),"%"),
                               r_final = paste0(round(cumsum(D)[nrow(dados_out_of_sample)]*100,4),"%"),
                               Q = Q, 
                               lamda_max = lamda_max,
                               eigen_up = paste0(length(g1)," (",round(length(g1)/(ncol(dados_in_sample) - 1)*100,2),"%)"),
                               eigen_v_up = paste0(round(sum(g1)/(sum(g1)+sum(g2))*100,2),"%"),
                               eigen_max = g1[1]/lamda_max,
                               eigen_v_max = paste0(round(sum(g1[1])/(sum(g1)+sum(g2))*100,2),"%"),
                               st_dev=sd(B),
                               down_dev=DownsideDeviation(D[,1],MAR = RF_mean),
                               sharpe_ratio = sr_month ,
                               p_sharp_ratio_sp = 1-pnorm(sr_month-sr_sp_month,0,sd_sr_q),
                               p_sharp_ratio_0 = 1-pnorm(sr_month,0,sd_sr_q),
                               rmt_improvement =((mean(D[,1]-C)+1)^(252)-1)*100,
                               sortino_ratio = sor ,
                               p_sort_ratio_sp = 1-pnorm(sor-sor_sp,0,sd_sr_q2),
                               p_sort_ratio_0 = 1-pnorm(sor,0,sd_sr_q2),
                               auto_cov_ar1 = auto_cov_out,
                               p_sharp_ratio_rmt =  1-pnorm(sr_month-sr_month_0,0,sd_sr_q),
                               p_sort_ratio_rmt =  1-pnorm(sor-sor_sp,0,sd_sr_q2))
  
  
  
  
  vec_parms <- rbind(vec_parms_temp,vec_parms)
  
  
  dados_out_of_sample$POLY4 <- cumsum(C)
  dados_out_of_sample$POLY4_RMT <- cumsum(D[,1])
  
  
# }

### PLOTS ----
  ggplot() + 
    theme_classic() +
    # theme(legend.position = right) +
    geom_line(data = dados_out_of_sample, aes(x = DT, y = pearson, color = "black"),size=0.7) +
    geom_line(data = dados_out_of_sample, aes(x = DT, y = RMT, color = "blue"),size=0.7) +
    geom_line(data = dados_out_of_sample, aes(x = DT, y = MCD, color = "red1"), color = "red1",size=0.7) +
    geom_line(data = dados_out_of_sample, aes(x = DT, y = MCD_RMT, color = "red4"), color = "red4",size=0.7) +
    geom_line(data = dados_out_of_sample, aes(x = DT, y = WMCD, color = "gold"), color = "gold",size=0.7) +
    geom_line(data = dados_out_of_sample, aes(x = DT, y = WMCD_RMT, color = "orange"), color = "orange",size=0.7) +
    # geom_line(data = dados_out_of_sample, aes(x = DT, y = STAHEL_DONOHO), color = "seagreen",size=0.7) +
    # geom_line(data = dados_out_of_sample, aes(x = DT, y = STAHEL_DONOHO_RMTDONOHO), color = "springgreen",size=0.7) +
    geom_line(data = dados_out_of_sample, aes(x = DT, y = OGK, color = "skyblue1"), color = "skyblue1",size=0.7) +
    geom_line(data = dados_out_of_sample, aes(x = DT, y = OGK_RMT, color = "turquoise1"), color = "turquoise1",size=0.7) +
    geom_line(data = dados_out_of_sample, aes(x = DT, y = POLY2, color = "purple1"), color = "purple1",size=0.7) +
    geom_line(data = dados_out_of_sample, aes(x = DT, y = POLY2_RMT, color = "purple4"), color = "purple4",size=0.7) +
    geom_line(data = dados_out_of_sample, aes(x = DT, y = POLY3, color = "olivedrab1"), color = "olivedrab1",size=0.7) +
    geom_line(data = dados_out_of_sample, aes(x = DT, y = POLY3_RMT, color = "olivedrab4"), color = "olivedrab4",size=0.7) +
    geom_line(data = dados_out_of_sample, aes(x = DT, y = POLY4, color = "khaki1"), color = "khaki1",size=0.7) +
    geom_line(data = dados_out_of_sample, aes(x = DT, y = POLY4_RMT, color = "khaki4"), color = "khaki4",size=0.7) +
    geom_line(data = dados_out_of_sample, aes(x = DT, y = GAUSS, color = "violetred1"), color = "violetred1",size=0.7) +
    geom_line(data = dados_out_of_sample, aes(x = DT, y = GAUSS_RMT, color = "violetred4"), color = "violetred4",size=0.7) +
    scale_y_continuous(labels = scales::percent) +
    xlab('Time') +
    ylab('Cumulative Return') +
    scale_colour_manual(name = "Legend",
                        values = c("black"="black",
                                   "blue"="blue","red1"="red1","red4"="red4","gold"="gold",
                                                  "orange"="orange","skyblue1"="skyblue1","turquoise1"="turquoise1",
                                                  "purple1"="purple1","purple4"="purple4","olivedrab1"="olivedrab1",
                                                  "olivedrab4"="olivedrab4","khaki1"="khaki1","khaki4"="khaki4",
                                                  "violetred1"="violetred1","violetred4"="violetred4"),
                        labels = c("Pearson","Pearson filtered","MCD","MCD filtered",
                                   "WMCD","WMCD filtered","OGK","OGK filtered",
                                   "POLY2","POLY2 filtered","POLY3","POLY3 filtered",
                                   "POLY4","POLY4 filtered","GAUSS","GAUSS filtered")) + 
    theme(legend.title=element_text(size=10), 
          legend.text=element_text(size=2))
  
  ggsave(width = 10,height = 8,paste0("program/results/GLOBAL_",i,"_CR_",nomes_arquivos[i],"_frequency.pdf"), device = "pdf")
  
  ggplot() + 
    theme_classic() +
    geom_line(data = dados_out_of_sample, aes(x = DT, y = RMT-pearson), color = "black") +
    geom_line(data = dados_out_of_sample, aes(x = DT, y = MCD_RMT-MCD), color = "red1") +
    geom_line(data = dados_out_of_sample, aes(x = DT, y = WMCD_RMT-WMCD), color = "gold") +
    # geom_line(data = dados_out_of_sample, aes(x = DT, y = STAHEL_DONOHO_RMTDONOHO), color = "springgreen") +
    # geom_line(data = dados_out_of_sample, aes(x = DT, y = STAHEL_DONOHO), color = "seagreen") +
    geom_line(data = dados_out_of_sample, aes(x = DT, y = OGK_RMT-OGK), color = "skyblue1") +
    scale_y_continuous(labels = scales::percent) +
    xlab('Time') +
    ylab('Cumulative Return') +
    geom_line(data = dados_out_of_sample, aes(x = DT, y = POLY2_RMT-POLY2), color = "purple1") +
    geom_line(data = dados_out_of_sample, aes(x = DT, y = POLY3_RMT-POLY3), color = "olivedrab1") +
    geom_line(data = dados_out_of_sample, aes(x = DT, y = POLY4_RMT-POLY4), color = "khaki1") + 
    geom_line(data = dados_out_of_sample, aes(x = DT, y = GAUSS_RMT-GAUSS), color = "violetred1") +
    scale_colour_manual(name = "Legend",values = c("black","red1","gold","skyblue1","purple1","olivedrab1","khaki1","violetred1"),
                        labels = c("Pearson","MCD","WMCD","OGK","POLY2","POLY3","POLY4","GAUSS"))
  
  
  ggsave(width = 10,height = 8,paste0("program/results/GLOBAL_",i,"_CRDIFF_",nomes_arquivos[i],"_frequency.pdf"), device = "pdf")
  
  
}

saveRDS(vec_parms,"KPCA_rottinas_kernel.rds")



saveRDS(vec_parms,"program/results/GLOBAL_KRMT_gaussian.rds")

vec_parms <- as.data.table(vec_parms)
vec_parms <- vec_parms[order(data,cov_method)]

vec_parms <- vec_parms %>%
  mutate(cov_method = ifelse(cov_method==1,"Pearson",
                             ifelse(cov_method==2,"RMT",
                                    ifelse(cov_method==3,"Kernel","KRMT"))))

# vec_parms <- vec_parms %>%
#   mutate(cov_method = ifelse(cov_method=="Pearson",1,
#                              ifelse(cov_method=="RMT",2,
#                                     ifelse(cov_method=="Kernel",3,4))))

nomes_vars<- names(vec_parms)
aprint <- cbind(vec_parms %>% select(nomes_vars[c(3,2)]),
                vec_parms[,"X"]-vec_parms[,"X_pos_n"],
                vec_parms %>% select(nomes_vars[c(12:15,19:25)]))

print(xtable(aprint),include.rownames = F)



print(xtable(cbind(vec_parms %>% select(nomes_vars[c(3,2)]),
                   vec_parms %>% select(nomes_vars[c(12:14,16,19:22)]))),include.rownames = F)

#################################################
# names(vec_parms)


# str_sub(as.character(vec_parms$r_final[1]),1,-2)
# round(vec_parms$st_dev[1]*100,6)
# str_sub(as.character(vec_parms$st_dev[1]),1,-2)

# xtable(vec_parms[vec_parms$data=="CAC_FRANCA",quallos],digits = 6)

quallos <- colnames(vec_parms)[c(1,2,15,22,23,24,28,18,19,20,21)]
ppp <- vec_parms %>% mutate(r_final=str_sub(as.character(r_final),1,-2),st_dev=round(st_dev*100,6),down_dev=round(down_dev*100,6))
xtable(ppp[ppp$data=="CAC_FRANCA",quallos],digits = 4)
xtable(ppp[ppp$data=="DAX_ALEMANHA",quallos],digits = 4)
xtable(ppp[ppp$data=="FTSE_REINO_UNIDO",quallos],digits = 4)
xtable(ppp[ppp$data=="IBOV",quallos],digits = 4)
xtable(ppp[ppp$data=="NASDAQ",quallos],digits = 4)
xtable(ppp[ppp$data=="NIKKEI_JAPAO",quallos],digits = 4)
xtable(ppp[ppp$data=="SHCOMP_CHINA",quallos],digits = 4)


##  1 - PEARSON
##  2 - RMT
##  3 - KERNEL
##  4 - KRMT

### ROBUST 3 ----

library(robust)

# covRob(dd,estim = 'mcd')
# covRob(dd,estim = 'weighted')
# covRob(dd,estim = 'donostah')
# covRob(dd,estim = 'pairwiseGK')
covariance <- covRob(dados_in_sample %>% select(var_modelo),estim = 'donostah')$cov
mpearson  <- cov2cor(covRob(dados_in_sample %>% select(var_modelo),estim = 'donostah')$cov)
vec_mean <- colMeans(dados_in_sample %>% select(var_modelo))



# Run the eff.frontier function based on no short and 50% alloc. restrictions
eff <- eff.frontier(covariance,vec_mean, short="no", max.allocation=NULL,
                    risk.premium.up=1, risk.increment=.1)

# Find the optimal portfolio
eff.optimal.point <- eff %>% filter(sharpe==max(eff$sharpe)) 

# # Conferindo valores
wm <- eff.optimal.point %>% select(var_modelo)
# sum(round(wm,10))


A <-  as.matrix(data.matrix(dados_out_of_sample %>% select(var_modelo)),
                ncol = ncol(dados_out_of_sample)) %*%  t(as.matrix(round(wm,10)))


dados_out_of_sample <- dados_out_of_sample %>%
  mutate(rs = cumsum(r))



## Não vamos comparar shape ratio de diferentes frequências
q <- 1


### Risk Free
RF_mean <- mean(dados_out_of_sample$RF,na.rm=T)


# Sp500
auto_cov_sp <- acf(dados_out_of_sample$r,lag=1)$acf[2]
n_q_sp <- q^(1/2)*(1 + (2*auto_cov_sp)/(1- auto_cov_sp)*(1 - (1-auto_cov_sp^q)/(q*(1-auto_cov_sp))))^(-1/2)

sr_sp <- (mean(dados_out_of_sample$rs)-RF_mean)/sd(dados_out_of_sample$rs)
sr_sp_month <-  n_q_sp*sr_sp 

# Carteira
auto_cov_out <- acf(A,lag=1)$acf[2]
n_q <- q^(1/2)*(1 + (2*auto_cov_out)/(1- auto_cov_out)*(1 - (1-auto_cov_out^q)/(q*(1-auto_cov_out))))^(-1/2)
sr <- (mean(A)-RF_mean)/sd(A)

sr_month <- n_q*sr 

sr_month_0 <- sr_month 
VIID_q <- (n_q^2)*(1+1/2*(sr^2))
VGMM_add <- ifelse(q==1,0,q*(sr)^2*(sum((1- (1:(q-1))/q)^2)))
##eq26 lo2002
VGMM_q <- VIID_q   + VGMM_add 

sd_sr_q <- (VGMM_q/nrow(dados_out_of_sample))^(1/2)



sor <- n_q*(mean(A)-RF_mean)/DownsideDeviation(A,MAR = RF_mean)
sor_sp <- (mean(dados_out_of_sample$rs)-RF_mean)/DownsideDeviation(dados_out_of_sample$rs,MAR = RF_mean)


VIID_q2 <- (n_q^2)*(1+1/2*(sor^2))
VGMM_add2 <- ifelse(q==1,0,q*(sor)^2*(sum((1- (1:(q-1))/q)^2)))
##eq26 lo2002
VGMM_q2 <- VIID_q2   + VGMM_add2 

sd_sr_q2 <- (VGMM_q2/nrow(dados_out_of_sample))^(1/2)



vec_parms_temp <- data.frame(data = nomes_arquivos[i],
                             cov_method = "STAHEL-DONOHO",
                             initial_in = min(dados_in_sample$DT),
                             end_in = max(dados_in_sample$DT),
                             initial_out = min(dados_out_of_sample$DT),
                             end_out = max(dados_out_of_sample$DT),
                             N_in = nrow(dados_in_sample),
                             N_out = nrow(dados_out_of_sample),
                             X = length(var_modelo),
                             X_pos_n = length(which(wm>0.00)),
                             X_pos_valid_n = length(which(wm>0.01)),
                             X_pos_valid_sum = paste0(round(sum(wm[(which(wm>0.01))])*100,2),"%"),
                             w_max = paste0(round(max(wm[(which(wm>0.01))])*100,2),"%"),
                             w_min = paste0(round(min(wm[(which(wm>0.01))])*100,2),"%"),
                             r_final = paste0(round(cumsum(A)[nrow(dados_out_of_sample)]*100,4),"%"),
                             Q = NA, 
                             lamda_max = NA,
                             eigen_up = NA,
                             eigen_v_up = NA,
                             eigen_max = NA,
                             eigen_v_max = NA,
                             st_dev=sd(A),
                             down_dev=DownsideDeviation(A,MAR = RF_mean),
                             sharpe_ratio = sr_month ,
                             p_sharp_ratio_sp = 1-pnorm(sr_month-sr_sp_month,0,sd_sr_q),
                             p_sharp_ratio_0 = 1-pnorm(sr_month,0,sd_sr_q),
                             rmt_improvement = NA ,
                             sortino_ratio = sor ,
                             p_sort_ratio_sp = 1-pnorm(sor-sor_sp,0,sd_sr_q2),
                             # rmt_improvement_so = NA ,
                             p_sort_ratio_0 = 1-pnorm(sor,0,sd_sr_q2),
                             auto_cov_ar1 = auto_cov_out,
                             p_sharp_ratio_rmt = NA,
                             p_sort_ratio_rmt = NA)



vec_parms <- rbind(vec_parms,vec_parms_temp)

### RMT ROBUST 3 ----
## ref: J. Bouchaud, M. Potters, Theory of Financial Risks—From Statistical Physics to Risk Management, Cambridge University Press, UK, 2000.
Q <- nrow(dados_in_sample)/(length(var_modelo))
## falta entender sigma!!!
eigen0 <- eigen(mpearson)
# eigen0 <- eigen(covariance)
lamda_max <- (1 + 1/Q + (1/Q)^.5)


g1 <- eigen0$values[eigen0$values>=lamda_max]
g2 <- eigen0$values[eigen0$values<lamda_max]
eigen1<- diag(c(g1,rep(mean(eigen0$values),length(g2))))
cor1 <- eigen0$vectors %*% eigen1 %*% t(eigen0$vectors)
covariance_rmt <-  cor1
for(ii in 1:ncol(cor1)){
  for(jj in 1:nrow(cor1)){
    covariance_rmt[ii,jj] <- cor1[ii,jj] * (covariance[ii,ii]*covariance[jj,jj])^.5 
  }
}

# Run the eff.frontier function based on no short and 50% alloc. restrictions
eff <- eff.frontier(covariance_rmt,vec_mean, short="no", max.allocation=NULL,
                    risk.premium.up=1, risk.increment=.1)

# Find the optimal portfolio
eff.optimal.point <- eff %>% filter(sharpe==max(eff$sharpe)) 

# Conferindo valores
wm <- eff.optimal.point %>% select(var_modelo)

B <-  as.matrix(data.matrix(dados_out_of_sample %>% select(var_modelo)),
                ncol = ncol(dados_out_of_sample)) %*%
  t(as.matrix(round(wm,10)))



# Carteira
auto_cov_out <- acf(B,lag=1)$acf[2]
n_q <- q^(1/2)*(1 + (2*auto_cov_out)/(1- auto_cov_out)*(1 - (1-auto_cov_out^q)/(q*(1-auto_cov_out))))^(-1/2)
sr <- (mean(B)-RF_mean)/sd(B)

sr_month <- n_q*sr 
VIID_q <- (n_q^2)*(1+1/2*(sr^2))
VGMM_add <- ifelse(q==1,0,q*(sr)^2*(sum((1- (1:(q-1))/q)^2)))
##eq26 lo2002
VGMM_q <- VIID_q   + VGMM_add 

sd_sr_q <- (VGMM_q/nrow(dados_out_of_sample))^(1/2)



sor <- n_q*(mean(B)-RF_mean)/DownsideDeviation(B,MAR = RF_mean)
sor_sp <- (mean(dados_out_of_sample$rs)-RF_mean)/DownsideDeviation(dados_out_of_sample$rs,MAR = RF_mean)


VIID_q2 <- (n_q^2)*(1+1/2*(sor^2))
VGMM_add2 <- ifelse(q==1,0,q*(sor)^2*(sum((1- (1:(q-1))/q)^2)))
##eq26 lo2002
VGMM_q2 <- VIID_q2   + VGMM_add2 

sd_sr_q2 <- (VGMM_q2/nrow(dados_out_of_sample))^(1/2)



vec_parms_temp <- data.frame(data = nomes_arquivos[i],
                             cov_method = "STAHEL-DONOHO RMT",
                             initial_in = min(dados_in_sample$DT),
                             end_in = max(dados_in_sample$DT),
                             initial_out = min(dados_out_of_sample$DT),
                             end_out = max(dados_out_of_sample$DT),
                             N_in = nrow(dados_in_sample),
                             N_out = nrow(dados_out_of_sample),
                             X = length(var_modelo),
                             X_pos_n = length(which(wm>0.00)),
                             X_pos_valid_n = length(which(wm>0.01)),
                             X_pos_valid_sum = paste0(round(sum(wm[(which(wm>0.01))])*100,2),"%"),
                             w_max = paste0(round(max(wm[(which(wm>0.01))])*100,2),"%"),
                             w_min = paste0(round(min(wm[(which(wm>0.01))])*100,2),"%"),
                             r_final = paste0(round(cumsum(B)[nrow(dados_out_of_sample)]*100,4),"%"),
                             Q = Q, 
                             lamda_max = lamda_max,
                             eigen_up = paste0(length(g1)," (",round(length(g1)/(ncol(dados_in_sample) - 1)*100,2),"%)"),
                             eigen_v_up = paste0(round(sum(g1)/(sum(g1)+sum(g2))*100,2),"%"),
                             eigen_max = g1[1]/lamda_max,
                             eigen_v_max = paste0(round(sum(g1[1])/(sum(g1)+sum(g2))*100,2),"%"),
                             st_dev=sd(B),
                             down_dev=DownsideDeviation(B,MAR = RF_mean),
                             sharpe_ratio = sr_month ,
                             p_sharp_ratio_sp = 1-pnorm(sr_month-sr_sp_month,0,sd_sr_q),
                             p_sharp_ratio_0 = 1-pnorm(sr_month,0,sd_sr_q),
                             rmt_improvement =((mean(B-A)+1)^(252)-1)*100,
                             sortino_ratio = sor ,
                             p_sort_ratio_sp = 1-pnorm(sor-sor_sp,0,sd_sr_q2),
                             p_sort_ratio_0 = 1-pnorm(sor,0,sd_sr_q2),
                             auto_cov_ar1 = auto_cov_out,
                             p_sharp_ratio_rmt =  1-pnorm(sr_month-sr_month_0,0,sd_sr_q),
                             p_sort_ratio_rmt =  1-pnorm(sor-sor_sp,0,sd_sr_q2))


vec_parms <- rbind(vec_parms,vec_parms_temp)


## saving results
dados_out_of_sample$STAHEL_DONOHO <- cumsum(A)
dados_out_of_sample$STAHEL_DONOHO_RMT <- cumsum(B)










### KERNEL poly MAVIKSSON  ----


poly1 <- polydot(degree=1)
K12 <- kernelMatrix(poly1,covariance)



K12_cor <- K12

vec_sd <- diag(K12)^.5
for(ii in 1:ncol(K12)){
  for(jj in 1:nrow(K12)){
    K12_cor[ii,jj] <- K12[ii,jj] / (vec_sd[ii]*vec_sd[jj]) 
  }
}


covariance_krmt_poly <-  K12_cor
for(ii in 1:ncol(cor1)){
  for(jj in 1:nrow(cor1)){
    covariance_krmt_poly[ii,jj] <- K12_cor[ii,jj] * (covariance[ii,ii]*covariance[jj,jj])^.5 
  }
}



# Run the eff.frontier function based on no short and 50% alloc. restrictions
eff <- eff.frontier(covariance_krmt_poly,vec_mean, short="no", max.allocation=NULL,
                    risk.premium.up=1, risk.increment=.1)

# Find the optimal portfolio
eff.optimal.point <- eff %>% filter(sharpe==max(eff$sharpe)) 


# Conferindo valores
wm <- eff.optimal.point %>% select(var_modelo)

E <-  as.matrix(data.matrix(dados_out_of_sample %>% select(var_modelo)),
                ncol = ncol(dados_out_of_sample)) %*%
  t(as.matrix(round(wm,10)))

dados_out_of_sample$KERNEL_POLY <- cumsum(E)

# Carteira
auto_cov_out <- acf(E,lag=1)$acf[2]
n_q <- q^(1/2)*(1 + (2*auto_cov_out)/(1- auto_cov_out)*(1 - (1-auto_cov_out^q)/(q*(1-auto_cov_out))))^(-1/2)
sr <- (mean(E)-RF_mean)/sd(E)

sr_month <- n_q*sr 
sr_month_0 <- sr_month
VIID_q <- (n_q^2)*(1+1/2*(sr^2))
VGMM_add <- ifelse(q==1,0,q*(sr)^2*(sum((1- (1:(q-1))/q)^2)))
##eq26 lo2002
VGMM_q <- VIID_q   + VGMM_add 

sd_sr_q <- (VGMM_q/nrow(dados_out_of_sample))^(1/2)



vec_parms_temp <- data.frame(data = nomes_arquivos[i],
                             cov_method = "KERNEL POLY",
                             initial_in = min(dados_in_sample$DT),
                             end_in = max(dados_in_sample$DT),
                             initial_out = min(dados_out_of_sample$DT),
                             end_out = max(dados_out_of_sample$DT),
                             N_in = nrow(dados_in_sample),
                             N_out = nrow(dados_out_of_sample),
                             X = length(var_modelo),
                             X_pos_n = length(which(wm>0.00)),
                             X_pos_valid_n = length(which(wm>0.01)),
                             X_pos_valid_sum = paste0(round(sum(wm[(which(wm>0.01))])*100,2),"%"),
                             w_max = paste0(round(max(wm[(which(wm>0.01))])*100,2),"%"),
                             w_min = paste0(round(min(wm[(which(wm>0.01))])*100,2),"%"),
                             r_final = paste0(round(cumsum(E)[nrow(dados_out_of_sample)]*100,2),"%"),
                             Q = Q, 
                             lamda_max = lamda_max,
                             eigen_up = NA,
                             eigen_v_up = NA,
                             eigen_max = NA,
                             eigen_v_max = NA,
                             sharpe_ratio = sr_month ,
                             p_sharp_ratio_sp = 1-pnorm(sr_month-sr_sp_month,0,sd_sr_q),
                             rmt_improvement =NA,
                             p_sharp_ratio_0 = 1-pnorm(sr_month,0,sd_sr_q),
                             auto_cov_ar1 = auto_cov_out,
                             p_sharp_ratio_rmt =  NA)


# vec_parms <- rbind(vec_parms_temp,vec_parms)


###KRMT POLY




eigen0 <- eigen(K12_cor)
g1 <- eigen0$values[eigen0$values>=lamda_max]
g2 <- eigen0$values[eigen0$values<lamda_max]
eigen1<- diag(c(g1,rep(mean(g2),length(g2))))
cor1 <- eigen0$vectors %*% eigen1 %*% t(eigen0$vectors)
covariance_krmt_poly <-  cor1
for(ii in 1:ncol(cor1)){
  for(jj in 1:nrow(cor1)){
    covariance_krmt_poly[ii,jj] <- cor1[ii,jj] * (covariance[ii,ii]*covariance[jj,jj])^.5 
  }
}


# Run the eff.frontier function based on no short and 50% alloc. restrictions
eff <- eff.frontier(covariance_krmt_poly,vec_mean, short="no", max.allocation=NULL,
                    risk.premium.up=1, risk.increment=.1)

# Find the optimal portfolio
eff.optimal.point <- eff %>% filter(sharpe==max(eff$sharpe)) 


# Conferindo valores
wm <- eff.optimal.point %>% select(var_modelo)

FF <-  as.matrix(data.matrix(dados_out_of_sample %>% select(var_modelo)),
                 ncol = ncol(dados_out_of_sample)) %*%
  t(as.matrix(round(wm,10)))

dados_out_of_sample$KRMT_POLY <- cumsum(FF)


# Carteira
auto_cov_out <- acf(FF,lag=1)$acf[2]
n_q <- q^(1/2)*(1 + (2*auto_cov_out)/(1- auto_cov_out)*(1 - (1-auto_cov_out^q)/(q*(1-auto_cov_out))))^(-1/2)
sr <- (mean(FF)-RF_mean)/sd(FF)

sr_month <- n_q*sr 
VIID_q <- (n_q^2)*(1+1/2*(sr^2))
VGMM_add <- ifelse(q==1,0,q*(sr)^2*(sum((1- (1:(q-1))/q)^2)))
##eq26 lo2002
VGMM_q <- VIID_q   + VGMM_add 

sd_sr_q <- (VGMM_q/nrow(dados_out_of_sample))^(1/2)





vec_parms_temp <- data.frame(data = nomes_arquivos[i],
                             cov_method = "KRMT POLY",
                             initial_in = min(dados_in_sample$DT),
                             end_in = max(dados_in_sample$DT),
                             initial_out = min(dados_out_of_sample$DT),
                             end_out = max(dados_out_of_sample$DT),
                             N_in = nrow(dados_in_sample),
                             N_out = nrow(dados_out_of_sample),
                             X = length(var_modelo),
                             X_pos_n = length(which(wm>0.00)),
                             X_pos_valid_n = length(which(wm>0.01)),
                             X_pos_valid_sum = paste0(round(sum(wm[(which(wm>0.01))])*100,2),"%"),
                             w_max = paste0(round(max(wm[(which(wm>0.01))])*100,2),"%"),
                             w_min = paste0(round(min(wm[(which(wm>0.01))])*100,2),"%"),
                             r_final = paste0(round(cumsum(FF)[nrow(dados_out_of_sample)]*100,2),"%"),
                             Q = Q, 
                             lamda_max = lamda_max,
                             eigen_up = paste0(length(g1)," (",round(length(g1)/(ncol(dados_in_sample) - 1)*100,2),"%)"),
                             eigen_v_up = paste0(round(sum(g1)/(sum(g1)+sum(g2))*100,2),"%"),
                             eigen_max = g1[1]/lamda_max,
                             eigen_v_max = paste0(round(sum(g1[1])/(sum(g1)+sum(g2))*100,2),"%"),
                             sharpe_ratio = sr_month ,
                             p_sharp_ratio_sp = 1-pnorm(sr_month-sr_sp_month,0,sd_sr_q),
                             rmt_improvement =(((mean(FF-E)+1)^(252)-1)*100),
                             p_sharp_ratio_0 = 1-pnorm(sr_month,0,sd_sr_q),
                             auto_cov_ar1 = auto_cov_out,
                             p_sharp_ratio_rmt =  1-pnorm(sr_month-sr_month_0,0,sd_sr_q))


# vec_parms <- rbind(vec_parms_temp,vec_parms)

