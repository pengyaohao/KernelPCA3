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


rm(list = ls());gc()
arquivo <- "data/bolsas_mundo.xlsx"
planilhas <- readxl::excel_sheets(arquivo)

detach(package:plyr)

i <- planilhas[9]
for (i in planilhas[-4]) {
  base <- read_excel(arquivo,
                     sheet = i,
                     skip = 1)
  names(base)[1] <- "DT"
  names(base) <- gsub("(.*)\\sUW\\sEquity","\\1",names(base))
  nomes_colunas_transform <- setdiff(names(base),"DT")
  setDT(base)[, (nomes_colunas_transform):= lapply(.SD,  function(x) as.numeric(gsub(",",".",x))),
               .SDcols=nomes_colunas_transform]
  
  
  base <- base  %>% 
    gather(key = "var",value = "valor",-DT) %>% data.table()
  
  
  base <- base[order(var,DT)]
  base[,r:=log(valor)-lag(log(valor)),by="var"]
  base[,m:=mean(r,na.rm=T),by = "var"]
  base[,v:=sd(r,na.rm=T),by = "var"]
  base[,value_z:=(r-m)/v]
  
  
  ## pulo do gato roubado ----
  base_RF <- base %>% select(DT,r) %>% dplyr::group_by(DT) %>%
    summarise(r = quantile(r,
                                 p = .45 ,
                                 na.rm=T)) %>%
    mutate(var = "RF")
  
  base_R <- base %>% select(DT,r) %>% dplyr::group_by(DT) %>%
    summarise(r = quantile(r,
                           p = .5 ,
                           na.rm=T)) %>%
    mutate(var = "r")
  
  base <- bind_rows(base,base_RF,base_R)
  
  base <-base %>% filter(!is.na(r)) %>%
    select(-value_z,-valor,-m,-v) %>%
    spread(key = var,value = r) %>% 
    filter(!is.na(DT))     %>%
    as.data.table() 
  
  
  dados_faltantes <- t(base[, lapply(.SD, function(x) sum(is.na(x)))])
  dados_faltantes <- cbind(data.table(dados_faltantes),
                           rownames(dados_faltantes)) %>% data.table()
  dados_faltantes <- dados_faltantes[order(V1)]

  vars_sem_na <- dados_faltantes$V2[dados_faltantes$V1==0]
  length(vars_sem_na)
  base <- base %>%
    select(vars_sem_na)
  
  saveRDS(base,paste0("data//global//",gsub("\\s","_",i),".rds"))

}


