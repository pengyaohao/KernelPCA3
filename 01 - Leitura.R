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
arquivo1 <- "data/BOLSAS2.xlsx"
planilhas1 <- readxl::excel_sheets(arquivo1)

arquivo2 <- "data/BOLSAS2 - Copy.xlsx"
planilhas2 <- readxl::excel_sheets(arquivo2)


arquivo <- bind_rows(data.frame(planilha = planilhas1,arquivo = arquivo1),
                     data.frame(planilha = planilhas2,arquivo = arquivo2))

arquivo <- arquivo %>%
  filter( !grepl("sheet|label|indices",planilha,ignore.case = T) ) %>%
  mutate(indice = ifelse(planilha=="Brasil", "IBOV",
                         ifelse(planilha=="China_shanghai", "SZSMEC",
                                ifelse(planilha=="Franca", "CAC",
                                       ifelse(planilha=="Alemanha", "SPX",
                                              ifelse(planilha=="Japao", "NKY",
                                                     ifelse(planilha=="Holanda", "AEX",
                                                            ifelse(planilha=="UK", "UKX","NDX"))))))))



base_indice <- read_excel("data/BOLSAS2 - Copy.xlsx",
                          sheet = "INDICES")
names(base_indice)[1] <- "DT"

detach(package:plyr)
j <- 2
for (j in (1:nrow(arquivo))[-c(1,3,5)] ) {
  i <- arquivo$planilha[j]
  base <- read_excel(arquivo$arquivo[j],
                     sheet = i)
  
  names(base)[1] <- "DT"
  names(base) <- gsub("(.*)\\s.*\\sEquity","\\1",names(base))
  
  
  base_indice_temp <- base_indice
  
  names(base_indice_temp)[grep(arquivo$indice[j],
                               names(base_indice_temp))] <- "r"
  
  base_indice_temp <- base_indice_temp %>% select(DT,r)
  
  base <- base %>% 
    left_join(base_indice_temp,
              by = "DT")
  
  
  nomes_colunas_transform <- setdiff(names(base),"DT")
  setDT(base)[, (nomes_colunas_transform):= lapply(.SD,  
                                                   function(x) as.numeric(gsub(",",".",x))),
              .SDcols=nomes_colunas_transform]
  
  
  
  
  ### adicionando risk free 
  base_rf <- readRDS("data/rf.rds")
  base_rf <- base_rf %>% 
    mutate(DT = as.POSIXct(DT, format="%Y-%m-%d") ,
           RFe  = treasury_10) %>%
    select(DT,RFe)
  
  
  base <- base %>%
    left_join(base_rf,
              by = "DT") %>% as.data.table()
  
  base[,RF:=ifelse(is.na(RFe),
                   lag(RFe),
                   RFe)]
  
  base <- base %>% 
    select(-RFe)
  ### corrigindo datas inexistentes na base americana
  
  
  base <- base  %>% 
    gather(key = "var",value = "valor",-DT) %>% data.table()
  
  base <- base[order(var,DT)]
  base[,ra:=log(valor)-lag(log(valor)),by="var"]
  base[,m:=mean(ra,na.rm=T),by = "var"]
  base[,v:=sd(ra,na.rm=T),by = "var"]
  base[,value_z:=(ra-m)/v]
  
  
  base <-base %>% filter(!is.na(ra)) %>%
    select(-value_z,-valor,-m,-v) %>%
    spread(key = var,value = ra) %>% 
    filter(!is.na(DT))     %>%
    as.data.table() 
  
  
  dados_faltantes <- t(base[, lapply(.SD, function(x) sum(is.na(x)))])
  dados_faltantes <- cbind(data.table(dados_faltantes),
                           rownames(dados_faltantes)) %>% data.table()
  dados_faltantes <- dados_faltantes[order(V1)]
  
  ## menos de 10% de missing
  vars_sem_na <- dados_faltantes$V2[dados_faltantes$V1<(nrow(base)/10)]
  length(vars_sem_na)
  # base <- base %>%
  #   select(vars_sem_na)
  
  saveRDS(base,paste0("data//global//",gsub("\\s","_NOVO_",i),".rds"))
  rm(base)
}

base$RF[1]
