# https://www.treasury.gov/resource-center/data-chart-center/interest-rates/Pages/TextView.aspx?data=longtermrateYear&year=2000

library(XML)
library(RCurl)
library(rlist)
library(dplyr)
library(tidyverse)
# base <- NULL
i <- 2003
for(i in 2000:2018){
  link <- paste0("https://www.treasury.gov/resource-center/data-chart-center/interest-rates/Pages/TextView.aspx?data=longtermrateYear&year=",
                 i)
  
  if(file.exists(paste("rf_",i,".rds"))){
    next
  }else{
    print(i)  
  theurl <- getURL(link,.opts = list(ssl.verifypeer = FALSE) )
  tables <- readHTMLTable(theurl)
  
  base_temp <- tables[[1]]
  base_temp <- base_temp %>%
    filter(grepl("^\\d\\d.\\d\\d.\\d\\d",V1))
  
  base_temp$ano <- i
  base_temp <- base_temp %>% select(V1,V2,ano)
  base_temp$indice <- rep(c("treasury_10","treasury_20"),
      each = nrow(base_temp)/2)
  base_temp <- base_temp %>%
    spread(key = indice,
           value = V2)
  # base <- bind_rows(base,base_temp)
  
  saveRDS(base_temp,paste("rf_",i,".rds"))

  
}
  
}

## jutando todos os arquivos 

arquivos <- list.files(path = getwd(),
           "rf_",
           full.names = T)

library(data.table)
base <- rbindlist(lapply(arquivos, readRDS))
base <- base %>%
  mutate(DT = as.Date(V1,format = "%m/%d/%y"),
         treasury_10 = as.numeric(as.character(treasury_10)),
         treasury_20 = as.numeric(as.character(treasury_20))) %>%
  select(-V1)

saveRDS(base,paste("rf.rds"))
