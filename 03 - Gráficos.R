## carregando bases
library(data.table)
library(dplyr)
library(tidyverse)

arquivos <- list.files(pattern = "retornos_.*rds",
                       full.names = T)

i <- 1
for(i in 1:length(arquivos)){
base <- readRDS(arquivos[i])

base <- base %>% gather(key = "variable" ,
                        value = "valor",
                        -DT)
# gráfico KPCA ----
base$type_line <- ifelse(grepl("_RMT",base$variable),"dotted","solid")
base$variable2 <- gsub('_RMT',"",base$variable)

# base <- base %>% select(variable2,DT,valor,type_line) %>% 
#   spread(key = "type_line" ,
#                         value = "valor")


g <- ggplot(data = base, 
       aes(x = DT, y = valor, color = variable2))+
  geom_line(aes(linetype = type_line),
            size=.9) + 
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  theme(legend.position = "bottom") +
  xlab('Time') +
  ylab('Cumulative Return') +
  scale_colour_manual(name = "",
                      values = c("pearson"="red4","RMT"="turquoise1",
                                 "GAUSS"="green","MCD"="blue",
                                 "OGK"="red1","WMCD"="purple1",
                                 "POLY2"="orange","POLY3"="khaki4",
                                 "POLY4"="skyblue1"),
                      labels = c("Pearson","RMT",
                                 "Gauss","MCD",
                                 "OGK","WMCD",
                                 "Poly2","Poly3",
                                 "Poly4"))+
  scale_linetype_manual(name = "",
                        values = c("solid"="solid",
                                   "dotted"="dotted"),
                        labels = c("Filtered","No-filtered"))
ggsave(paste0("graficos/grafico_",gsub(".*retornos_(.*).rds","\\1",arquivos[i]),".jpg"),
       device = "jpg",width = 12, height = 8)
ggsave(paste0("graficos/grafico_",gsub(".*retornos_(.*).rds","\\1",arquivos[i]),".pdf"),
       device = "pdf",width = 12, height = 8)
}
