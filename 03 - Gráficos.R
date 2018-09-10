## carregando bases
library(data.table)
library(dplyr)
library(tidyverse)

arquivos <- list.files(pattern = "retornos_.*rds",
                       full.names = T)

i <- 1
for(i in 1:length(arquivos)){
base <- readRDS(arquivos[i])
names(base)[grep("^RMT$",names(base))] <- "pearson_RMT"
base <- base %>% gather(key = "variable" ,
                        value = "valor",
                        -DT)


# gráfico KPCA ----
base$type_line <- ifelse(grepl("_RMT",base$variable),"dotted","solid")
base$variable2 <- gsub('_RMT',"",base$variable)

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
                      values = c("GAUSS"="green",
                                 "MCD"="blue",
                                 "OGK"="red1",
                                 "pearson"="red4",
                                 "POLY2"="orange",
                                 "POLY3"="pink",
                                 "POLY4"="skyblue1",
                                 "WMCD"="purple1"),
                      labels = c("Gauss","MCD",
                                 "OGK","Pearson",
                                 "Poly2","Poly3",
                                 "Poly4",
                                 "WMCD"))+
  scale_linetype_manual(name = "",
                        values = c("solid"="solid",
                                   "dotted"="dotted"),
                        labels = c("Filtered","Non-filtered"))
ggsave(paste0("graficos/grafico_",gsub(".*retornos_(.*).rds","\\1",arquivos[i]),".jpg"),
       device = "jpg",width = 12, height = 8)
ggsave(paste0("graficos/grafico_",gsub(".*retornos_(.*).rds","\\1",arquivos[i]),".pdf"),
       device = "pdf",width = 12, height = 8)



base <- base %>% select(DT,variable2,type_line,valor) %>% 
  spread (key = "type_line" , value = "valor")

g <- ggplot(data = base, 
            aes(x = DT, y = (dotted-solid), color = variable2))+
  geom_line(size=.9) + 
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  theme(legend.position = "bottom") +
  xlab('Time') +
  ylab('Cumulative return on filtered approach') +
  scale_colour_manual(name = "",
                      values = c("GAUSS"="green",
                                 "MCD"="blue",
                                 "OGK"="red1",
                                 "pearson"="red4",
                                 "POLY2"="orange",
                                 "POLY3"="pink",
                                 "POLY4"="skyblue1",
                                 "WMCD"="purple1"),
                      labels = c("Gauss","MCD",
                                 "OGK","Pearson",
                                 "Poly2","Poly3",
                                 "Poly4",
                                 "WMCD"))
ggsave(paste0("graficos/relativo_grafico_",gsub(".*retornos_(.*).rds","\\1",arquivos[i]),".jpg"),
       device = "jpg",width = 12, height = 8)
ggsave(paste0("graficos/relativo_grafico_",gsub(".*retornos_(.*).rds","\\1",arquivos[i]),".pdf"),
       device = "pdf",width = 12, height = 8)

g <- ggplot(data = base %>% filter(!grepl("POLY",variable2)), 
            aes(x = DT, y = (dotted-solid), color = variable2))+
  geom_line(size=.9) + 
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  theme(legend.position = "bottom") +
  xlab('Time') +
  ylab('Cumulative Return') +
  scale_colour_manual(name = "",
                      values = c("GAUSS"="green",
                                 "MCD"="blue",
                                 "OGK"="red1",
                                 "pearson"="red4",
                                 "WMCD"="purple1"),
                      labels = c("Gauss","MCD",
                                 "OGK","Pearson",
                                 "WMCD"))
ggsave(paste0("graficos/relativo_sem_poly_grafico_",gsub(".*retornos_(.*).rds","\\1",arquivos[i]),".jpg"),
       device = "jpg",width = 12, height = 8)
ggsave(paste0("graficos/relativo_sem_poly_grafico_",gsub(".*retornos_(.*).rds","\\1",arquivos[i]),".pdf"),
       device = "pdf",width = 12, height = 8)

print(i)
}
