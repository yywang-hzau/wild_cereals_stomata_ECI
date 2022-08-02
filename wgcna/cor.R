library("corrplot")
library(Hmisc)
library(pheatmap)
library(WGCNA)
library(tidyverse)
library(reshape2)
library(ggplot2)

setwd('wgcna/as')

data<-read.csv('merge_trait-as.txt',sep='\t',
               header = TRUE,row.names = 1,check.names = F)

data2<-read.csv('MEs_col.txt',sep='\t',
               header = TRUE,row.names = 1,check.names = F)
nSamples = nrow(data)
corr<-cor(data2,data)
corr<-as.matrix((corr))
corp = corPvalueStudent(corr, nSamples)

corr<-reshape2::melt(corr)
corp<-reshape2::melt(corp)
cor1<-merge(corr,corp,by=c('Var1','Var2'))
cor1<-cor1 %>%
  mutate(text=case_when(
    value.y <0.05 & value.y >0.01 ~paste(round(value.x,2),'\n*'),
    value.y < 0.01 ~ paste(round(value.x,2),'\n**')
  ))

ggplot(cor1, aes(Var2, Var1)) + 
  geom_tile(aes(fill = value.x), colour = "grey60", size = 0)+
  scale_fill_gradient2(low = "forestgreen",mid = "white",high = "firebrick4",breaks=c(-0.6,0,0.6))+
  geom_text(aes(label=text),col ="black",size = 4) +
  theme_minimal() + 
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0, size = 12, face = "bold"), # 调整x轴文字，字体加粗
        axis.text.y = element_text(size = 14, face = "bold"),
        plot.title=element_text(size=14,face='bold',hjust=0.5)) +
  labs(fill='Correlation',title='Module-trait relationships of African Slope')+
  scale_x_discrete(position = "top")

