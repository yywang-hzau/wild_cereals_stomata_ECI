setwd('')
library(ggplot2)
library(ggthemes)
library(Hmisc)
library(tidyverse)
library(agricolae)
library(car)
library(reshape2)
library(ggpubr)
library(dbplyr)
library(rstatix)
library(ggpmisc)

data<-read.csv('snp_number.txt',header=TRUE,sep='\t',check.names = FALSE)
data['snp2']<-data$snp/1000000
data$group<-factor(data$group,level=c('African Slope','European Slope'))
#compaired <- list(c('WT','null'),c('WT','R1'),c('WT','R2'),c('WT','R3'))
editing_type1<-ggboxplot(data, x = "group", y = "snp2",
                         width=0.5,
                         add = c("mean_se"),
                         color = 'group',
                         ylab=expression('The number of SNP'~("10"^"6")),
                         xlab="")+
  geom_jitter(size=0.5,alpha=0.5)+
  facet_grid(species~.,scales = "free_y")+
  stat_compare_means(aes(group=group),
                     label = "p.signif",
                     hide.ns=TRUE,
                     method = 't.test'
  )+theme(legend.title=element_blank())+
  theme(legend.text = element_text(size=14))+
  theme(legend.background = element_rect())+
  scale_color_brewer(palette = 'Set1')+
  theme(panel.background=element_rect(fill='transparent', color='gray'))+
  theme(text=element_text(family = 'sans',size=14))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
editing_type1
