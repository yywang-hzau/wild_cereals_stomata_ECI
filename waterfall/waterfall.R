#GenVisR waterfall
library(tidyverse)
library(reshape2)
library(ggpubr)
library(dbplyr)
library(dplyr)
library(rstatix)
library(RColorBrewer)
library(ggsignif)
library(cowplot)
library(data.table)
library(ggstatsplot)

library(GenVisR)

setwd('')
snp<-read.csv('hswaterfall_4.txt',sep='\t',header = TRUE)
snp<-data.table(snp)
snp2<-melt.data.table(snp,id.vars = c('Variant_Classification','Hugo_Symbol'), 
                           variable.name = "Tumor_Sample_Barcode",value.name = 'bool')
snp3<-snp2[snp2$bool==1]


# 
gene<-read.csv('8as31es.txt')
sample<-read.csv('hssample.txt',sep='\t')
sampOrder=sample$sample
geneOrder=gene$Gene
set.seed(383)
snp3<-data.frame(snp3)
waterfall(snp3,mainXlabel=TRUE,mainRecurCutoff = 0.1,
          sampOrder = sampOrder)
 