#snp
#BiocManager::install('trackViewer')

library(trackViewer)
library(Gviz)
library(rtracklayer)
library(RColorBrewer)
library(ggplot2)
library(ggthemes)
library(reshape2)
library(stringr)
#snp位点
setwd('/')
#同一个基因两个坡snp

essnppos<-read.csv('HORVU6Hr1G066050-snp-espos.txt',sep='\t',header = TRUE)
assnppos<-read.csv('HORVU6Hr1G066050-snp-aspos.txt',sep='\t',header = TRUE)
#注释信息
genegff<-read.csv('HORVU6Hr1G066050-asgff.txt',sep='\t',header=TRUE)
mycolors<-brewer.pal(9,"Set1") 

SNP <- essnppos$pos
sample.gr <- GRanges("chr4H", IRanges(SNP, width=1, names=paste0("snp-", essnppos$name)))
sample.gr$color <- mycolors[3]
sample.gr$border <- "gray30"
features <- GRanges("chr4H", IRanges(genegff$start, 
                                    width=genegff$width,
                                    names=genegff$type))

features$fill <- mycolors[c(2,9,5)]
features$height <- c(0.03,0.05, 0.03)
lolliplot(sample.gr, features)

#多个图
SNP2 <- assnppos$pos
sample2.gr <- GRanges("chr4H", IRanges(SNP2, width=1, names=paste0("snp-", assnppos$name)))
sample2.gr$color <- mycolors[1]
sample2.gr$border <- "gray30"

features2 <-GRanges("chr4H", IRanges(genegff$start, 
                                     width=genegff$width
                                     ))

features2$fill <- mycolors[c(2,9,5)]
features2$height <- c(0.03, 0.05, 0.03)
lolliplot(list(ES=sample.gr, AS=sample2.gr), 
          list(x=features, y=features2))

#snp density
hv_karyotype<-read.csv('hv.txt',sep='\t',header = TRUE)
gene<-read.csv('gene.txt',header=T,sep='\t',stringsAsFactors = F)

assnp<-read.csv('assnp.snpden.txt',sep='\t',header = TRUE)
assnp['End']=assnp$Start+999999
assnp<-assnp[,-4]
essnp<-read.csv('essnp.snpden.txt',sep='\t',header = TRUE)
essnp['End']=essnp$Start+999999

install.packages('RIdeogram')
library(RIdeogram)
data(karyotype_ternary_comparison, package="RIdeogram")
data(synteny_ternary_comparison, package="RIdeogram")

ideogram(karyotype = hv_karyotype, 
         overlaid = essnp,
         label=gene,
         label_type = 'marker',width=200,
         colorset1=c("forestgreen", "yellow", "red"),
          output = "es55chromosome.svg")
convertSVG("es55chromosome.svg", device = "tiff",width=10,height=8)


#共线性图
data(karyotype_ternary_comparison, package="RIdeogram")

data(synteny_ternary_comparison, package="RIdeogram")
data(synteny_ternary_comparison_graident, package="RIdeogram")

karyotype_ternary_comparison2<-read.table('karyotype_ternary_comparison_hvtd.txt',sep='\t',header = TRUE)
synteny_ternary_comparison2<-read.table('synteny_ternary_comparison_hvtd.txt',sep='\t',header = TRUE)
highlight<-read.table('highlight.txt',sep='\t',header = TRUE)

ideogram(karyotype = karyotype_ternary_comparison2, 
         label=highlight,
         label_type = 'marker',width=200,
         synteny = synteny_ternary_comparison2)
convertSVG("chromosome.svg", device = "png")
library("rsvg")
rsvg_pdf("chromosome.svg", "chromosome.pdf")

