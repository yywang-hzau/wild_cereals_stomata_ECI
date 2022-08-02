setwd('')
pi1<-read.csv('pop1_vs_pop2.windowed.weir.fst',sep='\t')
pi2<-pi1[which(pi1$N_VARIANTS>5 & pi1$WEIGHTED_FST>0),]
library(ggplot2)
p1<-ggplot(data=pi2,aes(x=BIN_START/1000000,y=WEIGHTED_FST))+
  geom_line(color='blue')+scale_fill_gradient(low="orange", high="purple")+
  facet_grid(~ CHROM,scales = 'free_x',space = 'free')+ylim(0,0.15)+
  theme_usual+labs(y=expression('F'[st]),x="Chromosome(Mb)")
p1

data<-read.csv('nonsyn_barplot_2.txt',
               sep='\t')
data2<-data[which(data$type=='A-G'),]

ggbarplot(data, x = "treat", y = "value", 
          add = c("mean_se"),
          fill = "vector",palette = 'npg',
          position = position_dodge(0.7),
          #ylab=expression(A ("μmol·m"^"-2"~"·s"^"-1")),a
          #ylab=expression(gs ("mol·m"^"-2"~"·s"^"-1")),gs
          #ylab=expression('C'[i]~('ppm')),ci
          #ylab=expression('VPD'~("kPa")),vpd
          #ylab=expression("T"[leaf]~(""^"o"*"C")),t
          #ylab=expression('Transpiration rate'~('mmolH'[2]*'O'~"m"^"-2"~"s"^"-1")),TR
          ylab=expression('Relative expression value'),
          xlab="",
          facet.by = 'gene')+
  stat_compare_means(aes(group = group), 
                     label = "p.signif",hide.ns=TRUE)+
  rotate_x_text(angle = 0)+
  theme(strip.text= element_text(size = 14,face='italic'))+
  theme(legend.position="right")+
  theme(legend.title  = element_blank())+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  theme(text=element_text(family = 'sans',size=14))