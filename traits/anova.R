#author: Wyy

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
library(RColorBrewer)
library(ggsignif)
library(cowplot)

ToothGrowth
summary(ToothGrowth)
table(ToothGrowth$supp,ToothGrowth$dose)
aggregate(ToothGrowth$len,by=list(ToothGrowth$supp,ToothGrowth$dose),FUN=mean)
class(ToothGrowth$dose)
ToothGrowth$dose<-factor(ToothGrowth$dose)
class(ToothGrowth$dose)
fit<-aov(len~supp*dose,data = ToothGrowth)
summary(fit)

library(HH)
interaction2wt(len~dose*supp,data = ToothGrowth)


setwd('')
data<-read.csv('photo.txt',sep='\t',check.names = F,header = T)
summary(data)
data$species<-factor(data$species)
data$slope<-factor(data$slope)
data$treat<-factor(data$treat)
list<-unique(data$group)
data$group<-factor(data$group,levels = list)
data$color<-factor(data$color)


bartlett.test(Photo~slope,data=data)
bartlett.test(Photo~treat,data=data)
bartlett.test(Photo~species,data=data)
Photo<-aov(Photo~slope*treat+species,data=data)
summary(Photo)

#a
model<-aov(Photo~group, data=data)
#vector
out <- LSD.test(model,"group", p.adj="none")
#out$group
out$groups
sig<-out$groups
sig['group']=row.names(sig)
sig<-sig[,-1]
mean<-out$means
mean['group']=row.names(mean)
mean<-mean[,-1]
add<-merge(sig,mean,by='group')
data2<-merge(data,add,by='group')
mycolors<-brewer.pal(5,"Set1")
p<-ggplot(data2, aes(x = group, y = Photo, fill = color)) + 
  geom_boxplot(alpha=0.8) +
  geom_jitter(size=0.5,alpha=0.4) +
  geom_vline(xintercept = 4.5,linetype=3)+
  geom_vline(xintercept = 8.5,linetype=3)+
  labs(y=expression(A~('u'*'mol·m'^"-2"~"·s"^"-1")),title='')+
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    legend.title  = element_blank(),
    legend.position="",
    plot.title = element_text(hjust = 0.5),
    text=element_text(family = 'sans',size=14),
    panel.background=element_rect(fill='transparent', color='black'),
    axis.text=element_text(color='black'),
    axis.text.x=element_text(hjust=0.5),
    axis.title.x = element_blank(),
    plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),
    panel.spacing=unit(c(0.1,0.1,0.1,0.1), "cm")
  ) +
  scale_fill_manual(values =mycolors[c(2,1,3,5)])+
  geom_text(aes(label = groups ,y= Max, x = group),
            vjust = -0.3,size = 5,color='black')+
  scale_x_discrete(breaks=NULL)
p
ggsave('aaa.pdf',p,width=8,height=6)

#gs
model<-aov(gs~group, data=data)

out <- LSD.test(model,"group", p.adj="none")
out$groups
sig<-out$groups
sig['group']=row.names(sig)
sig<-sig[,-1]
mean<-out$means
mean['group']=row.names(mean)
mean<-mean[,-1]
add<-merge(sig,mean,by='group')
data2<-merge(data,add,by='group')

p2<-ggplot(data2, aes(x = group, y = gs, fill = color)) + 
  geom_boxplot() +
  geom_jitter(size=0.5,alpha=0.4) +
  geom_vline(xintercept = 4.5,linetype=3)+
  geom_vline(xintercept = 8.5,linetype=3)+
  labs(y=expression('g'[s]~("mol·m"^"-2"~"·s"^"-1")),title='')+
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    legend.title  = element_blank(),
    legend.position="",
    plot.title = element_text(hjust = 0.5),
    text=element_text(family = 'sans',size=14),
    panel.background=element_rect(fill='transparent', color='black'),
    axis.text=element_text(color='black'),
    axis.title.x=element_blank(),
    axis.text.x=element_text(hjust=0.5),
    axis.ticks.x=element_blank(),
    plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),
    panel.spacing=unit(c(0.1,0.1,0.1,0.1), "cm")
  ) +
  scale_fill_manual(values =mycolors[c(2,1,3,5)])+
  geom_text(aes(label = groups ,y= Max, x = group),
            vjust = -0.3,size = 5,color='black')+
  scale_x_discrete(breaks=NULL)
p2
ggsave('gs.pdf',p2,width=8,height=6)

#WUE
model<-aov(wue~group, data=data)
#vector
out <- LSD.test(model,"group", p.adj="none")

out$groups
sig<-out$groups
sig['group']=row.names(sig)
sig<-sig[,-1]
mean<-out$means
mean['group']=row.names(mean)
mean<-mean[,-1]
add<-merge(sig,mean,by='group')
data2<-merge(data,add,by='group')
data2$species<-factor(data2$species,levels = c('Brachypodium','H.spontaneum',
                                             'T.dicoccoides'))

p3<-ggplot(data2, aes(x = group, y = wue, fill = color)) + 
  geom_boxplot() +
  geom_jitter(size=0.5,alpha=0.4) +
  geom_vline(xintercept = 4.5,linetype=3)+
  geom_vline(xintercept = 8.5,linetype=3)+
  labs(y=expression(WUEi),title='')+
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    legend.title  = element_blank(),
    legend.position="",
    plot.title = element_text(hjust = 0.5),
    text=element_text(family = 'sans',size=14),
    panel.background=element_rect(fill='transparent', color='black'),
    axis.text=element_text(color='black'),
    axis.title.x=element_blank(),
    axis.text.x=element_text(hjust=0.5),
    axis.ticks.x=element_blank(),
    plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),
    panel.spacing=unit(c(0.1,0.1,0.1,0.1), "cm")
  ) +
  scale_fill_manual(values =mycolors[c(2,1,3,5)])+
  geom_text(aes(label = groups ,y= Max, x = group),
            vjust = -0.3,size = 5,color='black')+
  scale_x_discrete(breaks=NULL)
p3
ggsave('wue.pdf',p3,width=8,height=6)

#Ci
model<-aov(ci~group, data=data)
#vector
out <- LSD.test(model,"group", p.adj="none")

out$groups
sig<-out$groups
sig['group']=row.names(sig)
sig<-sig[,-1]
mean<-out$means
mean['group']=row.names(mean)
mean<-mean[,-1]
add<-merge(sig,mean,by='group')
data2<-merge(data,add,by='group')

p4<-ggplot(data2, aes(x = group, y = ci, fill = color)) + 
  geom_boxplot() +
  geom_jitter(size=0.5,alpha=0.4) +
  geom_vline(xintercept = 4.5,linetype=3)+
  geom_vline(xintercept = 8.5,linetype=3)+
  labs(y=expression('C'[i]~('ppm')),title='')+
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    legend.title  = element_blank(),
    legend.position="",
    plot.title = element_text(hjust = 0.5),
    text=element_text(family = 'sans',size=14),
    panel.background=element_rect(fill='transparent', color='black'),
    axis.text=element_text(color='black'),
    axis.title.x=element_blank(),
    axis.text.x=element_text(hjust=0.5),
    axis.ticks.x=element_blank(),
    plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),
    panel.spacing=unit(c(0.1,0.1,0.1,0.1), "cm")
  ) +
  scale_fill_manual(values =mycolors[c(2,1,3,5)])+
  geom_text(aes(label = groups ,y= Max, x = group),
            vjust = -0.3,size = 5,color='black')+
  scale_x_discrete(breaks=NULL)
p4
ggsave('ci.pdf',p4,width=8,height=6)

#Tr
model<-aov(tr~group, data=data)
#vector
out <- LSD.test(model,"group", p.adj="none")
#out$group
out$groups
sig<-out$groups
sig['group']=row.names(sig)
sig<-sig[,-1]
mean<-out$means
mean['group']=row.names(mean)
mean<-mean[,-1]
add<-merge(sig,mean,by='group')
data2<-merge(data,add,by='group')

p5<-ggplot(data2, aes(x = group, y = tr, fill = color)) + 
  geom_boxplot() +
  geom_jitter(size=0.5,alpha=0.4) +
  geom_vline(xintercept = 4.5,linetype=3)+
  geom_vline(xintercept = 8.5,linetype=3)+
  labs(y=expression('Transpiration rate'~('mmolH'[2]*'O'~"·m"^"-2"~"·s"^"-1")),title='')+
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    legend.title  = element_blank(),
    legend.position="",
    plot.title = element_text(hjust = 0.5),
    text=element_text(family = 'sans',size=14),
    panel.background=element_rect(fill='transparent', color='black'),
    axis.text=element_text(color='black'),
    axis.title.x=element_blank(),
    axis.text.x=element_text(hjust=0.5),
    axis.ticks.x=element_blank(),
    plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),
    panel.spacing=unit(c(0.1,0.1,0.1,0.1), "cm")
  ) +
  scale_fill_manual(values =mycolors[c(2,1,3,5)])+
  geom_text(aes(label = groups ,y= Max, x = group),
            vjust = -0.3,size = 5,color='black')+
  scale_x_discrete(breaks=NULL)
p5
ggsave('tr.pdf',p5,width=8,height=6)

#VPD
model<-aov(vpdl~group, data=data)
#vector
out <- LSD.test(model,"group", p.adj="none")#结果显示：标记字母法out$group
out$groups
sig<-out$groups
sig['group']=row.names(sig)
sig<-sig[,-1]
mean<-out$means
mean['group']=row.names(mean)
mean<-mean[,-1]
add<-merge(sig,mean,by='group')
data2<-merge(data,add,by='group')

p6<-ggplot(data2, aes(x = group, y = vpdl, fill = color)) + 
  geom_boxplot() +
  geom_jitter(size=0.5,alpha=0.4) +
  geom_vline(xintercept = 4.5,linetype=3)+
  geom_vline(xintercept = 8.5,linetype=3)+
  labs(y=expression('VPD'~("kPa")),title='')+
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    legend.title  = element_blank(),
    legend.position="",
    plot.title = element_text(hjust = 0.5),
    text=element_text(family = 'sans',size=14),
    panel.background=element_rect(fill='transparent', color='black'),
    axis.text=element_text(color='black'),
    axis.title.x=element_blank(),
    axis.text.x=element_text(hjust=0.5),
    axis.ticks.x=element_blank(),
    plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),
    panel.spacing=unit(c(0.1,0.1,0.1,0.1), "cm")
  ) +
  scale_fill_manual(values =mycolors[c(2,1,3,5)])+
  geom_text(aes(label = groups ,y= Max, x = group),
            vjust = -0.3,size = 5,color='black')+
  scale_x_discrete(breaks=NULL)
p6
ggsave('vpd.pdf',p6,width=8,height=6)


#AL
model<-aov(AL~group, data=data)
#vector
out <- LSD.test(model,"group", p.adj="none")

out$groups
sig<-out$groups
sig['group']=row.names(sig)
sig<-sig[,-1]
mean<-out$means
mean['group']=row.names(mean)
mean<-mean[,-1]
add<-merge(sig,mean,by='group')
data2<-merge(data,add,by='group')

s1<-ggplot(data2, aes(x = group, y = AL, fill = color)) + 
  geom_boxplot() +
  geom_jitter(size=1) +
  geom_vline(xintercept = 4.5,linetype=3)+
  geom_vline(xintercept = 8.5,linetype=3)+
  labs(y=expression("Aperture length"~("um")),title='')+
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    legend.title  = element_blank(),
    legend.position="",
    plot.title = element_text(hjust = 0.5),
    text=element_text(family = 'sans',size=14),
    panel.background=element_rect(fill='transparent', color='black'),
    axis.text=element_text(color='black'),
    axis.title.x=element_blank(),
    axis.text.x=element_text(hjust=0.5),
    axis.ticks.x=element_blank(),
    plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),
    panel.spacing=unit(c(0.1,0.1,0.1,0.1), "cm")
  ) +
  scale_fill_manual(values =mycolors[c(2,1,3,5)])+
  geom_text(aes(label = groups ,y= Max, x = group),
            vjust = -0.3,size = 5,color='black')+
  scale_x_discrete(breaks=NULL)
s1
ggsave('al.pdf',s1,width=8,height=6)

#AW
data<-read.csv('aw.txt',sep='\t',check.names = F,header = T)

Dat <- data

Dat[, 2] <- as.numeric(as.character(Dat[, 2]))

outlier_limup <-
  3 * IQR(Dat[, 2], na.rm = TRUE) + quantile(Dat[, 2], 3 / 4, na.rm = TRUE, names = FALSE)
#max, Q3+k(Q3-Q1)

outlier_limdown <-
  quantile(Dat[, 2], 1 / 4, na.rm = TRUE, names = FALSE) - 3 * IQR(Dat [, 2], na.rm = TRUE) 
#min, Q1-k(Q3-Q1)

Dat[Dat[, 2] >= outlier_limup |
      Dat[, 2] <= outlier_limdown, 2] = ""

#replace
Dat <- Dat[!Dat[, 2] >= outlier_limup & ! Dat[, 2] <= outlier_limdown, ]
Dat[, 2] <- as.numeric(as.character(Dat[, 2]))
data<-Dat

summary(data)
data$species<-factor(data$species)
data$slope<-factor(data$slope)
data$treat<-factor(data$treat)

data$group<-factor(data$group,levels = list)
model<-aov(aw~group, data=data)

out <- LSD.test(model,"group", p.adj="none")

out$groups
sig<-out$groups
sig['group']=row.names(sig)
sig<-sig[,-1]
mean<-out$means
mean['group']=row.names(mean)
mean<-mean[,-1]
add<-merge(sig,mean,by='group')
data2<-merge(data,add,by='group')

s2<-ggplot(data2, aes(x = group, y = aw, fill = color)) + 
  geom_boxplot() +
  geom_jitter(size=0.1,alpha=0.3) +
  geom_vline(xintercept = 4.5,linetype=3)+
  geom_vline(xintercept = 8.5,linetype=3)+
  labs(y=expression("Aperture width"~("um")),title='')+
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    legend.title  = element_blank(),
    legend.position="right",
    plot.title = element_text(hjust = 0.5),
    text=element_text(family = 'sans',size=14),
    panel.background=element_rect(fill='transparent', color='black'),
    axis.text=element_text(color='black'),
    axis.title.x=element_blank(),
    axis.text.x=element_text(hjust=0.5),
    axis.ticks.x=element_blank(),
    plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),
    panel.spacing=unit(c(0.1,0.1,0.1,0.1), "cm")
  ) +
  scale_fill_manual(values =mycolors[c(2,1,3,5)])+
  geom_text(aes(label = groups ,y= Max, x = group),
            vjust = -0.3,size = 5,color='black')+
  scale_x_discrete(breaks=NULL)
s2
ggsave('aw.pdf',s2,width=8,height=6)

#AW/AL
model<-aov(AW/AL~group, data=data)

out <- LSD.test(model,"group", p.adj="none")
out$groups
sig<-out$groups
sig['group']=row.names(sig)
sig<-sig[,-1]
mean<-out$means
mean['group']=row.names(mean)
mean<-mean[,-1]
add<-merge(sig,mean,by='group')
data2<-merge(data,add,by='group')

s3<-ggplot(data2, aes(x = group, y = AW/AL, fill = color)) + 
  geom_boxplot() +
  geom_jitter(size=1) +
  geom_vline(xintercept = 4.5,linetype=3)+
  geom_vline(xintercept = 8.5,linetype=3)+
  labs(y=expression("Aperture Width/Apertyre Length"),title='')+
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    legend.title  = element_blank(),
    legend.position="",
    plot.title = element_text(hjust = 0.5),
    text=element_text(family = 'sans',size=14),
    panel.background=element_rect(fill='transparent', color='black'),
    axis.text=element_text(color='black'),
    axis.title.x=element_blank(),
    axis.text.x=element_text(hjust=0.5),
    axis.ticks.x=element_blank(),
    plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),
    panel.spacing=unit(c(0.1,0.1,0.1,0.1), "cm")
  ) +
  scale_fill_manual(values =mycolors[c(2,1,3,5)])+
  geom_text(aes(label = groups ,y= Max, x = group),
            vjust = -0.3,size = 5,color='black')+
  scale_x_discrete(breaks=NULL)
s3
ggsave('AWal.pdf',s3,width=8,height=6)

#GCV
model<-aov(GCV~group, data=data)
#vector
out <- LSD.test(model,"group", p.adj="none")
out$groups
sig<-out$groups
sig['group']=row.names(sig)
sig<-sig[,-1]
mean<-out$means
mean['group']=row.names(mean)
mean<-mean[,-1]
add<-merge(sig,mean,by='group')
data2<-merge(data,add,by='group')

s4<-ggplot(data2, aes(x = group, y = GCV, fill = color)) + 
  geom_boxplot() +
  geom_jitter(size=1) +
  geom_vline(xintercept = 4.5,linetype=3)+
  geom_vline(xintercept = 8.5,linetype=3)+
  labs(y=expression(GCV~("um"^"3")),title='')+
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    legend.title  = element_blank(),
    legend.position="",
    plot.title = element_text(hjust = 0.5),
    text=element_text(family = 'sans',size=14),
    panel.background=element_rect(fill='transparent', color='black'),
    axis.text=element_text(color='black'),
    axis.title.x=element_blank(),
    axis.text.x=element_text(hjust=0.5),
    axis.ticks.x=element_blank(),
    plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),
    panel.spacing=unit(c(0.1,0.1,0.1,0.1), "cm")
  ) +
  scale_fill_manual(values =mycolors[c(2,1,3,5)])+
  geom_text(aes(label = groups ,y= Max, x = group),
            vjust = -0.3,size = 5,color='black')+
  scale_x_discrete(breaks=NULL)
s4
ggsave('GCV.pdf',s4,width=8,height=6)

#SCV
model<-aov(SCV~group, data=data)
#vector
out <- LSD.test(model,"group", p.adj="none")
#out$group
out$groups
sig<-out$groups
sig['group']=row.names(sig)
sig<-sig[,-1]
mean<-out$means
mean['group']=row.names(mean)
mean<-mean[,-1]
add<-merge(sig,mean,by='group')
data2<-merge(data,add,by='group')

s5<-ggplot(data2, aes(x = group, y = SCV, fill = color)) + 
  geom_boxplot() +
  geom_jitter(size=1) +
  geom_vline(xintercept = 4.5,linetype=3)+
  geom_vline(xintercept = 8.5,linetype=3)+
  labs(y=expression(SCV~("um"^"3")),title='')+
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    legend.title  = element_blank(),
    legend.position="",
    plot.title = element_text(hjust = 0.5),
    text=element_text(family = 'sans',size=14),
    panel.background=element_rect(fill='transparent', color='black'),
    axis.text=element_text(color='black'),
    axis.title.x=element_blank(),
    axis.text.x=element_text(hjust=0.5),
    axis.ticks.x=element_blank(),
    plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),
    panel.spacing=unit(c(0.1,0.1,0.1,0.1), "cm")
  ) +
  scale_fill_manual(values =mycolors[c(2,1,3,5)])+
  geom_text(aes(label = groups ,y= Max, x = group),
            vjust = -0.3,size = 5,color='black')+
  scale_x_discrete(breaks=NULL)
s5
ggsave('SCV.pdf',s5,width=8,height=6)


data<-read.csv('part1.txt',sep='\t',check.names = FALSE)
#SA
model<-aov(SA~group, data=data)

out <- LSD.test(model,"group", p.adj="none")
out$groups
sig<-out$groups
sig['group']=row.names(sig)
sig<-sig[,-1]
mean<-out$means
mean['group']=row.names(mean)
mean<-mean[,-1]
add<-merge(sig,mean,by='group')
data2<-merge(data,add,by='group')

s6<-ggplot(data2, aes(x = group, y = SA, fill = color)) + 
  geom_boxplot() +
  geom_jitter(size=1) +
  geom_vline(xintercept = 4.5,linetype=3)+
  geom_vline(xintercept = 8.5,linetype=3)+
  labs(y=expression("Aperture area"~("um"^"2")),title='')+
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    legend.title  = element_blank(),
    legend.position="",
    plot.title = element_text(hjust = 0.5),
    text=element_text(family = 'sans',size=14),
    panel.background=element_rect(fill='transparent', color='black'),
    axis.text=element_text(color='black'),
    axis.title.x=element_blank(),
    axis.text.x=element_text(hjust=0.5),
    axis.ticks.x=element_blank(),
    plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),
    panel.spacing=unit(c(0.1,0.1,0.1,0.1), "cm")
  ) +
  scale_fill_manual(values =mycolors[c(2,1,3,5)])+
  geom_text(aes(label = groups ,y= Max, x = group),
            vjust = -0.3,size = 5,color='black')+
  scale_x_discrete(breaks=NULL)
s6
ggsave('sa.pdf',s6,width=8,height=6)

#SI
model<-aov(SI~group, data=data)

out <- LSD.test(model,"group", p.adj="none")
out$groups
sig<-out$groups
sig['group']=row.names(sig)
sig<-sig[,-1]
mean<-out$means
mean['group']=row.names(mean)
mean<-mean[,-1]
add<-merge(sig,mean,by='group')
data2<-merge(data,add,by='group')

s7<-ggplot(data2, aes(x = group, y = SI, fill = color)) + 
  geom_boxplot() +
  geom_jitter(size=1) +
  geom_vline(xintercept = 4.5,linetype=3)+
  geom_vline(xintercept = 8.5,linetype=3)+
  labs(y=expression("Stomatal index"),title='')+
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    legend.title  = element_blank(),
    legend.position="",
    plot.title = element_text(hjust = 0.5),
    text=element_text(family = 'sans',size=14),
    panel.background=element_rect(fill='transparent', color='black'),
    axis.text=element_text(color='black'),
    axis.title.x=element_blank(),
    axis.text.x=element_text(hjust=0.5),
    axis.ticks.x=element_blank(),
    plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),
    panel.spacing=unit(c(0.1,0.1,0.1,0.1), "cm")
  ) +
  scale_fill_manual(values =mycolors[c(2,1,3,5)])+
  geom_text(aes(label = groups ,y= Max, x = group),
            vjust = -0.3,size = 5,color='black')+
  scale_x_discrete(breaks=NULL)
s7
ggsave('sI.pdf',s7,width=8,height=6)

#SD
model<-aov(SD~group, data=data)
#vector
out <- LSD.test(model,"group", p.adj="none")

out$groups
sig<-out$groups
sig['group']=row.names(sig)
sig<-sig[,-1]
mean<-out$means
mean['group']=row.names(mean)
mean<-mean[,-1]
add<-merge(sig,mean,by='group')
data2<-merge(data,add,by='group')

s8<-ggplot(data2, aes(x = group, y = SD, fill = color)) + 
  geom_boxplot() +
  geom_jitter(size=1) +
  geom_vline(xintercept = 4.5,linetype=3)+
  geom_vline(xintercept = 8.5,linetype=3)+
  labs(y=expression("Stomatal density"),title='')+
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    legend.title  = element_blank(),
    legend.position="",
    plot.title = element_text(hjust = 0.5),
    text=element_text(family = 'sans',size=14),
    panel.background=element_rect(fill='transparent', color='black'),
    axis.text=element_text(color='black'),
    axis.title.x=element_blank(),
    axis.text.x=element_text(hjust=0.5),
    axis.ticks.x=element_blank(),
    plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),
    panel.spacing=unit(c(0.1,0.1,0.1,0.1), "cm")
  ) +
  scale_fill_manual(values =mycolors[c(2,1,3,5)])+
  geom_text(aes(label = groups ,y= Max, x = group),
            vjust = -0.3,size = 5,color='black')+
  scale_x_discrete(breaks=NULL)
s8
ggsave('sd.pdf',s8,width=8,height=6)

#DW
model<-aov(DW~group, data=data)
#vector
out <- LSD.test(model,"group", p.adj="none")
out$groups
sig<-out$groups
sig['group']=row.names(sig)
sig<-sig[,-1]
mean<-out$means
mean['group']=row.names(mean)
mean<-mean[,-1]
add<-merge(sig,mean,by='group')
data2<-merge(data,add,by='group')

data2$species<-factor(data2$species,levels = c())
d<-ggplot(data2, aes(x = group, y = DW, fill = color)) + 
  geom_boxplot() +
  geom_jitter(size=1) +
  geom_vline(xintercept = 4.5,linetype=3)+
  geom_vline(xintercept = 8.5,linetype=3)+
  labs(y=expression('Dry weight'~('g/plant')),title='')+
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    legend.title  = element_blank(),
    legend.position="",
    plot.title = element_text(hjust = 0.5),
    text=element_text(family = 'sans',size=14),
    panel.background=element_rect(fill='transparent', color='black'),
    axis.text=element_text(color='black'),
    axis.title.x=element_blank(),
    axis.text.x=element_text(hjust=0.5),
    axis.ticks.x=element_blank(),
    plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),
    panel.spacing=unit(c(0.1,0.1,0.1,0.1), "cm")
  ) +
  scale_fill_manual(values =mycolors[c(2,1,3,5)])+
  geom_text(aes(label = groups ,y= Max, x = group),
            vjust = -0.3,size = 5,color='black')+
  scale_x_discrete(breaks=NULL)
d
ggsave('dw.pdf',d,width=8,height=6)

photo<-plot_grid(p,p2,p5,s2,labels=c('A','B','C','D'))
photo
ggsave('photo.pdf',photo,width=18,height=10)
stomata<-plot_grid(p3,p4,s4,s1,s7,s5,labels=c('A','B','C','D','E','F'))
stomata
ggsave('stomata.pdf',stomata,width=18,height=10)