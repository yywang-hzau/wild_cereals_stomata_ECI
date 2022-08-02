#pinpis 
library(ggplot2)
setwd('pinpis/')
data<-read.csv('pinpis_merge.txt',sep='\t',header=TRUE)

theme_usual<-theme(
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  legend.title  = element_blank(),
  legend.position="right",
  plot.title = element_text(hjust = 0.5),
  text=element_text(family = 'sans',size=14),
  panel.background=element_rect(fill='transparent', color='black'),
  axis.text=element_text(color='black'),
  axis.text.x=element_text(angle=0,size=14,hjust=0.5))

p1<-ggplot(data,aes(x=pin.pis_bh,y=pin.pis_bs))+
  geom_jitter(alpha=0.5)+
  geom_point(data = dat2,aes(x=pin.pis_bh,y=pin.pis_bs),color='blue',alpha=0.5)+
  geom_point(data = dat3,aes(x=pin.pis_bh,y=pin.pis_bs),color='blue',alpha=0.5)+
  geom_point(data = dat4,aes(x=pin.pis_bh,y=pin.pis_bs),color='red',alpha=0.5)+
  geom_hline(yintercept = 1,linetype=2)+
  geom_vline(xintercept = 1,linetype=2)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),breaks = seq(0,10,1))+
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)),breaks = seq(0,10,1))+
  theme_usual+
  labs(x=expression(paste(pi,'N/',pi,'S in Bh homologs')),
       y=expression(paste(pi,'N/',pi,'S in Bs homologs')),
       title='')
  library(clusterProfiler)
library(org.At.tair.db)


#示例id
data2<-read.csv('pin_pis_forgo.txt',sep='\t',header = TRUE)

go2<-data2[data2$ns_bh<(1)&data2$ns_bs>(1),]
go3<-data2[data2$ns_bh>(1)&data2$ns_bs<(1),]
go4<-data2[data2$ns_bh>(1)&data2$ns_bs>(1),]

x<-data2$atid

eg=bitr(x,fromType = "TAIR" ,toType =c("ENTREZID","ENZYME" ,"SYMBOL"),OrgDb = 'org.At.tair.db')
write.table(eg,file = "pinpis_go.txt",sep = "\t",quote = FALSE, append = FALSE, na = "NA")

#write.table(eg,file = "rhythmnodeid.txt",sep = "\t",quote = FALSE, append = FALSE, na = "NA",row.names = FALSE)
index<-duplicated(eg$ENTREZID)
genelist<-eg[!index,]$ENTREZID
genelist2<-eg[!index,]$TAIR

#go
go <- enrichGO(genelist, OrgDb = 'org.At.tair.db', ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
               qvalueCutoff = 0.2,keyType = 'ENTREZID')

#go 图
p1<-dotplot(go, title = 'GO terms for Bh piN/piS of each sub-class', 
            showCategory = 10, color = 'p.adjust', split='ONTOLOGY')+ 
  facet_grid(ONTOLOGY~.,scale="free")
p1

write.table(go@result,file = "bsmorethan1-go.txt",sep = "\t",quote = FALSE, append = FALSE, na = "NA")
ggsave('bhpinpis.pdf',p1,width=10,height=6,dpi=600)
#kegg
kegg <- enrichKEGG(genelist2, organism = 'ath', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
p2<-dotplot(kegg, showCategory=30,title = 'KEGG for Bh piN/piS genes ',)
p2
write.table(kegg@result,file = "bhpiNpiSkegg.txt",sep = "\t",quote = FALSE, append = FALSE, na = "NA")
ggsave('bhpiNpiS.pdf',p2,width=10,height=6,dpi=600)

data<-read.csv('piNpiS_plot.txt',sep='\t',check.names = F)

p4<-
  ggplot(data,aes(x=species,y=`pin/pis`,color=species))+
  geom_boxplot()+
  geom_jitter(size=0.5,alpha=0.6) +
  labs(y=expression("πN/πS"),title='')+
  theme_usual+
  scale_color_manual(values=brewer.pal(5,"Set1"))+
  facet_grid(.~group)+
  stat_compare_means(aes(group=species),
                     label = "p.signif",
                     hide.ns=FALSE,
                     method = 't.test'
  )
p4
