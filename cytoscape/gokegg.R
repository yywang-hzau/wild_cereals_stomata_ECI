
library(DESeq2)
library(ggplot2)
library(dplyr)
library(pheatmap)

#DEseq2
#输入数据，两个表格，一个是表达量，一个是样本信息
setwd('')
mRNA_exprSet<-read.csv('tdtas.txt',sep='\t',header = TRUE,row.names = 1)
condition_table<-read.csv('condition_tas.txt',sep='\t',header = TRUE)
condition_table <- condition_table[, c('sample', 'sample_type')]
rownames(condition_table) <- condition_table$sample
condition_table[,2] <- as.factor(condition_table[,2])
#指定哪一组作为control
condition_table$sample_type <- relevel(condition_table$sample_type, ref = 'Td-A')
#过滤
qualified_genes <- c()
for (genes_in_sheet in rownames(mRNA_exprSet)) {
  qualification <- mRNA_exprSet[genes_in_sheet,] <= 10
  if (sum(qualification) < 0.8*length(mRNA_exprSet)) {
    qualified_genes <- append(qualified_genes, genes_in_sheet)
  }
}
mRNA_expr_for_DESeq <- mRNA_exprSet[qualified_genes,]

#封装函数
dds <- DESeqDataSetFromMatrix(mRNA_expr_for_DESeq, colData = condition_table, design = ~ `sample_type`)
dds_DE <- DESeq(dds)
#以表格的形式，按照自己的标准分别输出所有、上调和下调的差异基因。
res_DE <- results(dds_DE, alpha=0.05,contrast=c("sample_type","Td-A","Td-B"))
resOrdered <- res_DE[order(res_DE$padj),]
resSig <- subset(resOrdered,padj < 0.05)
resSigAll <- subset(resSig, (log2FoldChange < (-1)|log2FoldChange > 1))
resSigUp <- subset(resSig,log2FoldChange >1)
resSigDown <- subset(resSig,log2FoldChange <(-1))
#export the datasheets.
dir.create('tdtasSig_genes')
setwd('tdtasSig_genes')
write.table(resSigAll, 'tdtasDE_genes_all.txt', quote = F, sep = '\t')
write.table(resSigUp, 'tdtasDE_genes_up.txt', quote = F, sep = '\t')
write.table(resSigDown, 'tdtasDE_genes_down.txt', quote = F, sep = '\t')

#variance stablizing transformation
#使用vst函数对数据进行variance stablizing transformation，并获得转换后的表达矩阵。vst转换可以使表达矩阵适用于可视化、聚类以及机器学习等操作分析。
rld <- vst(dds_DE, blind = FALSE)
expr_vst <- assay(rld)
write.table(expr_vst, 'texpr_vst.txt', quote = F, sep = '\t')

#visualization of data normalization.
par(cex = 0.7)
par(mar=c(1,3,1,1))  
n.sample=ncol(expr_vst)
cols <- rainbow(n.sample*1.2)
par(mfrow = c(2, 1))
boxplot(expr_vst, col = cols,main="expression value",las=2)

#pca
plotPCA(rld, intgroup = 'sample_type')

#ggplot2火山
for_volcano <- data.frame('log2FoldChange'=res_DE$log2FoldChange, 
                          'padj'=res_DE$padj)

for_volcano$threshold = as.factor(ifelse(for_volcano$padj < 0.05 & abs(for_volcano$log2FoldChange) >= 1, 
                                     ifelse(for_volcano$log2FoldChange> 1 ,'Up','Down'),'NoSignifi'))
ggplot(data = for_volcano, aes(x =log2FoldChange , y = -log10(padj), colour=threshold)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8)+
  labs(x="log2(fold change)",y="-log10 (p-value)",title="Volcano plot for Td") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank())

#注释用到的库
library(clusterProfiler)
library(org.At.tair.db)
library(ggplot2)
keytypes(org.At.tair.db)
setwd('')
#示例id
test<-read.csv('rhythm default  node.csv',header = TRUE,sep=',')
x<-test$ncbi_gene_id
eg=bitr(x,fromType = "ENTREZID" ,toType =c("TAIR","ENZYME" ,"SYMBOL"),OrgDb = 'org.At.tair.db')
write.table(eg,file = "rhythmnodeid.txt",sep = "\t",quote = FALSE, append = FALSE, na = "NA",row.names = FALSE)
index<-duplicated(eg$ENTREZID)
genelist<-eg[!index,]$ENTREZID
genelist2<-eg[!index,]$TAIR
#go
go <- enrichGO(genelist, OrgDb = 'org.At.tair.db', ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
               qvalueCutoff = 0.2,keyType = 'ENTREZID')
head(go)
dim(go)
dim(go[go$ONTOLOGY=='BP',])
dim(go[go$ONTOLOGY=='CC',])
dim(go[go$ONTOLOGY=='MF',])

#go 图
p1<-dotplot(go, title = 'Top10 GO terms for Bs hub genes of each sub-class', 
            showCategory = 10, color = 'p.adjust', split='ONTOLOGY')+ 
  facet_grid(ONTOLOGY~.,scale="free")
p1
write.table(go@result,file = "bshub-go.txt",sep = "\t",quote = FALSE, append = FALSE, na = "NA")
ggsave('bshub.pdf',p1,width=10,height=6,dpi=600)
#kegg
kegg <- enrichKEGG(genelist2, organism = 'ath', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
p2<-dotplot(kegg, showCategory=30,title = 'KEGG for Bs hub genes ',)
p2
write.table(kegg@result,file = "bshubkegg.txt",sep = "\t",quote = FALSE, append = FALSE, na = "NA")
ggsave('bshub-kegg.pdf',p2,width=10,height=6,dpi=600)


library(clusterProfiler)
setwd('')
data<-read.csv('Forgo.csv',sep=',',check.names = FALSE)
xx <- compareCluster(colnames.dataExpr.~moduleColors, data=data, fun = enricher,
                     TERM2GENE=wp[,c("wpid", "gene")], TERM2NAME=wp[,c("wpid", "name")])


