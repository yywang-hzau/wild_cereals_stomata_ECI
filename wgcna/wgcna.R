#WGCNA
#导入各种包，确定文件夹
setwd('')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)
BiocManager::install("Cairo",library.dynam.unload(libpath ='/opt/X11/lib/libXrender.1.dylib' ))
BiocManager::install("impute")

library(WGCNA)
library(reshape2) 
library(stringr)
library(affy)
library(limma)
library(ggplot2)
library(pheatmap)
library(flashClust)
library(RColorBrewer)
library(cluster)
library(Cairo)
options(stringsAsFactors = FALSE)
disableWGCNAThreads()

type = "unsigned"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)

#读入转录组数据的格式
RNAseq_voom <- read.csv('merge_hvbsbhtdes.txt',header=TRUE,sep='\t',row.names = 1,check.names=F)
head(RNAseq_voom)
WGCNA_matrix = t(RNAseq_voom[order(apply(RNAseq_voom,1,mad), decreasing = T)[1:10000],])
datExpr0 <- WGCNA_matrix  ## top 10000 mad genes
dataExpr <- datExpr0 
gsg = goodSamplesGenes(dataExpr, verbose = 3)
gsg$allOK
#一般都会ok

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
dim(dataExpr)
GeneName=dimnames(dataExpr)[[2]]
#看一下样品是否偏离很厉害
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

#读入表型数据
datTraits = read.csv("merge_trait-es.txt",
                     sep='\t',header=TRUE,row.names = 1)
dim(datTraits)
table(rownames(datTraits)==rownames(dataExpr))
save(dataExpr, datTraits, file="esSamplesAndTraits.RData")
#Standard gene screening based on marginal correlation

#软阈值带图
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                        networkType=type, verbose=5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")

power = sft$powerEstimate
power=22

#构建网络
mergingThresh = 0.25
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = paste0('esexpr.txt', ".tom"),
                       verbose = 3)
table(net$colors)

#展示模块
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#模块之间相关性,特征向量
MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
write.table(MEs_col ,'MEs_col.txt',sep='\t',row.names = TRUE)

plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "networkConstruction-auto.RData")

#基因网络图
#load('bhSamplesAndTraits.RData')
#load("networkConstruction-auto.RData")
load('esexpr.txt.tom-block.1.RData')
TOM <- as.matrix(TOM)
dissTOM = 1-TOM
hierTOM = hclust(as.dist(dissTOM),method="average");
k <- softConnectivity(datE=dataExpr,power=power) 

sizeGrWindow(10, 5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main="Check Scale free topology\n")

#tom plot
nSelect = 1000
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window#sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^power;
diag(plotDiss) = NA;
TOMplot(plotDiss,
        selectTree,
        selectColors, 
        main = "Network heatmap plot, selected genes")
rm(plotDiss)
rm(selectTOM)

#模块和性状相关性
modTraitCor = cor(MEs_col, datTraits, use = "p")
write.table(modTraitCor,'trait_cor.txt',sep='\t',row.names = TRUE)
modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2),
                   "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

par(mar = c(5, 6, 4, 4)+0.1); 
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(datTraits), 
               yLabels = names(MEs_col), cex.lab = 0.9,  yColorWidth=0.01, 
               xColorWidth = 0.03,
               ySymbols = colnames(MEs_col), colorLabels = TRUE, colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))

# 画出样本聚类图（上）与样本性状热图（下）：
traitColors = numbers2colors(datTraits, signed = TRUE,centered=TRUE);
plotDendroAndColors(sampleTree, 
                    traitColors, 
                    groupLabels = names(datTraits), 
                    rowTextAlignment = "right-justified",
                    addTextGuide = TRUE ,
                    hang = 0.03,
                    dendroLabels = NULL,
                    addGuide = FALSE, 
                    guideHang = 0.05,
                    main = "Sample dendrogram and trait heatmap") 

#Finding modules that relate to a clinical trait
#性状和样品
sizeGrWindow(8,9)
colorlevels=c("turquoise", 'yellow','brown','green','blue')
par(mfrow=c(as.integer(0.5+length(colorlevels)/2),2)) 
par(mar = c(4,5,3,1))
for (n in c(1:length(colorlevels)))
{
  which.module=colorlevels[[n]];
  plotMat(t(scale(dataExpr[,moduleColors==which.module ]) ),nrgcols=30,rlabels=T,
          clabels=T,rcols=which.module, title=which.module )
}
#上调下调情况,有典型趋势的模块和特征向量展示
sizeGrWindow(8,7);
which.module="brown"
ME=MEs_col[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(dataExpr[,moduleColors==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2, 
        ylab="eigengene expression",xlab="array sample")

sizeGrWindow(8,7);
which.module="yellow"
ME=MEs_col[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(dataExpr[,moduleColors==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2, 
        ylab="eigengene expression",xlab="array sample")


HubGenes <- chooseTopHubInEachModule(dataExpr,moduleColors)
write.table(HubGenes,'hubgenes.csv',row.names = F, sep = ",")
#hub genes,GS1 datKME筛选
#cor_expr_trait
#correlation: q<0.2作为noiseindicator
for (i in c(1:16) )
{
  y=datTraits[,i]
  GS1 <- as.numeric(WGCNA::cor(y,dataExpr,use="p",method="pearson"))
  p.Standard=corPvalueFisher(GS1, nSamples =length(y))
  p.Standard2=p.Standard
  p.Standard2[is.na(p.Standard)]=1
  q.Standard=qvalue(p.Standard2)$qvalues
  StandardGeneScreeningResults=data.frame(GeneName,PearsonCorrelation=GS1,
                                          p.Standard, q.Standard,moduleColors)
  write.csv(StandardGeneScreeningResults,paste('trait',i,'bhexpr.csv', sep="_"))
}

#图解释：For the green and the brown module we observe that intramodular hub genes tend to have high gene significance. The opposite is true in the turquoise module.
ADJ1=abs(cor(dataExpr,use="p"))^power
Alldegrees1=intramodularConnectivity(ADJ1, moduleColors)
write.table(Alldegrees1,'alldegrees1.txt',sep='\t')
datKME=signedKME(dataExpr, MEs_col, outputColumnName="MM.")
###导出datKME
write.table(datKME,'datkme.txt',sep='\t')
#MM.purple,MM.blue等,针对每个module intramodular connectivity
sizeGrWindow(9,6)
par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2))) 
par(mar = c(4,5,3,1))
y=datTraits[,7]
GS1 <- as.numeric(WGCNA::cor(y,dataExpr,use="p",method="pearson"))
GeneSignificance <- abs(GS1)
ModuleSignificance <- tapply(GeneSignificance,moduleColors,mean,na.rm=T)
for (n in c(1:length(colorlevels)))
{
  whichmodule=colorlevels[[n]];
  restrict1 = (moduleColors==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     GeneSignificance[restrict1], col=moduleColors[restrict1],
                     main=whichmodule,
                     xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
}


#展示挑选出的，模块中和性状强关联的基因,基因的拓扑重叠热图，颜色越深表示拓扑重叠度越高
y=datTraits[,15]
GS1 <- as.numeric(WGCNA::cor(y,dataExpr,use="p",method="pearson"))
FilterGenes= abs(GS1)> .2 & abs(datKME$MM.turquoise)>.8 
table(FilterGenes)
dimnames(data.frame(dataExpr))[[2]][FilterGenes]
trait_hubGenes_spe <- colnames(dataExpr)[FilterGenes] 
plotNetworkHeatmap(dataExpr,
                   plotGenes = trait_hubGenes_spe,
                   networkType = "unsigned",
                   useTOM = TRUE,
                   power=power,
                   main="unsigned correlations of 1th trait in 'darkgreen' module")

#explore the relationship between the module membership measures (e.g. MM.turquoise) and intramodular connectivity
sizeGrWindow(9,6)
par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2))) 
par(mar = c(4,5,3,1))
for (n in c(1:length(colorlevels)))
{
  whichcolor=colorlevels[[n]]
  restrictGenes=(moduleColors==whichcolor);
  verboseScatterplot(Alldegrees1$kWithin[restrictGenes],
                     (datKME[restrictGenes, paste("MM.", whichcolor, sep="")])^power, 
                     col=whichcolor,
                     xlab="Intramodular Connectivity",
                     ylab="(Module Membership)^power")
}

# Output of eigengene information:
for (i in c(1:16) )
{
  y=datTraits[,i]
  datMEy = data.frame(y, MEs_col)
  eigengeneSignificance = cor(datMEy, y);
  eigengeneSignificance[1,1] = (1+max(eigengeneSignificance[-1, 1]))/2 
  eigengeneSignificance.pvalue = corPvalueStudent(eigengeneSignificance, nSamples = length(y)) 
  namesME=names(datMEy)
  out1=data.frame(t(data.frame(eigengeneSignificance,
                               eigengeneSignificance.pvalue, namesME, t(datMEy))))
  dimnames(out1)[[1]][1]="EigengeneSignificance"
  dimnames(out1)[[1]][2]="EigengeneSignificancePvalue"
  dimnames(out1)[[1]][3]="ModuleEigengeneName"
  dimnames(out1)[[1]][-c(1:3)]=dimnames(dataExpr)[[1]]
  write.table(out1, file=paste('trait',i,"ultsNetworkScreening.csv",sep='_'), row.names=TRUE, col.names = TRUE, sep=",")
  # here we output the module eigengenes and trait y without eigengene significances
  datTraitsME=data.frame(rownames(datTraits),datMEy)
  dimnames(datTraitsME)[[2]][2:length(namesME)]=paste("Trait",
                                                      dimnames(datTraits)[[2]][2:length(namesME)],
                                                      sep=".")
  write.table(datTraitsME, file=paste('trait',i,"without_eigengenesignificances.csv",sep='_'), row.names=F,sep=",")
}

#Output file for gene ontology analysis
datGS.Traits = data.frame(cor(dataExpr, datTraits, use = "p"))
names(datGS.Traits) = paste("cor", names(datGS.Traits), sep = ".")
datOutput = data.frame(colnames(dataExpr), moduleColors, datKME, datGS.Traits)
write.table(datOutput, "Forgo.csv", row.names = F, sep = ",")

#导出cyto需要的数据
probes = colnames(dataExpr) 
dimnames(TOM) <- list(probes, probes) 
cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = paste('bscyt',".edges.txt",sep=""), 
                               nodeFile = paste('bscyt', ".nodes.txt", sep=""), 
                               weighted = TRUE, threshold = 0.1,
                               nodeNames = probes, nodeAttr = moduleColors)

#NEO, One promising approach to this problem involves the integration of genetics and gene expression. SNP



