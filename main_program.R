rm(list = ls())
setwd("E:\\Subject\\data")

mir <- read.table("COAD_miRNA_exp.txt",sep = "\t",header = T,stringsAsFactors = F)
gene <- read.table("COAD_gene_exp.txt",sep = "\t",header = T,stringsAsFactors = F)
gene_id <- as.data.frame(rownames(gene))
mir_id <- as.data.frame(rownames(mir))

#miRNA数据预处理
mir_a <- mir[rowMeans(mir)!=0,] #去除行均值为0的miRNA
mir_tumor <- grep("Tumor.*",colnames(mir_a)) #提取疾病样本
mir_control <- grep("Normal.*",colnames(mir_a)) #提取正常样本
m_exp_c <- cbind(mir_id,mir[,mir_control]) #提取行名与数据矩阵合并
m_exp_t<- cbind(mir_id,mir[,mir_tumor]) #提取行名与数据矩阵合并

m_num1 <- apply(m_exp_c[,-1],1,function(x){sum(x==0)}) #计算一行中0的个数
m_data_con <- m_exp_c[m_num1<length(mir_control)*0.3,] #保留行（0的数量少于30%）

m_num2 <- apply(m_exp_t[,-1],1,function(x){sum(x==0)}) #计算一行中0的个数
m_data_tum <- m_exp_t[m_num2<length(mir_tumor)*0.3,] #保留行（0的数量少于30%）

mir_b <- merge(m_data_tum,m_data_con,by="rownames(mir)") #将疾病和正常样本再合并
rownames(mir_b) <- mir_b$`rownames(mir)`
mir_c <- mir_b[,-1]
#m_exp_a <- cbind(rownames(mir_a),mir_a) #提取行名与数据矩阵合并
#colnames(m_exp_a)[1] <- "id" 

#mRNA数据预处理
gene_tumor <- grep("Tumor.*",colnames(gene)) #提取疾病样本
gene_control <- grep("Normal.*",colnames(gene)) #提取正常样本
g_exp_c <- cbind(gene_id,gene[,gene_control]) #提取行名与数据矩阵合并
g_exp_t<- cbind(gene_id,gene[,gene_tumor]) #提取行名与数据矩阵合并

g_num1 <- apply(g_exp_c[,-1],1,function(x){sum(x==0)}) #计算一行中0的个数
g_data_con <- g_exp_c[g_num1<length(gene_control)*0.3,] #保留行（0的数量少于30%）

g_num2 <- apply(g_exp_t[,-1],1,function(x){sum(x==0)}) #计算一行中0的个数
g_data_tum <- g_exp_t[g_num2<length(gene_tumor)*0.3,] #保留行（0的数量少于30%）

gene_a <- merge(g_data_tum,g_data_con,by="rownames(gene)") #将疾病和正常样本再合并
ensem_id <- as.character(gene_a$`rownames(gene)`) #提取gene名称

library(clusterProfiler) #加载R包
entre_id <- bitr(ensem_id, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db") #将ENSEMBL转为ENTREZID
colnames(gene_a)[1] <- "ENSEMBL"
g_exp_b <- merge(gene_a,entre_id,by="ENSEMBL") #将数据，gene名称合并
g_exp_b <- g_exp_b[,c(length(colnames(g_exp_b)),2:(length(colnames(g_exp_b))-1))] #去除ENSEMBL列
geneidfactor <- factor(g_exp_b$ENTREZID)  ##将基因ID转化为因子
g_exp_d <- as.data.frame(apply(g_exp_b[,-1],2,function(x) tapply(x,geneidfactor,mean))) #将多个ENSEMBL名称对应一个ENTREZID取均值

#选取同样的样本
sid <- intersect(colnames(g_exp_d),colnames(mir_a)) #提取相同列名
mir_s <- mir_c[,sid]
gene_s <- g_exp_d[,sid]


#找出差异的miRNA
m_t <- grep("Tumor.*",colnames(mir_s))
m_c <- grep("Normal.*",colnames(mir_s))
m_p <- apply(mir_s,1,function(x) t.test(x[m_t],x[m_c],alternative = "two.sided")$p.value)
m_FDR <- p.adjust(m_p, method="BH",length(m_p))
m_FC <- apply(mir_s,1,function(x) log2(mean(x[m_t])/mean(x[m_c])))
mirR <- row.names(mir_s)[m_FDR<0.01 & (m_FC> 1 | m_FC< (-1))]

all_mir <- cbind.data.frame(m_p,m_FDR,m_FC) 
colnames(all_mir) <- c("p","FDR","FC")

dif_mir <- cbind.data.frame(m_p[mirR],m_FDR[mirR],m_FC[mirR])
colnames(dif_mir) <- c("p","FDR","FC")

mirR_Up <- all_mir[(all_mir$FDR < 0.01 & (all_mir$FC > 1)),]
mirR_Down <- all_mir[(all_mir$FDR < 0.01 & (all_mir$FC < (-1))),]


#火山图
pdf(file="miRNA_vol.pdf")
xMax <- max(all_mir$FC)+1
yMax <- max(-log10(all_mir$FDR))+1
plot(all_mir$FC,-log10(all_mir$FDR), ylab="-log10(FDR)",xlab="log2FC",main="Volcano", xlim=c(-xMax,xMax),ylim=c(0,yMax),xaxs="i",pch=20, cex=0.4)
points(mirR_Up$FC,-log10(mirR_Up$FDR),pch=20, col="red",cex=0.4) #上调gene标红
points(mirR_Down$FC,-log10(mirR_Down$FDR), pch=20, col="green",cex=0.4) #下调gene标绿
abline(v=1,lty=2,lwd=1)
abline(v=-1,lty=2,lwd=1)
abline(h=2,lty=2,lwd=1)
dev.off()

#heatmap
library(gplots)
Cluster_Method<-c( "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
myclust <- function(x){hclust(x,method=Cluster_Method[1])}
heatmapData_m <- as.matrix(mir_s[rownames(dif_mir),])
hmExp_m <- log10(heatmapData_m+0.001)
pdf(file="miRNA_heatmap.pdf",width=60,height=30)
par(oma=c(10,3,3,7)) #确定边界
heatmap.2(hmExp_m,col='greenred',key=T,keysize=0.8,cexCol=0.8,cexRow=0.5,trace="none",hclustfun=myclust)
dev.off()

#找出差异的gene
g_t <- grep("Tumor.*",colnames(gene_s))
g_c <- grep("Normal.*",colnames(gene_s))
g_p <- apply(gene_s,1,function(x) t.test(x[g_t],x[g_c],alternative = "two.sided")$p.value)
g_FDR <- p.adjust(g_p, method="BH",length(g_p))
g_FC <- apply(gene_s,1,function(x) log2(mean(x[g_t])/mean(x[g_c])))
gen <- row.names(gene_s)[g_FDR<0.01 & (g_FC>1 | g_FC<(-1))]

all_gene <- cbind.data.frame(g_p,g_FDR,g_FC)
colnames(all_gene) <- c("p","FDR","FC")

dif_gene <- cbind.data.frame(g_p[gen],g_FDR[gen],g_FC[gen])
colnames(dif_gene) <- c("p","FDR","FC")

gene_Up <- all_gene[(all_gene$FDR < 0.01 & (all_gene$FC > 1)),]
gene_Down <- all_gene[(all_gene$FDR < 0.01 & (all_gene$FC < (-1))),]

#火山图
pdf(file="gene_vol.pdf")
xMax <- max(all_gene$FC)+1
yMax <- max(-log10(all_gene$FDR))+1
plot(all_gene$FC,-log10(all_gene$FDR), ylab="-log10(FDR)",xlab="log2FC",main="Volcano", xlim=c(-xMax,xMax),ylim=c(0,yMax),xaxs="i",pch=20, cex=0.4)
points(gene_Up$FC,-log10(gene_Up$FDR),pch=20, col="red",cex=0.4) #上调gene标红
points(gene_Down$FC,-log10(gene_Down$FDR), pch=20, col="green",cex=0.4) #下调gene标绿
abline(v=1,lty=2,lwd=1)
abline(v=-1,lty=2,lwd=1)
abline(h=2,lty=2,lwd=1)
dev.off()

#heatmap
library(gplots)
Cluster_Method<-c( "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
myclust <- function(x){hclust(x,method=Cluster_Method[1])}
heatmapData_g <- as.matrix(gene_s[rownames(dif_gene),])
hmExp_g <- log10(heatmapData_g+0.001)
pdf(file="gene_heatmap.pdf",width=60,height=30)
par(oma=c(10,3,3,7)) #确定边界
heatmap.2(hmExp_g,col='greenred',key=T,keysize=0.8,cexCol=0.8,cexRow=0.5,trace="none",hclustfun=myclust)
dev.off()

#创建WGCNA所需数据
library(WGCNA)

mir_diff <- mir_s[mirR,]
gene_diff <- gene_s[gen,]
exp <- rbind.data.frame(mir_diff,gene_diff,stringsAsFactors = F)
datExpr <- t(exp)

subname <- sapply(colnames(exp),function(x) strsplit(x,".",fixed = T)[[1]][1])
datTraits  <-  data.frame(gsm=names(exp),subtype=subname)

powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="green")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#选择beta值
best_beta=sft$powerEstimate

net = blockwiseModules(datExpr, power = sft$powerEstimate,
                       maxBlockSize = 6000,TOMType = "unsigned", 
                       minModuleSize = 30,reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "module-TOM",
                       verbose = 3)

mergedColors = labels2colors(net$colors)
table(mergedColors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
#首先针对样本做个系统聚类树
datExpr_tree<-hclust(dist(datExpr), method = "average")
sample_colors <- numbers2colors(as.numeric(factor(datTraits$subtype)), 
                                colors = c("grey","blue","red","green"),signed = FALSE)
## 这个给样品添加对应颜色的代码需要自行修改以适应自己的数据分析项目。
#  sample_colors <- numbers2colors( datTraits ,signed = FALSE)
## 如果样品有多种分类情况，而且 datTraits 里面都是分类信息，那么可以直接用上面代码，
##当然，这样给的颜色不明显，意义不大。
#构造10个样品的系统聚类树及性状热图
par(mar = c(1,4,3,1),cex=0.8)
plotDendroAndColors(datExpr_tree, sample_colors,
                    groupLabels = colnames(sample),
                    cex.dendroLabels = 0.8,
                    marAll = c(1, 4, 3, 1),
                    cex.rowText = 0.01,
                    main = "Sample dendrogram and trait heatmap")

design <- model.matrix(~0+ datTraits$subtype)
colnames(design) <- levels(datTraits$subtype)
moduleColors <- labels2colors(net$colors)
# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <-  orderMEs(MEs0); ##不同颜色的模块的ME值矩阵(样本vs模块)
moduleTraitCor = cor(MEs, design , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#提取指定模块的基因名
# Select module
module <- "turquoise"
# Select module probes
probes <- colnames(datExpr) ## 我们例子里面的probe就是基因名
inModule <- (moduleColors==module)
modProbes <- probes[inModule]

#基因富集
ge_id <- modProbes[-c(1:95)]
kk <- enrichKEGG(ge_id,organism = 'hsa',pvalueCutoff = 0.05,minGSSize = 1,use_internal_data =FALSE)
barplot(kk,showCategory = 10, title="The KEGG enrichment analysis")
dotplot(kk,showCategory = 25, title="The KEGG enrichment analysis")

#gene、miRNA富集
recon_p <- read.table("reconstruct_pathway.txt",sep="\t",header=T,stringsAsFactors = F)
id_all <- read.table("risk_pathway_background.txt",sep="\t",header=T,stringsAsFactors = F)
l <- length(id_all$ID)
d <- length(modProbes)

fj <- function(x) {
  a <- length(intersect(modProbes,unlist(strsplit(x[3],"\\|"))))
  b <- length(unlist(strsplit(x[3],"\\|")))
  k <- phyper(a-1,b,l-b,d,lower.tail = F)
  return (k)
}
fj_p <- apply(recon_p,1,function(x) fj(x))
risk_p <- recon_p[fj_p<0.05,]
write.table(risk_p,"risk_pathway.txt",sep = "\t",quote = F,row.names = F)

#模块的导出
# Recalculate topological overlap
TOM <- TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate)
# Select the corresponding Topological Overlap
modTOM <- TOM[inModule, inModule]
dimnames(modTOM) <- list(modProbes, modProbes)

diff_value <- rbind.data.frame(dif_mir,dif_gene,stringsAsFactors = F)
diff_m_value <- diff_value[modProbes,]



test_path <- risk_p[1:2,]
id_diff <- row.names(diff_m_value)
aa <- intersect(id_diff,unlist(strsplit(test_path$Gene[1],"\\|")))
bb <- intersect(id_diff,unlist(strsplit(test_path$Gene[2],"\\|")))
dbb <- data.frame("0","0",stringsAsFactors = F)
colnames(dbb) <- c("id1","id2")
for(i in 1:length(aa)){
  cc <- setdiff(bb,aa[i])
  daa <- data.frame(aa[i],cc)
  colnames(daa) <- c("id1","id2")
  dbb <- rbind.data.frame(dbb,daa,stringsAsFactors = F)
}
dcc <- dbb[-1,]

xg <- apply(dcc,1,function(x) modTOM[x[1],x[2]])
FDR1 <- diff_m_value[dcc$id1,2]
FDR2 <- diff_m_value[dcc$id2,2]

res_d <- cbind.data.frame(dcc,FDR1,FDR2,xg)

cs <- apply(res_d,1,function(x) ((-log(as.numeric(x[3])))+(-log(as.numeric(x[4]))))*as.numeric(x[5]))
res <- cbind.data.frame(res_d,cs)
CS <- sum(res$cs)


