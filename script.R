####session10####
library(GEOquery)
library(limma)
library(pheatmap)
library(umap)
library(ggplot2)
library(reshape2)
library(plyr)
library(gplots)

setwd("~/Advanced_Bioinformatics1/")
res<- read.delim("table.txt")
head(res)
aml.up<- subset(res,logFC > 1 & adj.P.Val < 0.05)
dim(aml.up)
aml.up.gene <- unique(aml.up$Gene.symbol)
length(aml.up.gene)
head(aml.up.gene)
####session11####
setRepositories()
#1 2
#install.packages(c("GEOquery","limma","gplots","pheatmap","Biobase","umap"))

setwd("C:/Users/m/Documents/Advanced_Bioinformatics1")
series <- "GSE9476"
gset <- getGEO(series,GSEMatrix =TRUE, AnnotGPL=TRUE,destdir ="data/")#data is adress#

platform <- "GPL96"
if (length(gset) > 1) idx <- grep(platform , attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

gr <- c("CD34",rep("BM",10),rep("CD34",7),rep("AML",26),rep("PB",10),rep("CD34",10))
length(gr)
ex <- exprs(gset) #make matrix#
dim(ex)
max(ex)
min(ex)
#qualify control#
pdf("result/boxplot.pdf",width = 64)
boxplot(ex)
dev.off()
# ex <- normalizeQuantiles(ex)
# exprs(gset) <- ex
pdf("result/corheatmap1.pdf",width=15, height = 15)
pheatmap(cor(ex),labels_row = gr, labels_col = gr)
dev.off()
#cor(ex)
# excor<- cor(ex)
# excor[1:5,1:5]

pdf("result/corheatmap.pdf",width=15, height = 15)
pheatmap(cor(ex))
dev.off()

#library(gplots)
pdf("result/corheatmap2.pdf",width=15, height = 15)
pheatmap(cor(ex),labels_row = gr, labels_col = gr,color = greenred(256))
dev.off()

pdf("result/corheatmap3.pdf",width=15, height = 15)
pheatmap(cor(ex),labels_row = gr, labels_col = gr,color = greenred(256),border_color = NA)
dev.off()


pdf("result/corheatmap4.pdf",width=15, height = 15)
pheatmap(cor(ex),labels_row = gr, labels_col = gr,color = bluered(256),border_color = NA)
dev.off()
####principal component analysis####
pc<- prcomp(ex)
pdf("result/pc.pdf")
plot(pc$x[,1:2])
dev.off()
 
#names(pc)

#for genes#
ex.scale <- t(scale(t(ex),scale = F))
mean(ex.scale[1,])
pc<- prcomp(ex.scale)
pdf("result/pc_scale.pdf")
plot(pc$x[,1:2])
dev.off()

#for samples#
#pc$r equal pc$rotation
#library(ggplot2)
#library(reshape2)
#library(plyr)
pcr <- data.frame(pc$r[,1:3],group=gr)

pdf("result/pca_samples2.pdf")
ggplot(pcr,aes(PC1,PC2,color=group))+geom_point(size=3) + theme_bw()
dev.off()



head(aml.up.gene)

#x<-c("A","B","A","C","A")
#x
#factor(x)
#y <- factor(x)
#as.numeric(factor(x))
#levels(factor(x))
#factor(x)
gr <- factor(gr)
gset$description <- gr
design <- model.matrix(~ description +0,gset)
colnames(design) <- levels(gr)  
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(AML-CD34,levels = design)
fit2 <- contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit2,0.01)
tT<-topTable(fit2,adjust="fdr",sort.by="B",number=Inf)
tT <- subset(tT,select=c("Gene.symbol","Gene.ID","adj.P.Val","logFC"))
write.table(tT,"result/AML_CD34.txt",row.names =F,sep ="t",quote = F)
#####session12####
aml.up <- subset(tT,logFC > 1 & adj.P.Val < 0.05)
aml.up.genes <- unique(aml.up$Gene.symbol)
#aml.up.genes <- sub("///.*","",aml.up.genes)
aml.up.genes <- unique(as.character(strsplit2(aml.up.genes,"///")))
write.table(aml.up.genes,file = "result/AML_CD34_up.txt",quote = F,row.names = F,col.names = F)


aml.down <- subset(tT,logFC < -1 & adj.P.Val < 0.05)
#aml.down.genes <- sub("///.*","",aml.down.genes)
aml.down.genes <- unique(as.character(strsplit2(aml.down.genes,"///")))
write.table(aml.down.genes,file = "result/AML_CD34_down.txt",quote = F,row.names = F,col.names = F)

