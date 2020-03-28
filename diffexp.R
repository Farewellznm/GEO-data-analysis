options(stringsAsFactors = F)
rm(list = ls())
#读入数据
tumor_exp<-read.csv("data/HNSC/HNSC_tumor_exprssion_matrix.csv",row.names = 1)
tumor_exp<-log10(tumor_exp+0.0001)
annotation<-read.csv("data/HNSC/HNSC_annotation.csv")
dim(tumor_exp)
dim(annotation)
#differential
library(limma)
#使用limma差异分析，对于FPKM，log（FPKM+0.00001); eByes, add Trend = T
table(annotation$AJCC)
class<-c(rep("C",61),rep("T",56))
design<-model.matrix(~factor(class))
colnames(design)<-c("C","T")
fit<-lmFit(tumor_exp,design)
fit2<-eBayes(fit,trend = T)
allDiff=topTable(fit2,adjust='fdr',coef=2,number=200000)
write.csv(allDiff,file="data/HNSC/allDiff.csv")
#差异基因
diffLab<-allDiff[with(allDiff, ((logFC>1 |logFC<(-1)) & adj.P.Val<0.05)),]
write.csv(diffLab,file="data/HNSC/diffEXp.csv")
