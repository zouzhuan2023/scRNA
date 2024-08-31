library(Seurat)
library("limma")  
library(dplyr)
inputFile="symbol.txt"    

rt=read.table(inputFile, header=T, sep="\t", check.names=F)    
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

v=voom(data, plot=F, save.plot=F)
out=v$E
out=rbind(ID=colnames(out),out)
write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)    
source("ssGSEA18.CIBERSORT.R")
results=CIBERSORT("ref.txt", "uniq.symbol.txt", perm=1000, QN=TRUE)
file_path <- "symbol.txt" 
data <- read.table(file_path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
slpi_expression <- data["SLPI", ]
slpi_expression <- as.numeric(slpi_expression)
print(slpi_expression)
slpi_median_expression <- median(slpi_expression, na.rm = TRUE)
high_group <- colnames(data)[slpi_expression > slpi_median_expression]
low_group <- colnames(data)[slpi_expression <= slpi_median_expression]
new_data <- data.frame(
  id = colnames(data),
  SLPI = slpi_expression,
  Risk = ifelse(colnames(data) %in% high_group, "high", "low")
)
write.table(new_data, file="risk.train.txt", sep="\t", row.names=FALSE, quote=FALSE)

####################################
riskFile="risk.train.txt"           
immFile="CIBERSORT-Results.txt"     
pFilter=0.05                        
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<pFilter,]
data=as.matrix(immune[,1:(ncol(immune)-3)]) 

group=sapply(strsplit(row.names(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[group==0,]
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(data))
data=avereps(data)


risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(risk))
rt=cbind(data[sameSample,,drop=F], risk[sameSample,"Risk",drop=F])
rt=rt[order(rt$Risk, decreasing=T),]
conNum=nrow(rt[rt$Risk=="low",])

treatNum=nrow(rt[rt$Risk=="high",])
data=t(rt[,-ncol(rt)])## 
pdf("barplot.pdf", height=10, width=18)
col=rainbow(nrow(data), s=0.7, v=0.7)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(data,col=col,yaxt="n",ylab="Relative Percent",xaxt="n",cex.lab=1.8)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
par(srt=0,xpd=T)
rect(xleft = a1[1], ybottom = -0.01, xright = a1[conNum], ytop = -0.06,col="green")
text(a1[conNum]/2,-0.035,"Low risk",cex=2)
rect(xleft = a1[conNum], ybottom = -0.01, xright =a1[length(a1)] , ytop = -0.06,col="red")
text((a1[length(a1)]+a1[conNum])/2,-0.035,"High risk",cex=2)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.3)
dev.off()

data=rt
data=melt(data, id.vars=c("Risk"))
colnames(data)=c("Risk", "Immune", "Expression")

group=levels(factor(data$Risk))
data$Risk=factor(data$Risk, levels=c("low","high"))
bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", fill="Risk",
                  xlab="",
                  ylab="Fraction",
                  legend.title="Risk",
                  width=0.8,
                  palette=bioCol)+
  rotate_x_text(50)+
  stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")

pdf(file="immune.diff.pdf", width=8, height=6)
print(boxplot)
dev.off()
data = rt
data = melt(data, id.vars = c("Risk"))
colnames(data) = c("Risk", "Immune", "Expression")
data = data[data$Immune == "NK cells activated", ]
data$Risk = factor(data$Risk, levels = c("low", "high"))
bioCol = c("#0066FF", "#FF0000")
boxplot = ggboxplot(data, x = "Immune", y = "Expression", fill = "Risk",
                    xlab = "",
                    ylab = "Fraction",
                    legend.title = "Risk",
                    width = 0.5,
                    palette = bioCol) +
  rotate_x_text(0) +
  stat_compare_means(aes(group = Risk), symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")), label = "p.signif")
pdf(file = "immune.diff.filtered.pdf", width = 4, height = 4)


