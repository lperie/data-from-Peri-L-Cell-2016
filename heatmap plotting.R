#heatmap function
mydistM = function(z) dist(z,'eu') 
mydistK = function(z) as.dist(cor(z,method='pearson'))
myhclu = function(d) hclust(d,method='complete')
# make sure that your input is a matrix, else it won't work
library("gplots")

#heatmap with adjustable scale
d <- asinh(x1[,-1])
d <- d[which(rowSums(d)>0),]
pairs.breaks <- seq(round(min(d)), round(max(d))+1, by=0.5)
mycol <- colorpanel(n=length(pairs.breaks)-1,low="black",mid="green",high="red")
myhm = function(y) heatmap.2(y, distfun=mydistM, hclustfun=myhclu,scale="none", breaks=pairs.breaks, col=mycol, labCol = NULL, cexCol=0.8)
m <-as.matrix(d)
x <- myhm(m)

#heatmap with fixed scale
pdf("heatmap .pdf")
d <- asinh(x1[,-1])
d <- d[which(rowSums(d)>0),]
pairs.breaks <- seq(0, 12, by=0.5)
mycol <- colorpanel(n=length(pairs.breaks)-1,low="black",mid="green",high="red")
myhm = function(y) heatmap.2(y, distfun=mydistM, hclustfun=myhclu,scale="none", breaks=pairs.breaks, col=mycol, labCol = NULL, cexCol=0.9, density.info="none", trace="none", 
                             key = FALSE, keysize = 1, key.title=NA, # no title
                             key.xlab=NA)  # no xlab)
                             m <-as.matrix(d)
                             x <- myhm(m)
                             
                             dev.off()

#heatmap sans dendrogramme
d <- em[,-c(1,3)]
d <- d[which(rowSums(d)>0),]
#myhm = function(y) heatmap.2(y,dendrogram="none", density.info="none")
myhm = function(y) heatmap(y, distfun=mydistM, hclustfun=myhclu,scale="none", col=colorRampPalette(c("black","green","red"))(1500), labCol = NULL, cexCol=0.8)
m <-as.matrix(d)
x <- myhm(m)
                             