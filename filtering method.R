# D is the raw data where empty columns has been deleted as well as no tag and no index columns and indexes replaced by sample names
#load D
library("plyr", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
d <- rbind(names(D)[-1], colSums(D[,-1]))
d <- t(d)
d2 <- cbind(d, ldply(strsplit(as.character(d[,1]), "_"), identity) ) # split the identifier column
names(d2) <- c("id","sum", "exp", "mouse", "prog","type", "rep") # rename all columns

#mean reads filtering
#create table to check mean reads over replicate
da <- d2[which(d2$rep == 'a'),]#prend les lignes du replicat a 
da <- da[,!colnames(da)=="rep"]#elimine la colonne rep
da <- da[,!colnames(da)=="id"]
row.names(da) <- NULL #élimine la première colonne avec les numeros de ligne
dimnames(da) [[2]][1]<- "vara"#renommer la colonne var en vara

db <- d2[which(d2$rep == 'b'),]#idem avec b
row.names(db) <- NULL
db <- db[,!colnames(db)=="rep"]
db <- db[,!colnames(db)=="id"]
dimnames(db) [[2]][1]<- "varb"

dab <- merge(da,db,all=T)#fusionner da et db
dab[is.na(dab)]<-0
write.table( dab, "Desktop/LPX mean filtering a vs b.txt", sep="\t", row.names=F )

#compute mean of vara and varb
#if mean<1000 sample filtered out, done by hand on the D matrix.

#transform mean filtered data
#reload data without mean filtered samples as D
#norm and 100k
rownames(D) <- D[,1]
data <- D[,-1]
norm.data <- apply(data,2, function(x) (x/sum(x))*100000)
write.table(norm.data, "Desktop/LPX-meanfilt norm 100K.txt",row.names=TRUE,col.names=NA, sep="\t")

#transform wild to long
library("reshape2", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
d2 <- melt(norm.data, id.vars="rawnames") # convert from wide to long format
d2 <- cbind( d2, ldply(strsplit(as.character(d2[,2]), "_"), identity) ) # split the identifier column
names(d2) <- c("tag", "id","var", "exp","mouse","prog", "type", "rep") # rename all columns

#delete the zeros
d3<- d2[which(d2$var>0),]

#make matrix with replicates per samples
#a vs b
da <- d3[which(d3$rep == 'a'),]#prend les lignes du replicat a 
da <- da[,!colnames(da)=="rep"]#elimine la colonne rep
da <- da[,!colnames(da)=="id"]
row.names(da) <- NULL #élimine la première colonne avec les numeros de ligne
dimnames(da) [[2]][2]<- "vara"#renommer la colonne var en vara

db <- d3[which(d3$rep == 'b'),]#idem avec b
row.names(db) <- NULL
db <- db[,!colnames(db)=="rep"]
db <- db[,!colnames(db)=="id"]
dimnames(db) [[2]][2]<- "varb"

dab <- merge(da,db,all=T)#fusionner da et db
dab[is.na(dab)]<-0
write.table( dab, "Desktop/LPX -meanfilt norm 100K a vs b.txt", sep="\t", row.names=F )

#plot self/self, replicate against each other
library("ggplot2", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
#dev.copy(png,'selfself exp1 before cor filtering.png')
#x <- dab[which(dab$exp=="LP28"),]
qplot(asinh(vara), asinh(varb), data=dab) + facet_wrap(~prog~mouse~type)
#dev.off()
#ggsave(file="self e7 w4 before corr filtering.pdf")

#correlation filtering
#correlation on self/self
x <- ddply(dab, c("exp","prog","mouse","type"), summarize, cor=cor(vara,varb,use="na"))
x[is.na(x)]<-0
hist(x$cor, breaks=20)
write.table( x, "Desktop/LPX -meanfilt norm 100K a vs b corr linear result.txt", sep="\t", row.names=F )
y <- merge(dab,x, all=TRUE)
z<- y[which(y$cor>0.8),]
Z<- y[which(y$cor<0.8),]
write.table( z, "Desktop/LPX norm 100K-meanfilt a vs b-corrfilt.txt", sep="\t", row.names=F )
write.table( Z, "Desktop/LPX  norm 100K-meanfilt a vs b corfilt.txt", sep="\t", row.names=F )

#plot new self/self after corr filtering
qplot(asinh(vara), asinh(varb), data=z) + facet_wrap(~exp~mouse~type)

#eliminate bc that are not present in one half of a sample
z1 <- z[which(z$vara>0 & z$varb>0),]
filt3 <- z[which((z$vara>0 & z$varb==0) | (z$vara==0 & z$varb>0)),]
write.table(filt3, "Desktop/LPXnorm 100K a vs b-ratiofilt-corrfilt  abfilt.txt",quote=F, sep="\t", row.names=F )
write.table(z1, "Desktop/LPXnorm 100K a vs b-ratiofilt-corrfilt-abfilt.txt",quote=F, sep="\t", row.names=F )
qplot(asinh(vara), asinh(varb), data=z1) + facet_wrap(~prog~mouse~type)

#transform format to wild format of z1
z1[,8] <-paste(z1$exp,z1$prog,z1$mouse,z1$type, sep="_")
x1 <- z1[,-c(1,2,3,4)]
xa <- x1[,-3]
names(xa) <- c("tag","var", "id")
xa$id <-paste(xa$id,"a", sep="_")
xb <- x1[,-2]
names(xb) <- c("tag","var", "id")
xb$id <-paste(xb$id,"b", sep="_")
x1<- rbind(xa,xb)

x1 <- reshape(x1,direction="wide", timevar="id", idvar="tag" )
x1[is.na(x1)] <- 0
write.table( x1, "Desktop/LPX-ratiofilt norm 100K-meanfilt-corrfilt-abfilt wide.txt",quote=F, sep="\t", row.names=F )

#div table per mouse
z <- x1[,grep("m4", colnames(x1))]
z<- cbind(x1[,1],z)
dimnames(z)[[2]][1]<- "tag"
z <- z[which(rowSums(z[,-1])>0),]
write.table(z, "Desktop/LPX mX arcsin.txt",quote=F, sep="\t", row.names=F )

#averaging
y <- m7[,-1]
x <- numeric(0)
ave <-numeric(0)
for (i in seq(1, ncol(y-1), by = 2)) {
  x <- y[, c(i,i+1)]
  ave <- cbind(ave, (x[,1]+x[,2])/2 ) #same than before but average this time
}
ave <- cbind(rownames(y), ave)
x <- gsub("_a","",colnames(y[,seq(1, ncol(y), by = 2)]))#keep only the names of one of the column
colnames(ave) <- c("tag",x)
write.table(ave, "Desktop/LPX mX ab filt average linear.txt",quote=F, sep="\t", row.names=F )

#arcsin tranform for heatmap plotting
x <- x1[,-1]
y <- asinh(x)
y <- cbind(x1[,1],y)
dimnames(y)[[2]][1]<- "tag"
write.table(y, "Desktop/LPX-ratiofilt norm 100K-allfilt arcsin.txt",quote=F, sep="\t", row.names=F )

