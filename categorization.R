# here we will classify the barcodes (bc) for their contribution in erythrocytes (E) and myeloid cells (M)

#data
#the data should be organized as follow: first the read counts renormalized per column, then the data renormalized per row in excel keeping only E and M and serve as an input for this code

#classification with a threshold of 0
bc <- numeric(0) # create empty vector
pur <- numeric(0)

z <- ave[,grep("_3_", colnames(ave))] # select one mouse, here 3 but you should change it for every mouse and run the code again
z <- z[which(rowSums(z)>0),] # eliminate rows that are all zeros 
M <- z[which(z[,2]>0 & z[,1]==0),] # bc producing only M
E <- z[which(z[,2]==0 & z[,1]>0),] # bc producing only E
ME <- z[which(z[,2]>0 & z[,1]>0),] # bc producing both E and M
bc <- rbind(bc, cbind(nrow(M),nrow(E), nrow(ME), nrow(x[which(z[,2]>0 | z[,1]>0),]))) # count the number of bc in each categories
pur <- cbind(pur, rbind(colSums(M)/1000,colSums(E)/1000,colSums(ME)/1000)) # count the total number of reads per categories, the value to take into account is the one for the data renormalized per column

names(bc) <-c("M","E","M+E","tot") # name the column 
write.table(pur, "LP32 d6 E vs M threshold 0 output.txt",quote=F, sep="\t", row.names=F ) # record the output
write.table(bc, "LP32 d6 E vs M threshold 0.txt",quote=F, sep="\t", row.names=F )

#for different threshold
bc <- numeric(0)
pur <- numeric(0)

z <- ave[,grep("_1_", colnames(ave))]# select one mouse, here 3 but you should change it for every mouse and run the code again
z <- z[which(rowSums(z)>0),] # eliminate rows that are all zeros 
T=0.1 # define the threshold
M <- z[which(z[,2]/(z[,1]+z[,2])*100>T & z[,1]/(z[,1]+z[,2])*100<T),] # bc producing only M
E <- z[which(z[,2]/(z[,1]+z[,2])*100<T & z[,1]/(z[,1]+z[,2])*100>T),] # bc producing only E
ME <- z[which(z[,2]/(z[,1]+z[,2])*100>T & z[,1]/(z[,1]+z[,2])*100>T),]# bc producing both E and M
bc <- rbind(bc, cbind(nrow(M),nrow(E), nrow(ME)))# count the number of bc in each categories
pur <- cbind(pur, rbind(colSums(M)/1000,colSums(E)/1000,colSums(ME)/1000))# count the total number of reads per categories, the value to take into account is the one for the data renormalized per column

write.table(pur, "LP32 d6 E vs M output threshold 0.1.txt",quote=F, sep="\t", row.names=F )# record the output
write.table(bc, "LP32 d6 E vs M threshold 0.1.txt",quote=F, sep="\t", row.names=F )

#plotting is done in excel. 
