
##plotting pca for uoa_brahman_1
setwd("D:/maulana/GSE_second_paper/Revision/pca/5")
library(ggplot2)
#full_uoa
data = read.table("ready_combined_pca.eigenvec",header=T, sep="\t")
gp <- ggplot(data,aes(x=data$PC1, y=data$PC2, group=Breed, color=Breed, shape=Breed)) +
  scale_shape_manual(values=1:nlevels(data$Breed)) +
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(x="Component 1", y=" Component 2", title = "full_uoa")+
  geom_point(size=3)
gp
data$Breed

data[-c(119:129),]

#introgressed_regions
data1 = read.table("subset_ready_combined_pca.eigenvec",header=T, sep="\t")
gp1 <- ggplot(data1,aes(x=data1$PC1, y=data1$PC2, group=Breed, color=Breed, shape=Breed)) +
  scale_shape_manual(values=1:nlevels(data1$Breed)) +
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(x="Component 1", y=" Component 2", title = "introgressed_reg")+
  geom_point(size=3)
gp1
#non-introgressed_regions
data2 = read.table("subset1_ready_combined_pca.eigenvec",header=T, sep="\t")
gp2 <- ggplot(data2,aes(x=data2$PC1, y=data2$PC2, group=Breed, color=Breed, shape=Breed)) +
  scale_shape_manual(values=1:nlevels(data$Breed)) +
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(x="Component 1", y=" Component 2", title = "non-introgressed_reg")+
  geom_point(size=3)
gp2

#ars_full
data3 = read.table("D:/maulana/GSE_second_paper/Revision/pca/5/combined_filtered_extra_pca.eigenvec",header=T, sep="\t")
gp3 <- ggplot(data3,aes(x=data3$PC1, y=data3$PC2, group=Breed, color=Breed, shape=Breed)) +
  scale_shape_manual(values=1:nlevels(data$Breed)) +
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(x="Component 1", y=" Component 2", title = "ars_full")+
  geom_point(size=3)
gp3


##reading admixture analysis
#full uoa
setwd("D:/maulana/admixture_linux-1.3.0")
tbl=read.table("ready_combined_sr.3.Q")
data$Breed
barplot(t(as.matrix(tbl)), col=rainbow(3), xlab="Breed", ylab="K = 3",border=NA) #names.arg = data$Breed, 
#introgressed_uoa 
tbl1=read.table("subset_ready_combined.3.Q")
data$Breed
barplot(t(as.matrix(tbl1)), col=rainbow(3), xlab="Breed", ylab="K = 3",border=NA) #names.arg = data$Breed, 
#non-introgressed_uoa 
tbl2=read.table("subset1_ready_combined.3.Q")
data$Breed
barplot(t(as.matrix(tbl2)), col=rainbow(3), xlab="Breed", ylab="K = 3",border=NA) #names.arg = data$Breed, 
#full ars
tbl3=read.table("combined_filtered_extra.3.Q")
data$Breed
barplot(t(as.matrix(tbl3)), col=rainbow(3), xlab="Breed", ylab="K = 3",border=NA) #names.arg = data$Breed, 


text(seq(1.5, end_point, by = 2), par("usr")[3]-0.25, 
     srt = 60, adj = 1, xpd = TRUE,
     labels = paste(tbl$Individual), cex=0.65)

library(reshape2)
setwd("D:/maulana/admixture_linux-1.3.0")
tbl=read.table("ready_combined_sr.3.Q")
tbl$Individual = rownames(tbl)
tbl$Breed = data$Breed
tbl.m = melt(tbl, id.vars="Individual")
ggplot(tbl.m, aes(x = Individual, y =value, fill = variable)) +
  geom_bar(stat="identity")
library(lattice)
barchart(V1+V2+V3~Individual, data=tbl)


tbl = tbl[order(tbl$Individual),]
barplot(t(as.matrix(tbl)), col=rainbow(3), 
        xlab="Breed", ylab="Ancestry",border=NA)

x      = runif(8)
gender = factor(c("male","female","male","female","male","female","male","female"))
group  = c(0,0,1,1,2,2,3,3)
df     = data.frame(x,gender,group)

ggplot(df,aes(x=group,y=x,fill=gender)) + 
  geom_bar(stat="identity",position="dodge") + 
  scale_x_continuous("",breaks=c(0:3),
                     labels=c('G1','G2','G3','G4'))

###the tutorial
# library
library(ggplot2)
# create a dataset
specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
value <- abs(rnorm(12 , 0 , 15))
data <- data.frame(specie,condition,value)
# Stacked + percent
ggplot(data, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="fill", stat="identity")
##tweak the tutorial
specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3), rep("sorgho" , 3)  )
condition <- rep(c("normal" , "stress" , "Nitrogen") , 5)
value <- abs(rnorm(15 , 0 , 15))
data <- data.frame(specie,condition,value)
ggplot(data, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="fill", stat="identity")

#tweak ars
setwd("D:/maulana/GSE_second_paper/Revision/admixture")
tbl3=read.table("combined_filtered_extra.3.Q", na.strings = "NA", as.is = T)
library(reshape2)
tbl.m = melt(tbl3, id.vars="Individual")
ggplot(tbl.m, aes(x = Individual, y =value, fill = variable)) +
  geom_bar(stat="identity")

#tweak uoa
setwd("D:/maulana/GSE_second_paper/Revision/admixture")
tbl1=read.table("ready_combined_sr.3.Q")
library(reshape2)
tbl.m = melt(tbl3, id.vars="Individual")
ggplot(tbl.m, aes(x = Individual, y =value, fill = variable)) +
  geom_bar(stat="identity")

##Simple t-test for comparisons of alignment against ars_ucd and uoa_brahman
alignment = read.table("D:/maulana/GSE_second_paper/Revision/Additional_file_1.txt", sep = '\t', header = TRUE, na.strings = "NA")
t.test(as.integer(alignment$retain_reads_ars), as.integer(alignment$clean_reads_uoa))

##Simple t-test for number of SNVs retained in  alignment against ars_ucd and uoa_brahman
library(dplyr)
alignment = read.table("D:/maulana/Analysis/Plink/total_snps.csv", sep = ',', header = TRUE, na.strings = "NA")
taurus = filter(alignment, Type=="Btt")
t.test(as.integer(taurus$ARS_UCD1.2), as.integer(taurus$UOA_Brahman1))
indicus = filter(alignment, Type=="Bti")
t.test(as.integer(indicus$ARS_UCD1.2), as.integer(indicus$UOA_Brahman1))
