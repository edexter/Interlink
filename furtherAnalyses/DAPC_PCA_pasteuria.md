# DAPC clustering for Pasteuria

This script uses DAPC to find the optimal number of clusters for pasteuria and then output the best fit group memberships.

````R
#Load required packages
library(pegas)
library(adegenet)
library(ggplot2)

#Load and format VCF
df <- read.vcf("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/lineages/data/filtered.vcf", which.loci = 1:1e5)
host <- row.names(df) #Extract sample IDs

#Load and format heterozygosity data (needed to identify multiple infections)
het <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/lineages/data/filtered.het")
colnames(het)[1] <- "ID"

#Calculate proportion of variant sites that are heterozygous
het$poly <- (het$N_SITES-het$O.HOM.)/het$N_SITES 

# Filter multiple infections and format for DAPC. Filter cutoff based on iterative
# trials and visual inspection of sequences in IGV 
df <- df[het$poly <= 0.15,]

################################################################################
#DAPC clustering
################################################################################

#Format data for DAPC
df <- loci2genind(df)

#Set seed because output can differ depending on (random) starting sample
set.seed(1)

#Determine "optimal" number of groups to use
grp <- find.clusters(df, max.n.clust=10,n.pca = 999, stat ="BIC",
                     criterion = "goodfit", choose.n.clust = FALSE)
#Group with DAPC
dapc1 <- dapc(df, grp$grp,n.pca= 50, n.da= 50, var.contrib = TRUE, scale = FALSE)

#Check PCA plots
dfPCA <- dapc1$tab
dfPCA$group <- dapc1$grp
ggplot(dfPCA, aes(x=`PCA-pc.1`,y=`PCA-pc.2`))+ geom_point(size = 3, shape=21,aes(fill=group))
ggplot(dfPCA, aes(x=`PCA-pc.1`,y=`PCA-pc.3`))+ geom_point(size = 3, shape=21,aes(fill=group))
ggplot(dfPCA, aes(x=`PCA-pc.1`,y=`PCA-pc.4`))+ geom_point(size = 3, shape=21,aes(fill=group))

#Extract DAPC groups
temp <- data.frame(names(dapc1$grp),dapc1$grp)
colnames(temp) <- c("ID", "group")

#Rename DAPC groups, including multiple infection samples
groups<-merge(x = het, y = temp, by = "ID", all = TRUE)
groups$group <- as.character(groups$group)
groups$group[is.na(groups$group)] = "Mixed"
groups$group <- ifelse(groups$group==1,"Gamma",groups$group)
groups$group <- ifelse(groups$group==2,"Alpha-3",groups$group)
groups$group <- ifelse(groups$group==3,"Beta",groups$group)
groups$group <- ifelse(groups$group==4,"Alpha-2",groups$group)
groups$group <- ifelse(groups$group==5,"Alpha-1",groups$group)

#Manually correct any samples that IGV inspection shows are multiple infections
#but do not meet heterozygosity cutoff (only one - just at border of cutoff)
groups$group[groups$ID=="F005"] <- "Mixed"

#Format data for export
temp <- groups[,c(1,7)]

#Export DAPC groups as a text file for later use
write.table(temp,"C:/Users/ericd/Dropbox/Eric Work/Ebert lab/interlink3/results/dapc/PastDAPCgroups.txt",
quote = FALSE, row.names = FALSE, sep = "\t")
````

 

## DAPC and PCA plots

This creates a plot showing the major PCA axes with DAPC groups overlaid by point color. This script calls DAPC again, but only to generate the PCA.

````R
#################################################################################
#Load required packages
################################################################################

library(pegas)
library(reshape)
library(adegenet)
library(ggplot2)

################################################################################
#Load and format data
################################################################################

#Load and format VCF
df <- read.vcf("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/lineages/data/filtered.vcf", which.loci = 1:1e5)
host <- row.names(df) #Extract sample IDs

#Load PCL files
pcl <- read.csv("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/lineages/data/c1_pcl_annotated.csv")

#Load and format heterozygosity data
het <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/lineages/data/filtered.het")
colnames(het)[1] <- "ID"
het$poly <- (het$N_SITES-het$O.HOM.)/het$N_SITES #proportion of observed sites that are heterozygous

#Load DAPC groups
DAPCgroups <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink3/results/dapc/PastDAPCgroups.txt")

#Extract variant names from VCF
varNames <- colnames(df)
varNames <- as.numeric(sub("chr1_","",varNames))

################################################################################
#DAPC clustering (Only used for the PCA output)
################################################################################

#Format for DAPC
df <- loci2genind(df)

#Set seed because output can differ depending on (random) starting sample
set.seed(1)

#Determine number of groups to use
grp <- find.clusters(df, median.n.clust=10, n.pca = 50, n.clust = 2)

#Group with DAPC
dapc1 <- dapc(df, grp$grp,n.pca= 50, n.da= 5, var.contrib = TRUE, scale = FALSE)

#Extract and format PC loadings
loadPCA <- data.frame(dapc1$pca.loadings)
loadPCA$positions <- row.names(loadPCA)
loadPCA$positions <- sub("chr1_","",loadPCA$positions)
loadPCA$positions <- gsub("\\..*","",loadPCA$positions)
loadPCA$positions <- as.numeric(loadPCA$positions)

#Get percent variance for top PCs
dapc1$pca.eig[1] / sum(dapc1$pca.eig)
dapc1$pca.eig[2] / sum(dapc1$pca.eig)
dapc1$pca.eig[3] / sum(dapc1$pca.eig)
dapc1$pca.eig[4] / sum(dapc1$pca.eig)
dapc1$pca.eig[5] / sum(dapc1$pca.eig)

################################################################################
#Nice PCA plots
################################################################################
#Format data frame
dfPCA<-data.frame(dapc1$tab)
dfPCA$group<- factor(DAPCgroups$group)
dfPCA$het <- het$poly
dfPCA <- dfPCA[order(-dfPCA$het),]

#PC1 vs PC2
p1<-ggplot(dfPCA, aes(x=PCA.pc.1,y=PCA.pc.2))+
  geom_point(size = 4,alpha=0.7, shape=21,aes(fill=as.factor(group)))+
  theme_bw()+
  scale_fill_manual(labels=c("Alpha-1", "Alpha-2","Alpha-3","Beta","Gamma","Multiple infections"),values = c("#440154FF" ,"#1f968b","#95d840","#fde725","#31688EFF","white"),
                    name = "Pasteuria cluster")+
  #theme(legend.justification = c(1, 0), legend.position = c(1, 0))+
  theme(panel.grid.major = element_blank() ,panel.grid.minor = element_blank())+
  xlab("PC 1 (48% of variance)")+ylab("PC 2 (27% of variance)")+
  theme(legend.position = "none")+
  theme(aspect.ratio=1)
p1
ggsave("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/interlink3/figures/PCA_past_1v2.png",plot = p1, width = 3, height = 3)

#PC1 vs PC3
p2<-ggplot(dfPCA, aes(x=PCA.pc.1,y=PCA.pc.3))+
  geom_point(size = 4,alpha=0.7, shape=21,aes(fill=as.factor(group)))+
  theme_bw()+
  scale_fill_manual(labels=c("Alpha-1", "Alpha-2","Beta","Gamma","Multiple infections"),values = c("#440154FF" ,"#1f968b","#95d840","#fde725","#31688EFF","white"),
                    name = "Pasteuria cluster")+
  #theme(legend.justification = c(1, 0), legend.position = c(1, 0))+
  theme(panel.grid.major = element_blank() ,panel.grid.minor = element_blank())+
  xlab("PC 1 (48% of variance)")+ylab("PC 3 (6.5% of variance)")+
  theme(legend.position = "none")+
  theme(aspect.ratio=1)
p2
ggsave("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/interlink3/figures/PCA_past_1v3.png",plot = p2, width = 3, height = 3)
````