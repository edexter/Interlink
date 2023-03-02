# DAPC and PCA analysis - daphnia

Here I use DAPC to examine the population structure among the daphnia samples. The starting data is the quality filtered daphnia VCF, which is further filtered to MAF > 5, maximum alleles =2, and thinned to retain only 1 variant every 10 KB. The thinning is required because DAPC is too computationally expensive to run on the full data set and we don't want linked variants to influence the results. DAPC is then performed on this filtered VCF and the results are plotted in R.

````bash
#Prepare the VCF file for DAPC
/scicore/home/ebertd/dexter0000/myPlink/plink2 --vcf vcfsDaphnia/merged_daphnia_filtered.vcf --out daphnia_DAPC --allow-extra-chr --maf 0.05 --max-alleles 2 --bp-space 10000 --export vcf

/scicore/home/ebertd/dexter0000/myPlink/plink2 --pfile daphnia_DAPC --export vcf --out daphnia_DAPC --allow-extra-chr
````

This section of R code estimates the number of clusters using DAPC, performs PCA, calculates the percent variance explained by each PC, and makes a PCA plot.

````R
################################################################################
#Load packages and data
################################################################################

#Load required packages
library(pegas)
library(adegenet)
library(ggplot2)

#Load data and format data
vcf<-read.vcf("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink3/results/dapc_daphnia/daphnia_DAPC.vcf", which.loci = 1:10^6)
df<-loci2genind(vcf)

################################################################################
#DAPC
################################################################################

#Set random seed
set.seed(1)

#Find optimal number of groups
grp <- find.clusters(df, max.n.clust=10,n.pca = 999, stat ="BIC",
                     criterion = "goodfit", choose.n.clust = FALSE)
#Print the number of groups
grp$size

################################################################################
#PCA
################################################################################
#Only one group is the best fit, but we can manually specify 2 groups
#so that DAPC will run and we can look at the PCA output
grp <- find.clusters(df, max.n.clust=10,n.pca = 999, stat ="BIC",
                     criterion = "goodfit", n.clust = 2)
dapc1 <- dapc(df, grp$grp,n.pca= 999, n.da= 999, var.contrib = TRUE, 
              scale = TRUE, center = TRUE)

#Assemble dataframe
dfPCA<-data.frame(dapc1$tab)
dfPCA$group<- factor(dapc1$grp)

#Get percent variance for each PC
dapc1$pca.eig[1] / sum(dapc1$pca.eig)
dapc1$pca.eig[2] / sum(dapc1$pca.eig)
dapc1$pca.eig[3] / sum(dapc1$pca.eig)
dapc1$pca.eig[4] / sum(dapc1$pca.eig)
dapc1$pca.eig[5] / sum(dapc1$pca.eig)
dapc1$pca.eig[6] / sum(dapc1$pca.eig)
dapc1$pca.eig[7] / sum(dapc1$pca.eig)
dapc1$pca.eig[8] / sum(dapc1$pca.eig)
dapc1$pca.eig[9] / sum(dapc1$pca.eig)
dapc1$pca.eig[10] / sum(dapc1$pca.eig)

#Plot PCA
p1 <- ggplot(dfPCA[], aes(x=PCA.pc.1,y=PCA.pc.2))+
  geom_point(size = 4,alpha=0.7, shape=21, fill = "#619CFF")+
  theme_bw()+
  theme(panel.grid.major = element_blank() ,panel.grid.minor = element_blank())+
  xlab("PC 1 (3.5% of variance)")+ylab("PC 2 (2.8% of variance)")+
  #theme(legend.position = "none")+
  theme(aspect.ratio=1,legend.background = element_blank(),legend.title = element_blank())
ggsave(plot = p1, "C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink3/figures/pca_daph_1v2.png", width = 3, height = 3)

p2 <- ggplot(dfPCA[], aes(x=PCA.pc.1,y=PCA.pc.3))+
  geom_point(size = 4,alpha=0.7, shape=21, fill = "#619CFF")+
  theme_bw()+
  theme(panel.grid.major = element_blank() ,panel.grid.minor = element_blank())+
  xlab("PC 1 (3.5% of variance)")+ylab("PC 3 (2.6% of variance)")+
  #theme(legend.position = "none")+
  theme(aspect.ratio=1,legend.background = element_blank(),legend.title = element_blank())
ggsave(plot = p2, "C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink3/figures/pca_daph_1v3.png", width = 3, height = 3)
````

