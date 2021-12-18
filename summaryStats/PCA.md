# PCA and plots for both pasteuria and daphnia

````
#Note that the input files undergo considerable filtering and thinning before conducting PCA. This is documented in the main PLINK script.

#Daphnia PCA
scicore/home/ebertd/dexter0000/plinkDev/plink2 --pca --bfile daphnia.geno.mind.king.maf.hwe --no-pheno --out daphnia --allow-extra-chr

#Pasteuria PCA
/scicore/home/ebertd/dexter0000/plinkDev/plink2 --pca --bfile pasteuria.geno.mind.king.maf.hwe --no-pheno --out pasteuria --allow-extra-chr
````

Import the files into R for analysis and plotting

````R
################################################################################
#Load required data and packages
################################################################################

#Load packages
library(stringr)
library(ggplot2)
library(corrplot)

#Load the daphnia eigenvectors
daphnia <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink2/results/pca/daphnia.eigenvec", header=FALSE, comment.char="#")

#Load the pasteuria eigenvectors
past <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink2/results/pca/pasteuria.eigenvec", header=FALSE, comment.char="#")

#Load the phenotypes
pheno <- read.csv("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink/Data/stick tests/filtered_resistotypes_corrected.csv")

#Load DAPC haplotype groups
haplotype_groups <- read.csv("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink2/results/haplotypes/haplotype_groups.csv")

################################################################################
#Format data
################################################################################

#Reformat the sample names so they match in all data frames
pheno$host<- str_replace(pheno$host, "CH-H-2019-", "")

#Format the data for input to the corrplot package
df <- data.frame(daphnia[,3:12],past[,3:12])

colnames(df)<-c("H1","H2","H3","H4","H5","H6","H7","H8","H9","H10",
                "P1","P2","P3","P4","P5","P6","P7","P8","P9","P10")

#Find the intersection of the PCA and phenotype dataframes
#Note that this is smaller than either, because not all sequenced lines survived
#long enough to be phenotyped and not all phenotyped lines passed all the stages of
#extraction, sequencing, and filtering
daphMerged<-merge(x=daphnia, y=pheno, by.x= "V2", by.y = "host")
pastMerged<-merge(x=past, y=pheno, by.x= "V2", by.y = "host")

#Make a single resistotype variable
daphMerged$type<-as.factor(paste(daphMerged$C1,daphMerged$C19,daphMerged$P15,daphMerged$P20,daphMerged$P21,sep=""))
pastMerged$type<-as.factor(paste(pastMerged$C1,pastMerged$C19,pastMerged$P15,pastMerged$P20,pastMerged$P21,sep=""))

################################################################################
#Analysis
################################################################################

#Test correlation between daphnia and pasteuria PCs
df <-cor(df)
testRes = cor.mtest(df, conf.level = 0.95)

#Export the resistoype counts
temp2<-data.frame(table(daphMerged$type))
write.csv(temp2,"C:/Users/ericd/Downloads/temp.csv",quote=FALSE)

################################################################################
#Plots
################################################################################

png(height=1200, width=1200,pointsize=20, file="C:/Users/ericd/Downloads/haplotype_PCAcorrelations.png")
corrplot.mixed(df,lower = "number",upper="square")
dev.off()

png(height=1200, width=1200,pointsize=20, file="C:/Users/ericd/Downloads/haplotype_PCAcorrelations_2.png")
corrplot(df, p.mat = testRes$p, method = 'circle', type = 'lower', insig='blank',
         addCoef.col ='black', number.cex = 0.8, diag=FALSE)
dev.off()

#Plot daphnia PCA by resistotypes
ggplot(daphMerged, aes(x=V3,y=V4))+ geom_point(aes(color=type))
ggplot(daphMerged, aes(x=V3,y=V4))+ geom_point(aes(color=C1))
ggplot(daphMerged, aes(x=V3,y=V4))+ geom_point(aes(color=C19))
ggplot(daphMerged, aes(x=V3,y=V4))+ geom_point(aes(color=P15))
ggplot(daphMerged, aes(x=V3,y=V4))+ geom_point(aes(color=P20))
ggplot(daphMerged, aes(x=V3,y=V4))+ geom_point(aes(color=P21))

#Plot pasteuria PCA by resistotypes
ggplot(pastMerged, aes(x=V3,y=V4))+ geom_point(aes(color=type))
ggplot(pastMerged, aes(x=V3,y=V4))+ geom_point(aes(color=C1))
ggplot(pastMerged, aes(x=V3,y=V4))+ geom_point(aes(color=C19))
ggplot(pastMerged, aes(x=V3,y=V4))+ geom_point(aes(color=P15))
ggplot(pastMerged, aes(x=V3,y=V4))+ geom_point(aes(color=P20))
ggplot(pastMerged, aes(x=V3,y=V4))+ geom_point(aes(color=P21))

#Plot Daphnia PCs by P15 with nicer graphics
p<-ggplot(daphMerged, aes(x=V3,y=V4))+
  geom_point(size = 2, aes(color=P15),alpha=0.5)+
  geom_point(size = 2, data = subset(daphMerged, P15 == "R"),alpha=0.9,color="#F8766D")+
  theme_bw() +
  theme(panel.grid.major = element_blank() ,panel.grid.minor = element_blank())+
  xlab("PC 1 (4% of variance)")+ylab("PC 2 (3% of variance)")+ggtitle("PCA of host genotypes")+
  #guides(color=guide_legend(title="P15 isolate"))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_discrete(name = "Host susceptibility\n   to P15 isolate", labels = c("Resistant", "Susceptible"))
ggsave("C:/Users/ericd/Downloads/daph_pca_P15.png",plot = p, width = 7, height = 4)
p

#Plot Pasteuria PCs by P15 with nicer graphics
p<-ggplot(pastMerged, aes(x=V3,y=V4))+
  geom_point(size = 2, aes(color=P15),alpha=0.5)+
  geom_point(size = 2, data = subset(pastMerged, P15 == "R"),alpha=0.9,color="#F8766D")+
  theme_bw() +
  theme(panel.grid.major = element_blank() ,panel.grid.minor = element_blank())+
  xlab("PC 1 (47% of variance)")+ylab("PC 2 (35% of variance)")+ggtitle("PCA of parasite genotypes")+
  #guides(color=guide_legend(title="P15 isolate"))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_discrete(name = "Host susceptibility\n   to P15 isolate", labels = c("Resistant", "Susceptible"))
ggsave("C:/Users/ericd/Downloads/past_pca_P15.png",plot = p, width = 7, height = 4)
p

#Plot haplotype groups by daphnia PCA
p<-ggplot(daphnia, aes(x=V3,y=V4))+
  geom_point(size = 4,alpha=0.7, shape=21,aes(fill=as.factor(haplotype_groups$pastClusAll)))+
  theme_bw()+
  scale_fill_manual(labels=c("Cluster 1", "Cluster 2","Cluster 3","Multiple infections"),values = c("#440154FF" ,"#21908CFF","#FDE725FF","white"),
                      name = "Pasteuria cluster")+
  theme(legend.justification = c(1, 0), legend.position = c(1, 0))+
  theme(panel.grid.major = element_blank() ,panel.grid.minor = element_blank())+
  xlab("PC 1 (4.1% of variance)")+ylab("PC 2 (3.1% of variance)")
  
p
ggsave("C:/Users/ericd/Downloads/daph_pca_P15.png",plot = p, width = 8, height = 6)

#############################################################################
#End
#############################################################################
````

