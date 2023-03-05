# Admixture analysis for pasteuria

First an input VCF file needs to be specified that contains contain only bi-allelic markers for single infection samples. The input file are then converted into plink PED format for admixture analysis.

````bash
#Convert to bed
/scicore/home/ebertd/dexter0000/plinkDev/plink2 --vcf past_single_infection.vcf --make-bed --out past_single_infection --allow-extra-chr

#Filter missing genotypes
/scicore/home/ebertd/dexter0000/plinkDev/plink2 --bfile past_single_infection --allow-extra-chr --geno 0.50 --make-bed --out past_single_infection_geno

#Choosing the best value for K using cross vaidation
for K in 1 2 3 4 5 6 7 8 9 10; do admixture --cv past_single_infection_geno.bed $K | tee log${K}.out; done

#Check the cross-validation results
grep -h CV log*.out
````



## Admixture plots in R

````R
#Load required packages
library(reshape)
library(ggplot2)
library(ggthemes)
library(viridis)
library(patchwork)
library(RColorBrewer)
library(gridExtra)

################################################################################
#Elbow plot for admixture
K <- c(1,2,3,4,5,6,7,8,9,10)
CV <- c(0.70634,0.49691,0.17480,0.16413,0.14236,0.13225,0.11778,0.11441,0.11820,0.10828)
df2<-data.frame(K,CV)

p1 <- ggplot(df2, aes(x=K,y=CV))+
  geom_point()+
  geom_line()+
  theme_bw() +
  theme(panel.grid.major = element_blank() ,panel.grid.minor = element_blank())+
  xlab("Number of groups (K)")+ylab("Cross-validation value (CV)")+
  scale_x_continuous(breaks = seq(0, 10, by = 1))+
  theme(plot.title = element_text(hjust = 0.5))

ggsave("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/interlink3/figures/admixture_Kval.png",plot = p1, width = 3, height = 3)

################################################################################
#K = 3

#Load sample IDs
ID <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/markers/data/past_single_infection.fam", header=FALSE)
ID <- ID[,2]

#Load admixture results
Q <- read.table("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/markers/data/past_single_infection_geno.3.Q", quote="\"", comment.char="")

#Load and subset haplotype data
hap<- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink2/results/haplotypes/haplotype_groups3.txt")
hap<-hap[hap$host %in% ID,]

#Merge dataframe and melt
df<-data.frame(ID,hap$pastClusAll.1,Q)
dfL <- melt(df, id=c("ID","hap.pastClusAll.1"),measure.vars = c("V1","V2","V3"))

#Admixture plot

kplot3 <-  ggplot(dfL, aes(factor(ID), value, fill = factor(variable))) +
  geom_col(color = "gray", size = 0.1)+
  facet_grid(~hap.pastClusAll.1, switch = "x", scales = "free", space = "free")+ 
  theme_minimal() + labs(x = "Individuals", title = "K=3", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_gdocs(guide = FALSE)+
  scale_fill_brewer(palette = "Set1")+
  theme(legend.position = "none")

################################################################################
#K = 4

#Load sample IDs
ID <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/markers/data/past_single_infection.fam", header=FALSE)
ID <- ID[,2]

#Load admixture results
Q <- read.table("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/markers/data/past_single_infection_geno.4.Q", quote="\"", comment.char="")

#Load and subset haplotype data
hap<- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink2/results/haplotypes/haplotype_groups3.txt")
hap<-hap[hap$host %in% ID,]

#Merge dataframe and melt
df<-data.frame(ID,hap$pastClusAll.1,Q)
dfL <- melt(df, id=c("ID","hap.pastClusAll.1"),measure.vars = c("V1","V2","V3","V4"))

#Admixture plot

kplot4 <-  ggplot(dfL, aes(factor(ID), value, fill = factor(variable))) +
  geom_col(color = "gray", size = 0.1)+
  facet_grid(~hap.pastClusAll.1, switch = "x", scales = "free", space = "free")+ 
  theme_minimal() + labs(x = "Individuals", title = "K=4", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_gdocs(guide = FALSE)+
  scale_fill_brewer(palette = "Set1")+
  theme(legend.position = "none")

################################################################################
#K = 5

#Load sample IDs
ID <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/markers/data/past_single_infection.fam", header=FALSE)
ID <- ID[,2]

#Load admixture results
Q <- read.table("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/markers/data/past_single_infection_geno.5.Q", quote="\"", comment.char="")

#Load and subset haplotype data
hap<- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink2/results/haplotypes/haplotype_groups3.txt")
hap<-hap[hap$host %in% ID,]

#Merge dataframe and melt
df<-data.frame(ID,hap$pastClusAll.1,Q)
dfL <- melt(df, id=c("ID","hap.pastClusAll.1"),measure.vars = c("V1","V2","V3","V4","V5"))

#Admixture plot

kplot5 <-  ggplot(dfL, aes(factor(ID), value, fill = factor(variable))) +
  geom_col(color = "gray", size = 0.1)+
  facet_grid(~hap.pastClusAll.1, switch = "x", scales = "free", space = "free")+ 
  theme_minimal() + labs(x = "Individuals", title = "K=5", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_gdocs(guide = FALSE)+
  scale_fill_brewer(palette = "Set1")+
  theme(legend.position = "none")

################################################################################
#K = 6

#Load sample IDs
ID <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/markers/data/past_single_infection.fam", header=FALSE)
ID <- ID[,2]

#Load admixture results
Q <- read.table("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/markers/data/past_single_infection_geno.6.Q", quote="\"", comment.char="")

#Load and subset haplotype data
hap<- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink2/results/haplotypes/haplotype_groups3.txt")
hap<-hap[hap$host %in% ID,]

#Merge dataframe and melt
df<-data.frame(ID,hap$pastClusAll.1,Q)
dfL <- melt(df, id=c("ID","hap.pastClusAll.1"),measure.vars = c("V1","V2","V3","V4","V5","V6"))

#Admixture plot

kplot6 <-  ggplot(dfL, aes(factor(ID), value, fill = factor(variable))) +
  geom_col(color = "gray", size = 0.1)+
  facet_grid(~hap.pastClusAll.1, switch = "x", scales = "free", space = "free")+ 
  theme_minimal() + labs(x = "Individuals", title = "K=6", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_gdocs(guide = FALSE)+
  scale_fill_brewer(palette = "Set1")+
  theme(legend.position = "none")

################################################################################
#K = 7

#Load sample IDs
ID <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/markers/data/past_single_infection.fam", header=FALSE)
ID <- ID[,2]

#Load admixture results
Q <- read.table("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/markers/data/past_single_infection_geno.7.Q", quote="\"", comment.char="")

#Load and subset haplotype data
hap<- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink2/results/haplotypes/haplotype_groups3.txt")
hap<-hap[hap$host %in% ID,]

#Merge dataframe and melt
df<-data.frame(ID,hap$pastClusAll.1,Q)
dfL <- melt(df, id=c("ID","hap.pastClusAll.1"),measure.vars = c("V1","V2","V3","V4","V5","V6","V7"))

#Admixture plot

kplot7 <-  ggplot(dfL, aes(factor(ID), value, fill = factor(variable))) +
  geom_col(color = "gray", size = 0.1)+
  facet_grid(~hap.pastClusAll.1, switch = "x", scales = "free", space = "free")+ 
  theme_minimal() + labs(x = "Individuals", title = "K=7", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_gdocs(guide = FALSE)+
  scale_fill_brewer(palette = "Set1")+
  theme(legend.position = "none")

#export 1000 by 1000
kplot3 + kplot4 + kplot5 + kplot6 + kplot7 + plot_layout(ncol = 1)

#Report most probably group membership for each sample
################################################################################  
#Load and subset data
DAPCpast<-readRDS("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink2/results/haplotypes/confidence/DAPCpast.RDS")
PCA<-DAPCpast$tab
PCA$ID<-row.names(PCA)
PCA<-PCA[PCA$ID %in% ID,]
df<-data.frame(df,PCA)  

#Jittered pie chart view
library(scatterpie)

#Get fill colors from admixture plo
g <- ggplot_build(kplot3)
myCol<-unique(g$data[[1]]["fill"])

df3<-df

p2 <- ggplot() + geom_scatterpie(size=0.1,pie_scale=1,aes(x=PCA.pc.1,y=PCA.pc.3), data=df3,
                           cols=c("V1","V2","V3")) + coord_equal()+
  theme_bw() +
  theme(panel.grid.major = element_blank() ,panel.grid.minor = element_blank())+
  scale_fill_brewer(palette = "Set1",
                    name = "Dominant \nAdmixture \ngroup", labels = c("1 (Gamma)", "2 (Beta)", "3 (Alpha-3)","4 (Alpha-2)", "5 (Alpha-1)", "6 (Gamma)", "7 (Beta)"))+
  xlab("PC 1") +ylab("PC 3")+
  theme(legend.position = "none")
p2

ggsave("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/interlink3/figures/admixture_PCA_3groups.png",plot = p2, width = 3, height = 3)
````
