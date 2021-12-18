# Heterozygosity calculations for pasteuria genome

````bash
#Load required modules
Module load VCFtools

#Calculate per sample heterozygosity
vcftools --vcf pasteuria_plink.vcf --het
````

Examine and plot the data in R

````R
################################################################################
#Load required data and packages
###############################################################################

het <- read.csv("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink2/results/heterozygosity/summary.csv")

haplotype_groups <- read.csv("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink2/results/haplotypes/haplotype_groups.csv")

################################################################################
#Format data
################################################################################

colnames(het)[1]<-"host"
merged<-merge(haplotype_groups,het,by = "host" )

################################################################################
#Plots
################################################################################

#Histogram of polymorphism values
png(height=700, width=1000,pointsize=20, file="C:/Users/ericd/Downloads/past_hetero_hist.png")
hist(merged$N.polymorphic,breaks=50,main = "",xlab="Proportion of polymorphic sites")
abline(v=0.15,lty=2,col="red")
text(0.32, 50, "Double or triple infected")
text(0.085, 50, "Mostly single infected")
dev.off()

#Polymorphism by DAPC pasteuria cluster
boxplot(merged$N.polymorphic~merged$pastClusAll)
````