## Extract resistance loci regions from the full VCF file

````bash
#Request interactive job
srun --nodes=1 --cpus-per-task=4 --mem=8G --pty bash

#Load module
module load VCFtools

#Navigate to project folder
cd interlink

#ABC locus
vcftools --vcf vcfsDaphnia/merged_daphnia_filtered.vcf --chr "000011F|quiver" --from-bp 2170278  --to-bp 2357097 --out ABC_locus.vcf --recode --recode-INFO-all --maf 0.10 --max-alleles 2 --keep /scicore/home/ebertd/dexter0000/interlink/plink/plinkRun32/pheno_sample_include.txt

#D locus (narrow window on peak)
vcftools --vcf vcfsDaphnia/merged_daphnia_filtered.vcf --chr "000018F|quiver" --from-bp 1801000  --to-bp 1804000 --out D_locus_narrow.vcf --recode --recode-INFO-all --maf 0.10 --max-alleles 2 --keep /scicore/home/ebertd/dexter0000/interlink/plink/plinkRun32/pheno_sample_include.txt

#E locus (As defined by Ameline 2021)
vcftools --vcf vcfsDaphnia/merged_daphnia_filtered.vcf --chr "000067F|quiver" --from-bp 233835  --to-bp 388722 --out E_locus_camille.vcf --recode --recode-INFO-all --maf 0.10 --max-alleles 2 --keep /scicore/home/ebertd/dexter0000/interlink/plink/plinkRun32/pheno_sample_include.txt

````



## Run DAPC on resistance loci regions

````R
#Load required packages
library(pegas)
library(reshape)
library(adegenet)
library(ggplot2)

#Load required data
lineages <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink3/results/dapc/PastDAPCgroups.txt")

################################################################################
#Generate and validate ABC haplotype clusters
################################################################################

#ABC supergene

#Load data
vcf<-read.vcf("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink3/results/resistance_loci_sequences/ABC_locus.vcf", which.loci = 1:1e5)
pheno <- read.csv("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink/Data/stick tests/filtered_resistotypes_corrected.csv", header=TRUE)

#Format input files for DAPC
vcf$host<-row.names(vcf) #Format sample IDs
df<-vcf
df<-df[,-ncol(df)]
df2<-loci2genind(df)

#Set random seed for reproducibility
set.seed(1)

#Determine number of groups to use
grp <-  find.clusters(df2, max.n.clust=10,n.pca = 999, stat ="BIC",
                      criterion = "goodfit", n.clust = 3)

#Group with DAPC (Use all PCs)
dapc1 <- dapc(df2, grp$grp,n.pca= 3, n.da= 2, var.contrib = TRUE, scale = FALSE)

#Extract group membership
ABCclus<-dapc1$grp

#Subset to just the phenotyped samples to verify haplotypes
host<-names(ABCclus)
ABC<-data.frame(host,ABCclus)
merged<-merge(pheno,ABC,by = "host")
ABCsubset<-merged$ABCclus

#Make resistotype table
resist<-paste(merged$C1,merged$C19,sep="")
table(ABCsubset,resist)

#Make lineage-resistotype table
df2 <- merge(lineages,merged, by.x = "ID", by.y = "host")
table(df2$group,merged$C19)

#D locus
#Load data
vcf<-read.vcf("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink3/results/resistance_loci_sequences/D_locus_narrow.vcf", which.loci = 1:1e5)
pheno <- read.csv("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink/Data/stick tests/filtered_resistotypes_corrected.csv", header=TRUE)
vcf$host<-row.names(vcf)
df<-vcf
df<-df[,-ncol(df)]
df2<-loci2genind(df)

#Set random seed for reproducibility
set.seed(1)

#Determine number of groups to use
grp <-  find.clusters(df2, max.n.clust=10,n.pca = 999, stat ="BIC",
                      criterion = "goodfit", n.clust = 3)

#Group with DAPC (Use all PCs)
dapc1 <- dapc(df2, grp$grp,n.pca= 999, n.da= 999, var.contrib = TRUE, scale = FALSE)


#Subset to just the phenotyped samples to verify haplotypes
Dclus<-dapc1$grp

#Create phenotype subset
host<-names(Dclus)
D<-data.frame(host,Dclus)
merged<-merge(pheno,D,by = "host")
Dsubset<-merged$Dclus

#Make resistotype table
resist<-paste(merged$P15,merged$P21,sep="")
table(Dsubset,resist)

#Create binary groups
Dbinary<-ifelse(Dclus==2,"Both","Sus")

#E locus
#Load data
vcf<-read.vcf("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink3/results/resistance_loci_sequences/E_locus_camille.vcf", which.loci = 1:1e5)
pheno <- read.csv("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink/Data/stick tests/filtered_resistotypes_corrected.csv", header=TRUE)

#Format
vcf$host<-row.names(vcf)
df<-vcf
df<-df[,-ncol(df)]
df2<-loci2genind(df)

#Set random seed for reproducibility
set.seed(1)

#Find groups and run DAPC (10 PC, 3 Group, 2 LD)
grp <-  find.clusters(df2, max.n.clust=10,n.pca = 999, stat ="BIC",
                      criterion = "goodfit", n.clust = 3)

#Group with DAPC (Use all PCs)
dapc1 <- dapc(df2, grp$grp,n.pca= 999, n.da= 999, var.contrib = TRUE, scale = FALSE)

Eclus<-dapc1$grp

#Subset to just the phenotyped samples to verify haplotypes
host<-names(Eclus)
E<-data.frame(host,Eclus)
merged<-merge(pheno,E,by = "host")
Esubset<-merged$Eclus

#Make table
resist<-paste(merged$C1,merged$C19,merged$P20,sep="")
table(Esubset,resist)
Ebinary<-ifelse(Eclus==1,"Both","Res")

#Merge everything together
binaryAll<-paste(ABCbinary,Dbinary,Ebinary,sep = "/")

res <- data.frame(vcf$host,ABCbinary,Dbinary,Ebinary,binaryAll)
colnames(res)[1] <- "ID"

#Save results to disk
write.table(res, "C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink3/results/resistance_loci_sequences/resistance_loci_clusters.txt",
            quote = FALSE, row.names = FALSE)
````

