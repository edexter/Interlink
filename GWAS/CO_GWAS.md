# CO-GWAS

This bash script filters and formats the necessary input files for co-GWAS analysis in PLINK. It also aggregates and formats the resulting output files. Note that it calls a separate R script, which is given in the second code block. Note that some filter steps do not do anything or are set at extremely permissive values. They are kept in place to maintain consistent code across different versions of the model.

````bash
#!/bin/bash

#SBATCH --job-name=plinkRun    			#Job name
#SBATCH --cpus-per-task=20           	#Number of cores reserved
#SBATCH --mem-per-cpu=2G            	#Memory reserved per core.
                                     	#Total memory reserved: 32GB

#SBATCH --time=168:00:00            	#Maximum time the job will run
#SBATCH --qos=1week                    	#The job queue (time based)

#This is the stdout file
#SBATCH --output=/scicore/home/ebertd/dexter0000/interlink/logfiles/plinkRun.out

#Things remember when running on sciCore:
#This job runs from the current working directory
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.

#Load required modules
module load R/4.0.0-foss-2018b
module load VCFtools
module load GATK
################################################################################
#Preliminary steps
################################################################################
#Name the run
RUN=plinkRun

#Navigate to the project directory
cd /scicore/home/ebertd/dexter0000/interlink

#Create the run directory
mkdir plink/"$RUN"

#Create the results directory
mkdir plink/"$RUN"/results

################################################################################
#File preparation
################################################################################

#Create PLINK files
/scicore/home/ebertd/dexter0000/plinkDev/plink2 --vcf vcfsDaphnia/daphnia_plink.vcf --make-bed --out plink/"$RUN"/daphnia --allow-extra-chr

/scicore/home/ebertd/dexter0000/plinkDev/plink2 --vcf vcfsPasteuria/pasteuria_plink.vcf --make-bed --out plink/"$RUN"/pasteuria --allow-extra-chr

#Filter missing genotypes
/scicore/home/ebertd/dexter0000/plinkDev/plink2 --bfile plink/"$RUN"/daphnia --allow-extra-chr --geno 0.90 --make-bed --out plink/"$RUN"/daphnia.geno

/scicore/home/ebertd/dexter0000/plinkDev/plink2 --bfile plink/"$RUN"/pasteuria --allow-extra-chr --geno 0.90 --make-bed --out plink/"$RUN"/pasteuria.geno

#Filter missing individuals
/scicore/home/ebertd/dexter0000/plinkDev/plink2 --bfile plink/"$RUN"/daphnia.geno --allow-extra-chr --mind 0.05 --make-bed --out plink/"$RUN"/daphnia.geno.mind

#Filter missing individuals
/scicore/home/ebertd/dexter0000/plinkDev/plink2 --bfile plink/"$RUN"/pasteuria.geno --allow-extra-chr --mind 0.75 --make-bed --out plink/"$RUN"/pasteuria.geno.mind

#Filter daphnia samples based on daphnia relatedness
/scicore/home/ebertd/dexter0000/plinkDev/plink2 --bfile plink/"$RUN"/daphnia.geno.mind \
       --king-cutoff 0.350 \
       --make-bed \
       --out plink/"$RUN"/daphnia.geno.mind.king --allow-extra-chr

#Filter pasteuria samples to match above daphnia samples	   
/scicore/home/ebertd/dexter0000/plinkDev/plink2 --bfile plink/"$RUN"/pasteuria.geno.mind \
	--keep plink/"$RUN"/daphnia.geno.mind.king.king.cutoff.in.id \
	--make-bed \
	--out plink/"$RUN"/pasteuria.geno.mind.king --allow-extra-chr    

#Filter based on MAF
/scicore/home/ebertd/dexter0000/plinkDev/plink2 --bfile plink/"$RUN"/daphnia.geno.mind.king --allow-extra-chr --maf 0.10 --mac 10 --make-bed --out plink/"$RUN"/daphnia.geno.mind.king.maf.LD

/scicore/home/ebertd/dexter0000/plinkDev/plink2 --bfile plink/"$RUN"/pasteuria.geno.mind.king --allow-extra-chr --maf 0.10 --mac 10 --make-bed --out plink/"$RUN"/pasteuria.geno.mind.king.maf.LD

#Filter based on LD
/scicore/home/ebertd/dexter0000/plinkDev/plink2 --bfile plink/"$RUN"/daphnia.geno.mind.king.maf.LD \
     --indep-pairwise 1kb 1 0.9 \
     --out plink/"$RUN"/daphnia --allow-extra-chr

/scicore/home/ebertd/dexter0000/plinkDev/plink2 --bfile plink/"$RUN"/daphnia.geno.mind.king.maf.LD --make-bed --out plink/"$RUN"/daphnia.geno.mind.king.maf --allow-extra-chr

/scicore/home/ebertd/dexter0000/plinkDev/plink2 --bfile plink/"$RUN"/pasteuria.geno.mind.king.maf.LD \
     --indep-pairwise 1kb 1 0.9 \
     --out plink/"$RUN"/pasteuria --allow-extra-chr

/scicore/home/ebertd/dexter0000/plinkDev/plink2 --bfile plink/"$RUN"/pasteuria.geno.mind.king.maf.LD --make-bed --out plink/"$RUN"/pasteuria.geno.mind.king.maf --allow-extra-chr

#Filter based on HWE
/scicore/home/ebertd/dexter0000/plinkDev/plink2 --bfile plink/"$RUN"/daphnia.geno.mind.king.maf \
     --make-bed \
     --out plink/"$RUN"/daphnia.geno.mind.king.maf.hwe --allow-extra-chr

/scicore/home/ebertd/dexter0000/plinkDev/plink2 --bfile plink/"$RUN"/pasteuria.geno.mind.king.maf \
     --make-bed \
     --out plink/"$RUN"/pasteuria.geno.mind.king.maf.hwe --allow-extra-chr

################################################################################
#Phenotype file
################################################################################
#Create list of samples that match previous steps
awk 'FNR> 1 {print $2}' plink/"$RUN"/daphnia.geno.mind.king.king.cutoff.in.id > plink/"$RUN"/pheno_sample_include.txt

#Create list of variants that match previous steps
/scicore/home/ebertd/dexter0000/plinkDev/plink2 --write-snplist --bfile plink/"$RUN"/pasteuria.geno.mind.king.maf.hwe --out plink/"$RUN"/past_variants --allow-extra-chr

#The variants list requires some formatting
awk  '{gsub("_","\t",$0); print;}' plink/"$RUN"/past_variants.snplist > plink/"$RUN"/pheno_variants_include.txt

#filter VCF based on MAF and sample list
module purge
module load VCFtools
vcftools --vcf vcfsPasteuria/pasteuria_plink.vcf \
--keep plink/"$RUN"/pheno_sample_include.txt \
--positions  plink/"$RUN"/pheno_variants_include.txt \
--out plink/"$RUN"/pasteuria_plink_phenotypes.vcf \
--recode --recode-INFO-all

#Get reads supporting each allele from VCF files
module load GATK
gatk VariantsToTable -V plink/"$RUN"/pasteuria_plink_phenotypes.vcf.recode.vcf \
-F ID -F REF -F ALT -GF AD -O plink/"$RUN"/pastAD \
-R /scicore/home/ebertd/dexter0000/pasteuria_ref/Pramosa_C1_genome.gapfilled.13122016.fasta

#Get read depth for each site/sample
gatk VariantsToTable -V plink/"$RUN"/pasteuria_plink_phenotypes.vcf.recode.vcf \
-F ID -F REF -F ALT -GF DP -O plink/"$RUN"/pastDP \
-R /scicore/home/ebertd/dexter0000/pasteuria_ref/Pramosa_C1_genome.gapfilled.13122016.fasta

#Get genotype calls from Pasteuria VCF
gatk VariantsToTable -V plink/"$RUN"/pasteuria_plink_phenotypes.vcf.recode.vcf \
-F ID -F REF -F ALT -GF GT -O plink/"$RUN"/pastGT \
-R /scicore/home/ebertd/dexter0000/pasteuria_ref/Pramosa_C1_genome.gapfilled.13122016.fasta

#Rscript to format phenotype (parasite genotype) file. Too complicated for bash scripting.
module purge
module load R/4.0.0-foss-2018b
Rscript scripts/r_scripts/binary_pheno.R $RUN

################################################################################
#Covariates file
################################################################################
#NEW COVARIATES
/scicore/home/ebertd/dexter0000/plinkDev/plink2 --bfile plink/"$RUN"/pasteuria.geno.mind.king.maf.LD \
     --indep-pairwise 1kb 1 0.2 \
     --out plink/"$RUN"/pasteuria --allow-extra-chr

/scicore/home/ebertd/dexter0000/plinkDev/plink2 --bfile plink/"$RUN"/pasteuria.geno.mind.king.maf.LD --make-bed --out plink/"$RUN"/pasteuria.geno.mind.king.maf.covs --allow-extra-chr --extract plink/"$RUN"/pasteuria.prune.in

#Pasteuria PCA. Keep only samples present in daphnia PCA
/scicore/home/ebertd/dexter0000/plinkDev/plink2 --pca --bfile plink/"$RUN"/pasteuria.geno.mind.king.maf.covs --no-pheno --out plink/"$RUN"/pasteuria --allow-extra-chr

#Append label to pasteuria PCs so they don't duplicate daphnia PCs
sed -i \
-e '1s/PC1/PC1past/' \
-e '1s/PC2/PC2past/' \
-e '1s/PC3/PC3past/' \
-e '1s/PC4/PC4past/' \
-e '1s/PC5/PC5past/' \
-e '1s/PC6/PC6past/' \
-e '1s/PC7/PC7past/' \
-e '1s/PC8/PC8past/' \
-e '1s/PC9/PC9past/' \
-e '1s/PC10/PC10past/' \
plink/"$RUN"/pasteuria.eigenvec

#Make new daphnia eigenvectors
/scicore/home/ebertd/dexter0000/plinkDev/plink2 --bfile plink/"$RUN"/daphnia.geno.mind.king.maf.LD \
     --indep-pairwise 1000kb 1 0.2 \
     --out plink/"$RUN"/daphnia --allow-extra-chr

/scicore/home/ebertd/dexter0000/plinkDev/plink2 --bfile plink/"$RUN"/daphnia.geno.mind.king.maf.LD --extract plink/"$RUN"/daphnia.prune.in --make-bed --out plink/"$RUN"/daphnia.geno.mind.king.maf --allow-extra-chr

/scicore/home/ebertd/dexter0000/plinkDev/plink2 --pca --bfile plink/"$RUN"/daphnia.geno.mind.king.maf --no-pheno --out plink/"$RUN"/daphnia --allow-extra-chr

#Merge covariate files
paste plink/"$RUN"/pasteuria.eigenvec plink/"$RUN"/daphnia.eigenvec > plink/"$RUN"/combinedPC.eigenvec
cut -f13,14 --complement plink/"$RUN"/combinedPC.eigenvec > plink/"$RUN"/combinedPC.eigenvec

################################################################################
#Association analysis with additive logistic model
################################################################################

/scicore/home/ebertd/dexter0000/plinkDev/plink2 --bfile plink/"$RUN"/daphnia.geno.mind.king.maf.hwe --logistic firth-fallback hide-covar log10 --pheno plink/"$RUN"/binaryMerge.phe --out plink/"$RUN"/results/output --allow-extra-chr --no-pheno --pfilter 0.1e-5 --1 --covar plink/"$RUN"/combinedPC.eigenvec --covar-col-nums 13,14,15,16,17,18,19,20,21,22

#The output needs to be cleaned up a bit before R import.First append the locus name to each ignificant association
cd plink/"$RUN"/results

for i in output.chr1*; do awk '{print FILENAME"\t"$0}' $i > $i.bk; mv $i.bk $i; done

#Combine all the significant assocations into a single file
cat output.chr1* > results.txt

#Delete repeated headers from each file
sed '1!{/#CHROM/d;}' results.txt > results2.txt

#PLINK outputs contains a variable number of whitespaces between columns so redundant white spaces need to be stripped out before downstream #processing.
cat results2.txt | tr -s ' ' > results_"$RUN".txt
````



# R script for phenotype file

The above bash script calls this R script (binary_pheno.R) in order to generate and format the phenotype (parasite genotype) file.

````R
#Import arguments from shell
args <- commandArgs()
print(args[6])
RUN<-(args[6])

#Load packages
library(readr)
library(stringr)
library(reshape)
##Allele freqs histogram (optional for interative work)
#Requires objects from the main script to be loaded into the workspace
#mdata <- melt(data=refFreq[,c(1,4:ncol(refFreq))], id.vars = c("ID"))

#hist(mdata$value,breaks=30,xlim=c(0,1),xlab="Freq. of ref allele",ylab="Number of observations",main="Distribution of ref allele freqs.")

# Dosage phenotypes -------------------------------------------------------
#This calculates the number of reference alleles in the pool (max 2)
#This is scaled x 10 (eg, 0,10,20) otherwise plink will interpret the
#input as binary case/controls with 0 meaning missing data
fileGT<-paste("/scicore/home/ebertd/dexter0000/interlink/plink/",RUN,"/pastGT",sep="")
pastGT <- data.frame(read_delim(fileGT, "\t", escape_double = FALSE, trim_ws = TRUE))

#Convert to data frame
pastGT <- data.frame(pastGT)

#Substitute | for / in data.frame to make entries uniform
for (i in 1:ncol(pastGT)){
  pastGT[,i] <- sub("\\|", "/", pastGT[,i])
}

#Convert genotypes calls to the number of reference alleles (max 2)
refCount<-pastGT
for (j in 4:ncol(pastGT)){
  for (i in 1:length(pastGT$ID)) {
    refCount[i,j]<-sum(ifelse(unlist(strsplit(pastGT[i,j],split="/"))== pastGT$REF[i],10,0))
  }
}

#Annotate NA values with -9 for plink input
for (j in 4:ncol(refCount)){
  for (i in 1:length(refCount$ID)) {
    refCount[i,j]<-ifelse(pastGT[i,j] == "." | pastGT[i,j] == "./.", refCount[i,j] <- -9, refCount[i,j])
  }
}

#Recoding these as binary presence/absence variables
#########################################################################################
#Binary phenotypes
ref<-refCount
ref[ref==10 | ref==20]<-"1"
ref[ref==0 | ref==-9]<-"0"
ref$ID<-paste(ref$ID,"_R",sep="")

alt<-refCount
alt[alt==10 | alt==0]<-"1"
alt[alt==20 | alt==-9]<-"0"
alt$ID<-paste(alt$ID,"_A",sep="")

null<-refCount
null[null==10 | null==20]<-"0"
null[null=="-9"]<-"1"
null$ID<-paste(null$ID,"_N",sep="")

#Merge together different alleles
binary<-rbind(ref,alt,null)

#Filter out constant alleles (always present)
index<-4:ncol(binary)

binary[,index]<-data.frame(apply(binary[,index], 2, as.numeric))
binary<-binary[apply(binary[,index], 1, var, na.rm=TRUE) != 0,]

#Filter out alleles present in too low/high counts
freqs<-rowSums(binary[,-c(1:3)])/ncol(binary[,-c(1:3)])
binary<-binary[which(freqs >= 0.1 & freqs <= 0.9),]

#The final output needs to be formatted for PLINK input
binary <-binary[,-c(2,3)]
binary<-t(binary)
colnames(binary)<-binary[1,]
binary<-binary[-1,]
IID<-rownames(binary)
FID<-rep(0,length(IID))
binary<-cbind(FID,IID,binary)

binary[,2]<-str_replace(binary[,2], ".GT", "")

outFile<-paste("/scicore/home/ebertd/dexter0000/interlink/plink/",RUN,"/binaryMerge.phe",sep="")

write.table(binary,outFile,row.names = FALSE, sep = " ",quote =FALSE)
````





