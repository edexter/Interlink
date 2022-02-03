# Pipeline step 5 - GVCF to VCF

This script calls variants from a GVCF file to produce a VCF output. Variants are quality filtered separately by type (SNP, INDEL, MIXED) and the filtered variants are merged back together. A reference genome and index of input GVCF file names, which should be the output from step 4, is required. Summary stats are calculated for each step of filtering. This step is performed separately, but identically, for daphnia and pasteuria.

## Daphnia

````bash
#!/bin/bash

#SBATCH --job-name=GVCF2VCF_daph			    #Job name
#SBATCH --cpus-per-task=4                       #Number of cores reserved
#SBATCH --mem-per-cpu=32G                       #Memory reserved per core.
                                                #Total memory reserved: 128GB

#SBATCH --time=168:00:00                         #Maximum time the job will run
#SBATCH --qos=1week                              #The job queue (time based)

#This is the stdout file
#SBATCH --output=/scicore/home/ebertd/dexter0000/interlink/logfiles/GVCF2VCF_daph_out_%A

#This is the stderr file
#SBATCH --error=/scicore/home/ebertd/dexter0000/interlink/logfiles/GVCF2VCF_daph_err_%A

#Things remember when running on sciCore:
#This job runs from the current working directory
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.

################################################################################
#Define variables

#Reference genome
REF="/scicore/home/ebertd/dexter0000/daphnia_ref/Xinb3_ref.fasta"

#Sample index
INDEXFILE=/scicore/home/ebertd/dexter0000/interlink/scripts/gvcf.list.daphnia

################################################################################
#load required modules

module load GATK
module load picard

################################################################################
#Navigate to starting directory

cd /scicore/home/ebertd/dexter0000/interlink

################################################################################
#Merge GVCFs (less memory may cause crashes)
gatk --java-options "-Xmx64G" CombineGVCFs \
   -R "$REF" \
   -V scripts/gvcf.daphnia.list \
   -O gvcfsDaphnia/daphCohort.vcf.gz
   
# Call variants
gatk --java-options "-Xmx64G" GenotypeGVCFs \
    -R "$REF" \
    -V gvcfsDaphnia/daphCohort.vcf.gz \
    -O vcfsDaphnia/daphnia_unfiltered.vcf

################################################################################
#Variant filtering using GATK best practices for UnifiedGenotyper. First split 
#VCF into INDEL-only, SNP-only, and MIXED-only (SNP+INDEL) VCFs. The reason for 
#this is that different filters should be applied to INDELS and SNPs. MIXED 
#records get the same filters and INDEL only. These three files are recombined 
#after filtering.

gatk SelectVariants \
    -V vcfsDaphnia/daphnia_unfiltered.vcf \
    -select-type SNP \
    -O vcfsDaphnia/merged_daphnia_SNP.vcf

gatk SelectVariants \
    -V vcfsDaphnia/daphnia_unfiltered.vcf \
    -select-type INDEL \
    -O vcfsDaphnia/merged_daphnia_INDEL.vcf

gatk SelectVariants \
    -V vcfsDaphnia/daphnia_unfiltered.vcf \
    -select-type MIXED \
    -O vcfsDaphnia/merged_daphnia_MIXED.vcf

#Filter SNP-only VCF
gatk VariantFiltration \
    -V vcfsDaphnia/merged_daphnia_SNP.vcf \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
	-O vcfsDaphnia/merged_daphnia_SNP_filtered.vcf

#Remove SNPs that fail filter
gatk SelectVariants \
    -V vcfsDaphnia/merged_daphnia_SNP_filtered.vcf \
    --exclude-filtered \
    -O vcfsDaphnia/merged_daphnia_SNP_filtered_purged.vcf

#Filter INDEL-only VCF
gatk VariantFiltration \
	-V vcfsDaphnia/merged_daphnia_INDEL.vcf \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O vcfsDaphnia/merged_daphnia_INDEL_filtered.vcf

#Remove INDELS that fail filter
gatk SelectVariants \
    -V vcfsDaphnia/merged_daphnia_INDEL_filtered.vcf \
    --exclude-filtered \
    -O vcfsDaphnia/merged_daphnia_INDEL_filtered_purged.vcf

#Filter MIXED-only VCF
gatk VariantFiltration \
	-V vcfsDaphnia/merged_daphnia_MIXED.vcf \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O vcfsDaphnia/merged_daphnia_MIXED_filtered.vcf

#Remove MIXED that fail filter
gatk SelectVariants \
    -V vcfsDaphnia/merged_daphnia_MIXED_filtered.vcf \
    --exclude-filtered \
    -O vcfsDaphnia/merged_daphnia_MIXED_filtered_purged.vcf

#Merge the filtered INDEL and SNP VCF files
java -jar $EBROOTPICARD/picard.jar MergeVcfs \
	I=vcfsDaphnia/merged_daphnia_SNP_filtered_purged.vcf \
	I=vcfsDaphnia/merged_daphnia_INDEL_filtered_purged.vcf \
	I=vcfsDaphnia/merged_daphnia_MIXED_filtered_purged.vcf \
	O=vcfsDaphnia/merged_daphnia_filtered.vcf


###############################################################################
#Calculate pre and post filtering stats
gatk CountVariants -V vcfsDaphnia/daphnia_unfiltered.vcf > vcfsDaphnia/ALLpre.txt
gatk CountVariants -V vcfsDaphnia/merged_daphnia_filtered.vcf > vcfsDaphnia/ALLpost.txt
gatk CountVariants -V vcfsDaphnia/merged_daphnia_SNP.vcf > vcfsDaphnia/SNPpre.txt
gatk CountVariants -V vcfsDaphnia/merged_daphnia_SNP_filtered_purged.vcf > vcfsDaphnia/SNPpost.txt
gatk CountVariants -V vcfsDaphnia/merged_daphnia_INDEL.vcf > vcfsDaphnia/INDELpre.txt
gatk CountVariants -V vcfsDaphnia/merged_daphnia_INDEL_filtered_purged.vcf > vcfsDaphnia/INDELpost.txt
gatk CountVariants -V vcfsDaphnia/merged_daphnia_MIXED.vcf > vcfsDaphnia/MIXEDpre.txt
gatk CountVariants -V vcfsDaphnia/merged_daphnia_MIXED_filtered_purged.vcf > vcfsDaphnia/MIXEDpost.txt
````



## Pasteuria

````
#!/bin/bash

#SBATCH --job-name=GVCF2VCF_past			    #Job name
#SBATCH --cpus-per-task=4                       #Number of cores reserved
#SBATCH --mem-per-cpu=8G                        #Memory reserved per core.
                                                #Total memory reserved: 32GB

#SBATCH --time=24:00:00                         #Maximum time the job will run
#SBATCH --qos=1day                              #The job queue (time based)

#This is the stdout file
#SBATCH --output=/scicore/home/ebertd/dexter0000/interlink/logfiles/GVCF2VCF_past

#This is the stderr file
#SBATCH --error=/scicore/home/ebertd/dexter0000/interlink/logfiles/GVCF2VCF_past

#Things remember when running on sciCore:
#This job runs from the current working directory
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.

################################################################################
#Define variables

#Reference genome
REF="/scicore/home/ebertd/dexter0000/pasteuria_ref/Pramosa_C1_genome.gapfilled.13122016.fasta"

#Sample index
INDEXFILE=/scicore/home/ebertd/dexter0000/interlink/scripts/gvcf.list.pasteuria

################################################################################
#load required modules

module load GATK
module load picard

################################################################################
#Navigate to starting directory

cd /scicore/home/ebertd/dexter0000/interlink

################################################################################
#Merge GVCFs
gatk CombineGVCFs \
   -R "$REF" \
   --variant scripts/gvcf.pasteuria.list \
   -O gvcfsPasteuria/pastCohort.vcf.gz
   
# Call variants
gatk GenotypeGVCFs \
    -R "$REF" \
    -V gvcfsPasteuria/pastCohort.vcf.gz \
    -O vcfsPasteuria/pasteuria_unfiltered.vcf

################################################################################
#Variant filtering using GATK best practices for UnifiedGenotyper. First split 
#VCF into INDEL-only, SNP-only, and MIXED-only (SNP+INDEL) VCFs. The reason for 
#this is that different filters should be applied to INDELS and SNPs. MIXED 
#records get the same filters and INDEL only. These three files are recombined 
#after filtering.

gatk SelectVariants \
    -V vcfsPasteuria/pasteuria_unfiltered.vcf \
    -select-type SNP \
    -O vcfsPasteuria/merged_pasteuria_SNP.vcf

gatk SelectVariants \
    -V vcfsPasteuria/pasteuria_unfiltered.vcf \
    -select-type INDEL \
    -O vcfsPasteuria/merged_pasteuria_INDEL.vcf

gatk SelectVariants \
    -V vcfsPasteuria/pasteuria_unfiltered.vcf \
    -select-type MIXED \
    -O vcfsPasteuria/merged_pasteuria_MIXED.vcf

#Filter SNP-only VCF
gatk VariantFiltration \
    -V vcfsPasteuria/merged_pasteuria_SNP.vcf \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
	-O vcfsPasteuria/merged_pasteuria_SNP_filtered.vcf

#Remove SNPs that fail filter
gatk SelectVariants \
    -V vcfsPasteuria/merged_pasteuria_SNP_filtered.vcf \
    --exclude-filtered \
    -O vcfsPasteuria/merged_pasteuria_SNP_filtered_purged.vcf

#Filter INDEL-only VCF
gatk VariantFiltration \
	-V vcfsPasteuria/merged_pasteuria_INDEL.vcf \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O vcfsPasteuria/merged_pasteuria_INDEL_filtered.vcf

#Remove INDELS that fail filter
gatk SelectVariants \
    -V vcfsPasteuria/merged_pasteuria_INDEL_filtered.vcf \
    --exclude-filtered \
    -O vcfsPasteuria/merged_pasteuria_INDEL_filtered_purged.vcf

#Filter MIXED-only VCF
gatk VariantFiltration \
	-V vcfsPasteuria/merged_pasteuria_MIXED.vcf \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O vcfsPasteuria/merged_pasteuria_MIXED_filtered.vcf

#Remove MIXED that fail filter
gatk SelectVariants \
    -V vcfsPasteuria/merged_pasteuria_MIXED_filtered.vcf \
    --exclude-filtered \
    -O vcfsPasteuria/merged_pasteuria_MIXED_filtered_purged.vcf

#Merge the filtered INDEL and SNP VCF files
java -jar $EBROOTPICARD/picard.jar MergeVcfs \
	I=vcfsPasteuria/merged_pasteuria_SNP_filtered_purged.vcf \
	I=vcfsPasteuria/merged_pasteuria_INDEL_filtered_purged.vcf \
	I=vcfsPasteuria/merged_pasteuria_MIXED_filtered_purged.vcf \
	O=vcfsPasteuria/merged_pasteuria_filtered.vcf


###############################################################################
#Calculate pre and post filtering stats
gatk CountVariants -V vcfsPasteuria/pasteuria_unfiltered.vcf > vcfsPasteuria/ALLpre.txt
gatk CountVariants -V vcfsPasteuria/merged_pasteuria_filtered.vcf > vcfsPasteuria/ALLpost.txt
gatk CountVariants -V vcfsPasteuria/merged_pasteuria_SNP.vcf > vcfsPasteuria/SNPpre.txt
gatk CountVariants -V vcfsPasteuria/merged_pasteuria_SNP_filtered_purged.vcf > vcfsPasteuria/SNPpost.txt
gatk CountVariants -V vcfsPasteuria/merged_pasteuria_INDEL.vcf > vcfsPasteuria/INDELpre.txt
gatk CountVariants -V vcfsPasteuria/merged_pasteuria_INDEL_filtered_purged.vcf > vcfsPasteuria/INDELpost.txt
gatk CountVariants -V vcfsPasteuria/merged_pasteuria_MIXED.vcf > vcfsPasteuria/MIXEDpre.txt
gatk CountVariants -V vcfsPasteuria/merged_pasteuria_MIXED_filtered_purged.vcf > vcfsPasteuria/MIXEDpost.txt
````