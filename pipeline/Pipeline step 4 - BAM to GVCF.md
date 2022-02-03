# Pipeline step 4 - BAM to GVCF

This script creates a GVCF file by marking variants jointly across all BAM files using GATK HAPLOTYPE CALLER. Requires a reference genome and an index of input files, which should be the BAM files created in step 3 of the pipeline. The step is performed separately, but identically, for Pasteuria.

## Daphnia

````bash
#!/bin/bash

#SBATCH --job-name=BAM_2_GVCF_daphnia    		#Job name
#SBATCH --cpus-per-task=4                       #Number of cores reserved
#SBATCH --mem-per-cpu=8G                        #Memory reserved per core.
                                                #Total memory reserved: 32GB

#SBATCH --time=24:00:00                         #Maximum time the job will run
#SBATCH --qos=1day                              #The job queue (time based)

#This is the stdout file
#SBATCH --output=/scicore/home/ebertd/dexter0000/interlink/logfiles/daphGATKout_%A_%a

#This is the stderr file
#SBATCH --error=/scicore/home/ebertd/dexter0000/interlink/logfiles/daphGATKerr_%A_%a

#Specifies an array of jobs from 1-258
#SBATCH --array=1-258

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
INDEXFILE=/scicore/home/ebertd/dexter0000/interlink/scripts/index

#Sample i from index
SAMP=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXFILE)

################################################################################
#load required modules

module load GATK

################################################################################
#Navigate to starting directory

cd /scicore/home/ebertd/dexter0000/interlink

################################################################################
#Variant calling using GATK HAPLOTYPE CALLER to call SNPs jointly from many BAM
#files. GATK recommends the product of the options -nct and -nt should approx. 
#equal the number of processors requested.

gatk --java-options "-Xmx32G" HaplotypeCaller \
-R "$REF" \
-I bamsDaphnia/"$SAMP"_daphnia_Rm_rdup_indelrealigner.bam \
-ploidy 2 \
-O gvcfsDaphnia/"$SAMP".vcf \
-ERC GVCF

echo "daphnia gvcf finished for "$SAMP""
````



## Pasteuria

````bash
#!/bin/bash

#SBATCH --job-name=BAM_2_GVCF_pasteuria    		#Job name
#SBATCH --cpus-per-task=4                       #Number of cores reserved
#SBATCH --mem-per-cpu=8G                        #Memory reserved per core.
                                                #Total memory reserved: 32GB

#SBATCH --time=24:00:00                         #Maximum time the job will run
#SBATCH --qos=1day                              #The job queue (time based)

#This is the stdout file
#SBATCH --output=/scicore/home/ebertd/dexter0000/interlink/logfiles/pastGATKout_%A_%a

#This is the stderr file
#SBATCH --error=/scicore/home/ebertd/dexter0000/interlink/logfiles/pastGATKerr_%A_%a

#Specifies an array of jobs from 1-258
#SBATCH --array=1-258

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
INDEXFILE=/scicore/home/ebertd/dexter0000/interlink/scripts/index

#Sample i from index
SAMP=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXFILE)

################################################################################
#load required modules

module load GATK

################################################################################
#Navigate to starting directory

cd /scicore/home/ebertd/dexter0000/interlink

################################################################################
#Variant calling using GATK HAPLOTYPE CALLER to call SNPs jointly from many BAM
#files. GATK recommends the product of the options -nct and -nt should approx. 
#equal the number of processors requested.

gatk --java-options "-Xmx32G" HaplotypeCaller \
-R "$REF" \
-I bamsPasteuria/"$SAMP"_pasteuria_Rm_rdup_indelrealigner.bam \
-ploidy 2 \
-O gvcfsPasteuria/"$SAMP".vcf \
-ERC GVCF
````