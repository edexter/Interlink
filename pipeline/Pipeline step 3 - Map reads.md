# Pipeline step 3 - Map reads to reference

This script maps the trimmed reads to a reference genome using BWA-MEM. It then converts the SAM output to BAM format, calculates mapping statistics with SAMTOOLS, appends read group data to the alignments, removes duplicate reads, performs indel realignment with GATK, removes unmapped reads, and indexes the BAM. Requires a reference genome and index of input files, which should be the output from step 2 of the pipeline. This entire process is repeated separately, but identically, relative to the pasteuria genome.

````bash
#!/bin/bash

#SBATCH --job-name=bam_preparation_daphnia		#Job name
#SBATCH --cpus-per-task=8	                  	#Number of cores reserved
#SBATCH --mem-per-cpu=4G              			#Memory reserved per core.
												#Total memory reserved: 32GB
#SBATCH --time=24:00:00	        			#Maximum time the job will run
#SBATCH --qos=1day	           					#The job queue (time based)

#This is the stdout file
#SBATCH --output=/scicore/home/ebertd/dexter0000/interlink/logfiles/daphFq2BamOut_%A_%a

#This is the stderr file
#SBATCH --error=/scicore/home/ebertd/dexter0000/interlink/logfiles/daphFq2BamErr_%A_%a

#Specifies an array of jobs from 1-8 with 8 max simultaneous
#SBATCH --array=1-258%25

#Things remember when running on sciCore:
#This job runs from the current working directory
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.
#Make sure reference is indexed before running (only needs to be done once)
#bwa index "$REF"

###############################################################################
#Define variables

#Sample index
INDEXFILE=/scicore/home/ebertd/dexter0000/interlink/scripts/index

#Sample i from index
SAMP=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXFILE)

REF="/scicore/home/ebertd/dexter0000/daphnia_ref/Xinb3_ref.fasta"

################################################################################

#load required modules
#Due to some workarounds with java licensing issues, GATK must be loaded first
module load GATK/3.7-0-Java-1.8.0_92
module load bwakit
module load SAMtools/1.9-foss-2018b
module load picard

################################################################################
#Part 1

cd /scicore/home/ebertd/dexter0000/interlink/bamsDaphnia

#Verify that the reference genome has already been indexed
if [ -s "$REF".fai ]; then
        echo "The reference genome has already been indexed"
else
        echo "ERROR: Index the reference genome before running the script"
        
fi

#Start a timer to record how long the pipeline takes (wall time)
start_time=`date +%s`

#Map reads to reference genome, convert SAM output to BAM, and sort BAM.
#It's best to combine these steps to avoid writing large intermediate files
#to disk
echo "mapping reads of "$SAMP""
bwa mem -t 8 -M "$REF" /scicore/home/ebertd/dexter0000/interlink/reads/"$SAMP"_paired_R1_trimmed.fq.gz \
/scicore/home/ebertd/dexter0000/interlink/reads/"$SAMP"_paired_R2_trimmed.fq.gz \
| samtools sort -@ 8 -o "$SAMP"_daphnia_sorted.bam

# print mapping statistics to screen
echo "Mapping stats for "$SAMP"_daphnia are:"
samtools flagstat "$SAMP"_daphnia_sorted.bam | tee "$SAMP"_daphnia_flagstat.txt

#Append read group metadata
echo "Adding readgroup to "$SAMP"_daphnia_sorted.bam"
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
I="$SAMP"_daphnia_sorted.bam \
O="$SAMP"_daphnia_R.bam \
RGID="$SAMP"  \
RGLB="$SAMP"  \
RGPL=illumina  \
RGPU="$SAMP"  \
RGSM="$SAMP"

#Index BAM
echo ""$SAMP"_daphnia_R.bam file is being indexed"
samtools index -@ 8 -b "$SAMP"_daphnia_R.bam
	
#remove duplicates
echo "duplicates are being removed"
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
REMOVE_DUPLICATES=TRUE \
I="$SAMP"_daphnia_R.bam \
O="$SAMP"_daphnia_R_rdup.bam \
M="$SAMP"_daphnia_duplicate_metrics.txt
	
#Index new BAM
echo ""$SAMP"_daphnia_R.bam file is being indexed"
samtools index -@ 8 -b "$SAMP"_daphnia_R_rdup.bam

#remove intermediate files
if [ -s "$SAMP"_daphnia_R_rdup.bam ]; then
	echo "removing intermediate files. ALMOST DONE!"
	rm "$SAMP"_daphnia_sorted.bam
	rm "$SAMP"_daphnia_R.bam
	rm "$SAMP"_daphnia_R.bam.bai
else 
	echo "Something went wrong! the sorted bam file has size 0" 
fi

################################################################################
#Part 2: Prepre BAM files for variant calling
################################################################################

#Create sequence dictionary for reference genome (needed for indel realignment)
#Only needs to be performed once for the reference
if [ -s "$REF".dict ]; then
	echo "The reference genome dictonary already exists"
else
	echo "Building the reference dictionary"
	java -jar $EBROOTPICARD/picard.jar  CreateSequenceDictionary \
	REFERENCE= "$REF" \
	OUTPUT= /scicore/home/ebertd/dexter0000/daphnia_ref/Xinb3_ref.dict
fi

#Create preliminary set of known sites for indel realignment
echo "preliminary vcf for indel realignment is being created"
java -Xmx32G -jar $EBROOTGATK/GenomeAnalysisTK.jar -nt 8 \
	-T UnifiedGenotyper \
	-R "$REF"    \
        -I "$SAMP"_daphnia_R_rdup.bam \
        -glm INDEL \
        -o "$SAMP"_INDEL.vcf \

	#Create targets for indel realignment
echo "Creating target intervals list using RealignerTargetCreator"
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -nt 8 \
	-T RealignerTargetCreator \
	-R "$REF" \
	-I "$SAMP"_daphnia_R_rdup.bam \
	-o "$SAMP"_daphnia_realignertargetcreator.intervals \
	-known "$SAMP"_INDEL.vcf

#This step could be multithreaded. This is a major bottleneck, also I might need
#to request more RAM
echo "Realigning reads using IndelRealigner"
java -Xmx32G -jar $EBROOTGATK/GenomeAnalysisTK.jar \
	-T IndelRealigner \
	-R "$REF"    \
	-targetIntervals "$SAMP"_daphnia_realignertargetcreator.intervals \
	-I "$SAMP"_daphnia_R_rdup.bam \
	-o "$SAMP"_daphnia_R_rdup_indelrealigner.bam \
	-known "$SAMP"_INDEL.vcf

#Remove unmapped reads from bam
echo "removing unmapped reads of "$SAMP" bam file"
samtools view -b -F 4 -@ 8 "$SAMP"_daphnia_R_rdup_indelrealigner.bam \
> "$SAMP"_daphnia_Rm_rdup_indelrealigner.bam

#Index new bam
echo "indexing new "$SAMP" bam file"
samtools index -@ 8 "$SAMP"_daphnia_Rm_rdup_indelrealigner.bam

#remove intermediate files
if [ -s "$SAMP"_daphnia_Rm_rdup_indelrealigner.bam ]; then
	echo "removing intermediate bam files."
	rm "$SAMP"_daphnia_R_rdup_indelrealigner.bam
	rm "$SAMP"_daphnia_R_rdup_indelrealigner.bai
	rm "$SAMP"_daphnia_R_rdup.bam
	rm "$SAMP"_daphnia_R_rdup.bam.bai
	rm "$SAMP"_daphnia_realignertargetcreator.intervals
	rm "$SAMP"_INDEL.vcf
	rm "$SAMP"_INDEL.vcf.idx
else
	echo "ERROR: the sorted bam file has size 0"
fi
````



## Pasteuria

````bash
#!/bin/bash

#SBATCH --job-name=bam_preparation_pasteuria	#Job name
#SBATCH --cpus-per-task=16	                  	#Number of cores reserved
#SBATCH --mem-per-cpu=2G              			#Memory reserved per core.
												#Total memory reserved: 32GB
#SBATCH --time=24:00:00	        				#Maximum time the job will run
#SBATCH --qos=1day	           					#The job queue (time based)

#This is the stdout file
#SBATCH --output=/scicore/home/ebertd/dexter0000/interlink/logfiles/pastFq2BamOut_%A_%a

#This is the stderr file
#SBATCH --error=/scicore/home/ebertd/dexter0000/interlink/logfiles/pastFq2BamErr_%A_%a

#Specifies an array of jobs from 1-128 with 50 max simultaneous
#SBATCH --array=1-258%50

#Things remember when running on sciCore:
#This job runs from the current working directory
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.
#Make sure reference is indexed before running (only needs to be done once)
#bwa index "$REF"

################################################################################
#Define variables

#Sample index
INDEXFILE=/scicore/home/ebertd/dexter0000/interlink/scripts/index

#Sample i from index
SAMP=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXFILE)

REF="/scicore/home/ebertd/dexter0000/pasteuria_ref/Pramosa_C1_genome.gapfilled.13122016.fasta"

################################################################################

#load required modules
#Due to some workarounds with java licensing issues, GATK must be loaded first
module load GATK/3.7-0-Java-1.8.0_92
module load bwakit
module load SAMtools/1.9-foss-2018b
module load picard

################################################################################
#Part 1

cd /scicore/home/ebertd/dexter0000/interlink/bamsPasteuria

#Verify that the reference genome has already been indexed
if [ -s "$REF".fai ]; then
        echo "The reference genome has already been indexed"
else
        echo "ERROR: Index the reference genome before running the script"
        
fi

#Start a timer to record how long the pipeline takes (wall time)
start_time=`date +%s`

#Map reads to reference genome, convert SAM output to BAM, and sort BAM.
#It's best to combine these steps to avoid writing large intermediate files
#to disk
echo "mapping reads of "$SAMP""
bwa mem -t 16 -M "$REF" /scicore/home/ebertd/dexter0000/interlink/reads/"$SAMP"_paired_R1_trimmed.fq.gz \
/scicore/home/ebertd/dexter0000/interlink/reads/"$SAMP"_paired_R2_trimmed.fq.gz \
| samtools sort -@ 16 -o "$SAMP"_pasteuria_sorted.bam

# print mapping statistics to screen
echo "Mapping stats for "$SAMP"_pasteuria are:"
samtools flagstat "$SAMP"_pasteuria_sorted.bam | tee "$SAMP"_pasteuria_flagstat.txt

#Append read group metadata
echo "Adding readgroup to "$SAMP"_pasteuria_sorted.bam"
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
I="$SAMP"_pasteuria_sorted.bam \
O="$SAMP"_pasteuria_R.bam \
RGID="$SAMP"  \
RGLB="$SAMP"  \
RGPL=illumina  \
RGPU="$SAMP"  \
RGSM="$SAMP"

#Index BAM
echo ""$SAMP"_pasteuria_R.bam file is being indexed"
samtools index -@ 16 -b "$SAMP"_pasteuria_R.bam
	
#remove duplicates
echo "duplicates are being removed"
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
REMOVE_DUPLICATES=TRUE \
I="$SAMP"_pasteuria_R.bam \
O="$SAMP"_pasteuria_R_rdup.bam \
M="$SAMP"_pasteuria_duplicate_metrics.txt
	
#Index new BAM
echo ""$SAMP"_pasteuria_R.bam file is being indexed"
samtools index -@ 16 -b "$SAMP"_pasteuria_R_rdup.bam

#remove intermediate files
if [ -s "$SAMP"_pasteuria_R_rdup.bam ]; then
	echo "removing intermediate files. ALMOST DONE!"
	rm "$SAMP"_pasteuria_sorted.bam
	rm "$SAMP"_pasteuria_R.bam
	rm "$SAMP"_pasteuria_R.bam.bai
else 
	echo "Something went wrong! the sorted bam file has size 0" 
fi

################################################################################
#Part 2: Prepare BAM files for variant calling
################################################################################

#Create sequence dictionary for reference genome (needed for indel realignment)
#Only needs to be performed once for the reference
if [ -s "$REF".dict ]; then
	echo "The reference genome dictonary already exists"
else
	echo "Building the reference dictionary"
	java -jar $EBROOTPICARD/picard.jar  CreateSequenceDictionary \
	REFERENCE= "$REF" \
	OUTPUT= /scicore/home/ebertd/dexter0000/pasteuria_ref/Xinb3_ref.dict
fi

#Create preliminary set of known sites for indel realignment
echo "preliminary vcf for indel realignment is being created"
java -Xmx32G -jar $EBROOTGATK/GenomeAnalysisTK.jar -nt 16 \
	-T UnifiedGenotyper \
	-R "$REF"    \
        -I "$SAMP"_pasteuria_R_rdup.bam \
        -glm INDEL \
        -o "$SAMP"_INDEL.vcf \

#Create targets for indel realignment
echo "Creating target intervals list using RealignerTargetCreator"
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -nt 16 \
	-T RealignerTargetCreator \
	-R "$REF" \
	-I "$SAMP"_pasteuria_R_rdup.bam \
	-o "$SAMP"_pasteuria_realignertargetcreator.intervals \
	-known "$SAMP"_INDEL.vcf

#This step could be multithreaded. This is a major bottleneck, also I might need
#to request more RAM
echo "Realigning reads using IndelRealigner"
java -Xmx32G -jar $EBROOTGATK/GenomeAnalysisTK.jar \
	-T IndelRealigner \
	-R "$REF"    \
	-targetIntervals "$SAMP"_pasteuria_realignertargetcreator.intervals \
	-I "$SAMP"_pasteuria_R_rdup.bam \
	-o "$SAMP"_pasteuria_R_rdup_indelrealigner.bam \
	-known "$SAMP"_INDEL.vcf

#Remove unmapped reads from bam
echo "removing unmapped reads of "$SAMP" bam file"
samtools view -b -F 4 -@ 16 "$SAMP"_pasteuria_R_rdup_indelrealigner.bam \
> "$SAMP"_pasteuria_Rm_rdup_indelrealigner.bam

#Index new bam
echo "indexing new "$SAMP" bam file"
samtools index -@ 16 "$SAMP"_pasteuria_Rm_rdup_indelrealigner.bam

#remove intermediate files
if [ -s "$SAMP"_pasteuria_Rm_rdup_indelrealigner.bam ]; then
	echo "removing intermediate bam files."
	rm "$SAMP"_pasteuria_R_rdup_indelrealigner.bam
	rm "$SAMP"_pasteuria_R_rdup_indelrealigner.bai
	rm "$SAMP"_pasteuria_R_rdup.bam
	rm "$SAMP"_pasteuria_R_rdup.bam.bai
	rm "$SAMP"_pasteuria_realignertargetcreator.intervals
	rm "$SAMP"_INDEL.vcf
	rm "$SAMP"_INDEL.vcf.idx
else
	echo "ERROR: the sorted bam file has size 0"
fi
````

