# Pipeline step 2 - Read trimming

This script trims Illumina adaptors and performs light quality trimming using TRIMMOMATIC. It requires an index of the input files to be trimmed, which should be the output from the previous step of the pipeline.  

````bash
#!/bin/bash

#SBATCH --job-name=read_trimming				#Job name
#SBATCH --cpus-per-task=8	                  	#Number of cores reserved
#SBATCH --mem-per-cpu=2G              			#Memory reserved per core.
												#Total memory reserved: 16GB
#SBATCH --time=24:00:00	        				#Maximum time the job will run
#SBATCH --qos=1day           					#The job queue (time based)

#This is the stdout file
#SBATCH --output=/scicore/home/ebertd/dexter0000/interlink/logfiles/readTrim_%A_%a

#This is the stderr file
#SBATCH --error=/scicore/home/ebertd/dexter0000/interlink/logfiles/readTrim_%A_%a

#Specifies an array of jobs from 1-n with n max simultaneous
#SBATCH --array=1-258%25

#Things remember when running on sciCore:
#This job runs from the current working directory
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.

########################################################################
#Define variables

#Interlink IDs list
INDEXFILE=/scicore/home/ebertd/dexter0000/interlink/scripts/index

#Interlink ID for sample i
SAMP=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXFILE)

#Load required modules
module load Trimmomatic/0.39-Java-1.8

################################################################################
cd /scicore/home/ebertd/dexter0000/interlink/reads/

echo "trimming paired reads from "$SAMP""

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
-threads 8 \
${SAMP}_paired_R1.gz ${SAMP}_paired_R2.gz \
${SAMP}_paired_R1_trimmed.fq.gz ${SAMP}_unpaired_R1_trimmed.fq.gz \
${SAMP}_paired_R2_trimmed.fq.gz ${SAMP}_unpaired_R2_trimmed.fq.gz \
ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 \
SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36
````