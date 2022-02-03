# Pipeline step 1 - Read aggregation

This script splits raw sequencer reads according to sample ID. It merges across lanes. Left and right read files are created for each sample. It also runs BBTOOLS to make sure that read pairs are properly matched. Unpaired reads are redirected to a separate file. Requires an index file with the sequencer file names for the input and an index file with the corresponding desired sample names for the output.

````
#!/bin/bash

#SBATCH --job-name=read_renaming	#Job name
#SBATCH --cpus-per-task=8	        #Number of cores reserved
#SBATCH --mem-per-cpu=2G              	#Memory reserved per core.
					#Total memory reserved: 16GB
#SBATCH --time=24:00:00	        	#Maximum time the job will run
#SBATCH --qos=1day           		#The job queue (time based)

#This is the stdout file
#SBATCH --output=/scicore/home/ebertd/dexter0000/interlink/logfiles/readRenameOut_%A_%a

#This is the stderr file
#SBATCH --error=/scicore/home/ebertd/dexter0000/interlink/logfiles/readRenameErr_%A_%a

#Specifies an array of jobs from 1-8 with 8 max simultaneous
#SBATCH --array=1-258%50

#Things remember when running on sciCore:
#This job runs from the current working directory
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.


###############################################################################
# This bash script:
#	Aggregates together samples on the openBIS server that have the same ID
#	Renames them using a custom naming scheme
#	Makes sure read1 and read2 are ordered correctly
#	Splits single-end reads into a separate file
#
# Revised 23.09.2020
# author: Eric Dexter
################################################################################
#Define variables

#Interlink IDs list
INDEXFILE=/scicore/home/ebertd/dexter0000/interlink/scripts/index

#Open BIS IDs list
OPENBIS=/scicore/home/ebertd/dexter0000/interlink/scripts/openbisID

#Interlink ID for sample i
SAMP=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXFILE)

#Open BIS ID for sample i
BIS=$(sed -n ${SLURM_ARRAY_TASK_ID}p $OPENBIS)

#Load required modules
module load Java/8u212b03

################################################################################
cd /scicore/home/ebertd/dexter0000/interlink/reads/

echo "concatenating all reads named "$BIS" and renaming to "$SAMP""

#Find and concatenate
find /scicore/projects/openbis/userstore/duw_ebert/ -name "*${BIS}*R1*fastq.gz" \
-exec cat > /scicore/home/ebertd/dexter0000/interlink/reads/${SAMP}_R1.gz {} +

find /scicore/projects/openbis/userstore/duw_ebert/ -name "*${BIS}*R2*fastq.gz" \
-exec cat > /scicore/home/ebertd/dexter0000/interlink/reads/${SAMP}_R2.gz {} +

/scicore/home/ebertd/dexter0000/bbtools/bbmap/repair.sh in1=${SAMP}_R1.gz in2=${SAMP}_R2.gz out1=${SAMP}_paired_R1.gz out2=${SAMP}_paired_R2.gz outs=${SAMP}_unpaired.gz repair

rm ${SAMP}_R1.gz ${SAMP}_R2.gz
````
