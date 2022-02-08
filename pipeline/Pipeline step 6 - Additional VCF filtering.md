# Additional VCF filtering

Before further analyses, a number of additional filtering and formatting steps should be undertaken. This are very quick and can be performed interactively.



## Create unique variant IDs

Downstream analysis will require that each variant has a unique name. Using BCFtools each variant can be named based on the chromosome number and position.

````bash
srun --nodes=1 --cpus-per-task=4 --mem=16G --pty bash

module load BCFtools/1.10.2-foss-2018b

bcftools annotate --set-id +'%CHROM\_%POS' vcfsPasteuria/merged_pasteuria_filtered.vcf \

> vcfsPasteuria/merged_pasteuria_filtered_annotated.vcf

bcftools annotate --set-id +'%CHROM\_%POS' vcfsDaphnia/merged_daphnia_filtered.vcf \

> vcfsDaphnia/merged_daphnia_filtered_annotated.vcf
````



## Filter sites with excessive depth

Sites with excessive depth may represent mapping problems and should be filtered out. Excessive is a relative term, but a conservative rule of thumb is more than twice the mean read depth across all sites is potentially problematic. GATK VariantsToTable can be used calculate the total depth all each site or the mean depth at each site. Examination of the data shows that both metrics show near identical patterns. Mean depth allows for taking missing data into account and is a little more intuitive than total depth. It does require a few extra calculation in R, but nothing difficult. 

````bash
#Request an interactive compute node
srun --nodes=1 --cpus-per-task=4 --mem=16G --pty bash

module load GATK

gatk VariantsToTable \
     -V vcfsPasteuria/merged_pasteuria_filtered_annotated.vcf \
     -F CHROM -F POS -F DP -F AN -GF DP\
     -O stats/variant.depth.pasteuria

gatk VariantsToTable \
     -V vcfsDaphnia/merged_daphnia_filtered_annotated.vcf \
     -F CHROM -F POS -F DP -F AN -GF DP\
     -O stats/variant.depth.daphnia
````

The following R script uses these output files to calculate cutoff values defining excessive depth for both the host and the parasite.

````
#Request an interactive compute node
srun --nodes=1 --cpus-per-task=4 --mem=16G --pty bash

module load GATK

gatk VariantsToTable \
     -V vcfsPasteuria/merged_pasteuria_filtered_annotated.vcf \
     -F CHROM -F POS -F DP -F AN -GF DP\
     -O stats/variant.depth.pasteuria

gatk VariantsToTable \
     -V vcfsDaphnia/merged_daphnia_filtered_annotated.vcf \
     -F CHROM -F POS -F DP -F AN -GF DP\
     -O stats/variant.depth.daphnia
````

This cutoff is then passed to GATK for filtering. Note that the values here are entirely specific to a given dataset.

````
#Request an interactive compute node
srun --nodes=1 --cpus-per-task=4 --mem=16G --pty bash srun --nodes=1 --cpus-per-task=4 --mem=16G --pty bash 

module load GATK

#Filter variants with > 2X the mean depth
gatk VariantFiltration \
    -V vcfsDaphnia/merged_daphnia_filtered_annotated.vcf \
    -filter "DP > 14540" --filter-name "DP_2*mean" \
	-O vcfsDaphnia/merged_daphnia_filtered_annotated_maxDP.vcf

gatk VariantFiltration \
    -V vcfsPasteuria/merged_pasteuria_filtered_annotated.vcf \
    -filter "DP > 132890" --filter-name "DP_2*mean" \
	-O vcfsPasteuria/merged_pasteuria_filtered_annotated_maxDP.vcf

#Remove variants that fail filter
gatk SelectVariants \
    -V vcfsPasteuria/merged_pasteuria_filtered_annotated_maxDP.vcf \
    --exclude-filtered \
    -O vcfsPasteuria/merged_pasteuria_filtered_annotated_maxDP_purged.vcf

gatk SelectVariants \
    -V vcfsDaphnia/merged_daphnia_filtered_annotated_maxDP.vcf \
    --exclude-filtered \
    -O vcfsDaphnia/merged_daphnia_filtered_annotated_maxDP_purged.vcf
````



## Filter dual mapping regions

Approximately 0.1% of mapped reads aligned to both host and parasite genomes. If there are variants at these locations, this could be a source of false positives in the association test. This module identifies genomic ranges in the host and parasite where reads consistently dual-align and then filters variants that are located in those ranges. The process involves several blocks of code in both bash and R script.

````bash
#!/bin/bash

#SBATCH --job-name=dual_map            #Job name
#SBATCH --cpus-per-task=1              #Number of cores reserved
#SBATCH --mem-per-cpu=16G              #Memory reserved per core.
                                       #Total memory reserved: 16GB
#SBATCH --time=6:00:00                 #Maximum time the job will run
#SBATCH --qos=6hours                   #The job queue (time based)

#This is the stdout file
#SBATCH --output=/scicore/home/ebertd/dexter0000/interlink/logfiles/dual_map_%A_%a

#This is the stderr file
#SBATCH --error=/scicore/home/ebertd/dexter0000/interlink/logfiles/dual_map_%A_%a

#Specifies an array of jobs from 1-258
#SBATCH --array=1-258

#Things remember when running on sciCore:
#This job runs from the current working directory
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.

###############################################################################
#Define variables

#Sample index
INDEXFILE=/scicore/home/ebertd/dexter0000/interlink/scripts/index

#Sample i from index
SAMP=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXFILE)

################################################################################
#load required modules

module load SAMtools/1.10-foss-2018b
module load picard

################################################################################
#Navigate to starting directory

cd /scicore/home/ebertd/dexter0000/interlink

################################################################################

#Get read names from each pair of BAMs and find the union of the two lists
samtools view bamsDaphnia/"$SAMP"_daphnia_Rm_rdup_indelrealigner.bam | cut -f1 | sort | uniq > bamsDaphnia/readNames/"$SAMP"_read_names_daph.txt

samtools view bamsPasteuria/"$SAMP"_pasteuria_Rm_rdup_indelrealigner.bam | cut -f1 | sort | uniq > bamsPasteuria/readNames/"$SAMP"_read_names_past.txt

comm -12 <(sort bamsDaphnia/readNames/"$SAMP"_read_names_daph.txt) <(sort bamsPasteuria/readNames/"$SAMP"_read_names_past.txt) > dualMapped/"$SAMP"_dual_mapped.txt

#Create BAM files that contains only dual-mapped reads
java -jar $EBROOTPICARD/picard.jar FilterSamReads \
        I=bamsPasteuria/"$SAMP"_pasteuria_Rm_rdup_indelrealigner.bam \
        O=dualMapped/"$SAMP"_dualMapped_past.bam \
        READ_LIST_FILE=dualMapped/"$SAMP"_dual_mapped.txt \
        FILTER=includeReadList \
        WRITE_READS_FILES=false

echo "finished finding dual-mapped reads for "$SAMP""
dual_mapped_exclude.sh (END)
````



The above script creates BAM files that contain only the dual mapped reads. After all of the files are completed from the previous script (takes about 2 hours using ~250 nodes), the following code merges the reduced BAM files and calculate coverage across the whole set of samples.

````bash
#Request an interactive compute node, navigate to working directory, and load required modules
srun --nodes=1 --cpus-per-task=4 --mem=16G --pty bash
cd interlink
module load SAMtools/1.10-foss-2018b
module load picard

#merge pasteuria BAMs and merge daphnia BAMs
samtools merge -@ 4 dualMapped/merged_past.bam dualMapped/*_dualMapped_past.bam

samtools merge -@ 4 dualMapped/merged_daph.bam dualMapped/*_dualMapped_daph.bam

#Calculate coverage across the new BAM files
#The -m 0 flag is needed to override the 8000X coverage cap for the calculations
samtools depth -a -m 0 dualMapped/merged_past.bam \
-o dualMapped/dual_mapped_coverage_past.txt

samtools depth -m 0 dualMapped/merged_daph.bam \
-o dualMapped/dual_mapped_coverage_daph_smaller.txt

#Intermediate BAMs can be deleted once satisfied that this step ran properly
````



The output files can be examined in R and a set of loci to exclude can be generated

````R
#Load the coverage file for pasteuria
coverage <- read.delim("~/interlink/dualMapped/dual_mapped_coverage_past.txt", header=FALSE)

#Custom function to average across non-overlapping windows
n.colmeans = function(df, n = x){
  aggregate(x = df,
            by = list(gl(ceiling(nrow(df)/n), n)[1:nrow(df)]),
            FUN = mean)
}

#Aggregate over 100 bp windows (too many points to plot otherwise)
windows<-n.colmeans(coverage, 100)

#Histogram of coverages (limited X axis)
hist(coverage$V3,xlim=c(0,10^4),breaks=10000)

#Pick a cutoff value
cutoff<-mean(coverage$V3)+2*sd(coverage$V3)

#Plot with cutoff
plot(windows$V3~windows$V2,type="b",xlab="Position",ylab="Coverage",main="Dual-mapped reads on Pasteuria genome")
abline(h=cutoff,col="red")

#Plot with cutoff
plot(windows$V3~windows$V2,type="b",xlab="Position",ylab="Coverage",main="Dual-mapped reads on Pasteuria genome (zoomed)",ylim=c(0,10^5))
abline(h=cutoff,col="red")

#Get excluded loci
excluded<-coverage[which(coverage$V3 >= cutoff),]

#Find the positions with highest coverage
peak<-coverage[which(coverage$V3 >= 10^5),]

plot(peak$V3~peak$V2,type="p",xlab="Position",ylab="Coverage",main="Dual-mapped reads on Pasteuria genome")

#Make the output file for VCF filtering based on exlcuded loci
write.table(excluded[,1:2], "~/interlink/scripts/dual_mapped_exclude_past.txt",
            sep = "\t",row.names = FALSE, col.names = FALSE, quote = FALSE)

#############################################################################
#DAPHNIA
#############################################################################

#Load the linear positions
contigs <- read.csv("~/scripts/daphnia_contig_length.csv", stringsAsFactors=FALSE)

#Load the coverage file for pasteuria
coverage <- read.delim("~/interlink/dualMapped/dual_mapped_coverage_daph.txt", header=FALSE)

#Pick a cutoff value
cutoff<-mean(coverage$V3)+2*sd(coverage$V3)

#Get excluded loci
excluded<-coverage[which(coverage$V3 >= cutoff),]

#Proportion of loci excluded
length(excluded$V3)/length(coverage$V3)

#Make the output file for VCF filtering
write.table(excluded[,1:2], "~/interlink/scripts/dual_mapped_exclude_daph.txt",
            sep = "\t",row.names = FALSE, col.names = FALSE)
````

Finally, positions from the VCF file that are hotspots of dual-mapping are excluded from the VCF file

````bash
#Request interactive compute node
srun --nodes=1 --cpus-per-task=4 --mem=16G --pty bash

module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0

vcftools --vcf vcfsPasteuria/merged_pasteuria_filtered_annotated_maxDP_purged.vcf --out vcfsPasteuria/merged_pasteuria_filtered_annotated_maxDP_purged_dualmap --exclude-positions scripts/dual_mapped_exclude_past.txt --recode

vcftools --vcf vcfsDaphnia/merged_daphnia_filtered_annotated_maxDP_purged.vcf --out vcfsDaphnia/merged_daphnia_filtered_annotated_maxDP_purged_dualmap --exclude-positions scripts/dual_mapped_exclude_daph.txt --recode
````

