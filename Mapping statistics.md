# Mapping statistics

This script generates more extensive mapping statistics than the main bioinformatics pipeline produces, and organizes the statistics produces by the pipeline into an easier to use format. This code executes quickly and can be run interactively through the terminal.

````bash
module load SAMtools/1.10-foss-2018b

#########################################################################################
Pasteuria stats
#########################################################################################
#Extract the sample name from the filenames
ls *_pasteuria_flagstat.txt | awk -F'[_.]' '{print $1}'

#Iterate an awk function over a list of files to get # and % reads mapped
for i in *_pasteuria_flagstat.txt; do
awk -F'[(% ]' 'NR==5{print $1,$6}' $i
done

#Calculate coverage across BAMS
cat /scicore/home/ebertd/dexter0000/interlink/scripts/index | while read line 
do
   samtools coverage "$line"_pasteuria_Rm_rdup_indelrealigner.bam -o "$line"_pasteuria_coverage.txt
done

#Extract the sample name from the filenames
ls *_pasteuria_coverage.txt | awk -F'[_.]' '{print $1}'

#Iterate an awk function over a list of files to get coverage stats
for i in *_pasteuria_coverage.txt; do
awk 'NR==2{print $4,$5,$6,$7,$8,$9}' $i
done

#########################################################################################
Daphnia stats
#########################################################################################

#Extract the sample name from the filenames

ls *_daphnia_flagstat.txt | awk -F'[_.]' '{print $1}'

#Iterate an awk function over a list of files to get # and % reads mapped
for i in *_daphnia_flagstat.txt; do
awk -F'[(% ]' 'NR==5{print $1,$6}' $i
done

#Calculate coverage across BAMS
cat /scicore/home/ebertd/dexter0000/interlink/scripts/index | while read line 
do
   samtools coverage "$line"_daphnia_Rm_rdup_indelrealigner.bam -o "$line"_daphnia_coverage.txt
done

# Extract the sample name from the filenames
ls *_daphnia_coverage.txt | awk -F'[_.]' '{print $1}'

#Iterate an awk function over a list of files to get coverage stats
for i in *_daphnia_coverage.txt; do
awk 'NR==2{print $4,$6,$7}' $i
done

#Calculate coverage across all contigs (pasteuria has only a single contig)
for i in *_daphnia_coverage.txt; do
awk '{totalBases+=$3; coveredBases+=$5;} END {print coveredBases/totalBases;}' $i 
done

````

