# LD dissolution with distance code

This script calculates LD dissolution with distance using the pasteuria C1 reference genome

````bash
#Request interactive node
srun --nodes=1 --cpus-per-task=8 --mem=32G --pty bash

#Load VCFtools for file conversion
module load VCFtools

#Load older version of plink for code compatibility
module load PLINK/1.90b_170113

#Create plink input files
vcftools --vcf vcfsPasteuria/merged_pasteuria_filtered_annotated.vcf --plink --out LD/C1

#Convert plink files to binary
plink --file LD/C1 --recode --out LD/C1 --make-bed

#Calculate LD
plink --bfile LD/C1 --r2 --ld-window-r2 0 --ld-window 999999 --ld-window-kb 8000 --out LD/C1 --maf 0.05

#Format results for plotting
cat LD/C1.ld | sed 1,1d | awk -F " " 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($5-$2),$7}' | sort -k1,1n> LD/C1.ld.summary
````

This code block repeats the  same calculations, but using the pasteuria P21 reference genome

````bash
#Request interactive node
srun --nodes=1 --cpus-per-task=8 --mem=32G --pty bash

#Load VCFtools for file conversion
module load VCFtools

#Load older version of plink for code compatibility
module load PLINK/1.90b_170113

#Create plink input files
vcftools --vcf vcfsPasteuria_P21/merged_pasteuria_filtered_annotated.vcf --plink --out LD/P21

#Convert plink fils to binary
plink --file LD/P21 --recode --out LD/P21 --make-bed

#Calculate LD
plink --bfile LD/P21 --r2 --ld-window-r2 0 --ld-window 999999 --ld-window-kb 8000 --out LD/P21  --maf 0.05

#Format results for plotting
cat LD/P21.ld | sed 1,1d | awk -F " " 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($5-$2),$7}' | sort -k1,1n > LD/P21.ld.summary
````

This code block plots LD dissolution for the C1 and P21 pasteuria reference genomes in a single panel.

````R
#Load required packages
library(dplyr)
library(stringr)
library(ggplot2)

#Load data
dfr <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink3/results/LD_dissolution/C1.ld.summary",sep="",header=F,check.names=F,stringsAsFactors=F)
dfs <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink3/results/LD_dissolution/P21.ld.summary",sep="",header=F,check.names=F,stringsAsFactors=F)

#Format data
colnames(dfr) <- c("dist","rsq")
colnames(dfs) <- c("dist","rsq")

#Calculate stats
dfr$distc <- cut(dfr$dist,breaks=seq(from=min(dfr$dist)-1,to=max(dfr$dist)+1,by=100))
dfr1 <- dfr %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
dfr1 <- dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))

#Calculate stats
dfs$distc <- cut(dfs$dist,breaks=seq(from=min(dfs$dist)-1,to=max(dfs$dist)+1,by=100))
dfs1 <- dfs %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
dfs1 <- dfs1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))

df<-rbind(dfr1,dfs1)
df$genome<-rep(c("C1","P21"),c(nrow(dfr1),nrow(dfs1)))

#Check the range
range(df$distc)

#Plot data
ggplot()+
  geom_point(data=df,aes(x=start,y=mean,color=genome),size=0.5)+
  geom_line(data=df,aes(x=start,y=mean,color=as.factor(genome)))+
  labs(x="Distance (bp)",y=expression(LD~(r^{2})))+
  scale_x_continuous(breaks=c(2000,4000,6000,8000,10000),labels=c("2","4","6","8","10"))+
  theme_classic()+xlim(0,10000)+ylim(0.25,0.75)+
  labs(color = "Reference\n  genome")+
  theme(legend.position = c(0.9, 0.8))

ggsave("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink3/figures/LD_dissolution.png",width = 5, height = 3,)
````
