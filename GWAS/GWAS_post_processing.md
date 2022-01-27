# CO-GWAS post-processing and plotting

This script accepts CO-GWAS results from PLINK and performs some post-processing before producing Manhattan plots and circular interlinkage plots.

````R
################################################################################
#This R script performs post-processing of co-GWAS results from PLINK output
#and produces several types of plots
#Written by Eric Dexter
#Updated 27.01.2022
################################################################################

#Load required packages
library(readr)
library(stringr)
library(dplyr)
library(ggplot2)
library(shiny)
library(circlize)

################################################################################

#Load PLINK association test results
df <- read.csv("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink2/results/gwas_primary/plink_results.txt", sep="", stringsAsFactors=FALSE) 

#Load daphnia contig file
contigs <- read.csv("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink2/data/daphnia_genome/contig2chrom.csv", stringsAsFactors=FALSE)

#Load PCL files
pcl <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink2/data/pcl/PCL_locations_forEric_25012021.bed")

################################################################################
#Format data

#Merge contig file with results
df<-merge(contigs,df,by.x = "contig_quiver",by.y="X.CHROM")

#flip orientation if needed
df$POSflip<-ifelse(df$orientation=="flip",df$contig_length-df$POS,df$POS)

#Recalculate position in genome based on contig position and orientation within
#chromosomes
df$daphBase<-as.numeric(df$POSflip)+df$full_position

#Create numeric vector of pasteuria bases
df$pastBase<-gsub("output.chr1[_]([^.]+)[.].*", "\\1", df[,18])
df$pastBase<-gsub('.{2}$','',df$pastBase)
df$pastBase<-as.numeric(df$pastBase)

################################################################################
#Filter data

#Set a p-value filter to reduce file size
df<-df[df$LOG10_P>=7,]

#Thin results based on linked variation in host. Current window size is 10 kb.

positions<-sort(unique(df$daphBase))
singles<-NULL
for (i in 1:length(positions)){
  singles[i]<-ifelse(positions[i+1]<= (positions[i]+5000 ),1,0)
}

singles2<-NULL
for (i in 1:length(positions)){
  singles2[i]<-ifelse(positions[i-1]>= (positions[i]-5000 ),1,0)
}

combined<-singles+singles2

positions2<-positions[combined ==2]
df<-df[df$daphBase %in% positions2,]

#Thin results based on linked variation in parasites. Current window size is 10 kb.

positions<-sort(unique(df$pastBase))
singles<-NULL
for (i in 1:length(positions)){
  singles[i]<-ifelse(positions[i+1]<= (positions[i]+5000 ),1,0)
}

singles2<-NULL
for (i in 1:length(positions)){
  singles2[i]<-ifelse(positions[i-1]>= (positions[i]-5000 ),1,0)
}

combined2<-singles+singles2

positions2<-positions[combined2 ==2 ]
df<-df[df$pastBase %in% positions2,]

################################################################################
#Annotate data

#Attach PCL names
PCLname<-df$pastBase
for (i in 1:length(PCLname)) {
  for (j in 1:nrow(pcl)) {
    PCLname[i]<-ifelse(df$pastBase[i] >= pcl$start[j] & df$pastBase[i] <= pcl$end[j],
                       PCLname[i]<-as.character(pcl$name[j]),PCLname[i])
  }
}

#Determine in a variant is located (+/- 1 kb) within a PCL
PCL<-df$pastBase
for (i in 1:length(PCL)) {
  for (j in 1:nrow(pcl)) {
    PCL[i]<-ifelse(df$pastBase[i] >= pcl$start[j]-1000 & df$pastBase[i] <= pcl$end[j]+1000,
                   PCL[i]<-"PCL",PCL[i])
  }
}

df<-cbind(df,PCL,PCLname)
df$PCL<-ifelse(df$PCL == "PCL",1,0)
df$PCL<-as.factor(df$PCL)

################################################################################
#Manhattan plots

#Plot against daphnia genome
p<-ggplot(df, aes(x=daphBase, y=LOG10_P))+
  annotate("rect", xmin=13817783, xmax=31418508, ymin=7, ymax=Inf, alpha=0.25, fill="gray")+ #Chr2
  annotate("rect", xmin=46370085, xmax=57921315, ymin=7, ymax=Inf, alpha=0.25, fill="gray")+ #Chr4
  annotate("rect", xmin=72372247, xmax=82443146, ymin=7, ymax=Inf, alpha=0.25, fill="gray")+ #Chr6
  annotate("rect", xmin=95294318, xmax=106203609, ymin=7, ymax=Inf, alpha=0.25, fill="gray")+ #Chr8
  annotate("rect", xmin=115152047, xmax=125578570, ymin=7, ymax=Inf, alpha=0.25, fill="gray")+ #Chr10
  theme_classic()+ 
  scale_y_continuous(limits = c(7, 10), breaks = c(7,8,9,10))+
  theme(legend.position = "none")+
  xlab("Position in D. magna genome (MB)") +ylab("P-value (-log10)")+
  annotate("rect", xmin=57119482, xmax=57278672, ymin=7, ymax=Inf, alpha=0.5, fill="blue")+ #ABC locus
  annotate("rect", xmin=57259625, xmax=57290855, ymin=7, ymax=Inf, alpha=0.5, fill="blue")+ #F locus
  annotate("rect", xmin=64965787, xmax=65120674, ymin=7, ymax=Inf, alpha=0.5, fill="green")+ #E locus
  annotate("rect", xmin=90872554, xmax=91095554, ymin=7, ymax=Inf, alpha=0.5, fill="red")+ #D locus (enlarged to show on plot)
  geom_point()+
  theme(text = element_text(size = 16)) +
  scale_x_continuous(breaks=c(0,25000000,50000000,75000000,100000000,125000000), labels= c("0", "25", "50","75","100","125"))

ggsave("C:/Users/ericd/Downloads/manhattan_primary_model_daphnia.png",plot = p, width = 9, height = 3)

#Plot against Pasteuria genome
p<-ggplot(df, aes(x=pastBase, y=LOG10_P))+ 
  geom_point(aes(color=PCL))+
  scale_y_continuous(limits = c(7, 10), breaks = c(7,8,9,10))+
  theme_classic() +
  theme(legend.position = "none")+
  annotate("rect", xmin=7278, xmax=10679, ymin=7, ymax=Inf, alpha=0.25, fill="blue")+ #PCL 6-7-8
  annotate("rect", xmin=86443, xmax=90028, ymin=7, ymax=Inf, alpha=0.25, fill="blue")+ #PCL 23-22-21A
  annotate("rect", xmin=317101, xmax=319830, ymin=7, ymax=Inf, alpha=0.25, fill="blue")+ #PCL 21B-27-28
  annotate("rect", xmin=1199952, xmax=1203511, ymin=7, ymax=Inf, alpha=0.25, fill="blue")+ #PCL 15-16-17
  annotate("rect", xmin=1517919, xmax=1521372, ymin=7, ymax=Inf, alpha=0.25, fill="blue")+ #PCL 24-25-26
  annotate("rect", xmin=1571213, xmax=1574559, ymin=7, ymax=Inf, alpha=0.25, fill="blue")+ #PCL 9-10-11
  annotate("rect", xmin=1577719, xmax=1581168, ymin=7, ymax=Inf, alpha=0.25, fill="blue")+ #PCL 12-13-14
  annotate("rect", xmin=1586627, xmax=1589937, ymin=7, ymax=Inf, alpha=0.25, fill="blue")+ #PCL 36-37-38
  scale_color_manual(values=c("gray", "black"))+
  xlab("Position in P. ramosa genome (MB)") +ylab("P-value (-log10)")+
  theme(text = element_text(size = 16)) +
  scale_x_continuous(breaks=c(0,250000,500000,750000,1000000,1250000,1500000,1750000), labels= c("0","0.25","0.50","0.75","1.0","1.25","1.50","1.75"))

ggsave("C:/Users/ericd/Downloads/manhattan_primary_model_pasteuria.png",plot = p, width = 9, height = 3)

#Plot against single contigs
################################################################
#ABC and F locus
p<-ggplot(df[which(df$contig==11),], aes(x=POS, y=LOG10_P))+ theme_classic()+
  xlab("Position in D. magna genome contig 11 (MB)") +ylab("P-value (-log10)")+
  #geom_vline(xintercept=c(2141056,2343857))+ #Lucas TSP
  #annotate("rect", xmin=2329948, xmax=2361178, ymin=7, ymax=Inf, alpha=0.7, fill="yellow")+ #F locus
  annotate("rect", xmin=2189805, xmax=2348995, ymin=7, ymax=Inf, alpha=0.2, fill="blue")+ #ABC locus
  geom_point()+
  theme(text = element_text(size = 22))+
  scale_y_continuous(limits = c(7, 10), breaks = c(7,8,9,10))+
  scale_x_continuous(limits = c(2000000,3000000),
                     breaks = c(2000000,2200000,2400000,2600000,2800000,3000000),
                     labels = c(2.00,2.2,2.4,2.6,2.8,3.00 ))
ggsave("C:/Users/ericd/Downloads/manhattan_primary_contig11_thin_1MB.png",plot = p, width = 9, height = 3)

#D locus
p<-ggplot(df[which(df$contig==18),], aes(x=POS, y=LOG10_P)) + theme_classic()+
  xlab("Position in D. magna genome contig 18 (MB)") +ylab("P-value (-log10)")+
  annotate("rect", xmin=1801000, xmax=1804000, ymin=7, ymax=Inf, alpha=0.2, fill="red")+ #Dlocus
  geom_point()+
  theme(text = element_text(size = 22))+
  scale_y_continuous(limits = c(7,9), breaks = c(7,7.5,8,8.5,9))+
  scale_x_continuous(limits = c(1700000,1900000),
                     breaks = c(1700000,1750000,1800000,1850000,1900000),
                     labels = c(1.7,1.75,1.8,1.85,1.9 ))
ggsave("C:/Users/ericd/Downloads/manhattan_primary_contig18_thin_1MB.png",plot = p, width = 9, height = 3)

#E locus
p<-ggplot(df[which(df$chr == 5),], aes(x=daphBase, y=LOG10_P))+ geom_point() + theme_classic()+
  xlim(57921315,75984711)+ #entire left arm of chr 5
  xlab("Position on D. magna chromosome 5 left arm") +ylab("P-value (-log10)")+
  annotate("rect", xmin=64731952, xmax=65154682, ymin=7, ymax=Inf, alpha=0.25, fill="green") #Elocus (just contig 67)+
ggsave("C:/Users/ericd/Downloads/manhattan_primary_chrome5L_thin.png",plot = p, width = 9, height = 3)

################################################################################
#Interactive plots to query individual points
df2<-df[df$LOG10_P>=7.0,]
df3<-df2[df2$contig==11,]


ui <- basicPage(
  plotOutput("plot1", brush = "plot_brush"),
  verbatimTextOutput("info")
)

server <- function(input, output) {
  output$plot1 <- renderPlot({
    plot(df3$LOG10_P~df3$daphBase)
  })
  
  output$info <- renderPrint({
    # With base graphics, need to tell it what the x and y variables are.
    brushedPoints(df3, input$plot_brush, xvar = "daphBase", yvar = "LOG10_P")
  })
}

shinyApp(ui, server)

################################################################################
#Circular plots
################################################################################
#The results have to be thinned or plotting is impossible
results<-df[df$LOG10_P>=7.5,]

#Load PCL triplet locations
pcl2 <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink2/data/pcl/PCL_C1_triplets.txt")

#Order descending by p-value
results<-results[order(-results$LOG10_P),]

#Standardize genomic position to 0-1 scale
results$index<-1:length(results$pastBase)
results$daphBase<-results$daphBase/125578570
results$pastBase<-results$pastBase/1723598

#Recalculate PCL positions
pcl2$plotLines<-pcl2$start/1723598

#Create index
index<-nrow(results)
genome<-as.factor(c(rep("Host",index),rep("Parasite",index)))
x<-c(results$daphBase,results$pastBase)
y<-c(results$LOG10_P,results$LOG10_P)
test<-data.frame(genome,x,y)

#Initialize PDF (more efficient)
pdf(file="C:/Users/ericd/Downloads/circle_secondary_thin8.pdf",
    width=5, height=5)
# Step1: Initialise the chart giving factor and x-axis.
circos.clear()
circos.par(gap.degree=15,"clock.wise"=TRUE,start.degree=-8)
circos.initialize( factors=test$genome, x=test$x, xlim = c(0, 1) )
circos.trackPlotRegion(factors = test$genome, y = test$y, panel.fun = function(x, y) {
  circos.axis(major.at = c(seq(0,1,0.1)))
})
#Add lines for PCLs
circos.lines(pcl2$plotLines, rep(9.9,length(pcl2$plotLines)), type = "h", baseline = 7.9,col="darkgoldenrod1")

#Add points
circos.trackPoints(test$genome, test$x, test$y, col = "black", pch = 16, cex = 0.5) 

#Add lines
blues<-results[which(results$LOG10_P>=7.0 & results$contig_quiver== "000011F|quiver"),]
index<-nrow(blues)

host<-rep("Host",index)
parasite<-rep("Parasite",index)
hostPos<-blues$daphBase
parPos<-blues$pastBase

for (i in 1:index){
  circos.link(host[i], hostPos[i], parasite[i], parPos[i], h = 1,col=rgb(red = 0, green = 0, blue = 1.0, alpha = 0.1),lwd=1)
}

#Add lines
reds<-results[which(results$LOG10_P>=7.0 & results$contig_quiver== "000018F|quiver"),]
index<-nrow(reds)

host<-rep("Host",index)
parasite<-rep("Parasite",index)
hostPos<-reds$daphBase
parPos<-reds$pastBase

for (i in 1:index){
  circos.link(host[i], hostPos[i], parasite[i], parPos[i], h = 1,col=rgb(red = 1.0, green = 0, blue = 0.0, alpha = 0.1),lwd=1)
}

#Add lines
greens<-results[which(results$LOG10_P>=7.0 & results$chr == 5),]
index<-nrow(greens)

host<-rep("Host",index)
parasite<-rep("Parasite",index)
hostPos<-greens$daphBase
parPos<-greens$pastBase

for (i in 1:index){
  circos.link(host[i], hostPos[i], parasite[i], parPos[i], h = 1,col=rgb(red = 0.0, green = 1.0, blue = 0.0, alpha = 0.1),lwd=1)
}

grays<-results[which(results$LOG10_P>=7.0 & results$chr != 5 & results$contig_quiver!= "000011F|quiver" & results$contig_quiver!= "000018F|quiver" ),]
index<-nrow(grays)

host<-rep("Host",index)
parasite<-rep("Parasite",index)
hostPos<-grays$daphBase
parPos<-grays$pastBase

for (i in 1:index){
  circos.link(host[i], hostPos[i], parasite[i], parPos[i], h = 1,col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.1),lwd=1)
}

dev.off()
````
