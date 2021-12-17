# Pi calculations for P21 and C1

## C1 version

Calculate pi from the VCF file across a sliding window

````bash
#Load required module
Module load VCFtools

#Calculate pi across a sliding window
vcftools --vcf daphnia_plink.vcf --window-pi 100000 --window-pi-step 1000
````

Analyze and plot the data in R

````R
################################################################################
#Load required packages and data
################################################################################

#Load packages
library(shiny)
library(ggplot2)

#Load pi data
df <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink2/results/pi/daph100kb.windowed.pi")

#Load daphnia contig data
contigs <- read.csv("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink2/data/daphnia_genome/daphnia_contig_length.csv", stringsAsFactors=FALSE)

################################################################################
#Format data
################################################################################

#Merge contig file with results
df<-merge(contigs,df,by.x = "contig",by.y="CHROM")
df$daphBase<-as.numeric(df$BIN_START)+df$startPos

################################################################################
#Plots
################################################################################

#Plot all pi values across genome
plot(df$PI~df$daphBase)

#Plot pi values across contig 11
p<-ggplot(df[which(df$contig =="000011F|quiver"),],aes(x=BIN_START+50000,y=PI))+geom_point()+
  scale_color_manual(values=c("gray", "black"))+
  theme_classic()+
  xlab("Position in D. magna genome") +ylab("Nucleotide diversity (pi)")+
  geom_vline(xintercept=c(2141056,2343857))+ #Lucas TSP
  annotate("rect", xmin=2329948, xmax=2361178, ymin=0, ymax=Inf, alpha=0.2, fill="red")+ #F locus
  annotate("rect", xmin=2189805, xmax=2348995, ymin=0, ymax=Inf, alpha=0.2, fill="blue")#ABC locus
p
ggsave("C:/Users/ericd/Downloads/daphPi100kbABClocus.png",plot = p, width = 9, height = 3)

################################################################################
#Examine plot interactively
################################################################################
#Thin the data to 90th percentile to make plotting easier
cutoff<-quantile(df$PI, 0.90)
df2<-df[which(df$PI>=cutoff),]

#Examine points interactively
ui <- basicPage(
  plotOutput("plot1", brush = "plot_brush"),
  verbatimTextOutput("info")
)

server <- function(input, output) {
  output$plot1 <- renderPlot({
    plot(df2$PI~df2$daphBase)
  })
  
  output$info <- renderPrint({
    # With base graphics, need to tell it what the x and y variables are.
    brushedPoints(df2, input$plot_brush, xvar = "daphBase", yvar = "PI")
  })
}

shinyApp(ui, server)
################################################################################
#END
################################################################################
````

