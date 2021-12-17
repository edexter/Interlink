# Pi calculations for pasteuria C1 genome

Calculate pi from the VCF file across a sliding window

````bash
vcftools --vcf pasteuria_plink.vcf --window-pi 500 --window-pi-step 250
````

Analyze and plot the data in R

````R
################################################################################
#Load packages and data
################################################################################

#Load packages
library(shiny)
library(ggplot2)

#Load pi data
df <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink2/results/pi/C1.windowed.pi")

#Load PCL location file
pcl <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/Interlink2/data/pcl/PCL_locations_forEric_25012021.bed")

################################################################################
#Format data
################################################################################
#Attach PCL names
PCLname<-df$BIN_START
for (i in 1:length(PCLname)) {
  for (j in 1:nrow(pcl)) {
    PCLname[i]<-ifelse(df$BIN_START[i] >= pcl$start[j] & df$BIN_START[i] <= pcl$end[j],
                       PCLname[i]<-as.character(pcl$name[j]),PCLname[i])
  }
}

for (i in 1:length(PCLname)) {
  for (j in 1:nrow(pcl)) {
    PCLname[i]<-ifelse(df$BIN_END[i] >= pcl$start[j] & df$BIN_END[i] <= pcl$end[j],
                       PCLname[i]<-as.character(pcl$name[j]),PCLname[i])
  }
}

#Determine in a variant is located within a PCL
PCL1<-df$BIN_START
for (i in 1:length(PCL1)) {
  for (j in 1:nrow(pcl)) {
    PCL1[i]<-ifelse(df$BIN_START[i] >= pcl$start[j] & df$BIN_START[i] <= pcl$end[j],
                    PCL1[i]<-"PCL",PCL1[i])
  }
}

PCL2<-df$BIN_START
for (i in 1:length(PCL2)) {
  for (j in 1:nrow(pcl)) {
    PCL2[i]<-ifelse(df$BIN_END[i] >= pcl$start[j] & df$BIN_END[i] <= pcl$end[j],
                    PCL2[i]<-"PCL",PCL2[i])
  }
}
df2<-cbind(df,PCL1,PCL2,PCLname)
df2$PCL<-ifelse(df2$PCL1 == "PCL" | df2$PCL2 == "PCL",1,0)

################################################################################
#Plots
################################################################################

#Plot the data with the position of PCLs and PCL triplets shown
p<-ggplot(df2,aes(x=BIN_START,y=PI))+geom_point(aes(color=as.factor(PCL)))+
  scale_color_manual(values=c("gray", "black"))+
  theme_classic()+
  xlab("Position in P. ramosa genome") +ylab("Nucleotide diversity (pi)")+
  theme(legend.position = "none")+
  annotate("rect", xmin=7278, xmax=10679, ymin=0, ymax=Inf, alpha=0.5, fill="blue")+ #PCL [6-7-8]
  annotate("rect", xmin=86443, xmax=90028, ymin=0, ymax=Inf, alpha=0.5, fill="blue")+ #PCL [21-22-23]
  annotate("rect", xmin=317101, xmax=319830, ymin=0, ymax=Inf, alpha=0.5, fill="blue")+ #PCL [21B-27-28]
  annotate("rect", xmin=1199952, xmax=1203511, ymin=0, ymax=Inf, alpha=0.5, fill="blue")+ #PCL [15-16-17]
  annotate("rect", xmin=1517919, xmax=1521372, ymin=0, ymax=Inf, alpha=0.5, fill="blue")+ #PCL [24-25-26]
  annotate("rect", xmin=1571213, xmax=1574559, ymin=0, ymax=Inf, alpha=0.5, fill="blue")+ #PCL [9-10-11]
  annotate("rect", xmin=1577719, xmax=1581168, ymin=0, ymax=Inf, alpha=0.5, fill="blue")+ #PCL [12-13-14]
  annotate("rect", xmin=1586627, xmax=1589937, ymin=0, ymax=Inf, alpha=0.5, fill="blue") #PCL [36-37-38]

p  

ggsave("C:/Users/ericd/Downloads/C1pi.png",plot = p, width = 9, height = 3)

################################################################################
#Examine points interactively
################################################################################

ui <- basicPage(
  plotOutput("plot1", brush = "plot_brush"),
  verbatimTextOutput("info")
)

server <- function(input, output) {
  output$plot1 <- renderPlot({
    plot(df2$PI~df$BIN_START)
  })
  
  output$info <- renderPrint({
    # With base graphics, need to tell it what the x and y variables are.
    brushedPoints(df2, input$plot_brush, xvar = "BIN_START", yvar = "PI")
  })
}

shinyApp(ui, server)

################################################################################
#END
################################################################################
````
