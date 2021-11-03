###################
## LOAD PACKAGES ##
###################
setwd("D:/codes/GalaxyFrogs_metabarcoding")
## Load custom graphical functions stored in custom R-scripts folder in your working directory
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.14")
library("BiocManager")
BiocManager::install("phyloseq")

library(phyloseq)
## Load additional packages (mostly for nice graphs)
library(ggplot2)
library(plyr)
library(reshape2)
library(grid)
library(gridExtra)
library(scales)
library(ape)
library(vegan)
##############################
## IMPORT 16S AND GYRB DATA ##
##############################
## Import data from tabular files - These files can be created from the POIRIER_ET_AL_FOOD_MICROBIOME_DATASET_V2.xlsx
## available at figshare DOI:10.6084/m9.figshare.7083209
otuTable <- as.matrix(read.csv("otu_table_redlosses.tsv", sep = "\t", row.names = 1))
taxTable <- as.matrix(read.csv("tax_table_redlosses.tsv", sep = "\t", row.names = 1))
sampleData <- read.csv("sample_data_redlosses.tsv", sep = "\t", row.names = 1)
rdldata <- phyloseq(otu_table(otuTable, taxa_are_rows = TRUE), 
                    tax_table(taxTable), 
                    sample_data(sampleData))
#Subsetting the dataset between markers and type of samples
rdl16S<-subset_samples(rdldata, MarkerType %in% "16S")
rdlgyrB<-subset_samples(rdldata, MarkerType %in% "GYRB")
rdl_mockdata<-subset_samples(rdldata, EnvType %in% "6_Mock_data")
rdl_fooddata<-subset_samples(rdldata, EnvType %in% c("1_CodFillet","2_SalmonFillet", "3_BeefBurger","4_PoultrySausage","5_PorkSausage"))
## Subsetting the data at the scale of phyla
firmicutes.rdldata <- subset_taxa(rdldata, Phylum == "Firmicutes")
proteobacteria.rdldata <- subset_taxa(rdldata, Phylum == "Proteobacteria")
###############
##  FIGURE 1 ##
###############
## Rarefaction curves (no downsampling)
##ggrare function is called from the redlosses_richness.R.
pRare <- ggrare(rdldata, step = 100, color = "MarkerType", plot = FALSE)
p<- pRare + facet_wrap(~EnvType, ncol = 6) + theme_bw() + ylab("OTU richness") + xlab("Number of sequences")
#p <- p + annotate("text", label = c("","","","","","MC5"), size = 4, x = 50000, y = 200, col="#1dc6c9")
#p <- p + annotate("text", label = c("","","","","","MC5"), size = 4, x = 50000, y = 125, col="#f47167")
plot(p)
#figure 1 was saved as pdf using 11.98 x 4.53 inches portrait.
###############
##  FIGURE 2 ##
###############
#loading the gyrB-parE ratio data. Ratio.csv can be obtained from the POIRIER_ET_AL_FOODMICROBIOME_DATASET.xlsx dataset
## available at figshare DOI:10.6084/m9.figshare.7083209
data = read.csv("Ratio.tsv", header = TRUE, sep = "\t", row.names = 8,stringsAsFactors=FALSE)
data_barplot = data[,c(8:10)]
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE),  widths=c(4,1), heights=c(5,1))
par(fg = "NA", col.lab = "black", col.axis = "NA",xaxs="i", yaxs="i",mar=c(0,1,1,1),font=1)
plot(0, type = "n")
brackets(1.25, 0.27, 1.25, 1, lwd=0.5, type=4, h=0.025,col="black")
mtext("Firmicutes",side=1,line=-33.75,at=1.185,col="black",font=3)
brackets(1.25, -0.69, 1.25, 0.23, lwd=0.5, type=4, h=0.025,col="black")
mtext("Proteobacteria",side=1,line=-16.2,at=1.17,col="black",font=3)
brackets(1.25, -0.757, 1.25, -0.727, lwd=0.5, type=4, h=0.025,col="black")
mtext("Fusobacteria",side=1,line=-5.8,at=1.177,col="black",font=3)
brackets(1.25, -0.90, 1.25, -0.79, lwd=0.5, type=4, h=0.025,col="black")
mtext("Bacteroidetes",side=1,line=-3.7,at=1.175,col="black",font=3)
brackets(1.25, -0.999, 1.25, -0.935, lwd=0.5, type=4, h=0.025,col="black")
mtext("Actinobacteria",side=1,line=-1.25,at=1.173,col="black",font=3)
par(fg = "black", col.lab = "black", col.axis = "black",
    xaxs="i",yaxs="i",mar=c(0,3,1,1),mgp=c(2.5,0.1,0),tcl=-0.1,lwd=0.2)
espace=c(0.25,0.25,1.25,0.25,0.25,1.25,1.25,0.25,0.25,0.25,
         0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,
         0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,
         0.25,1.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,
         0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,
         0.25)
color<-c("#00d820","#ff9811","#ff000a")
barplot(height=as.matrix(t(data_barplot)*100),
        xaxs="i", yaxs="i",
        legend.text=NULL,beside=FALSE,horiz=T,col=color,lwd=0.5,
        cex.axis=0.85,cex.names=0.85,las=1,space=espace,font.axis=3,xaxt="n")
barplot(height=as.matrix(t(data_barplot)*100),
        xaxs="i", yaxs="i",
        legend.text=NULL,beside=FALSE,horiz=T,lwd=0.5,
        density = c(0,70,0),angle=c(0,50,0),col="black",
        cex.axis=0.85,cex.names=0.85,las=1,space=espace,xaxt="n",yaxt="n",add=T)
axis(1,cex.axis=0.8,lwd=0.5)
mtext(text = "Relative abundance (%)", side=1, line=1.5,cex=0.8,col="black")
par(fg = NA, col.lab = NA, col.axis = NA,xaxs="i", yaxs="i",mar=c(0,1,1,1),font=1)
plot(0, type = "n")
par(fg = NA, col.lab = NA, col.axis = NA,xaxs="i", yaxs="i",mar=c(0,1,1,1),font=1)
plot(0, type = "n")
legende = c(expression(paste(italic(gyrB)," amplicon")),
            "Uncertainty interval",
            expression(paste(italic(parE)," amplicon")))
legend(x=0.6,y=0.5,
       legend = legende,
       horiz=F,fill=color,density=c(NA, NA, NA),
       text.col ="black",bty = "n",adj=0,cex=1)
legend(x=0.6,y=0.5,
       legend = legende,
       horiz=F,fill=c("white","black",NA),density=c(0, 70,0),
       text.col ="black",bty = "n",adj=0,cex=1)
dev.off()
###############
##  FIGURE 4 ##
###############
################################################
## ORDINATION BRAY-CURTIS GENUS LEVEL FIGURE 4##
################################################
## 1 - Merging the data at the scale of genus ##
mergedTaxa_all <- tax_glom(rdldata, "Genus")
mergedTaxa_all <- prune_taxa(taxa_sums(mergedTaxa_all ) > 0, mergedTaxa_all)
mergedTaxa_all_mockdata<-subset_samples(mergedTaxa_all, EnvType %in% "6_Mock_data")
mergedTaxa_all_mockdata
mergedTaxa_all_food<-subset_samples(mergedTaxa_all, EnvType %in% c("1_CodFillet","2_SalmonFillet", "3_BeefBurger","4_PoultrySausage","5_PorkSausage"))
mergedTaxa_all_food
mockPalette <- c('#FF0000','#FF0000','#00CD00','#00CD00',
                 '#0000FF','#0000FF','#FFA500','#FFA500','#EE1289','#EE1289')
white.mockPalette <- c('white','white','white','white',
                       'white','white','white','white','white','white')
foodPalette <- c('#FF0000','#FF0000','#FF0000','#FF0000','#FF0000','#FF0000',
                 '#00CD00','#00CD00','#00CD00','#00CD00','#00CD00','#00CD00',
                 '#0000FF','#0000FF','#0000FF','#0000FF','#0000FF','#0000FF',
                 '#FFA500','#FFA500','#FFA500','#FFA500','#FFA500','#FFA500',
                 '#EE1289','#EE1289','#EE1289','#EE1289','#EE1289','#EE1289')
vjust.mock <- c(-2,2,-2,2,-2,2,2,2,-2,-2)
hjust.mock <- c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5)
mock.labels <-c("MC1_16S","MC1_GYRB","MC2_16S","MC2_GYRB","MC3_16S","MC3_GYRB","MC4_16S","MC4_GYRB","MC5_16S","MC5_GYRB")
#2 - plotting the PcoA
p1 <- plot_samples(mergedTaxa_all_mockdata, ordinate(mergedTaxa_all_mockdata, "MDS", "bray"), color = "SampleID")
p1<- p1 + theme_bw() + scale_color_manual(values = white.mockPalette) + theme(legend.position = "none")
p1 <- p1 + geom_point(size = 3.5, fill = mockPalette, shape=21, color = "black")
p1 <- p1 + geom_text(label=mock.labels, vjust=vjust.mock, hjust=hjust.mock, color = mockPalette)
p1<- p1 + theme(text = element_text(size=14)) + theme(axis.text = fontaxis)
p1<- p1 + scale_x_continuous(limits = c(-0.5,0.6))
plot(p1)
white.foodPalette <- c('white','white','white','white','white',
                       'white','white','white','white','white',
                       'white','white','white','white','white',
                       'white','white','white','white','white',
                       'white','white','white','white','white',
                       'white','white','white','white','white')
vjust.food <- c(1,2,0.5,2,0.5,-1.5,0.5,-1.0,2,-1.0,2,-1.0,0.5,-1.5,-1.5,-1.5,0.5,0.5,2,0.5,0.5,0.5,-1.5,2,2,-1.5,1,1.5,2,0.5)
hjust.food <- c(-0.25,0.5,-0.25,0.5,-0.25,0.0,-0.25,1.0,1.0,1.0,1.0,0.5,-0.25,0.5,0.5,0.5,-0.25,-0.25,0.5,1.25,-0.25,-0.25,1.0,1.0,0.0,0.0,1.25,-0.15,0.5,-0.25)
food.labels <-c("CF1_16S","CF1_GYRB","CF2_16S","CF2_GYRB","CF3_16S","CF3_GYRB",
                "SF1_16S","SF1_GYRB","SF2_16S","SF2_GYRB","SF3_16S","SF3_GYRB",
                "GB1_16S","GB1_GYRB","GB2_16S","GB2_GYRB","GB3_16S","GB3_GYRB",
                "CS1_16S","CS1_GYRB","CS2_16S","CS2_GYRB","CS3_16S","CS3_GYRB",
                "PS1_16S","PS1_GYRB","PS2_16S","PS2_GYRB","PS3_16S","PS3_GYRB")
p2 <- plot_samples(mergedTaxa_all_food, ordinate(mergedTaxa_all_food, "MDS", "bray"), color = "SampleID")
p2 <- p2 + theme_bw() + scale_color_manual(values = white.foodPalette) + theme(legend.position = "none")
p2 <- p2 + geom_point(size = 3.5, shape=21, fill=  foodPalette,color = "black")
p2 <- p2 + geom_text(label=food.labels, vjust=vjust.food, hjust=hjust.food, color = foodPalette)
p2<- p2 + theme(text = element_text(size=14)) + theme(axis.text = fontaxis)
p2<- p2 + scale_x_continuous(limits = c(-0.75,0.75)) + scale_y_continuous(limits = c(-0.45,0.50))
plot(p2)
#All plots were saved as 5.27 x 5.27 inches portrait.
###############
##  FIGURE 5 ##
###############
#############################################################
## Estimating Phylum and Genus richness in samples FIGURE 5##
#############################################################
whitepalette<-c("white","white","white","white","white","white","white","white","white","white")
blackpalette<-c("black","black","black","black","black","black","black","black","black","black")
phylum.palette<-c("#6cbe45","#f7981d","#003b7d","#9a37bf","#cb2027","white","white","white","white","white","white")
firmicutes.palette<-c("#003b7d","#3182bd","#6baed6","#bdd7e7","#fafb97","#c2e699","#78c679","#31a354","#006837","black")
proteobacteria.palette<-c("#800026","#cb2027","#f03b20","#fd8d3c","#fecc5c","#ffffb2","#9ecae1","#4291c6","#08519c","black")
fontaxis <- element_text(color = "black", size = 14)
#plotting
p1<- plot_composition(rdldata, "Kingdom", "Bacteria", "Phylum", fill = "Phylum") + facet_grid(~EnvType, scales = "free_x", space = "free_x")
p1<- p1 + theme_bw() + theme(axis.text.x = element_blank()) + theme(axis.text = fontaxis) 
p1 <- p1 + scale_fill_manual(values=phylum.palette) + scale_color_manual(values=blackpalette)
p1 <- p1 + theme(text = element_text(size=14))
plot(p1)
p2<-plot_composition(firmicutes.rdldata, "Phylum", "Firmicutes", "Genus", fill = "Genus") + facet_grid(~EnvType, scales = "free_x", space = "free_x")
p2<- p2 + theme_bw() + theme(axis.text.x = element_blank()) + theme(axis.text = fontaxis) 
p2<- p2 + scale_fill_manual(values=firmicutes.palette) + scale_color_manual(values=blackpalette) 
p2<- p2 + theme(text = element_text(size=14))
plot(p2)
p2<- p2 + theme(legend.margin = unit(1, "cm"))
p3<- plot_composition(proteobacteria.rdldata, "Phylum", "Proteobacteria", "Genus", fill = "Genus")  + facet_grid(~EnvType, scales = "free_x", space = "free_x")
p3<- p3 + theme_bw() + theme(axis.text.x = element_blank()) + theme(axis.text = fontaxis) 
p3<- p3 + scale_fill_manual(values=proteobacteria.palette) + scale_color_manual(values=blackpalette) 
p3<- p3 + theme(text = element_text(size=14))
plot(p3)
#optional : p3<- p3 + theme(legend.margin = unit(1, "cm"))
#All plots were saved as 4.98 x 17.27 inches landscape.
###############
##  FIGURE 6 ##
###############
#data qPCR_vs_16S_V4.csv and qPCR_vs_gyrB_V3.csv can be produced from loading from the POIRIER_ET_AL_FOODMICROBIOME_DATASET_V2.xlsx.
#available at figshare DOI:10.6084/m9.figshare.7083209
data_16S = read.csv(file = "qPCR_vs_16S_V4.csv",sep=";",header = T)
data_gyrb = read.csv(file = "qPCR_vs_gyrB_V3.csv",sep=";",header = T)
data_gyrb = data_gyrb[c(c(1:235),c(238:248)),]
#selection of data for which log reads>2 et log CFU >3
data_16S = subset(data_16S,data_16S[,9]>3)
data_16S = subset(data_16S,data_16S[,11]>2)
data_gyrb = subset(data_gyrb,data_gyrb[,10]>3)
data_gyrb = subset(data_gyrb,data_gyrb[,12]>2)
data_16S = subset(data_16S,data_16S[,6]!="Photobacterium")
data_gyrb = subset(data_gyrb,data_gyrb[,6]!="Photobacterium")
#Tiff
# tiff("C:/Users/sipoirier/Desktop/Simon Poirier/Articles/1-En pr?paration/7-16S gyrB/Figures/Figure6_V3.tiff",
#      res=600,width=4200,height=3200,pointsize=12)
# matrix layout
layout(matrix(c(1,2,3,4,5,5,5,5,6,6,6,6), 3, 4, byrow = T),  widths=c(1,1,1,1), heights=c(1.25,1.25,0.4))
#subseting files according to phyla
data_16S_firmicutes = subset(data_16S, data_16S[,2] == "Firmicutes")
data_16S_proteobacteria = subset(data_16S, data_16S[,2] == "Proteobacteria")
data_gyrb_firmicutes = subset(data_gyrb, data_gyrb[,2] == "Firmicutes")
data_gyrb_proteobacteria = subset(data_gyrb, data_gyrb[,2] == "Proteobacteria")
pink16S<-("#f47167")
bluegyrB<-("#1dc6c9")
# pink50_gyrB <- rgb(253, 227, 225, max = 255, alpha = 125, names = "pink50")
# blue50_16S <- rgb(204, 242, 243, max = 255, alpha = 125, names = "blue50")
########################################################################################
#16S FIRMICUTES
par(fg="black",col.lab="black",col.axis="black")
par(xaxs="i", yaxs="i",mar=c(2,3,1,1),mgp=c(2.5,0.25,0),tcl=-0.2,lwd=0.5)
plot(x=NULL,y=NULL,xlim=c(2,9),ylim=c(2,8.4),type="o",
     ylab="", xlab="",xaxt="n",yaxt="n",
     font.axis=1,cex.axis=1,pch=16,cex=0.8,cex.lab=1)
axis(lwd=0.5,side=1,at=c(seq(from=2,to=9,by=1)),labels=c(seq(from=2,to=9,by=1)))
axis(lwd=0.5,side=2,at=c(seq(from=2,to=8,by=1)),labels=c(seq(from=2,to=8,by=1)))
abline(h = seq(3, 8, by = 1), col="grey40",lty="dotted")
mtext(expression(paste("CFU/g (in log"[10],")")), side=2, line=1.25,cex=0.95)
mtext(expression(paste("16S reads (in log"[10],")")),side=1, line=2,cex=0.95)
mtext("within Firmicutes",side=1, line=3,cex=0.95)
points(data_16S_firmicutes[,9]~data_16S_firmicutes[,15],col= alpha(pink16S, 0.5),pch=16,cex=1.5)
points(data_16S_firmicutes[,9]~data_16S_firmicutes[,15],col= pink16S,pch=21,cex=1.5)
#Model
x=data_16S_firmicutes[,15]
y=data_16S_firmicutes[,9]
lm1 <- lm(y~x)
abline(lm(y~x),lwd=1.5)
summary(lm1)
mtext("r? = 0.64",side=1, line=-9.5,at=7.75,cex=1.1)
# mtext("Firmicutes",side=1, line=-21,at=4,cex=1.1,col="navyblue")

# par(fg = NA, col.lab = NA, bg=NA,col.axis = NA,xaxs="i", yaxs="i",mar=c(0,0,0,0))
# legend(x=1.5,y=8,cex=1.75,
#        legend = "Firmicutes",
#        ncol=1, col = "navyblue", pch=16,bty = "o", pt.cex = 1.5, text.col = "navyblue",adj=0, inset = c(0, -0.5))

########################################################################################
#gyrB FIRMICUTES
par(fg="black",col.lab="black",col.axis="black")
par(xaxs="i", yaxs="i",mar=c(2,3,1,1),mgp=c(2.5,0.25,0),tcl=-0.2,lwd=0.5)
plot(x=NULL,y=NULL,xlim=c(2,9),ylim=c(2,8.4),type="o",
     ylab="", xlab="",xaxt="n",yaxt="n",
     font.axis=1,cex.axis=1,pch=16,cex=0.8,cex.lab=1)
axis(lwd=0.5,side=1,at=c(seq(from=0,to=9,by=1)),labels=c(seq(from=0,to=9,by=1)))
axis(lwd=0.5,side=2,at=c(seq(from=0,to=8,by=1)),labels=c(seq(from=0,to=8,by=1)))
abline(h = seq(3, 8, by = 1), col="grey40",lty="dotted")
mtext(expression(paste("CFU/g (in log"[10],")")), side=2, line=1.25,cex=0.95)
mtext(expression(paste(italic(gyrB)," reads (in log"[10],")")),side=1, line=2,cex=0.95)
mtext("within Firmicutes",side=1, line=3,cex=0.95)
points(data_gyrb_firmicutes[,10]~data_gyrb_firmicutes[,16],col= alpha(bluegyrB, 0.5),pch=16,cex=1.5)
points(data_gyrb_firmicutes[,10]~data_gyrb_firmicutes[,16],col= bluegyrB,pch=21,cex=1.5)
#Model
x=data_gyrb_firmicutes[,16]
y=data_gyrb_firmicutes[,10]
lm1 <- lm(y~x)
abline(lm(y~x),lwd=1.5)
summary(lm1)
mtext("r? = 0.66",side=1, line=-9.5,at=7.75,cex=1.1)
# mtext("Firmicutes",side=1, line=-21,at=4,cex=1.1,col="navyblue")
##############################################################################################
#16S PROTEOBACTERIA
par(fg="black",col.lab="black",col.axis="black")
par(xaxs="i", yaxs="i",mar=c(2,3,1,1),mgp=c(2.5,0.25,0),tcl=-0.2,lwd=0.5)
plot(x=NULL,y=NULL,xlim=c(2,9),ylim=c(2,8.4),type="o",
     ylab="", xlab="",xaxt="n",yaxt="n",
     font.axis=1,cex.axis=1,pch=16,cex=0.8,cex.lab=1)
axis(lwd=0.5,side=1,at=c(seq(from=2,to=9,by=1)),labels=c(seq(from=2,to=9,by=1)))
axis(lwd=0.5,side=2,at=c(seq(from=2,to=8,by=1)),labels=c(seq(from=2,to=8,by=1)))
abline(h = seq(3, 8, by = 1), col="grey40",lty="dotted")
points(data_16S_proteobacteria[,9]~data_16S_proteobacteria[,15],col=alpha(pink16S, 0.5),pch=16,cex=1.5)
points(data_16S_proteobacteria[,9]~data_16S_proteobacteria[,15],col=pink16S,pch=21,cex=1.5)
mtext(expression(paste("CFU/g (in log"[10],")")), side=2, line=1.25,cex=0.95)
mtext(expression(paste("16S reads (in log"[10],")")),side=1, line=2,cex=0.95)
mtext("within Proteobacteria",side=1, line=3,cex=0.95)
#Model
x=data_16S_proteobacteria[,15]
y=data_16S_proteobacteria[,9]
lm1 <- lm(y~x)
abline(lm(y~x),lwd=1.5)
summary(lm1)
mtext("r? = 0.61",side=1, line=-9.5,at=7.75,cex=1.1)
# mtext("Proteobacteria",side=1, line=-21,at=4,cex=1.1,col="red")

# par(fg = NA, col.lab = NA, bg=NA,col.axis = NA,xaxs="i", yaxs="i",mar=c(0,0,0,0))
# legend(x=1.5,y=8,cex=1.75,
#        legend = "Proteobacteria",
#        ncol=1, col = "red", pch=16,bty = "o", pt.cex = 1.5, text.col = "red",adj=0, inset = c(0, -0.5))
##############################################################################################
#gyrB PROTEOBACTERIA
par(fg="black",col.lab="black",col.axis="black")
par(xaxs="i", yaxs="i",mar=c(2,3,1,1),mgp=c(2.5,0.25,0),tcl=-0.2,lwd=0.5)
plot(x=NULL,y=NULL,xlim=c(2,9),ylim=c(2,8.4),type="o",
     ylab="", xlab="",xaxt="n",yaxt="n",
     font.axis=1,cex.axis=1,pch=16,cex=0.8,cex.lab=1)
axis(lwd=0.5,side=1,at=c(seq(from=0,to=9,by=1)),labels=c(seq(from=0,to=9,by=1)))
axis(lwd=0.5,side=2,at=c(seq(from=0,to=8,by=1)),labels=c(seq(from=0,to=8,by=1)))
abline(h = seq(3, 8, by = 1), col="grey40",lty="dotted")
mtext(expression(paste("CFU/g (in log"[10],")")), side=2, line=1.25,cex=0.95)
mtext(expression(paste(italic(gyrB)," reads (in log"[10],")")),side=1, line=2,cex=0.95)
mtext("within Proteobacteria",side=1, line=3,cex=0.95)
points(data_gyrb_proteobacteria[,10]~data_gyrb_proteobacteria[,16],col=alpha(bluegyrB, 0.5),pch=16,cex=1.5)
points(data_gyrb_proteobacteria[,10]~data_gyrb_proteobacteria[,16],col=bluegyrB,pch=21,cex=1.5)
#Model
x=data_gyrb_proteobacteria[,16]
y=data_gyrb_proteobacteria[,10]
lm1 <- lm(y~x)
abline(lm(y~x),lwd=1.5)
summary(lm1)
mtext("r? = 0.78",side=1, line=-9.5,at=7.75,cex=1.1)
# mtext("Proteobacteria",side=1, line=-21,at=4,cex=1.1,col="red")
########################################################################################
#16S GENUS - FIRMICUTES
#subsetting files according to species
data_16S_Lactobacillus_group = subset(data_16S, data_16S[,6] == "Lactobacillus sakei group")
data_16S_Lactobacillus_algidus = subset(data_16S, data_16S[,6] == "Lactobacillus algidus")
data_16S_Lactococcus = subset(data_16S, data_16S[,6] == "Lactococcus")
data_16S_Leuconostoc = subset(data_16S, data_16S[,6] == "Leuconostoc")
data_16S_Carnobacterium = subset(data_16S, data_16S[,6] == "Carnobacterium")
data_16S_Brochothrix = subset(data_16S, data_16S[,6] == "Brochothrix")

data_gyrb_curvatus = subset(data_gyrb, data_gyrb[,7] == "curvatus")
data_gyrb_piscium = subset(data_gyrb, data_gyrb[,7] == "piscium")
data_gyrb_sakei = subset(data_gyrb, data_gyrb[,7] == "sakei")
data_gyrb_carnosum = subset(data_gyrb, data_gyrb[,7] == "carnosum")
data_gyrb_divergens = subset(data_gyrb, data_gyrb[,7] == "divergens")
data_gyrb_gelidum = subset(data_gyrb, data_gyrb[,7] == "gelidum")
data_gyrb_algidus = subset(data_gyrb, data_gyrb[,7] == "algidus")
data_gyrb_thermosphacta = subset(data_gyrb, data_gyrb[,7] == "thermosphacta")

data_16S_Pseudomonas = subset(data_16S, data_16S[,6] == "Pseudomonas")
data_16S_Morganella = subset(data_16S, data_16S[,6] == "Morganella")
data_16S_Serratia = subset(data_16S, data_16S[,6] == "Serratia")
data_16S_Hafnia = subset(data_16S, data_16S[,6] == "Hafnia")

data_gyrb_lundensis = subset(data_gyrb, data_gyrb[,7] == "lundensis")
data_gyrb_fragi = subset(data_gyrb, data_gyrb[,7] == "fragi")
data_gyrb_psychrotolerans = subset(data_gyrb, data_gyrb[,7] == "psychrotolerans")
data_gyrb_proteamaculans = subset(data_gyrb, data_gyrb[,7] == "proteamaculans")
data_gyrb_alvei = subset(data_gyrb, data_gyrb[,7] == "alvei")


par(fg="black",col.lab="black",col.axis="black")
par(xaxs="i", yaxs="i",mar=c(3,3,3,1),mgp=c(2.5,0.25,0),tcl=-0.2,lwd=0.5)
boxplot(x=NULL,y=NULL,xlim=c(0,32),ylim=c(-3,3),type="o",
        ylab="", xlab="",xaxt="n",yaxt="n",
        font.axis=1,cex.axis=1,cex.lab=1)
axis(lwd=0.5,side=2,at=c(seq(from=-4,to=4,by=1)),labels=c(seq(from=-4,to=4,by=1)))
mtext("Deviation from linear model", side=2, line=1.25,cex=0.95)
abline(h = seq(-3, 3, by = 1), col="grey40",lty="dotted")
abline(h=0, col= "red")
par(new=T)
col16S = alpha(pink16S, 0.5)
colgyrB = alpha(bluegyrB, 0.5)
couleurs = c(col16S,colgyrB,
             col16S,colgyrB,
             col16S,colgyrB,
             col16S,colgyrB,colgyrB,
             col16S,colgyrB,
             col16S,colgyrB,colgyrB,
             col16S,colgyrB,colgyrB,
             col16S,colgyrB,
             col16S,colgyrB,
             col16S,colgyrB)
espace=c(1,2,
         4,5,
         7,8,
         10,11,12,
         14,15,
         17,18,19,
         21,22,23,
         25,26,
         28,29,
         31,32)
boxplot(data_16S_Carnobacterium[,18],data_gyrb_divergens[,19],
        data_16S_Brochothrix[,18],data_gyrb_thermosphacta[,19],
        data_16S_Lactobacillus_algidus[,18],data_gyrb_algidus[,19],
        data_16S_Lactobacillus_group[,18],data_gyrb_sakei[,19],data_gyrb_curvatus[,19],
        data_16S_Lactococcus[,18],data_gyrb_piscium[,19],
        data_16S_Leuconostoc[,18],data_gyrb_gelidum[,19],data_gyrb_carnosum[,19],
        data_16S_Pseudomonas[,18],data_gyrb_lundensis[,19],data_gyrb_fragi[,19],
        data_16S_Serratia[,18],data_gyrb_proteamaculans[,19],
        data_16S_Hafnia[,18],data_gyrb_alvei[,19],
        data_16S_Morganella[,18],data_gyrb_psychrotolerans[,19],
        ylab="", xlab="",xaxt="n",yaxt="n",
        col=couleurs,
        at=espace,ylim=c(-3,3))
legende = c(expression(paste(italic("C.")," ",italic("divergens")," ")),
            expression(paste(italic("C.")," ",italic("divergens")," ")),
            expression(paste(italic("B.")," ",italic("thermosphacta")," ")),
            expression(paste(italic("B.")," ",italic("thermosphacta")," ")),          
            expression(paste(italic("L.")," ",italic("algidus")," ")),
            expression(paste(italic("L.")," ",italic("algidus")," ")),
            expression(paste(italic("L.")," ",italic(sakei)," clade ")),
            expression(paste(italic("L.")," ",italic("sakei")," ")),
            expression(paste(italic("L.")," ",italic("curvatus")," ")),
            expression(paste(italic("L.")," ",italic("piscium")," ")),
            expression(paste(italic("L.")," ",italic("piscium")," ")),
            expression(paste(italic("L.")," ",italic(gelidum)," clade ")),
            expression(paste(italic("L.")," ",italic("gelidum")," ")),
            expression(paste(italic("L.")," ",italic("carnosum")," ")),
            expression(paste(italic("P.")," ",italic(fragi)," clade ")),
            expression(paste(italic("P.")," ",italic("lundensis")," ")),
            expression(paste(italic("P.")," ",italic("fragi")," ")),
            expression(paste(italic("S.")," ",italic(proteamaculans)," clade ")),
            expression(paste(italic("S.")," ",italic("proteamaculans")," ")),
            expression(paste(italic("H.")," ",italic(halvei)," clade ")),
            expression(paste(italic("H.")," ",italic("alvei")," ")),
            expression(paste(italic("M.")," ",italic("psychrotolerans")," ")),
            expression(paste(italic("M.")," ",italic("psychrotolerans")," ")))
axis(side=1,las=2,srt=90,at=espace,
     labels=legende,cex=1.75)
###############
##  FIGURE 7 ##
###############
##################################
#ANALYSIS INTRA-SPECIES diversity#
##################################

my_16S_subset <- subset(otu_table(rdl16S), rownames(otu_table(rdl16S)) %in% c('16S_Cluster_1', '16S_Cluster_2','16S_Cluster_3','16S_Cluster_5','16S_Cluster_7','16S_Cluster_8','16S_Cluster_10','16S_Cluster_11','16S_Cluster_13','16S_Cluster_16','16S_Cluster_20','16S_Cluster_30'))
new_rdl16S <- merge_phyloseq(my_16S_subset, tax_table(rdl16S), sample_data(rdl16S))
new_rdl16S

sample.order.16S<-c("CF1_SSU","CF2_SSU","CF3_SSU",
                "SF1_SSU","SF2_SSU","SF3_SSU",
                    "GB1_SSU","GB2_SSU","GB3_SSU",
                    "CS1_SSU","CS2_SSU","CS3_SSU",
                    "PS1_SSU","PS2_SSU","PS3_SSU",
                "MC5_SSU","MC4_SSU","MC3_SSU","MC2_SSU","MC1_SSU")

fontaxis.otus <- element_text(color = "black", size = 10, hjust=1, face="italic")
fontaxis.samples <- element_text(color = "black", size = 10)
p<-plot_heatmap(new_rdl16S, taxa.label= "Species", taxa.order = "Species", sample.order = sample.order.16S)
p <- p + facet_grid(~EnvType, scales = "free", space = "free")
p <- p + scale_fill_gradient2(low = "#1a9850", mid = "#feff9a", high = "#d73027", 
                              na.value = "white", trans = log_trans(10), 
                              midpoint = log(100, base = 10))
p<- p + theme(axis.text.y = fontaxis.otus) + theme(axis.text.x = fontaxis.samples)
plot(p)


my_gyrB_subset <- subset(otu_table(rdlgyrB), rownames(otu_table(rdlgyrB)) %in% c('gyrB_Cluster_7','gyrB_Cluster_108','gyrB_Cluster_23','gyrB_Cluster_141','gyrB_Cluster_2','gyrB_Cluster_3','gyrB_Cluster_16','gyrB_Cluster_4','gyrB_Cluster_13','gyrB_Cluster_6','gyrB_Cluster_11','gyrB_Cluster_21','gyrB_Cluster_30','gyrB_Cluster_39','gyrB_Cluster_8','gyrB_Cluster_32','gyrB_Cluster_26','gyrB_Cluster_37','gyrB_Cluster_12','gyrB_Cluster_25','gyrB_Cluster_60','gyrB_Cluster_42','gyrB_Cluster_43','gyrB_Cluster_14','gyrB_Cluster_52','gyrB_Cluster_56','gyrB_Cluster_17','gyrB_Cluster_22','gyrB_Cluster_44','gyrB_Cluster_54','gyrB_Cluster_18','gyrB_Cluster_50','gyrB_Cluster_63','gyrB_Cluster_89','gyrB_Cluster_19','gyrB_Cluster_145','gyrB_Cluster_176','gyrB_Cluster_95','gyrB_Cluster_96','gyrB_Cluster_101','gyrB_Cluster_152','gyrB_Cluster_161','gyrB_Cluster_201','gyrB_Cluster_124'))
new_rdlgyrB <- merge_phyloseq(my_gyrB_subset, tax_table(rdlgyrB), sample_data(rdlgyrB))
new_rdlgyrB

sample.order.gyrB<-c("CF1_GYRB","CF2_GYRB","CF3_GYRB",
                    "SF1_GYRB","SF2_GYRB","SF3_GYRB",
                    "GB1_GYRB","GB2_GYRB","GB3_GYRB",
                    "CS1_GYRB","CS2_GYRB","CS3_GYRB",
                    "PS1_GYRB","PS2_GYRB","PS3_GYRB",
                    "MC1_GYRB","MC2_GYRB","MC3_GYRB","MC4_GYRB","MC5_GYRB")
otus.colors<-c("#003b7d","#003b7d","#003b7d","#003b7d","#003b7d","#003b7d",
               "#003b7d","#003b7d","#003b7d","#003b7d","#003b7d","#003b7d",
               "#003b7d","#003b7d","#003b7d","#003b7d","#003b7d","#003b7d",
               "#cb2027","#cb2027","#cb2027","#cb2027","#cb2027","#cb2027",
               "#cb2027","#cb2027","#cb2027","#cb2027","#cb2027","#cb2027",
               "#cb2027","#cb2027","#cb2027","#cb2027","#cb2027","#cb2027",
               "#cb2027","#cb2027","#cb2027","#cb2027","#cb2027","#cb2027","#cb2027","#cb2027")

fontaxis.otus <- element_text(color = otus.colors, size = 10, hjust=1, face="italic")
fontaxis.samples <- element_text(color = "black", size = 10)
p<-plot_heatmap(new_rdlgyrB, taxa.label= "Species", taxa.order = "Species", sample.order = sample.order.gyrB)
p <- p + facet_grid(~EnvType, scales = "free", space = "free")
p <- p + scale_fill_gradient2(low = "#1a9850", mid = "#feff9a", high = "#d73027", 
                              na.value = "white", trans = log_trans(10), 
                              midpoint = log(100, base = 10))
p<- p + theme(axis.text.y = fontaxis.otus) + theme(axis.text.x = fontaxis.samples)
plot(p)


####################
# to produce the final version, a new phyloseq object is done to merge 16S-gyrB samples and to introduce fake otus for blank spaces#
####################

#The new phyloseq dataset can be constructed from the POIRIER_ET_AL_FOODMICROBIOME_DATASET.xlsx and subsetting the corresponding OTUs

otuTable3 <- as.matrix(read.csv("otu_table_redlosses_v3.tsv", sep = "\t", row.names = 1))
taxTable3 <- as.matrix(read.csv("tax_table_redlosses_v3.tsv", sep = "\t", row.names = 1))
sampleData3 <- read.csv("sample_data_redlosses_v3.tsv", sep = "\t", row.names = 1)
rdldata3 <- phyloseq(otu_table(otuTable3, taxa_are_rows = TRUE), 
                    tax_table(taxTable3), 
                    sample_data(sampleData3))


sample.order.merged<-c("CF1","CF2","CF3",
                     "SF1","SF2","SF3",
                     "GB1","GB2","GB3",
                     "CS1","CS2","CS3",
                     "PS1","PS2","PS3",
                     "MC1","MC2","MC3","MC4","MC5")

# 16S and gyrB will be plotted separately and with order #

otu.order.merged<-c("16S_Cluster_20","16S_Cluster_5","16S_Cluster_11","16S_Cluster_16","16S_Cluster_13","16S_Cluster_30",
                    "16S_Cluster_2","16S_Cluster_33","16S_Cluster_22",
                    "16S_Cluster_7","16S_Cluster_19","16S_Cluster_1","16S_Cluster_3","16S_Cluster_8","16S_Cluster_10","16S_Cluster_6",
                    "gyrB_Cluster_f1","gyrB_Cluster_f2",
                    "gyrB_Cluster_95","gyrB_Cluster_96","gyrB_Cluster_201",
                    "gyrB_Cluster_12","gyrB_Cluster_43","gyrB_Cluster_25","gyrB_Cluster_60","gyrB_Cluster_42",
                    "gyrB_Cluster_14","gyrB_Cluster_52","gyrB_Cluster_56",
                    "gyrB_Cluster_18","gyrB_Cluster_50","gyrB_Cluster_89","gyrB_Cluster_123","gyrB_Cluster_63","gyrB_Cluster_19","gyrB_Cluster_145","gyrB_Cluster_176",
                    "gyrB_Cluster_2","gyrB_Cluster_3","gyrB_Cluster_16","gyrB_Cluster_30","gyrB_Cluster_59",
                    "gyrB_Cluster_6","gyrB_Cluster_11","gyrB_Cluster_21","gyrB_Cluster_39","gyrB_Cluster_27",
                    "gyrB_Cluster_4","gyrB_Cluster_13","gyrB_Cluster_23","gyrB_Cluster_7","gyrB_Cluster_141",
                    "gyrB_Cluster_26","gyrB_Cluster_37",
                    "gyrB_Cluster_17","gyrB_Cluster_22","gyrB_Cluster_44","gyrB_Cluster_54",
                    "gyrB_Cluster_8","gyrB_Cluster_32","gyrB_Cluster_46"
                    )

C1<-("white")
C2<-("black")
C3<-("red")

otus.colors<-c(C1,
               C3,C2,C2,
               C3,C2,C2,C2,C2,
               C3,C2,C2,
               C3,C2,C2,C2,C2,C3,C2,
               C2,C2,C3,C2,
               C2,C3,C2,C2,
               C3,C3,C2,C2,C2,
               C3,C3,
               C2,C3,C2,C2,
               C2,C3,
               C1,C1,
               C2,C2,C2,C2,C2,C2,
               C2,C2,
               C2,C2,C2,C2,C2)


fontaxis.otus <- element_text(color = "black", size = 12, hjust=1, face="italic")
fontaxis.samples <- element_text(color = "black", size = 15)
p<-plot_heatmap(rdldata3, taxa.label= "Species", taxa.order = rev(otu.order.merged), sample.order = sample.order.merged)
p <- p + facet_grid(~EnvType, scales = "free", space = "free")
p <- p + scale_fill_gradient2(low = "#1a9850", mid = "#feff9a", high = "#d73027", 
                              na.value = "white", trans = log_trans(10), 
                              midpoint = log(100, base = 10))
p<- p + theme(axis.text.y = fontaxis.otus) + theme(axis.text.x = fontaxis.samples)
p<- p+ theme(strip.text.x = element_text(size = 12))
p<- p+ theme(legend.text=element_text(size=14),legend.key.height = unit(1.0, "cm"), legend.key.width = unit(0.75, "cm"))
plot(p)
