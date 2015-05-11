library(ggplot2)
# library("plyr")
# library(reshape2)
# setwd("Z:././././04_projects/Quantitative_Modeling//Training//QCPA-trainingset//TestData")
########################## Kinetics of PTMs #########################################

TestData <- read.table("tp_data.csv", sep=',', header = T)

# This case has one protein

# %Phosphorylation for each pep will be calculated as abundance of P state over
# abundance of F state. %Phosphorylation of protein will be calculated as the 
# sum of the two states. Total abundance will be calculated by Q 

# % Phospho peptide
#Normalize heavy to light
TestData=TestData[order(TestData$Condition),]

timepoints <- unique(TestData$Condition)
data.p <- data.frame()
for (i in timepoints) { 
  norm.p <- TestData[which(TestData$Condition==i & TestData$Modification=="P"  & TestData$Label=="L"),]$AUC/
    TestData[which(TestData$Condition==i & TestData$Modification=="P" & TestData$Label=="H"),]$AUC  
  norm.np <- TestData[which(TestData$Condition==i & TestData$Modification=="F" & TestData$Label=="L"),]$AUC/
    TestData[which(TestData$Condition==i & TestData$Modification=="F" & TestData$Label=="H"),]$AUC
  x <- rbind(TestData[which(TestData$Condition==i & TestData$Modification=="P" & 
                                        TestData$Label=="L"),], 
                       TestData[which(TestData$Condition==i & TestData$Modification=="F" & 
                                        TestData$Label=="L"),])
  x$normAUC <- c(norm.p, norm.np)
  data.p <- rbind(data.p, x)
}
rm(x)
row.names(data.p) =NULL
data.p$Label = NULL
data.p$AUC = NULL
data.p$X = NULL

#Calculate percent phospho for each peptide

phospho <- data.frame()
for (i in timepoints){
  x <- data.p[which(data.p$Condition==i & data.p$Modification=="P"),]$normAUC/
    (data.p[which(data.p$Condition==i & data.p$Modification=="F"),]$normAUC+
       data.p[which(data.p$Condition==i & data.p$Modification=="P"),]$normAUC)
  y <- data.frame(data.p[which(data.p$Condition==i & data.p$Modification=="P"),])
  y$perc.p <- x
  phospho <- rbind(phospho, y)
                      
}
rm(x)
rm(y)
row.names(phospho) =NULL

#Calculate percent phospho for each protein
num.proteins <- length(unique(data.p$Protein))
proteins <- unique(data.p$Protein)

data.prot<- data.frame( )
for (i in 1:num.proteins){
  x <- data.p[grep(proteins[i],data.p$Protein),]
  data.prot <- rbind(c(proteins[i], (x[which(x$Modification=="P"),]$AUC/x[which(x$Modification=="F"),]$AUC) ))
  
}

# Abundance of protein (log-transformed)
norm.q <- log(TestData[which(TestData$Modification=="Q" & TestData$Label=="L"),]$AUC/
  TestData[which(TestData$Modification=="Q" & TestData$Label=="H"),]$AUC)
norm.q = data.frame(TestData[which(TestData$Modification=="Q" & TestData$Label=="L"),], norm.q)
names(norm.q)[length(norm.q)]= "Abundance"


runs <- length(unique(data.p$Run))
for (i in 1:num.proteins){
  abund <-data.frame(rep(unique(data.p$Protein)[i],length.out=(dim(norm.q)[1]/runs)),
                         sapply(timepoints, function(x) mean(norm.q[which(norm.q$Condition==x), ]$Abundance)))
  colnames(abund) <- c("Protein", "Abundance")
}


# Plot for each tp for each peptide

peptides <- length(unique(phospho$Peptide))

for (i in 1:peptides){
  
  y1 <- phospho[grep(unique(phospho$Peptide)[i], phospho$Peptide),]
  y1 <- y1[order(y1$Condition),]
  y1$Condition=as.factor(y1$Condition)
  y1$abundance = rep(abund$Abundance, each = runs)
  p<-ggplot(y1,aes(x=Condition,y=perc.p,group=Condition,fill=abundance))+
    
    geom_boxplot()+ ggtitle(paste("Percent Phosphorylation for", unique(data.p$Peptide)[i]))+ 
    theme(plot.title = element_text(lineheight=.8, face="bold")) +
    xlab("Timepoint")+ ylab("Phosphorylated Fraction")+ 
    scale_fill_continuous(name="log(Protein Abundance)")
  show(p)
  #ggsave(filename = "boxplot.jpg")
  #ggplot(y1,aes(x=Condition,y=perc.p,group=Condition, 
  #                color=Condition))+geom_jitter(position=position_jitter(width=0.1, height=0))
}




