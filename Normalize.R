library(ggplot2)
source("QC.R")

# setwd("Z:././04_projects/Quantitative_Modeling//ProteoModlr")
#data <- read.csv("Emandel_input", sep=',', header=T, na.strings=" ")
data <- read.csv("molten_filtered_ 2015-05-27_19_55_08", sep=',', header=T)
###################################################################################
################################ NORMALIZATION ####################################
# Normalyze takes in the data along with 3 options. If no options are chosen, data
# is returned unchanged. If the normalyzation step is carried out, the normalized
# intensity values replace the raw intensity values in the Intensit column.
# iso.norm takes in a character string corresponding to the reference isotopologue
# to which the other isotopologue will be normalized. internal.norm takes in a 
# vector of character strings corresponding to the peptide sequence of the 
# reference peptides to which other peptides are normalized. If more than one
# reference peptide is indicated, an arithmetic mean is take. ref.state takes in
# a character sting indicating the reference condition to which the intensity of the other conditions will be 
# normalyzed. 

Normalyze <- function(data, iso.norm = "", internal.norm = "", ref.state = "", mod="phospho", skyline=F, MaxQuant=T){

  data <- QC(data, mod, skyline, MaxQuant) # Send data to QC function 
  
  ########################## NORMALYZING BY AN ISOTOPOLOGUE ######################
  if(nchar(iso.norm) != 0){ 
    # Check to see if there are 2 labels 
    if(length(unique(data$Label))==1){
      stop(message("Only one isotopologue found. Data must contain two isotopologues 
                   in the Label column for this type of normalization. \n"))
    }
    ref = iso.norm #reference isotope
    
    if(any(data$Label==ref)==F){
      stop(message("Reference isotope not found in the Label column of data. 
                   Match spelling and case to format used in data. \n"))
    }
    
    exp = unique(data$Label)[-which(match(unique(data$Label),ref)==1)][1] #non-ref isotope
    
    #Normalize heavy to light
    data$Label = factor(data$Label, levels=c(exp, ref))
    data = data[order(data$Label),]
    #EntryID is a unique identifier of each isotopologue for each peptide
    data$EntryID = paste(data$Protein,data$Peptide,data$Modification,data$Condition,data$PatientID,data$Run,sep = "_")
    
    ratios=sapply(unique(data$EntryID), function(i) {
      ind=grep(i,data$EntryID)
      if(length(ind)!=2){
        stop(message(sprintf("EntryID %s missing a labeled or unlabeled state", ind)))
      }
      ratio=data$Intensity[ind[1]]/data$Intensity[ind[2]]
      })    
    
    # Cut data in half (leaving out one isotopologue) and store normalized intensities
    data=data[,(1:(0.5*length(data)))]
    data$Intensity = ratios
  }
  
  #################### NORMALYZING BY AN INTERNAL REFERENCE ######################
  
  if(length(internal.norm) != 0){
    
      # Check to see if reference peptides exist. Retain only those that do. 
      ref = sapply(internal.norm, function(i){
      if(any(as.character(data$Peptide)==i)){
        return(i)
      }
      else{
        cat(sprintf("Reference peptide %s not found in the dataset provided. \n", i))
      }
    })
    
    if(length(ref)==0){
      stop(message("No reference peptide found. \n"))
    }
    
    # Store the references in a separate data frame only if all intensities for each are non-zero
    # and plot the intensities over the conditions
    data$Condition = factor(data$Condition)
    data=data[order(data$Condition),]
    ref.df <- ldply(ref, function(i){
      if(length(which(data[grep(i, data$Peptide),]$Intensity==0))==0)
      ref.df=data[grep(i, data$Peptide),]
      ref.df = ref.df[order(ref.df$Condition),]
      ggplot(ref.df, aes(x=Condition,y=Intensity))+
        geom_boxplot()+ ggtitle(paste("Intensity of reference peptide ", i))+ 
        theme(plot.title = element_text(lineheight=.8, face="bold")) +
        xlab("Condition")+ ylab("Intensity")
      ggsave(file=sprintf("reference_%s.pdf", i))
      return(ref.df)
    })
    
    ref.gm.mean <- ldply(unique(ref.df$Condition), function(i){
       ref.gm.mean = exp(sum(log(ref.df$Intensity[ref.df$Condition==i])) / 
            length(ref.df$Intensity[ref.df$Condition==i]))
       return(ref.gm.mean)
    })
    ref.gm.mean = data.frame(t(ref.gm.mean), row.names=NULL)
    colnames(ref.gm.mean)=unique(ref.df$Condition)
    
    data=ldply(1:length(unique(data$Condition)), function(i){
      temp = subset(data, Condition==unique(data$Condition)[i])
      temp$Intensity = temp$Intensity/as.numeric(ref.gm.mean[i])
      return(temp)
    })
    
    
  }
    
  #################### NORMALYZING BY AN REFERENCE STATE ######################  
  
  if(length(ref.state) != 0){
    ref = ref.state
    levels = c(ref, as.character(unique(data$Condition)[which(unique(data$Condition)!=ref)]))
    
    # Check for existence and validity of reference state (i.e. only one such state)
    if(any(as.character(data$Condition)==ref)==F | length(ref)>1){
      stop(message(sprintf("Invalid or non-existant reference state %s. \n", ref)))
    }
  
    #EntryID is a unique identifier of each sample for each peptide
    data$EntryID = paste(data$Protein,data$Peptide,data$Modification,data$PatientID,data$Run, sep = "_")
    
    data$Condition = factor(data$Condition, levels=levels)
    data = data[order(data$Condition),]
    
    data.n = ldply(unique(data$EntryID), function(i){
      temp<- subset(data, EntryID==i)
      intensity = sapply(1:nrow(temp), function(k){
        intensity = temp$Intensity[k]/ temp$Intensity[1]
      })
      
      temp$Intensity=intensity
      
      return(temp) 
    }, .progress = "text") 
  }
  
  
  return (data)

}
  

