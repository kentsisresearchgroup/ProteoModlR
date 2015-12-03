
library(plyr)
library(ggplot2)
library(reshape2)
###################################################################################
################################ QUALITY CONTROL ##################################
###################################################################################
# Quality Control takes in a Skyline or a MaxQuant input file. Inputs in the form
# of a Skyline output need the following columns in the specified order even if 
# the entries are left empty for the particular experiment: Protein (list of UniProtID's),
# Peptide (list of peptide sequences), Modification (type of modification for each peptide 
# in UniMod format), Gene.Name (GO), Protein.name (name of protein), Condition (time or
# study group), PatientID (unique ID corresponding to biological replicate), Label (isotope
# label if applicable, otherwise label all as L), Run (technical replicate) and Intensity.
# For MaxQuant output, first 5 columns must be "Protein", "Peptide", "Modification", 
# "Gene.Name" and "Protein.Name" in the specified order. Intensities corresponding to 
# each peptide should be labled with a header containing one-letter desgnators in the 
# following format: BioReplicate_Condition_Level_Run. For example, intensity values for
# the 3rd heavy-labeled sample of patient A in the 1st timepoint will be in a column
# labeled with a header "A_1_H_3."
# Data is then filtered, with entries that have zero intensity across all replicates and 
# all conditions, along with those with missing (NA) intensities removed. Data is then
# classified based on whether exact or approximate stoichiometry and abundance calculations
# can be performed. 

QC <- function(data, mod){
  data$X=NULL
  # Data must be specified to be either in a MaxQuant or in a Skyline output format
#   if (sum(skyline,MaxQuant)>1){
#     stop(message("Data should be formatted in either the form of Skyline output or that of MaxQuant. \n"))
#   }
  
  #For data in skyline format
  # if(skyline==T){
    # Set column names to what is required for the remainder of the analysis
    if(ncol(data)!=12){
      stop(message("Data in the form of Skyline output should include 12 columns (even if empty) in the following order:Protein, Peptide, Modification, Gene.Name,Protein.Name, Condition, PatientID, Label, Run, Intensity. \n"))
    }
    colnames(data) <- c("Protein", "Peptide", "Pathway","Modification","Position", "Gene.Name",
                              "Protein.Name", "Condition", "PatientID", "Label", "Run", "Intensity")
    #Select one (first) uniprotID for each peptide
    if (sum(nchar(as.character(data$Protein))<6)>0){
      stop(message("Use standard protein names(e.g. Uniprot-ID). \n"))
    }
    data$Protein <- substr(data$Protein, 1, 6) #Take first protein identifier
    
    # Print and save the number of proteins and peptides before filtering 
    
    tracker = matrix(c(length(unique(data$Protein)),length(unique(data$Peptide))), nrow=2)
    
    # Checking for Duplicates: If identical samples have appear, sum intensitites
    data$EntryID = paste(data$Protein, data$Peptide, data$Modification, data$Condition, data$PatientID, data$Label, data$Run, sep = "_")
    
    data.dup <- ldply(1:length(unique(data$EntryID)), function(i){
      if(nrow(data[data$EntryID==unique(data$EntryID)[i],])>1){
        temp<- subset(data,EntryID==unique(data$EntryID)[i])
        temp$Intensity = sum(temp$Intensity)
        return(temp[1,])
      }
      else{
        return(data[data$EntryID==unique(data$EntryID)[i],])
      }
    }, .progress="text")
    
    #Cast into Wide format
#     data.dup$EntryID = paste(data.dup$PatientID,data.dup$Condition,data.dup$Run, data.dup$Label, sep = "_")
#     data.dup = data.dup[,!(colnames(data.dup) %in% c("PatientID","Run","Condition","Label"))]
#     data = dcast(data.dup,Protein+Peptide+Pathway+Position+Modification+Gene.Name+Protein.Name~EntryID, value.var="Intensity")
#   
  
  #For data in MaxQuant format
#   if(MaxQuant==T){
#     
#     #Change column names 
#     colnames(data)[1:7] <-c("Protein", "Peptide", "Pathway","Modification", "Position","Gene.Name", "Protein.Name")
#     if (sum(nchar(as.character(data$Protein))<6)>0){
#       stop(message("Use standard protein names(e.g. Uniprot-ID). \n"))
#     }
#     #Take first UniprotID if multiple ID's are provided 
#     data$Protein <- substr(data$Protein, 1, 6) #Take first protein identifier
#     # Print and save the number of proteins and peptides before filtering 
#     
#     tracker = matrix(c(length(unique(data$Protein)),length(unique(data$Peptide))), nrow=2)
#   }
  
  # Select modification to look at
  modifications<- c(mod, "unmodified") 
  
  data.mod<- data.dup[grep(paste(modifications,collapse="|"),data.dup$Modification, ignore.case = T),]
  
  tracker = cbind(tracker, c(length(unique(data.mod$Protein)),length(unique(data$Peptide))))

  # Report if all intensity values are zero and remove from further analysis
 
  data.mod.int <- ldply(unique(data.mod$Peptide), function(i){
    temp <- subset(data.mod,Peptide==i)
      if(sum(is.na(temp$Intensity)) | sum(temp$Intensity==0)<(nrow(temp))){
      return(temp)
    }
  }, .progress="text")
  
  tracker = cbind(tracker, c(length(unique(data.mod.int$Protein)),length(unique(data.mod.int$Peptide))))
  
  #Remove intensity values that are NA's
  # data.mod.int = na.omit(data.mod.int)
  
  # tracker = cbind(tracker, c(length(unique(data.mod.int$Protein)),length(unique(data.mod.int$Peptide))))
  
  write.csv(data.mod.int, "preclass")
  # Classify
  data.mod.int$Modification <- as.character(data.mod.int$Modification)
  data.mod.int$Abundance <- NA
  data.mod.int$Stoichiometry <- NA
  data.mod.int$EntryID = paste(data.mod.int$Protein, data.mod.int$Peptide, data.mod.int$Condition, data.mod.int$PatientID, data.mod.int$Label, data.mod.int$Run, sep = "_")
  
  classified.data <- ldply(unique(data.mod.int$EntryID),function(i){
    temp<- subset(data.mod.int,EntryID==i)
    # CASE 1: Single Peptide for the Protein
    if(nrow(temp)<2){
      # If it is alone and "unmodified" use for quantitation and exact abundance
      if(temp[,"Modification"]=="Unmodified"){
        temp[,"Modification"]="U"
        temp[,"Abundance"] = "Exact"
        return(as.data.frame(temp))
      }else{
        #Otherwise both abundance and stoichiometry will be approximate
        temp[,"Abundance"] = "Approx"
        temp[,"Stoichiometry"] = "Approx"
        return(as.data.frame(temp))
      }
    }
    #CASE 2: More than 2 peptides per protein
    if(nrow(temp)>=2){
      return(ldply(unique(temp$Peptide),function(j){
        #Subset into unique peptides
        temp2 <- subset(temp, Peptide==j)
        # If the unique sequence is unpaired 
        if(nrow(temp2)<2){
          # and if peptide is unmodified, use for quantitation and abundance is exact
          if(temp2[,"Modification"]=="Unmodified"){
            temp2[,"Modification"]="U"
            temp2[,"Abundance"] = "Exact"
            return(as.data.frame(temp2))
            #Otherwise both stoichiometry and abundance can be approximated 
          }else{
            temp2[,"Stoichiometry"] = "Approx"
            temp2[,"Abundance"] = "Approx"
            return(as.data.frame(temp2))
          }
          # If the unique sequence has a pair(s)
        }else{
          # If any of them are labled unmodified, use as U for exact stoichiometry
          if(any("Unmodified"==temp2$Modification)){
            temp2[which(temp2$Modification=="Unmodified"),"Modification"]="U"
            temp2[,"Stoichiometry"] = "Exact"
            temp3 <- subset(temp2, Modification!="U")
            temp3[1,]$Intensity = sum(temp3$Intensity)
            temp4 = rbind(temp2[which(temp2$Modification=="U"), ], temp3[1,])
            return(as.data.frame(temp4))
            # Otherwise, return as is (they are used as modified) and label stoichiometry
            # and abundance as approximate 
          } else {
            temp2[,"Stoichiometry"] = "Approx"
            temp2[,"Abundance"] = "Approx"
            temp2[1,]$Intensity = sum(temp2$Intensity)
            return(as.data.frame(temp2[1,]))
          }
        }
      }))
    }
  },.progress="text")

  # Anything remaining will be labeled as modified "M"
  classified.data[which(classified.data$Modification!="Q" & classified.data$Modification!="U"  ),"Modification"] <- "M"
  
  # Save file 
  t<- Sys.time()
  t=gsub(" ", "_", t)
  t=gsub(":", "_", t)
  write.csv(classified.data, paste("filtered_", t))
  
  tracker = cbind(tracker, c(length(unique(classified.data$Protein)),length(unique(classified.data$Peptide))))
  
  
  # Melt into MSsta/Skyline format and store
  
  # data.melt <- melt(classified.data, id.vars=c("Protein", "Peptide", "Pathway","Modification", "Position","Abundance", "Stoichiometry", "Gene.Name", "Protein.Name" ))
  

#   data.melt$PatientID <- substr( data.melt$variable,1,1)
#   data.melt$Condition <- substr( data.melt$variable,3,3)

  
#   if(nchar(as.character(data.melt$variable[1]))==5){
#     data.melt$Run <- substr(data.melt$variable,5,5)
#     
#   }else if(nchar(as.character(data.melt$variable[1]))==7){
#     data.melt$Label <- substr(data.melt$variable,7,7)
#   }else{
#     data.melt$Label <- "L"
#     data.melt$Run <- "1"
#   }
  
  
#   data.melt$Intensity <- data.melt$value
#   data.melt$value = NULL
#   data.melt$variable =NULL
#   
#   t<- Sys.time()
#   t=gsub(" ", "_", t)
#   t=gsub(":", "_", t)
#   write.csv(data.melt, paste("molten_filtered_", t))
  
  ################################TRACKING FILTERING STEPS##########################
  tracker=as.data.frame(t(tracker))
  tracker=cbind(c("Raw Data", "After selecting for PTM of interest", 
                  "After removing peptides with zero intensity values", 
                  "After classification and aggregation"), tracker)
  
  rownames(tracker)=NULL
  colnames(tracker)=c("process", "protein", "peptide")
  tracker$process=factor(tracker$process, levels=c("Raw Data", "After selecting for PTM of interest", 
                                                   "After removing peptides with zero intensity values", 
                                                   "After classification and aggregation"))
  
  tracker=melt(tracker, id.vars="process")
  
  ggplot(data=as.data.frame(tracker), aes(x=process, y=value, fill=variable)) + 
    geom_bar(stat="identity", position="dodge") +
    scale_fill_brewer(palette="Set1") +
    xlab("Filtering Step") + 
    ylab("count") +
    ggtitle("Number of Proteins and Peptides after Each Filtering Step") +
    theme(plot.title = element_text(size=12),
          axis.text.x = element_text(size=12,angle=90, hjust=1),
          axis.text.y = element_text(size=12), 
          legend.text=element_text(size=12), 
          axis.title = element_text(size=12))
  ggsave(file="Filtering_tracker.pdf")
  
  return(classified.data)
}
