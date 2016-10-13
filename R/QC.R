#' QC
#'
#' In this quality control function, the data is checked for the presence of all specified columns in the correct order. Data is then filtered to remove peptides without any valid (i.e. < 0) intensity measurement and classified based on the possibility to perform exact or approximate stoichiometry and abundance calculations according to the modification specified by the user (e.g. "phosphorylation" or "acetylation"). The filtered and classified data is exported as a CSV file into the working directory (file name), along with a bar graph summarizing the number of peptides and proteins remaining after each step of quality control.
#' @param data data frame with 12 columns as decribed in documentation.
#' @param mod If any chemoform of a given peptide bears the modification specified in mod="" (f.ex. mod="phosphorylation") all peptides with the same sequence and containing the modification will be classified as Modified, while all peptides with the same sequence and not containing the modification will be classified as Not-Modified. If none of the chemoforms of a peptide contains the specified modification, the peptide is classified Abundance.
#' @export
#' @return data frame of filtered data
#' @examples 
#' QC(testData,mod="Phosphorylation")

###################################################################################
################################ QUALITY CONTROL ##################################
###################################################################################
# Quality Control takes in a Skyline or a MaxQuant input file. Inputs need the following columns 
# in the specified order even if the entries are left empty for the particular experiment: "Protein", 
# "Peptide", "Pathway","Modification","Position", "Gene.Name", "Protein.Name", "Condition", "PatientID", 
# "Label", "Run", "Intensity".
# Data is filtered, with entries that have zero intensity across all replicates and 
# all conditions, along with those with missing (NA) intensities removed. Data is then
# classified based on whether exact or approximate stoichiometry and abundance calculations
# can be performed. 

QC <- function(data, mod){
  data$X=NULL
  
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
  
  data.dup <- ldply(unique(data$EntryID), function(i){
    if(nrow(data[data$EntryID==i,])>1){
      temp<- subset(data,EntryID==i)
      temp$Intensity = sum(temp$Intensity)
      return(temp[1,])
    }
    else{
      return(data[data$EntryID==i,])
    }
  }, .progress="text")
  
  
  # Select proteins with peptides containing the modification of interest or unmodified peptides
  
  data.mod <- ldply(unique(data.dup$Protein), function(i){
    temp<- subset(data.dup,Protein==i)
    if((sum(temp$Modification==" ")==length(temp$Modification)|any(grepl(temp$Modification, pattern = mod)))==TRUE){
      return(temp)
    }
  }, .progress="text")
  
  tracker = cbind(tracker, c(length(unique(data.mod$Protein)),length(unique(data.mod$Peptide))))
  
  # Report if all intensity values are zero and remove from further analysis
  
  data.mod.int <- ldply(unique(data.mod$Peptide), function(i){
    temp <- subset(data.mod,Peptide==i)
    if(sum(is.na(temp$Intensity)) | sum(temp$Intensity==0)<(nrow(temp))){
      return(temp)
    }
  }, .progress="text")
  
  tracker = cbind(tracker, c(length(unique(data.mod.int$Protein)),length(unique(data.mod.int$Peptide))))
  
  # write.csv(data.mod.int, "preclass")
  # Classify
  data.mod.int$Modification <- as.character(data.mod.int$Modification)
  data.mod.int$Abundance <- NA
  data.mod.int$Stoichiometry <- NA
  data.mod.int$EntryID = paste(data.mod.int$Protein,data.mod.int$Condition, data.mod.int$PatientID, data.mod.int$Label, data.mod.int$Run, sep = "_")
  
  classified.data <- ldply(unique(data.mod.int$EntryID),function(i){
    temp<- subset(data.mod.int,EntryID==i)
    # CASE 1: Single Peptide for the Protein
    if(nrow(temp)<2){
      # If it is alone and "unmodified" use for quantitation and exact abundance
      if(temp[,"Modification"]=="Unmodified"){
        temp[,"Abundance"] = "Exact"
        temp[,"Stoichiometry"] = "Zero"
        return(as.data.frame(temp))
      }else{
        #Otherwise both abundance and stoichiometry will be approximate
        temp[,"Abundance"] = "Approx"
        temp[,"Stoichiometry"] = "Hundred"
        return(as.data.frame(temp))
      }
    }
    #CASE 2: More than 2 peptides per protein
    if(nrow(temp)>=2){
      return(ldply(unique(temp$Peptide),function(j){
        #Subset into unique peptides
        temp2 <- subset(temp, Peptide==j)
          # If the peptide only exists in unmodified form 
        if(nrow(temp2)==1){
          if(temp2$Modification=="Unmodified"){
            temp2[,"Stoichiometry"]="Zero"
            temp2[,"Abundance"]="Exact"
            
          } else{
            temp2[,"Stoichiometry"]="Hundred"
            temp2[,"Abundance"]="Approx"
          }
          return(as.data.frame(temp2))
        } else if(any(temp2$Modification=="Unmodified")){
            temp2[,"Stoichiometry"] = "Exact"
            temp2[,"Abundance"] = "Exact"
            return(as.data.frame(temp2))
            # Otherwise, return as is (they are used as modified) and label stoichiometry
            # and abundance as approximate 
          } else {
            temp2[,"Stoichiometry"] = "Approx"
            temp2[,"Abundance"] = "Approx"
            # temp2[1,]$Intensity = sum(temp2$Intensity)
            return(as.data.frame(temp2))
          }
    
      }))
    }
  
    },.progress="text")
  
  # Anything remaining will be labeled as modified "M"
  # classified.data[which(classified.data$Modification!="Q" & classified.data$Modification!="U"  ),"Modification"] <- "M"
  
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
  tracker=cbind(c("Raw Data", "Selection of proteins with PTM of interest", 
                  "Removal of undetected proteins", 
                  "Classification"), tracker)
  
  rownames(tracker)=NULL
  colnames(tracker)=c("process", "protein", "peptide")
  tracker$process=factor(tracker$process, levels=c("Raw Data", "Selection of proteins with PTM of interest", 
                                                   "Removal of undetected proteins", 
                                                   "Classification"))
  
  tracker=melt(tracker, id.vars="process")
  
  ggplot(data=as.data.frame(tracker), aes(x=process, y=value, fill=variable)) + 
    geom_bar(stat="identity", position="dodge") +
    scale_fill_brewer(palette="Set1") +
    xlab(" ") + 
    ylab("Count") +
    ggtitle("Number of Proteins and Peptides after Each Filtering Step") +
    theme(plot.title = element_text(size=12),
          axis.text.x = element_text(size=12,angle=90, hjust=1),
          axis.text.y = element_text(size=12), 
          legend.text=element_text(size=12), 
          axis.title = element_text(size=12))
  ggsave(file="Filtering_tracker.pdf")
  
  return(classified.data)
}