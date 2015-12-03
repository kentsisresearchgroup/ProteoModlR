library(ggplot2)
library(reshape2)
library(plyr)




########################## ABUNDANCE AND STOICHIOMETRY ##################################
# Analyze.R calculates the approximate and/or exact site occupancy and/or abundance 
# of the peptides. The choice between exact and approximate calculations will be
# based on user input, while the selection of appropriate peptides is based on labeling
# done in QC.R. If a reference state is specified, stoichiometric ratios and abundances
# for all conditions are normalized to the reference state. 

Analyze <- function(data, stoich, abund, ref.state){
  data$X <- NULL
  data$log2FC.Abundance = NA
  data$Occupancy = NA 
  if(stoich=="Exact"){
    data$EntryID = paste(data$Protein,data$Peptide,data$Condition,data$PatientID,data$Run, sep = "_")
    temp= subset(data, Stoichiometry=="Exact")
    exact.mod = ldply(unique(temp$EntryID), function(i){
      temp2 = subset(temp, EntryID==i)
      stoich = temp2[which(temp2$Modification=="M"),]$Intensity / 
        (temp2[which(temp2$Modification=="M"),]$Intensity+
           temp2[which(temp2$Modification=="U"),]$Intensity)
      temp2[which(temp2$Modification=="M"),]$Occupancy = stoich
      return(temp2[which(temp2$Modification=="M"),])
    }, .progress='text')
    rm(temp)
    
    #Normalize to reference state
    if(length(ref.state) != 0){
      ref = ref.state
      levels = c(ref, as.character(unique(exact.mod$Condition)[which(unique(exact.mod$Condition)!=ref)]))
      
      # Check for existence and validity of reference state (i.e. only one such state)
      if(any(as.character(exact.mod$Condition)==ref)==F | length(ref)>1){
        stop(message(sprintf("Invalid or non-existant reference state %s. \n", ref)))
      }
      
      #EntryID is a unique identifier of each sample for each peptide
      exact.mod$EntryID = paste(exact.mod$Protein,exact.mod$Peptide,exact.mod$Pathway,exact.mod$Modification,exact.mod$PatientID,exact.mod$Run, sep = "_")
      
      exact.mod$Condition = factor(exact.mod$Condition, levels=levels)
      
      exact.mod = exact.mod[order(exact.mod$Condition),]
      
      exact.mod = ldply(unique(exact.mod$EntryID), function(i){
        temp<- subset(exact.mod, EntryID==i)
        norm.occ = sapply(1:nrow(temp), function(i){
          norm.occ = log2(temp$Occupancy[i]/ temp$Occupancy[1])
        })
        
        temp$norm.occ=norm.occ
        
        return(temp) 
      }, .progress = "text") 
    }
    
    # Save file 
    t<- Sys.time()
    t=gsub(" ", "_", t)
    t=gsub(":", "_", t)
    write.csv(exact.mod, paste("Exact_Mod_", t))
    
  }

  if(stoich=="Approx"){
    data$EntryID = paste(data$Protein,data$Condition,data$PatientID,data$Run, sep = "_")
    temp=subset(data[-which(data$Intensity==0),], (Abundance=="Exact" | Stoichiometry=="Approx"))
    approx.mod = ldply(unique(temp$EntryID),function(i){
      temp2<- subset(temp,EntryID==i)
      if(any(temp2$Abundance=="Exact")){
        stoich = temp2[which(temp2$Modification=="M"),]$Intensity / 
          (exp(sum(log(temp2[which(temp2$Modification=="Q"),]$Intensity))/
             nrow(temp2[which(temp2$Modification=="Q"),])))
        temp2[which(temp2$Modification=="M"),]$Occupancy = stoich
        return(as.data.frame(temp2[which(temp2$Modification=="M"),]))
      }else{
        temp2$Occupancy = temp2$Intensity
        return(as.data.frame(temp2))
      }
      
    }, .progress="text")
    rm(temp)
    
    if(length(ref.state) != 0){
      ref = ref.state
      levels = c(ref, as.character(unique(approx.mod$Condition)[which(unique(approx.mod$Condition)!=ref)]))
      
      # Check for existence and validity of reference state (i.e. only one such state)
      if(any(as.character(approx.mod$Condition)==ref)==F | length(ref)>1){
        stop(message(sprintf("Invalid or non-existant reference state %s. \n", ref)))
      }
      
      #EntryID is a unique identifier of each sample for each peptide
      approx.mod$EntryID = paste(approx.mod$Protein,approx.mod$Peptide,approx.mod$Modification,approx.mod$PatientID,approx.mod$Run, sep = "_")
      
      approx.mod$Condition = factor(approx.mod$Condition, levels=levels)
      
      approx.mod = approx.mod[order(approx.mod$Condition),]
      
      approx.mod = ldply(unique(approx.mod$EntryID), function(i){
        temp<- subset(approx.mod, EntryID==i)
        norm.occ = sapply(1:nrow(temp), function(i){
          norm.occ = log2(temp$Occupancy[i]/ temp$Occupancy[1])
        })
        
        temp$norm.occ=norm.occ
        
        return(temp) 
      }, .progress = "text") 
    }  
    # Save file 
    t<- Sys.time()
    t=gsub(" ", "_", t)
    t=gsub(":", "_", t)
    write.csv(approx.mod, paste("Approx_Mod_", t))
  }
  
  if(abund=="Exact"){
    exact.abund=data[which(data$Abundance=="Exact"),]
    exact.abund$Log2fc.Abundance=exact.abund$Intensity

    if(length(ref.state) != 0){
      ref = ref.state
      levels = c(ref, as.character(unique(exact.abund$Condition)[which(unique(exact.abund$Condition)!=ref)]))
      
      # Check for existence and validity of reference state (i.e. only one such state)
      if(any(as.character(exact.abund$Condition)==ref)==F | length(ref)>1){
        stop(message(sprintf("Invalid or non-existant reference state %s. \n", ref)))
      }
      
      #EntryID is a unique identifier of each sample for each peptide
      exact.abund$EntryID = paste(exact.abund$Protein,exact.abund$Peptide,exact.abund$Modification,exact.abund$PatientID,exact.abund$Run, sep = "_")
      
      exact.abund$Condition = factor(exact.abund$Condition, levels=levels)
      
      exact.abund = exact.abund[order(exact.abund$Condition),]
      
      exact.abund = ldply(unique(exact.abund$EntryID), function(i){
        temp<- subset(exact.abund, EntryID==i)
        norm.abund = sapply(1:nrow(temp), function(i){
        norm.abund = log2(temp$Log2fc.Abundance[i]/ temp$Log2fc.Abundance[1])
        })
        
        temp$Log2fc.Abundance=norm.abund
        
        return(as.data.frame(temp)) 
      }, .progress = "text") 
    }
    # Save file 
    t<- Sys.time()
    t=gsub(" ", "_", t)
    t=gsub(":", "_", t)
    write.csv(exact.abund, paste("Exact_Abund_", t))
  }
  
  if(abund=="Approx"){
    data$EntryID = paste(data$Protein,data$Condition,data$PatientID,data$Run, sep = "_")
    temp=data[which(data$Abundance=="Approx" | (data$Stoichiometry=="Exact")),]
    approx.abund = ldply(unique(temp$EntryID), function(i){
      temp2 = subset(temp, EntryID==i)
      if(any(temp2$Modification=="U")){
        temp2$Log2fc.Abundance = sum(temp2$Intensity)
        return(as.data.frame(temp2))
      }else{
        temp2$Log2fc.Abundance = temp2$Intensity
        return(as.data.frame(temp2))
      }
    }, .progress="text")

    
    if(length(ref.state) != 0){
      ref = ref.state
      levels = c(ref, as.character(unique(approx.abund$Condition)[which(unique(approx.abund$Condition)!=ref)]))
      
      # Check for existence and validity of reference state (i.e. only one such state)
      if(any(as.character(approx.abund$Condition)==ref)==F | length(ref)>1){
        stop(message(sprintf("Invalid or non-existant reference state %s. \n", ref)))
      }
      
      #EntryID is a unique identifier of each sample for each peptide
      approx.abund$EntryID = paste(approx.abund$Protein,approx.abund$Peptide,approx.abund$Modification,approx.abund$PatientID,approx.abund$Run, sep = "_")
      
      approx.abund$Condition = factor(approx.abund$Condition, levels=levels)
      
      approx.abund = approx.abund[order(approx.abund$Condition),]
      
      approx.abund = ldply(unique(approx.abund$EntryID), function(i){
        temp<- subset(approx.abund, EntryID==i)
        norm.abund = sapply(1:nrow(temp), function(i){
          norm.abund = log2(temp$Log2fc.Abundance[i]/ temp$Log2fc.Abundance[1])
        })
        
        temp$Log2fc.Abundance=norm.abund
        
        return(as.data.frame(temp)) 
      }, .progress = "text") 
    } 
    # Save file 
    t<- Sys.time()
    t=gsub(" ", "_", t)
    t=gsub(":", "_", t)
    write.csv(approx.abund, paste("Approx_Abund_", t))
  }
  
}
  
  
    
    
    
    
    
    
    
    
  