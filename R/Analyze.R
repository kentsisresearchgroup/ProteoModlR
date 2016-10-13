#' Analyze
#'
#' This function calculates peptides site occupancy and/or abundance.
#' @param data data frame with 12 columns as decribed in documentation. Output of Normalize function
#' @param mod If any chemoform of a given peptide bears the modification specified in mod="" (f.ex. mod="phosphorylation") all peptides with the same sequence and containing the modification will be classified as Modified, while all peptides with the same sequence and not containing the modification will be classified as Not-Modified. If none of the chemoforms of a peptide contains the specified modification, the peptide is classified Abundance.
#' @param stoich user's preference for calculating the site occupancy of peptides for which either calculations are possible. Possible arguments: "Exact" or "Approximate".
#' @param adund user's preference for exact or approximate abundance calculations. Possible arguments: "Exact" or "Approximate".
#' @param ref.state option to normalize to a reference condition (e.g. ref.state="Disease").
#' @export
#' @return data frame of analyzed data
#' @examples 
#' normalized <- Normalize(testData,tot.current = T,mod="Phosphorylation")
#' Analyze(data=normalized,mod="Phosphorylation",stoich="Exact")

########################## ABUNDANCE AND STOICHIOMETRY ##################################
# Analyze.R calculates the approximate and/or exact site occupancy and/or abundance 
# of the peptides. The choice between exact and approximate calculations will be
# based on user input, while the selection of appropriate peptides is based on labeling
# done in QC.R. If a reference state is specified, stoichiometric ratios and abundances
# for all conditions are normalized to the reference state. 

Analyze <- function(data, mod, stoich="", abund="", ref.state=NA){
  data$X <- NULL
  data$Abundance.Calc = NA
  data$Occupancy = NA 
  if(stoich=="Exact"){
    data$Modification <- as.character(data$Modification)
    temp= subset(data, Stoichiometry=="Exact")
    temp$EntryID <- paste(temp$Protein,temp$Peptide,temp$Condition, sep = "_")
    exact.mod = ldply(unique(temp$EntryID), function(i){
      temp2 = subset(temp, EntryID==i)
      stoich = sum(temp2[which(grepl(mod, temp2$Modification)),]$Intensity)/sum(temp2$Intensity)
      temp2$Occupancy <- stoich
      return(as.data.frame(temp2))
    }, .progress='text')
    rm(temp)
    
    #Normalize to reference state
    if(!is.na(ref.state)){
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
    temp= subset(data, Stoichiometry=="Approx")
    approx.mod = ldply(unique(temp$EntryID),function(i){
      temp2<- subset(temp,EntryID==i)
      if(any(temp2$Abundance=="Exact")){
        stoich = sum(temp2[which(grepl(mod, temp2$Modification)),]$Intensity)/
          (exp(sum(log(temp2[which(temp2$Modification=="Unmodified"),]$Intensity))/
             nrow(temp2[which(temp2$Modification=="Unmodified"),])))
        temp2$Occupancy = stoich
        return(as.data.frame(temp2))
      }else{
        temp2$Occupancy= 1
        return(as.data.frame(temp2))
      }
      
    }, .progress="text")
    rm(temp)
    
    if(!is.na(ref.state)){
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
    exact.abund=data[which(data$Stoichiometry=="Zero"),]
    exact.abund$Abundance.Calc=exact.abund$Intensity

    if(!is.na(ref.state)){
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
        norm.abund = log2(temp$Abundance.Calc[i]/ temp$Abundance.Calc[1])
        })
        
        temp$Abundance.Calc=norm.abund
        
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
    data$EntryID = paste(data$Protein,data$Peptide, data$Condition,data$PatientID,data$Run, sep = "_")
    temp=data[which(data$Abundance=="Approx"),]
    approx.abund = ldply(unique(temp$EntryID), function(i){
      temp2 = subset(temp, EntryID==i)
      temp2$Abundance.Calc= (exp(sum(log(temp2$Intensity))/
             nrow(temp2)))
        return(as.data.frame(temp2))
    }, .progress="text")

    
    if(!is.na(ref.state)){
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
          norm.abund = log2(temp$Abundance.Calc[i]/ temp$Abundance.Calc[1])
        })
        
        temp$Abundance.Calc=norm.abund
        
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
  
  
    
    
    
    
    
    
    
    
  