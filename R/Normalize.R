#' Normalization
#'
#' This function offers several modes of normalization to correct for different sources of error, and enables calculating differential protein expression and stoichiometry of post-tranlsational modification even from sparsely annotated datasets.
#' @param data data frame with 12 columns as decribed in documentation.
#' @param iso.norm takes a character string corresponding to the reference isotopologue to which the other isotopologue will be normalized (e.g. iso.norm = "H"). The string must match the corresponding annotation in the "Label" column. 
#' @param equim.iso.norm takes a character string corresponding to the reference isotopologue to which the other isotopologue will be normalized (e.g. iso.norm = "H"). The string must match the corresponding annotation in the "Label" column.
#' @param internal.norm takes in a vector of character strings corresponding to the peptide sequence of the reference peptides to which other peptides are normalized (e.g. internal.norm = "ATDVIVP"). If more than one reference peptide is indicated (e.g. internal.norm = c("ATDVIVP","AAATDVI")), a geometric mean of the corresponding intensity values is taken. This normalization equalizes the intensities for the specified peptides in each sample.
#' @param tot.current takes in a boolean input indicating whether or not normalization to total ion current is to be made (e.g. tot.current = T). This normalization equalizes the sum of all intensities in each sample, used as a proxy for total ion current.
#' @param mod If any chemoform of a given peptide bears the modification specified in mod="" (f.ex. mod="phosphorylation") all peptides with the same sequence and containing the modification will be classified as Modified, while all peptides with the same sequence and not containing the modification will be classified as Not-Modified. If none of the chemoforms of a peptide contains the specified modification, the peptide is classified Abundance.
#' @export
#' @return data frame of normalized data
#' @examples 
#' Normalize(testData,tot.current = T,mod="Phosphorylation")

###################################################################################
################################ NORMALIZATION ####################################
# Normalyze takes in the data along with 3 options. If no options are chosen, data
# is returned unchanged. If the normalyzation step is carried out, the normalized
# intensity values replace the raw intensity values in the Intensity column.
# iso.norm takes in a character string corresponding to the reference isotopologue
# to which the other isotopologue will be normalized. internal.norm takes in a 
# vector of character strings corresponding to the peptide sequence of the 
# reference peptides to which other peptides are normalized. If more than one
# reference peptide is indicated, an arithmetic mean is taken. tot.current takes in 
# a boolean input indicating whether or not normalization to total current should
# be made. 

Normalize <- function(data, iso.norm=NA , equim.iso.norm=NA, internal.norm=NA , tot.current=F, mod){

  data <- QC(data, mod) # Send data to QC function 
  
  ########################## NORMALYZING BY AN ISOTOPOLOGUE ######################
  if(sum(is.na(iso.norm)) == 0){ 
    # Check to see if there are 2 labels 
    if(length(unique(data$Label))==1){
      stop(message("Only one isotopologue found. Data must contain two isotopologues 
                   in the Label column for this type of normalization. \n"))
    }
    exp = as.character(unique(data$Label)[which(unique(data$Label)!=iso.norm)]) #non-ref isotope
    ref = iso.norm #reference isotope

    
    if(any(data$Label==ref)==F){
      stop(message("Reference isotope not found in the Label column of data. 
                   Match spelling and case to format used in data. \n"))
    }
    
    if(any(data$Label==exp)==F){
      stop(message("Non-reference isotope not found in the Label column of data. 
                   Match spelling and case to format used in data. \n"))
    }
    
    
    #Normalize to reference
    data$Label = factor(data$Label, levels=c(exp, ref))
    data = data[order(data$Label),]
    #EntryID is a unique identifier of each isotopologue for each peptide
    data$EntryID = paste(data$Protein,data$Peptide,data$Condition,sep = "_")

    medians=ldply(unique(data$EntryID), function(i) {
      ind=grep(i,data$EntryID)
      median.int = median(data[ind,]$Intensity)
      temp = data[ind,]
      temp$Intensity=median.int*temp$Intensity[1:(0.5*nrow(temp))]/temp$Intensity[(0.5*nrow(temp)+1):nrow(temp)]
      return(temp)
    }, .progress="text")
        

     
    # Cut data in half (leaving out one isotopologue) and store normalized intensities
    medians = medians[order(medians$Label),]
    data=medians[(1:(0.5*(dim(medians)[1]))),]

  }
  
  ########################## NORMALYZING BY AN EQUIMOLAR ISOTOPOLOGUE ######################
  if(sum(is.na(equim.iso.norm)) == 0){ 
    # Check to see if there are 2 labels 
    if(length(unique(data$Label))==1){
      stop(message("Only one isotopologue found. Data must contain two isotopologues 
                   in the Label column for this type of normalization. \n"))
    }
    exp = as.character(unique(data$Label)[which(unique(data$Label)!=equim.iso.norm)]) #non-ref isotope
    ref = equim.iso.norm #reference isotope
    
    
    if(any(data$Label==ref)==F){
      stop(message("Reference isotope not found in the Label column of data. 
                   Match spelling and case to format used in data. \n"))
    }
    
    if(any(data$Label==exp)==F){
      stop(message("Non-reference isotope not found in the Label column of data. 
                   Match spelling and case to format used in data. \n"))
    }
    
    
    #Normalize to reference
    #Median Intensity of reference
    data$Label = as.character(data$Label)
    median.int = median(data[data$Label==ref,]$Intensity)
    data$Label = factor(data$Label, levels=c(exp, ref))
    data = data[order(data$Label),]
    #EntryID is a unique identifier of each isotopologue for each peptide
    data$EntryID = paste(data$Protein,data$Peptide,data$Modification,data$Condition,sep = "_")

    ratios=sapply(unique(data$EntryID), function(i) {
      ind=grep(i,data$EntryID)
      if(length(ind)!=2){
        stop(message(sprintf("EntryID %s missing a labeled or unlabeled state", ind)))
      }
      ratio=median.int*data$Intensity[ind[1]]/data$Intensity[ind[2]]
      })    
    
    # Cut data in half (leaving out one isotopologue) and store normalized intensities
    data=data[(1:(0.5*(dim(data)[1]))),]
    data$Intensity = ratios
  }
  
  
  #################### NORMALYZING BY AN INTERNAL REFERENCE ######################
  
  if(sum(is.na(internal.norm)) == 0){
    
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
      # if(length(which(data[grep(i, data$Peptide),]$Intensity==0))==0)
      ref.df=data[which(data$Peptide==i),]
      ref.df = ref.df[order(ref.df$Condition),]
      if(length(unique(ref.df$Condition))<length(unique(data$Condition))){
        stop(message("Reference peptides must be measured under all conditions. \n"))
      }
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
    
  #################### NORMALYZING BY TOTAL ION CURRENT ######################  
  
  if(tot.current==T){
    data$Condition = as.character(data$Condition)
    tic=ldply(unique(data$Condition), function(i) {
      temp = data[which(data$Condition==i),]
      tic.condition=sum(temp$Intensity)
      return(tic.condition)
    }, .progress="text")
    median.tic = median(tic$V1)
    data.normalized=ldply(unique(data$Condition), function(i) {
      temp = data[which(data$Condition==i),]
      temp$Intensity*median.tic/tic[i,1]
      return(temp)
    }, .progress="text")
  data <- data.normalized
    
  }
  

  # Save file 
  t<- Sys.time()
  t=gsub(" ", "_", t)
  t=gsub(":", "_", t)
  write.csv(data, paste("Normalized_", t))
  
  return (data)

}
  

