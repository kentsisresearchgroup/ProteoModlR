# This function performs basic QC befor analysis

QC <- function(data){
  data$X <- NULL # reading in df gives an extra column. How to avoid?
  # Check for the correct number of columns. The input must be in the required order
  if(ncol(data)!=9){
    stop(message("Input does not contain the required number of columns. \n"))
  }
  # Set column names to what is required for the remainder of the analysis
  colnames(data) <- c("Protein", "Peptide", "Modification", "Position", "Condition",
                      "PatientID", "Label", "Run", "Intensity")
  # Report intensity values that are zero or missing (NA) and remove from further analysis
  if(sum(data$Intensity<=0)!=0) {
    cat(sprintf("%s-labeled peptide %s for patient %s condition %s run %s has 0 intensity: Peptide is removed from further analysis \n",
            data$Label[data$Intensity<=0], data$Peptide[data$Intensity<=0],data$PatientID[data$Intensity<=0],
            data$Condition[data$Intensity<=0],data$Run[data$Intensity<=0]))
    data <- data[-which(data$Intensity<=0),]
  }
  
  if(sum(is.na(data$Intensity))!=0) {
    cat(sprintf("%s-labeled peptide %s for patient %s condition %s run %s has NA as intensity: Peptide is removed from further analysis \n",
            data$Label[is.na(data$Intensity)], data$Peptide[is.na(data$Intensity)],data$PatientID[is.na(data$Intensity)],
            data$Condition[is.na(data$Intensity)],data$Run[is.na(data$Intensity)])) 
    data <- data[-which(is.na(data$Intensity)),]  
  }

  
  
  #Turn contents into factors 
  data$Protein <- factor(data$Protein)
  data$Peptide <- factor(data$Peptide)
  data$Modification <- factor(data$Modification)
  data$Condition <- factor(data$Condition)
  data$PatientID <- factor(data$PatientID)
  data$Label <- factor(data$Label)
  
  return(data)
}