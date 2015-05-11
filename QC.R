# data <- read.table("tp_data.csv", sep=',', header = T)

# This function performs basic QC for the analysis

QC <- function(data){
  data$X <- NULL # reading in df gives an extra column. How to avoid?
  if(ncol(data)!=9){
    stop(message("Input does not contain the required number of columns. \n"))
  }
  colnames(data) <- c("Protein", "Peptide", "Modification", "Position", "Condition",
                      "PatientID", "Label", "Run", "Intensity")
  if(sum(data$Intensity<=3000)!=0) {
    sprintf("Peptides %s for patients %s runs %s has 0 intensity: Peptide is removed from further analysis",
            data$Peptide[data$Intensity<=3000],data$PatientID[data$Intensity<=3000],data$Run[data$Intensity<=3000])
  }
  data <- data[-which(data$Intensity<=3000),]
  if(sum(is.na(data$Intensity))!=0) {
    sprintf("Peptides %s for patients %s runs %s has 0 intensity: Peptide is removed from further analysis",
            data$Peptide[data$Intensity<=3000],data$PatientID[data$Intensity<=3000],data$Run[data$Intensity<=3000])
    data <- data[-which(is.na(data$Intensity)),]
  }
  
  data$Protein <- factor(data$Protein)
  data$Peptide <- factor(data$Peptide)
  data$Modi <- factor(data$Peptide)
  
  
  return(data)
}