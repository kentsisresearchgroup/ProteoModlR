# STEP 1: Generate Skyline output

########################Generating Case-Control Data#############################

#Import File
#e.g. For LC-MS case-control, file must contain 10 columns with ProteinName,
#PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge, 
#IsotopeLabelType, Condition, BioReplicate, Run, Intensity 

RawData<-read.table("Input.csv", sep=',', header=T, na.strings="")
TestData<-RawData[3:9,]

# Remove rows with type D (specific to Paolo's dataset)
TestData = TestData[-which(TestData$Modification=="D"),]
#TestData = TestData[-which(TestData$Modification=="NC"),]


#Duplicate data for 2 conditions, keeping it as dataframe
TestData<-as.data.frame(lapply(TestData,function(x)rep(x,2)))

#Make first half chemo resistant (CR) and second half non-resistant (NR)
TestData$Condition[1:(nrow(TestData)/2)]="CR"
TestData$Condition[(nrow(TestData)/2+1):nrow(TestData)]="NR"


#Duplicate data for 2 isotopes and order based on condition
TestData<-as.data.frame(lapply(TestData,function(x)rep(x,2)))
TestData=TestData[order(TestData$Condition),]

#Make first half of each condition L and second half H 
TestData$Label[1:(nrow(TestData)/4)]="L"
TestData$Label[(nrow(TestData)/4+1):(nrow(TestData)/2)]="H"
TestData$Label[(nrow(TestData)/2+1):(0.75*nrow(TestData))]="L"
TestData$Label[(0.75*nrow(TestData)+1):nrow(TestData)]="H"

TestData$Condition <- factor(TestData$Condition, levels=unique(TestData$Condition))
TestData$Label <- factor(TestData$Label, levels=unique(TestData$Label))

# Generate intensity randomly between 10^3 to 10^10 with higher values
# in CR
TestData$AUC[1:(nrow(TestData)/2)] = runif((nrow(TestData)/2), min=1000, max(10000))
TestData$AUC[(nrow(TestData)/2+1):nrow(TestData)] = runif((nrow(TestData)/2), min=10, max(100))

#Creating replicates
Testdata.2 <- TestData
Testdata.2$PatientID[1:(nrow(Testdata.2)/2)] = 3
Testdata.2$PatientID[(nrow(Testdata.2)/2+1):(nrow(Testdata.2))] = 4
Testdata.2$PatientID <- factor(Testdata.2$PatientID, levels=unique(Testdata.2$PatientID))
Testdata.2$AUC = Testdata.2$AUC + rnorm(dim(Testdata.2)[1], mean=100, sd=20)

Testdata.3 <- TestData
Testdata.3$PatientID[1:(nrow(Testdata.3)/2)] = 5
Testdata.3$PatientID[(nrow(Testdata.3)/2+1):(nrow(Testdata.3))] = 6
Testdata.3$PatientID <- factor(Testdata.3$PatientID, levels=unique(Testdata.3$PatientID))
Testdata.3$AUC = Testdata.3$AUC + rnorm(dim(Testdata.3)[1], mean=100, sd=20)

TestData$PatientID[1:(nrow(TestData)/2)] = 1
TestData$PatientID[(nrow(TestData)/2+1):(nrow(TestData))] = 2
TestData$PatientID <- factor(TestData$PatientID, levels=unique(TestData$PatientID))
TestData=rbind(TestData, Testdata.2, Testdata.3)
TestData$Run = TestData$PatientID
rm(Testdata.2)
rm(Testdata.3)
rm(RawData)
rownames(TestData)=NULL

write.csv(TestData, "stat_data.csv")

########################Generating Timepoint Data#############################

#Import File
#e.g. For LC-MS case-control, file must contain 10 columns with ProteinName,
#PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge, 
#IsotopeLabelType, Condition, BioReplicate, Run, Intensity 

RawData<-read.table("Input.csv", sep=',', header=T, na.strings="")
TestData<-RawData[3:9,]

# Remove rows with type D (specific to Paolo's dataset)
TestData = TestData[-which(TestData$Modification=="D"),]
#TestData = TestData[-which(TestData$Modification=="NC"),]


#Duplicate data for 10 time points, keeping it as dataframe, with tp indicated by
#numbers
TestData<-as.data.frame(lapply(TestData,function(x)rep(x,10)))

for (i in 1:5){
  TestData$Condition[(10*(i-1)+1):(10*i)] = i  
}


#Make first half of each timepoint L and second half H 
for (i in seq(1,50,10)){
  TestData$Label[i:(4+i)] = "H"
  TestData$Label[(5+i):(9+i)] = "L"
}

TestData$Condition <- factor(TestData$Condition, levels=unique(TestData$Condition))
TestData$Label <- factor(TestData$Label, levels=unique(TestData$Label))

# Generate intensity randomly between 10^3 to 10^10 with higher values
# later time points
for (i in as.numeric(unique(TestData$Condition))){
  for (j in which(TestData$Condition==i)){
  TestData$AUC[j] = runif(1, min=1000, max(10000))
  + rnorm(1, (1000*i), 10)
  }
}

#Creating replicates
Testdata.2 <- TestData
Testdata.2$PatientID[1:nrow(Testdata.2)] = 2
Testdata.2$PatientID <- factor(Testdata.2$PatientID, levels=unique(Testdata.2$PatientID))
Testdata.2$AUC = Testdata.2$AUC + rnorm(dim(Testdata.2)[1], mean=1000, sd=100)

Testdata.3 <- TestData
Testdata.3$PatientID[1:nrow(Testdata.3)] = 3
Testdata.3$PatientID <- factor(Testdata.3$PatientID, levels=unique(Testdata.3$PatientID))
Testdata.3$AUC = Testdata.3$AUC + rnorm(dim(Testdata.3)[1], mean=1000, sd=100)

TestData$PatientID[1:nrow(TestData)] = 1
TestData$PatientID <- factor(TestData$PatientID, levels=unique(TestData$PatientID))
TestData=rbind(TestData, Testdata.2, Testdata.3)
TestData$Run = TestData$PatientID
rm(Testdata.2)
rm(Testdata.3)
rm(RawData)
rownames(TestData)=NULL


write.csv(TestData, "tp_data.csv")


