## Helper
library(dplyr)
Actionability<- read.csv("./Actionability Assessment.csv")
Functional<- read.csv("./Functional Category.csv")
Legend.label<- read.csv("./Legend.csv")

load("molecular_data.RData")

source("functions.R")

my.molecular <- molecular
##get rid of missing test results
my.molecular <- my.molecular[!is.na(my.molecular$Positive),] 


## Filter the Gene and Sample to be greater than 10
gene<-(unique(my.molecular$Biomarker))
important.gene<-c()
for (i in (1:length(gene))){

  if (sum(my.molecular$Biomarker==gene[i])>=10){
    important.gene<-c(important.gene,gene[i])
  }

}


sample<-unique(my.molecular$PatientID)
important.sample<-c()
for (i in (1:length(sample))){

  if (sum(my.molecular$PatientID==sample[i])>=10){
    important.sample<-c(important.sample,sample[i])
  }

}

my.molecular=my.molecular[which(my.molecular$PatientID%in%important.sample & my.molecular$Biomarker%in%important.gene),]


### Label
## Special: IHC, Microarray-Illumina 
## CISH: Amplified; Not Amplified
## CNV: Amplifiedl Amplification Not Detected/QNS (Quantity Not Sufficient); Amplified
## FA: MSI-High; Negative/Stable
## FISH: Amplified/Positive-Not Amplified (Positive);Not Amplified/Negative
## FISH-Mutational: Positive; Negative
## Fusion (0): Fusion Not Detected;
## IHC: Negative/Above Threshold/Below Threshold/Technical Issues(Test not performed);Positive
## IHC H-Score: Negative;Positive
## Microarray-Illumina (0): Over Expressed/No Change/Under Expressed
## Molecular: Wild-Type genotype/Mutated (c.35G>T,p.G12V); Mutated
## Next Gen Seq: Wild Type; Pathogenic (Mutated/G13D/Exon 20...)
## NGS Q2: Mutated(Pathogenic);Wild Type (Mutation Not Detected)/ No Result
## NGS Q3: Mutated (Pathogenic)/G12/G13/Exon11/Persumed Pathogenic; Wild Type/Variant not detected
## RT-PCR: High; Low
## Sanger SEQ: Absent/Wild Type;


## Count the number of Disease Categories for each patient ID
# my.molecular$Disease.Category<-as.character(my.molecular$Disease.Category)
# my.molecular$Disease.Category[is.na(my.molecular$Disease.Category)]<-FALSE

##PatientID represents the patient ID (one row per gene or protein tested, usually have multiple rows per patient)
##Test represents the type of test
##Biomarker is usually gene or protein name

##display just a subset of genes (otherwse takes forever to make and load plots)
# genes.of.interest=c("cMET", "EGFR", "ER", "ERCC1", "Her2/Neu", "MGMT", "MLH1", "MSH2", "MSH6", "PD-1", "PD-1 IHC", "PD-L1", "PD-L1 IHC", "PGP", "PMS2", "PR", "PTEN", "RPM1", "SPACR Monoclonal", "SPARC Polyclonal", "TLE3", "TOP2A", "TOPO1", "TS", "TUBB3")
# ##only show IHC markers
# my.molecular=my.molecular[which(my.molecular$Technology=="IHC" & my.molecular$Biomarker%in%genes.of.interest),]
# 
# ##also take a subset of samples
# set.seed(304819)
# sample.subset <- sample(my.molecular$PatientID, 200)
# 
# my.molecular <- filter(my.molecular, PatientID %in% sample.subset)
