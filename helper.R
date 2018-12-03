## Helper
library(dplyr)
Actionability<- read.csv("./Actionability Assessment.csv")
Functional<- read.csv("./Functional Category.csv")
Legend.label<- read.csv("./Legend.csv")

load("../molecular_data.RData")

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


