setwd("/Users/Aida/Documents/Bioinformatics_Master/2nd_TERM/AGB/project")
#USEFUL: To read the large files that we will be producing along all the process: 
# library(data.table)
# psi.tissue <- fread("psi.tissue.tsv")

#Merge table and data preprocessing
PSI <- read.table("gtex_heart_kidney_liver_lung_muscle_nerve_samples_formatted.psi.txt")
library(reshape2)
psi_long <- melt(as.matrix(PSI), value.name = "PSI")
colnames(psi_long)[1] <- "splicing_event"
colnames(psi_long)[2] <- "sample_id"
description <- read.table("gtex_Heart_Kidney_Liver_Lung_Muscle_Nerve_phenotype.txt", 
                          col.names = c('sample_id', 'subtissue', 'tissue', 'sample_type', 'gender', 'source'))
description$sample_id <- gsub("-", ".", description$sample_id, fixed = T)
sample_id.tissue <- description[,c("sample_id", "tissue", "gender")]
psi.tissue <- merge(psi_long, sample_id.tissue, by="sample_id")
psi.tissue <- psi.tissue[psi.tissue$tissue != "Kidney", ]
#Filtering NAs
notNA_tissue <- dcast(psi.tissue, splicing_event~tissue, value.var = "PSI", fun.aggregate = function(x) {length(x[is.na(x)])/length(x)})
#Other way to do it: dcast(psi.tissue, splicing_event~tissue, value.var = "PSI", fun.aggregate = function(x) {table(is.na(x))[1]/length(x)})
#To remove splicing events that contain more than 10% of NAs across at least one tissue.
notNA_tissue <- notNA_tissue[notNA_tissue$Muscle <= 0.1, ]
notNA_tissue <- notNA_tissue[notNA_tissue$Nerve <= 0.1, ]
notNA_tissue <- notNA_tissue[notNA_tissue$Heart <= 0.1, ]
notNA_tissue <- notNA_tissue[notNA_tissue$Liver <= 0.1, ]
notNA_tissue <- notNA_tissue[notNA_tissue$Lung <= 0.1, ]
SE.filtered<-notNA_tissue$splicing_event
psi.tissue.filtered<-psi.tissue[which(psi.tissue$splicing_event%in%SE.filtered),]

#In psi.tissue.filtered we have all the events that fulfill our restrictions (less than 10%NAs across all the tissue)

#To obtain our training_set
Heart_male <- psi.tissue.filtered[psi.tissue.filtered$tissue == "Heart" & psi.tissue.filtered$gender == "Male" , ]
random_Heart_male<-sample(unique(Heart_male$sample_id),50)
Heart_male_randomtable<-Heart_male[which(Heart_male$sample_id%in%random_Heart_male),]
unique(Heart_male_randomtable$sample_id) #prove 50 samples per Heart tissue

Heart_female <- psi.tissue.filtered[psi.tissue.filtered$tissue == "Heart" & psi.tissue.filtered$gender == "Female" , ]
random_Heart_female<-sample(unique(Heart_female$sample_id),50)
Heart_female_randomtable<-Heart_female[which(Heart_female$sample_id%in%random_Heart_female),]
unique(Heart_female_randomtable$sample_id) #prove 50 samples per Heart tissue

Liver_male <- psi.tissue.filtered[psi.tissue.filtered$tissue == "Liver" & psi.tissue.filtered$gender == "Male" , ]
random_Liver_male<-sample(unique(Liver_male$sample_id),50)
Liver_male_randomtable<-Liver_male[which(Liver_male$sample_id%in%random_Liver_male),]
unique(Liver_male_randomtable$sample_id) #prove 50 samples per Liver tissue

Liver_female <- psi.tissue.filtered[psi.tissue.filtered$tissue == "Liver" & psi.tissue.filtered$gender == "Female" , ]
random_Liver_female<-sample(unique(Liver_female$sample_id),33)
Liver_female_randomtable<-Liver_female[which(Liver_female$sample_id%in%random_Liver_female),]
unique(Liver_female_randomtable$sample_id) #prove 50 samples per Liver tissue

Lung_male <- psi.tissue.filtered[psi.tissue.filtered$tissue == "Lung" & psi.tissue.filtered$gender == "Male" , ]
random_Lung_male<-sample(unique(Lung_male$sample_id),67)
Lung_male_randomtable<-Lung_male[which(Lung_male$sample_id%in%random_Lung_male),]
unique(Lung_male_randomtable$sample_id) #prove 50 samples per Lung tissue

Lung_female <- psi.tissue.filtered[psi.tissue.filtered$tissue == "Lung" & psi.tissue.filtered$gender == "Female" , ]
random_Lung_female<-sample(unique(Lung_female$sample_id),50)
Lung_female_randomtable<-Lung_female[which(Lung_female$sample_id%in%random_Lung_female),]
unique(Lung_female_randomtable$sample_id) #prove 50 samples per Lung tissue

Nerve_male <- psi.tissue.filtered[psi.tissue.filtered$tissue == "Nerve" & psi.tissue.filtered$gender == "Male" , ]
random_Nerve_male<-sample(unique(Nerve_male$sample_id),50)
Nerve_male_randomtable<-Nerve_male[which(Nerve_male$sample_id%in%random_Nerve_male),]
unique(Nerve_male_randomtable$sample_id) #prove 50 samples per Nerve tissue

Nerve_female <- psi.tissue.filtered[psi.tissue.filtered$tissue == "Nerve" & psi.tissue.filtered$gender == "Female" , ]
random_Nerve_female<-sample(unique(Nerve_female$sample_id),50)
Nerve_female_randomtable<-Nerve_female[which(Nerve_female$sample_id%in%random_Nerve_female),]
unique(Nerve_female_randomtable$sample_id) #prove 50 samples per Nerve tissue

Muscle_male <- psi.tissue.filtered[psi.tissue.filtered$tissue == "Muscle" & psi.tissue.filtered$gender == "Male" , ]
random_Muscle_male<-sample(unique(Muscle_male$sample_id),50)
Muscle_male_randomtable<-Muscle_male[which(Muscle_male$sample_id%in%random_Muscle_male),]
unique(Muscle_male_randomtable$sample_id) #prove 50 samples per Muscle tissue

Muscle_female <- psi.tissue.filtered[psi.tissue.filtered$tissue == "Muscle" & psi.tissue.filtered$gender == "Female" , ]
random_Muscle_female<-sample(unique(Muscle_female$sample_id),50)
Muscle_female_randomtable<-Muscle_female[which(Muscle_female$sample_id%in%random_Muscle_female),]
unique(Muscle_female_randomtable$sample_id) #prove 50 samples per Muscle tissue

training_set <- rbind(Muscle_female_randomtable, Muscle_male_randomtable, Nerve_female_randomtable, Nerve_male_randomtable, Heart_female_randomtable, Heart_male_randomtable, Liver_female_randomtable, Liver_male_randomtable, Lung_female_randomtable, Lung_male_randomtable)
training_set.def <- na.omit(training_set)
training_set.def$PSI_Discrete[training_set.def$PSI < 0.5] <- 0
training_set.def$PSI_Discrete[training_set.def$PSI >= 0.5] <- 1

#To obtain testing set 
training_ids<-unique(training_set$sample_id)
testing_set <- subset(psi.tissue.filtered, !(psi.tissue.filtered$sample_id %in% training_ids))
testing_set.def <- na.omit(testing_set)
testing_set.def$PSI_Discrete[testing_set.def$PSI < 0.5] <- 0
testing_set.def$PSI_Discrete[testing_set.def$PSI >= 0.5] <- 1

# TRAINING SET INPUT (training_set.description)
install.packages("data.table")
library(data.table)
training_set.def <- fread("training_set.def.tsv")
training_samples<-unique(as.factor(training_set.def$sample_id))
training_set.description<-description[which(description$sample_id%in%training_samples),]
write.table(training_set.description, file = "training_set.description", sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)

# TESTING SET INPUT (testing_set.description)
install.packages("miceadds")
library(miceadds)
library(mice)
library(lattice)
library(miceadds)
load.Rdata(filename = "testing_set.def.RData", "testing_set.def") #testing_set.def object
testing_samples<-unique(as.factor(testing_set.def$sample_id))
testing_set.description<-description[which(description$sample_id%in%testing_samples),]
write.table(testing_set.description, file = "testing_set.description", sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)

# PSI INPUT with our filtered training and testing samples (psi.formatted)
training_set.def <- training_set.def[,-1] #delete 1st column (not necessary)
training_testing_set.def<-rbind(training_set.def, testing_set.def)
training_testing_set.def.long<-training_testing_set.def[,c("sample_id","splicing_event","PSI_Discrete")]

training_testing_samples<-append(as.vector(training_samples),as.vector(testing_samples)) #training_samples are already unique

training_testing_set.def.long<-training_testing_set.def.long[which(as.vector(training_testing_set.def.long$sample_id)%in%as.factor(training_testing_samples)),]
save(training_testing_set.def.long, file = "training_testing_set.def.long.RData")

library(tidyr)
training_testing_set.def.wide<-spread(training_testing_set.def.long, key=sample_id, value=PSI_Discrete)
write.table(training_testing_set.def.wide, file = "psi.formatted", sep = "\t", col.names = TRUE, quote = FALSE)

######################################## INPUT files of NaiveBayes.py #################################################### 
############################## names () are Python input files of NaiveBayes function ####################################

# training_set.description (training_file)
# testing_set.description (testing_file)
# psi.formatted (psi_file)

