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
#To obtain testing set 
training_ids<-unique(training_set$sample_id)
testing_set <- subset(psi.tissue.filtered, !(psi.tissue.filtered$sample_id %in% training_ids))

