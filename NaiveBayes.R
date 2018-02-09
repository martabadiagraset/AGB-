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
sample_id.tissue <- description[,c("sample_id", "tissue")]
psi.tissue <- merge(psi_long, sample_id.tissue, by="sample_id")

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
psi.tissue.filtered <- psi.tissue.filtered[psi.tissue.filtered$tissue != "Kidney", ]
#In psi.tissue.filtered we have all the events that fulfill our restrictions (less than 10%NAs across all the tissue

#To obtain our training_set
Heart <- psi.tissue.filtered[psi.tissue.filtered$tissue == "Heart", ]
random_Heart<-sample(unique(Heart$sample_id),100)
Heart_randomtable<-Heart[which(Heart$sample_id%in%random_Heart),]
unique(Heart_randomtable$sample_id) #prove 100 samples per Heart tissue

Lung <- psi.tissue.filtered[psi.tissue.filtered$tissue == "Lung", ]
random_Lung <- sample(unique(Lung$sample_id),100)
Lung_randomtable <- Lung[which(Lung$sample_id%in%random_Lung),]
unique(Lung_randomtable$sample_id) #prove 100 samples per Lung tissue

Muscle <- psi.tissue.filtered[psi.tissue.filtered$tissue == "Muscle", ]
random_Muscle <- sample(unique(Muscle$sample_id),100)
Muscle_randomtable <- Muscle[which(Muscle$sample_id%in%random_Muscle),]
unique(Muscle_randomtable$sample_id) #prove 100 samples per Muscle tissue

Nerve <- psi.tissue.filtered[psi.tissue.filtered$tissue == "Nerve", ]
random_Nerve <- sample(unique(Nerve$sample_id),100)
Nerve_randomtable <- Nerve[which(Nerve$sample_id%in%random_Nerve),]
unique(Nerve_randomtable$sample_id) #prove 100 samples per Nerve tissue

Liver <- psi.tissue.filtered[psi.tissue.filtered$tissue == "Liver", ]
random_Liver <- sample(unique(Liver$sample_id),100)
Liver_randomtable <- Liver[which(Liver$sample_id%in%random_Liver),]
unique(Liver_randomtable$sample_id) #prove 100 samples per Liver tissue


