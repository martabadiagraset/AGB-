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
