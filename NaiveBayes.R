setwd("/Users/Aida/Documents/Bioinformatics_Master/2nd_TERM/AGB/project")

#PSI<-read.table("data/PSI/gtex_heart_kidney_liver_lung_muscle_nerve_samples_formatted.psi")
psi_test <- read.table("data/PSI/test.psi")

library(reshape2)
psi_long <- melt(as.matrix(psi), value.name = "PSI")
colnames(psi_long)[1] <- "splicing_event"
colnames(psi_long)[2] <- "sample_id"

description <- read.table("data/description/gtex_Heart_Kidney_Liver_Lung_Muscle_Nerve_phenotype.txt", 
                        col.names = c('sample_id', 'subtissue', 'tissue', 'sample_type', 'gender', 'source'))

description$sample_id <- gsub("-", ".", description$sample_id, fixed = T)
sample_id.tissue <- description[,c("sample_id", "tissue")]

psi.tissue <- merge(psi_long, sample_id.tissue, by="sample_id")
