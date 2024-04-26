##############################################################
# Authors: 
#
# Alfredo Sanchez-Tojar (@ASanchez_Tojar)
# Profile: https://goo.gl/PmpPEB
# Bielefeld University
# Email: alfredo.tojar@gmail.com

##############################################################
# Description of script and Instructions
##############################################################

# Script first created the 23rd of July 2020

# This script is to import the results of a fulltext screening
# on indirect genetic effects conducted using rayyan and to
# generate the databases to proceed to the data extraction
#  phase.


##############################################################
# Packages needed
##############################################################

# load pacakges
pacman::p_load(dplyr,stringr,openxlsx)

# cleaning up
rm(list=ls())


##############################################################
# Functions needed
##############################################################

# none


##############################################################
# Cleaning data from rayyan
##############################################################

# Import the full list of studies for which title and abstract screening was performed
rayyan.output <- read.table("literature_review/title_and_abstract_screening/title-and-abstract_decisions_rayyan_studyID.csv",
                            header=TRUE,sep=",",quote="\"")

# Import the list of studies finally included after fulltext screening
included.after.fulltext <- read.xlsx("literature_review/fulltext_screening/Final_fulltext_screening_responses_including_conflict_resolution_google_form_data.xlsx",
                                     colNames=T,sheet = 3)

# Import the studies included after fulltext screening but only those using animal models
included.after.fulltext.animal.model <- read.table("literature_review/data_extraction/animal_model_papers_to_be_assigned_20200703.csv",
                                                   header=TRUE,sep=",",quote="\"")


##############################################################
# Creating a database with the animal model studies
##############################################################

# reducing to only choose relevant columns
rayyan.output.reduced <- select(rayyan.output,studyID,title,year,journal,volume,issue,pages,authors,url)

# reducing database to only those studies included after fulltext screening and that used animal models
data.extraction.animal.models.full <- rayyan.output.reduced[rayyan.output.reduced$studyID %in% included.after.fulltext.animal.model$StudyID,]

# write.csv(data.extraction.animal.models.full,
#           "literature_review/data_extraction/animal_model_papers_to_be_assigned_full.csv",row.names=FALSE)


# excluding those 6 references that we have already extracted as part of our trial-testing extraction method where all four observers extracted and discussed conflicts
already.extracted <- c("IGE0418","IGE0501","IGE0505","IGE0509","IGE0515","IGE0852")

data.extraction.animal.models.reduced <- data.extraction.animal.models.full[!(data.extraction.animal.models.full$studyID %in% already.extracted),]


##############################################################
# Assigning each reviewer a share of the extraction
##############################################################

# reshuffling references before assigning them to reviewers
set.seed(77)
data.extraction.animal.models.reduced.2 <- data.extraction.animal.models.reduced[sample(c(1:nrow(data.extraction.animal.models.reduced)),
                                                                                        nrow(data.extraction.animal.models.reduced)),]


# final templates
write.xlsx(data.extraction.animal.models.reduced.2[1:14,],
           "literature_review/data_extraction/data_extraction_screening_template_DF.xlsx",
           sheetName="Sheet1",col.names=TRUE, row.names=F,
           append=FALSE, showNA=TRUE, password=NULL)

write.xlsx(data.extraction.animal.models.reduced.2[15:28,],
           "literature_review/data_extraction/data_extraction_screening_template_MM.xlsx",
           sheetName="Sheet1",col.names=TRUE, row.names=F,
           append=FALSE, showNA=TRUE, password=NULL)

write.xlsx(data.extraction.animal.models.reduced.2[29:42,],
           "literature_review/data_extraction/data_extraction_screening_template_AST.xlsx",
           sheetName="Sheet1",col.names=TRUE, row.names=F,
           append=FALSE, showNA=TRUE, password=NULL)

write.xlsx(data.extraction.animal.models.reduced.2[43:57,],
           "literature_review/data_extraction/data_extraction_screening_template_FS.xlsx",
           sheetName="Sheet1",col.names=TRUE, row.names=F,
           append=FALSE, showNA=TRUE, password=NULL)


##############################################################
# Assigning each reviewer a share of the double-checking
##############################################################


# final templates
write.xlsx(data.extraction.animal.models.reduced.2[c(15:18,29:33,43:47),],
           "literature_review/data_extraction/double-checking/data_extraction_double-checking_template_DF.xlsx",
           sheetName="Sheet1",col.names=TRUE, row.names=F,
           append=FALSE, showNA=TRUE, password=NULL)

write.xlsx(data.extraction.animal.models.reduced.2[c(1:5,34:37,48:52),],
           "literature_review/data_extraction/double-checking/data_extraction_double-checking_template_MM.xlsx",
           sheetName="Sheet1",col.names=TRUE, row.names=F,
           append=FALSE, showNA=TRUE, password=NULL)

write.xlsx(data.extraction.animal.models.reduced.2[c(6:10,19:22,53:57),],
           "literature_review/data_extraction/double-checking/data_extraction_double-checking_template_AST.xlsx",
           sheetName="Sheet1",col.names=TRUE, row.names=F,
           append=FALSE, showNA=TRUE, password=NULL)

write.xlsx(data.extraction.animal.models.reduced.2[c(11:14,23:28,38:42),],
           "literature_review/data_extraction/double-checking/data_extraction_double-checking_template_FS.xlsx",
           sheetName="Sheet1",col.names=TRUE, row.names=F,
           append=FALSE, showNA=TRUE, password=NULL)
