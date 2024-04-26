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

# Script first created the 24th of Jan 2020

# This script is to import the results of a title-and-abstract
# screening on indirect genetic effects conducted using rayyan
# and to generate the database to proceed to the fulltext
# screening phase.


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

# We (four screeners: DF, MM, AST, FS) screened all 1058 references
# and include 199 for fulltext screening.

rayyan.output <- read.table("literature_review/title_and_abstract_screening/title-and-abstract_decisions_rayyan_studyID.csv",
                            header=TRUE,sep=",",quote="\"")


# title-and-abstract decision
rayyan.output$t.and.a_final_decision <- ifelse(str_detect(rayyan.output$notes, "Included"),"yes","no")
rayyan.output$t.and.a_final_decision <- as.factor(rayyan.output$t.and.a_final_decision)

# extracting exclusion reasons, and standardizing word use
t.and.a_exclusion_reason <- rayyan.output$notes %>% str_match("RAYYAN-EXCLUSION-REASONS: (\\X+)")

t.and.a_exclusion_reason <- str_to_lower(t.and.a_exclusion_reason[,2], locale = "en")

t.and.a_exclusion_reason <- str_replace(t.and.a_exclusion_reason, "maternal effects only", "maternal effects")
t.and.a_exclusion_reason <- str_replace(t.and.a_exclusion_reason, "no-social effects", "no social effects")
t.and.a_exclusion_reason <- str_replace(t.and.a_exclusion_reason, "off-topic", "off topic")
t.and.a_exclusion_reason <- str_replace(t.and.a_exclusion_reason, "maternal,", "maternal effects,")
t.and.a_exclusion_reason <- str_replace(t.and.a_exclusion_reason, "phenotypic,", "phenotypic level only,")
t.and.a_exclusion_reason <- str_replace(t.and.a_exclusion_reason, "phenotypic$", "phenotypic level only")

# adding this column to the database
rayyan.output$t.and.a_exclusion_reasons <- t.and.a_exclusion_reason


##############################################################
# Assigning fulltext screening
##############################################################

set.seed(77) #setting seed for reproducibility

# subset all those included for fulltext screening
rayyan.output.include <- rayyan.output[rayyan.output$t.and.a_final_decision=="yes",] # there are five included references with exclusion reasons assigned, this is because those were part of the conflicts, and we forgot to remove the exclusion labels in rayyan. Nothing to worry about.

# subset for testing decision tree
rayyan.output.include.trial <- rayyan.output.include[sample(c(1:nrow(rayyan.output.include)),20),]

# adding label of whether this is trial or remaining to make the process later simpler
rayyan.output.include.trial$type <- "trial"

# remaining subset for 
set.seed(77)
rayyan.output.include.remaining <- rayyan.output.include[-sample(c(1:nrow(rayyan.output.include)),20),]

# adding label of whether this is trial or remaining to make the process later simpler
rayyan.output.include.remaining$type <- "remaining"


# # adding two information variables
# rayyan.output.include.remaining$pdf_downloaded <- ""
# rayyan.output.include.remaining$paper_screened <- ""

##############################################################
# Trial screening assignment so that all 20 refs are 
# double-screened twice, and each observer overlaps with all 
# the other observers
rayyan.output.include.trial.FS <- rayyan.output.include.trial[c(1:5,13:17),]
rayyan.output.include.trial.MM <- rayyan.output.include.trial[c(3:7,11:12,18:20),]
rayyan.output.include.trial.DF <- rayyan.output.include.trial[c(6:15),]
rayyan.output.include.trial.AST <- rayyan.output.include.trial[c(1:2,8:10,16:20),]


##############################################################
# Final screening

# estimating the number of references to be double-screened if we aim at double-screening ca. 25% of all references
number.additional.refs <- round(round(nrow(rayyan.output.include.remaining)*0.25,0)/4,0)

# assigning 1/4 references randomly to each observer


list.studies <- sample(rayyan.output.include.remaining$studyID,length(rayyan.output.include.remaining$studyID))

DF.list <- as.character(list.studies[1:45])
MM.list <- as.character(list.studies[46:89])
AST.list <- as.character(list.studies[90:134])
FS.list <- as.character(list.studies[135:179])

# extracing additional double-screening studies
DF.additional <- sample(c(MM.list,AST.list,FS.list),number.additional.refs)
MM.additional <- sample(c(DF.list,AST.list,FS.list),number.additional.refs)
AST.additional <- sample(c(DF.list,MM.list,FS.list),number.additional.refs)
FS.additional <- sample(c(DF.list,MM.list,AST.list),number.additional.refs)

# adding them all
DF.fulltext.screening.list <- c(DF.list,DF.additional)
MM.fulltext.screening.list <- c(MM.list,MM.additional)
AST.fulltext.screening.list <- c(AST.list,AST.additional)
FS.fulltext.screening.list <- c(FS.list,FS.additional)


##############################################################
# Creating output
##############################################################

DF.full <- rbind(rayyan.output.include.trial.DF,
                 rayyan.output.include.remaining[rayyan.output.include.remaining$studyID %in% DF.fulltext.screening.list,])

MM.full <- rbind(rayyan.output.include.trial.MM,
                 rayyan.output.include.remaining[rayyan.output.include.remaining$studyID %in% MM.fulltext.screening.list,])

AST.full <- rbind(rayyan.output.include.trial.AST,
                 rayyan.output.include.remaining[rayyan.output.include.remaining$studyID %in% AST.fulltext.screening.list,])

FS.full <- rbind(rayyan.output.include.trial.FS,
                 rayyan.output.include.remaining[rayyan.output.include.remaining$studyID %in% FS.fulltext.screening.list,])


# final templates
write.xlsx(DF.full, "literature_review/fulltext_screening/fulltext_screening_template_DF.xlsx",
           sheetName="Sheet1",col.names=TRUE, row.names=F,
           append=FALSE, showNA=TRUE, password=NULL)
write.xlsx(MM.full, "literature_review/fulltext_screening/fulltext_screening_template_MM.xlsx",
           sheetName="Sheet1",col.names=TRUE, row.names=F,
           append=FALSE, showNA=TRUE, password=NULL)
write.xlsx(AST.full, "literature_review/fulltext_screening/fulltext_screening_template_AST.xlsx",
           sheetName="Sheet1",col.names=TRUE, row.names=F,
           append=FALSE, showNA=TRUE, password=NULL)
write.xlsx(FS.full, "literature_review/fulltext_screening/fulltext_screening_template_FS.xlsx",
           sheetName="Sheet1",col.names=TRUE, row.names=F,
           append=FALSE, showNA=TRUE, password=NULL)


# saving versions used for reproducibility purposes
sink("literature_review/fulltext_screening/fulltext_templates_Rpackages_session.txt")
sessionInfo()
sink()
