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

# Script first created the 16th of Aug 2019

# This script is to import the results of a systematic review
# on indirect genetic effects conducted in Web of Science and
# Scopus.


##############################################################
# Packages needed
##############################################################

# load pacakges
pacman::p_load(dplyr, revtools)

# cleaning up
rm(list=ls())


##############################################################
# Functions needed
##############################################################

# none


##############################################################
# Importing reference data
##############################################################

# The literature searches were conducted on the 8th of August
# 2019. Timespan = all years. Remaining settings = default.

# For more information about the search see:
# https://docs.google.com/document/d/1aJByb9gY79vtSHs_qZ9D1ZVv0zAoaxI-FT7XgGNsjpc/edit?usp=sharing


################################
# WEB OF SCIENCE: general search

# importing the .bib files
wos <- read_bibliography("literature_review/systematic_search/WoS_IGEs_search.bib")


# reducing fields to the minimum number of fields
# so that all databases have the same columns. Also, these fields
# are the important ones for the screening (see below).
reducing.fields.wos <- c("label","title","author","journal","issn","volume",
                         "number","pages","year","publisher","doi","abstract") #number to issue, 

wos.red <- wos[,reducing.fields.wos]


########################
# SCOPUS: general search

# importing the .bib files
scopus <- read_bibliography("literature_review/systematic_search/Scopus_IGEs_search.bib")

scopus$publisher <- "" #simply adding an empty column for the publisher since we did not extracted.


# reducing fields to the minimum number of fields
# so that all databases have the same columns. Also, these fields
# are the important ones for the screening (see below).
reducing.fields.scopus <- c("label","title","author","journal","issn","volume",
                            "number","pages","year","publisher","doi","abstract") #number to issue,

scopus.red <- scopus[,reducing.fields.scopus]


################################
# Snowbawling

# importing the .bib files
Bailey.scopus <- read_bibliography("literature_review/snowballing/Bailey_et_al_2017_cited_references_Scopus.bib")
Bijma.scopus <- read_bibliography("literature_review/snowballing/Bijma_and_Wade_2008_citing_references_Scopus.bib")
Ellen.scopus <- read_bibliography("literature_review/snowballing/Ellen_et_al_2014_cited_references_Scopus.bib")
Moore.scopus <- read_bibliography("literature_review/snowballing/Moore_et_al_1997_citing_references_Scopus.bib")
Wolf.scopus <- read_bibliography("literature_review/snowballing/Wolf_et_al_1998_citing_references_Scopus.bib")
#Bailey.wos <- read_bibliography("literature_review/snowballing/Bailey_et_al_2017_cited_references_WoS.txt") #not needed
Bijma.wos <- read_bibliography("literature_review/snowballing/Bijma_and_Wade_2008_citing_references_WoS.bib")
#Ellen.wos <- read_bibliography("literature_review/snowballing/Ellen_et_al_2014_cited_references_WoS.txt") #not needed
Moore.wos <- read_bibliography("literature_review/snowballing/Moore_et_al_1997_citing_references_WoS.bib")
Wolf.wos <- read_bibliography("literature_review/snowballing/Wolf_et_al_1998_citing_references_WoS.bib")


Bailey.scopus$publisher <- "" #simply adding an empty column for the publisher since we did not extracted.
Bijma.scopus$publisher <- "" #simply adding an empty column for the publisher since we did not extracted.
Ellen.scopus$publisher <- "" #simply adding an empty column for the publisher since we did not extracted.
Moore.scopus$publisher <- "" #simply adding an empty column for the publisher since we did not extracted.
Wolf.scopus$publisher <- "" #simply adding an empty column for the publisher since we did not extracted.


# reducing fields to the minimum number of fields
Bailey.scopus.red <- Bailey.scopus[,reducing.fields.scopus]
Bijma.scopus.red <- Bijma.scopus[,reducing.fields.scopus]
Ellen.scopus.red <- Ellen.scopus[,reducing.fields.scopus]
Moore.scopus.red <- Moore.scopus[,reducing.fields.scopus]
Wolf.scopus.red <- Wolf.scopus[,reducing.fields.scopus]
Bijma.wos.red <- Bijma.wos[,reducing.fields.scopus]
Moore.wos.red <- Moore.wos[,reducing.fields.scopus]
Wolf.wos.red <- Wolf.wos[,reducing.fields.scopus]


##############################################################
# Full reference data: before deduplication
##############################################################

wos.red$origin <- "wos.full"
scopus.red$origin <- "scopus.full"
Bailey.scopus.red$origin <- "Bailey.scopus"
Bijma.scopus.red$origin <- "Bijma.scopus"
Ellen.scopus.red$origin <- "Ellen.scopus"
Moore.scopus.red$origin <- "Moore.scopus"
Wolf.scopus.red$origin <- "Wolf.scopus"
Bijma.wos.red$origin <- "Bijma.wos"
Moore.wos.red$origin <- "Moore.wos"
Wolf.wos.red$origin <- "Wolf.wos"


# building the full reference list before deduplication
full.ref.data <- rbind(wos.red,
                       scopus.red,
                       Bailey.scopus.red,
                       Bijma.scopus.red,
                       Ellen.scopus.red,
                       Moore.scopus.red,
                       Wolf.scopus.red,
                       Bijma.wos.red,
                       Moore.wos.red,
                       Wolf.wos.red)

# write.csv(full.ref.data,"output_cleaning_reference_list/search_2/full_ref_data_search_2_bf_deduplication.csv",row.names=FALSE)
# full.ref.data <- read.table("output_cleaning_reference_list/search_2/full_ref_data_search_2_bf_deduplication.csv",
#                             header=T,sep=",")


##############################################################
# Full reference data: after deduplication
##############################################################

# # searching duplicates using revtools: 
# search.duplicated <- find_duplicates(data = full.ref.data,
#                                      match_variable = "title",
#                                      group_variable = NULL,
#                                      match_function = "fuzzdist",
#                                      method = "fuzz_m_ratio",
#                                      remove_punctuation = T,
#                                      threshold = 0.1)
# 
# # extracing duplicates
# screening.ref.data <- extract_unique_references(full.ref.data,search.duplicated)
# write.csv(screening.ref.data,"literature_review/search_unique_references_extracted.csv",row.names=FALSE)
# Note: we exported the file to manually find for duplicates in Excel and then 
# unfortunately forgot to modify the name of the file before importing it back
# Given that we provide the edited file, this should not be a problem for 
# the computational reproducibility of this code.
screening.ref.data <- read.table("literature_review/search_unique_references_extracted.csv",
                                 header=T,sep=",")


##############################################################
# Formatting data for RAYYAN QCRI
##############################################################

# choose only the fields needed for creating a .csv file importable by: https://rayyan.qcri.org

# example of a valid .csv file. The fields are the following:
# key,title,authors,journal,issn,volume,issue,pages,year,publisher,url,abstract
names.rayyan <- c("key","title","authors","journal","issn","volume","issue","pages","year","publisher","url","abstract")
names.rayyan
names(screening.ref.data)

# standardizing fields according to rayyan.example

# what's different between the two?
setdiff(names.rayyan,names(screening.ref.data))
setdiff(names(screening.ref.data),names.rayyan)


# excluding variables that are not needed
screening.ref.data$n_duplicates <- NULL
screening.ref.data$origin <- NULL


# what's different now?
setdiff(names.rayyan,names(screening.ref.data))
setdiff(names(screening.ref.data),names.rayyan)


screening.ref.data.rayyan <- plyr::rename(screening.ref.data, c("label"="key", "author"="authors", "number"="issue", "doi"="url"))
names(screening.ref.data.rayyan)

# reorder
screening.ref.data.rayyan <- screening.ref.data.rayyan[,names.rayyan]


# finding authors with missing initial(s) as that causes an error when importing into rayyan
table(grepl(",  ",screening.ref.data.rayyan$authors,fixed=T))

for(i in 1:nrow(screening.ref.data.rayyan)){
  
  if(grepl(",  ",screening.ref.data.rayyan$authors[i],fixed=T)){
    
    print(i)
  }
  
}

# manual fixes (i might differ for the same data)
# none


##############################################################
# Creating output
##############################################################

write.csv(screening.ref.data.rayyan[order(screening.ref.data.rayyan$title),],"literature_review/IGE_search_unique_references_rayyan.csv",row.names=FALSE)
#remember to manually remove the quotes for the column names only in the .csv file


# saving versions used for reproducibility purposes
sink("literature_review/deduplicating_Rpackages_session.txt")
sessionInfo()
sink()
