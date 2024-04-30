################################################################################
# Authors: 
#
# Alfredo Sanchez-Tojar (@ASanchez_Tojar)
# Profile: https://goo.gl/PmpPEB
# Bielefeld University
# Email: alfredo.tojar@gmail.com

# David Fisher
# Profile: https://evoetholab.com/
# University of Aberdeen
# Email: david.fisher@abdn.ac.uk

################################################################################
# Description of script and Instructions####
################################################################################

# Script first created the 10th of November 2021

# This script is to run a meta-analysis on Indirect Genetic
# Effects across species.

# Revised between the 5th and the 9th of December 2022 by David Fisher, Francesca 
# Santostefano, and Maria Moiron
# Completed between the 2nd and the 25th August 2023 by David Fisher
# Revised and finished in September 2023 by Alfredo Sanchez-Tojar
# Computational reproducibility checked by Maria Moiron in September 2023


################################################################################
# README (added here to make things easier, might delete from final version)
################################################################################

# Metadata with description of variables used for the analyses presented in the paper
# "xxx"
# F. Santostefano, D. Fisher, M. Moiron, A. Sanchez-Tojar
# 
# Note: some unused variables may get deleted, to be updated after we have the final dataset for the analyses
# 
# Record_id: Progressive row number
# Paper_id: unique identifier assigned to a paper
# Year: publication year
# Species_name: latin name of the species 
# Taxon: taxon to which the study species belongs. Categorical: insects/mammals/reptiles/amphibians/birds/fish
# Population: name of the population, combination of country and location or institution. Categorical: captive/semicaptive/wild 
# Population_type: captive (>5 generations in captivity) / semi-captive (<5 generations in captivity) / wild
# Study_type: Experimental: researchers set up the interactions; Observational: no manupulation of interactions. 
# Sex: Sex of the individuals. male/female/both/NA
# Age: Age the individuals. adult/juvenile/both/NA
# N_id_pedigree: # of total individuals in the pedigree 
# N_records: # of total records in the study
# N_id_w_records: # of individuals with records
# N_sires: # sires in the pedigree
# N_dams: # dams in the pedigree
# N_families: # families in the pedigree
# Depth_of_ped: # of max generations in the pedigree
# Trait_name: name of trait used by the study
# Trait_category: Type of trait: development/metabolism and physiology/morphology/reproduction/survival 
# Mean_standardized: whether the trait was mean standardized. yes/no	
# Trait_mean: trait mean value
# Variance_standardized_ whether the trait was variance standardized. yes/no	
# Trait_sd: whether the trait was mean standardized. yes/no	
# Teatment_group: Name of treatment or group if there are any	
# Fixed_eff_of_partner_trait: Does the model used to estiamte IGEs include partner traits as fixed effects? yes/no	
# Other_fixed_eff:  Does the model used to estimate IGEs include fixed effects? yes/no	 	
# Mean_group_size: group size of the group of individuals where IGEs are estimated. Mean if there are groups with different sizes.	
# Va: Addittive genetic variance from the main model estimating IGEs
# V_ige:Indirect genetic effects variance from the main model estimating IGEs	
# Vpe_foc: Permanent environmental effects variance (focal) from the main model estimating IGEs	
# Vpe_soc: Permanent environmental effects variance (social) from the main model estimating IGEs	
# V_other_1: Other variance component	from the main model estimating IGEs
# V_other_2: Other variance component	from the main model estimating IGEs
# V_other_3: Other variance component	from the main model estimating IGEs
# V_other_4:Other variance component	from the main model estimating IGEs
# V_other_5:Other variance component	from the main model estimating IGEs
# V_other_6:Other variance component from the main model estimating IGEs
# V_residual:Residual variance from the main model estimating IGEs
# Total_v_phen: Total phenotypic variance in the trait from the main model estimating IGEs
# H2: Trait heritability from the main model estimating IGEs as presented by the authors (VA/VP) 	
# Social_h2: VIGE/VP (equivalent of h2 but for IGEs) from the main model estimating IGEs as presented by the authors	
# R_a_ige: genetic correlation between VA and VIGE, from the main model estimating IGEs	
# R_pe_pe_soc: permanent environmental correlation between focal and social from the main model estimating IGEs	
# Cov_a_ige: genetic covariance between VA and VIGE, from the main model estimating IGEs		
# Cov_pe_pe_soc: permanent environmental covariance between focal and social from the main model estimating IGEs	
# V_tbv: Variance in total breeding value	
# T2: estimated as the variance in total breeding value / Total phen variance	
# Data_location: Page or table location of the data 	
# Screener_id: Initials of the screener who extracted the data	
# Notes: Comments from the first screener	
# Second_screener_id: Initials of the screener who extracted the data	
# Second_screener_notes:  Comments from the second screener		
# Changes_after_second_screening: Changes to data extraction after revision from the second screener	
# Data_transf: whether the data was transformed prior to the analyses, and if yes, what type of transformation. arcsin, non_transf, log, sqrt 	
# Trait_mean_transf: if the reported mean for transformed traits was before or after transformation	
# Trait_error_distrib: Error distribution for the modelled trait. : gaussian, binomial_logit, poisson_loglink, bernoulli_logit, binary_ordinal, NA 	
# Estimates_scale_nongauss. Category: latent, liability, link, NA.  If data was modelled assuming a non-Gaussian distribution, write down the scale of the reported estimates	
# Vp_group_formula: yes/no. Yes whenever vp is estimated not as a sum of the variance components but including also either group size or other components (e.g. cage mates) 
# Species_name.2: recoded species name	
# Population2: recoded population name	
# Group_id: paper ID recoded  	
# Second_screener_original_decision: decision of the second screener, whether to accept the values of the first screener o discuss/change 	
# missingSocial_h2: whether the paper is missing social_h2 (yes/no)	
# missingVP: whether the paper is missing VP (yes/no)		
# missingVIGE: whether the paper is missing VIGE (yes/no)	
# missingMean: whether the paper is missing the trait mean (yes/no)	
# missingVA: whether the paper is missing the trait VA (yes/no)		
# missingIdRecords: whether the paper is missing the # of records (yes/no)		
# Total_v_phen2: Total_v_phen2 is calculated as:  = 1 if trait is variance standardized, else copy value of reported VP from Total_v_phen	
# Trait_mean2: mean = 0 if trait is mean standardized, else copy value of reported mean,	
# Reported_values_flag: to flag studies where instead of using our calculated variables like ratios we should use what they report, it is usually non gaussian traits and specific cases we went through	
# Excluded_after_control: NA	
# social_h2: NA

#Some decisions made####

# We agreed to remove non-Gaussian traits, there aren???t many (9 I think) and they cause headaches for various reasons. Removing them removed the 2 big outliers on the funnel plots.
# 
# For VP, while summing all the reported variance components seems a rubbish way of back calculating VP, dividing Va by H2 actually seems quite good. Therefore, we think we should use the value for VP the authors reported if it is there, and use Va/H2 if it is not (new variable: ???Total_v_phen3??? which is used throughout moving forward)).
# 
# For H2 and Social H2, we are happy to use the author reported values if it is there, and do the calculations Va/VP and Vige/VP if not reported. I???ve moved the new back-calculation of VP to the start of the script so that we can use the new values for these calculations if necessary. Doing this adds some values but not a very large number (new variables: ???H2_2??? and ???Social_h2_2??? which are to be used for the analyses).
# 
# We looked at doing something similar for calculating the total heritability (T2), as we can 1) calculate the Total breeding values (TBV) for studies that do not report it and 2) do TBV/VP to get T2 for studies that do not report T2. However, this does not actually add that many values (110 -> 123), and the relationship between reported and calculated T2 (in studies where we can get both) isn???t as good as for the other variable we???ve back-calculated, so for this variable it seems better to stick to only author reported values. Perhaps this slight inconsistency will seem odd to readers though?
#   
#   An important variable we want to analyse is the DGE-IGE correlation. We can calculate this if authors do not report it, and then use either reported values or this one for analyses moving forward (new variable: ???R_a_ige_2??? which is to be used for the analyses). We think we can subject this to formal meta-analysis since it is bound between -1 and 1 but the function r.to.Zr did not work, perhaps due to the negative values?
#   
#   We also discussed what measure of evolvability we want to do and settled on ???I??? e.g. I_A= Va/trait_mean^2 and I_IGE= Vige/trait_mean^2.  We now calculate that for most values in the script, but not if the trait has been variance standardised and not if the trait had a very low mean (e.g. < 0.001) as this created huge values which just seem unrealistic.
#   
#   Finally, is there something you think we should be doing to in general assess how important IGEs are compared to just DGEs? Are we comparing means of Va vs Vige, h2 vs h2 social, h2 vs T2, I_A vs I_IGE, or something more than that?
#     
#     As a reminder, these are the questions we want to answer, with the names of the variable to use in the analysis from this updated R script:
#     
#     What is the proportion of the trait variance attributable to IGE?
#     Test: meta-analysis of h2 social (i.e., VIGE/VP)
#   
#   R variable name: social_h2_2
#   
#   Explore the moderators influence on IGE
#   Test: test for all moderators
#   
#   R variable name: multiple, described in the R code
#   
#   What is the relative importance of Direct versus Indirect genetic effects?
#     Test: comparison of Vdirect vs Vsocial; h2 direct vs h2social; and IA direct vs IA social
#   
#   R variable name: VA, V_IGE, H2_2, social_h2_2, I_A, I_IGE
#   
#   Do IGE change evolutionary trajectories?
#     Test: comparison of h2 and ??2TBV
#   R variable name: H2_2, T2
#   
#   Test: meta-analysis of genetic correlation between Direct versus Indirect genetic effects
#   R variable name: R_A_IGE_2

################################################################################
# Packages needed####
################################################################################

# load packages
pacman::p_load(ape,metafor,rotl,treebase,dplyr,stringr,ggplot2,tidyverse,
               #orchaRd, 
               devtools, patchwork, R.rsp, emmeans,
               huxtable,dfoptim,subplex, BB)# packages for rma.mv optimizers used below

#devtools::install_github("daniel1noble/orchaRd", force = TRUE)
library(orchaRd)

# cleaning up
rm(list=ls())


################################################################################
# Functions needed####
################################################################################

# function to convert r (ICC) to Zr
# Here, r correspond to our variance standardized variance estimates and the
# correlations
r.to.Zr <- function(r){
  Zr <- round(0.5*(log(1+r)-log(1-r)),3)
}

# function to obtain variance of Zr
VZr <- function(N){
  VZr <- 1/(N-3)
}

# function to back-transform r from Zr
Zr.to.r <- function(Zr){
  Zr <- round((exp(2*Zr)-1)/(exp(2*Zr)+1),3)
}


################################################################################
# Data preparation ####
################################################################################

# setp = to create a unique dataset where all relevant NAs are removed, to then 
# proceed to build the VCV and phylogenetic tree

# importing final dataset
dataset.IGE <- read.csv("data/dataset_final_after_cleaning_and_adding_author_contact_FS_MM.csv",
                        header=T)

nrow(dataset.IGE)
length(unique(dataset.IGE$Paper_id))

# removing humans from the dataset (only 1 row) and Eucalyptus (the only plant,
# 7 rows)
dataset.IGE <- as.data.frame(dataset.IGE 
                             %>% filter(Species_name.2!="Homo sapiens",
                                        Species_name.2!="Eucalyptus globulus"))

nrow(dataset.IGE)
length(unique(dataset.IGE$Paper_id))

# Include only traits that are Gaussian (this excludes 10 rows only)
dataset.IGE <- as.data.frame(dataset.IGE 
                             %>% filter(Trait_error_distrib == "gaussian"))

nrow(dataset.IGE)
length(unique(dataset.IGE$Paper_id))

# removing variables that are no longer useful at this stage to keep the dataset
# clean and more easy-to-use
summary(dataset.IGE$Trait_mean2==0)
dataset.IGE <- as.data.frame(dataset.IGE %>% 
                               select(-c(X,missingSocial_h2,missingVP,
                                         missingVIGE,missingMean,missingVA,
                                         missingIdRecords)))

# converting (character) variables to factors
names <- c('Record_id','Paper_id','Year','Species_name','Taxon',
           'Population_type','Study_type','Sex','Age','Trait_category',
           'Mean_standardized','Variance_standardized',
           'Fixed_eff_of_partner_trait','Other_fixed_eff','Screener_id',
           'Second_screener_id','Population2','Group_id')#,'Species_name.2'

dataset.IGE[,names] <- lapply(dataset.IGE[,names] , factor)


################################################################################
# Quick data exploration ####
################################################################################

# number of unique studies
length(unique(dataset.IGE$Paper_id))

# number of unique species
length(unique(dataset.IGE$Species_name.2))

# number of unique studies per species
as.data.frame(dataset.IGE %>% 
                group_by(Species_name.2) %>% 
                summarise(count = n_distinct(Paper_id)))

#Maria
a=as.data.frame(dataset.IGE %>% 
                  group_by(Species_name.2) %>% 
                  summarise(count = n_distinct(Paper_id)))
mean(a$count)
median(a$count)
max(a$count)

# quick triple-check for duplicated rows, which should not exist
table(duplicated(dataset.IGE[,-1]))


################################################################################
# Phylogeny information ####
################################################################################

# keep in mind that a different phylogenetic matrix will have to be created or
# subset for the different models if those models differ in which species are
# included due to missing information here and there, which can be a bit
# cumbersome.

# ##### loading the taxonomic data created on the 25th Aug 2023. It is important
# to load the created tree rather than creating a new one each time to ensure
# reproducibility. The reason being that information, format, etc is constantly 
# updated at the Open Tree of Life
load("data/taxa_Open_Tree_of_Life.RData") #taxa
# taxa <- tnrs_match_names(names = unique(dataset.IGE$Species_name.2))
# save(taxa, file = "data/taxa_Open_Tree_of_Life.RData")
# 
# # retrieving phylogenetic relationships among taxa in the form of a trimmed
# # sub-tree
# tree <- tol_induced_subtree(ott_ids =
#                                taxa[,"ott_id"],
#                              label_format = "name")
# 
# # following code is to exclude Eucalyptus, which is the only plant species, from
# # the phylo to see the consequences of doing so when computing the matrix of
# # phylogenetic relationships: as far as I can see, the matrix seems to work a
# # bit better when keeping Eucalyptus
# tree <- tol_induced_subtree(ott_ids = taxa[
#    taxa$search_string!="homo sapiens" &
#     taxa$search_string!="eucalyptus globulus","ott_id"],
#   label_format = "name")
# 
# 
# # we need to check for the existence of polytomies
# is.binary(tree) # yes, meaning there are no polytomies. We can proceed as it is
# 
#  
# # to confirm that our tree covers all the species we wanted it to include, and
# # make sure that the species names in our database match those in the tree, we
# # use the following code
# tree$tip.label <- gsub("_"," ", tree$tip.label)
# intersect(as.character(tree$tip.label), as.character(dataset.IGE$Species_name.2))
# setdiff(as.character(dataset.IGE$Species_name.2), as.character(tree$tip.label)) #listed in our database but not in the tree
# setdiff(as.character(tree$tip.label),as.character(dataset.IGE$Species_name.2)) # listed in the tree but not in our database
# 
# 
# # let's fix the following name in the tree and check again: all in order
# tree$tip.label[tree$tip.label == "Gadus morhua (in domain Bacteria)"] <- "Gadus morhua"
# intersect(as.character(tree$tip.label), as.character(dataset.IGE$Species_name.2))
# setdiff(as.character(dataset.IGE$Species_name.2), as.character(tree$tip.label)) #listed in our database but not in the tree
# setdiff(as.character(tree$tip.label),as.character(dataset.IGE$Species_name.2)) # listed in the tree but not in our database
#   
# # finally, save matrix for future analyses
# save(tree, file = "data/tree.Rdata")

# we can now load the saved tree
load("data/tree.Rdata") #tree


# # compute branch lengths of tree
# phylo_branch <- compute.brlen(tree, method = "Grafen", power = 1)
# 
# # check tree is ultrametric
# is.ultrametric(phylo_branch) # TRUE
# 
# # matrix to be included in the models
# phylo_cor <- vcv(phylo_branch, cor = T)
# 
# # finally, save matrix for future analyses
# save(phylo_cor, file = "data/phylo_cor.Rdata")

# we can now load the saved matrix
load("data/phylo_cor.Rdata") #phylo_cor

# # we can then plot the tree
# plot(tree, cex=1, label.offset =.1, no.margin = TRUE)


################################################################################
# Defining questions and variables ####
# QUESTIONS OF THE PROJECT AND VARIABLES WE NEED TO ANSWER THEM, including
# some further data exploration and calculation
################################################################################

### IMPORTANT NOTE: we decided to contact all the authors who didn't provide VP

# Additionally, calculating VP as Va / H2 seemed quite effective
# We add the calculation of VP up here so that the new VPs can be used to help
# fill in any missing Social H2s

################################################################################
# Calculating VP by dividing VA by H2, which leads to 56 VP values
dataset.IGE$Total_v_phen_calc <- dataset.IGE$Va / dataset.IGE$H2

with(dataset.IGE, plot(Total_v_phen2, Total_v_phen_calc))
abline(a = 0, b = 1) 
sum(is.na(dataset.IGE$Total_v_phen2))

# we only use the VP calculated value if the authors did not provide the value 
# otherwise

dataset.IGE$Total_v_phen3 = ifelse(is.na(dataset.IGE$Total_v_phen2),
                                   dataset.IGE$Total_v_phen_calc,
                                   dataset.IGE$Total_v_phen2)

sum(is.na(dataset.IGE$Total_v_phen3)) # There are now 39 missing VPs, meaning
# that we managed to add 13 VP values to the dataset by our calculation


################################################################################
# 1. A) What is the proportion of variance in traits explained by IGEs?

# OTHER  METAANALYSES ANALYZING VARIANCE/VP did:
# Moore et al. take m2 or they do the vM/VP ratio, did not mention lack of VP 
# Dochtermann gets directly h2, didn't have to do ratios to VP
# Holtmann et al. also gets directly repeatability (r or ICC) 

# We need the IGE/VP ratio 
# What is the number of effect sizes with NO V_ige?
sum(is.na(dataset.IGE$V_ige))

# What is the number of effect sizes with NO VP (i.e. Total_v_phen3)? [given that 
# Total_v_phen2 was calculated as 1 if trait was variance standardized, else 
# copied VP from Total_v_phen and then additional VP's were calculated and added
# for missing VP's]
sum(is.na(dataset.IGE$Total_v_phen3))

# What is the number of effect sizes with NO Social_h2?
sum(is.na(dataset.IGE$Social_h2)) 

# When Variance is standardized, Total_v_phen2 is always set to 1. Check:
summary(dataset.IGE[dataset.IGE$Variance_standardized=="yes","Total_v_phen2"])

# That means, we can proceed to estimate the ratio of variance explained by IGE as
dataset.IGE$Social_h2_calc <- dataset.IGE$V_ige/dataset.IGE$Total_v_phen3

with(dataset.IGE, plot(Social_h2, Social_h2_calc, ylim=c(0,0.15)))
abline(0,1) 

# Use social H2 if the authors reported it, or calculate it ourselves if they
# did not
dataset.IGE$Social_h2_2 <- ifelse(is.na(dataset.IGE$Social_h2), 
                                  dataset.IGE$Social_h2_calc,
                                  dataset.IGE$Social_h2)

hist(dataset.IGE$Social_h2_2) # None over 1.

################################################################################
# 1.	B) Are there differences between: (i.e. moderators for the meta-analysis  
# on point 1A)

# this checking was all done previously

# moderators
summary(dataset.IGE$Population_type)
summary(dataset.IGE$Study_type) #only 6 observational
summary(dataset.IGE$Sex)
summary(dataset.IGE$Age)
summary(dataset.IGE$Trait_category)
summary(dataset.IGE$Fixed_eff_of_partner_trait)
summary(dataset.IGE$Other_fixed_eff) #likely not to be used, not very informative, and really unbalanced


################################################################################
# 2. What's the relative and absolute importance of direct and indirect genetic 
# effects?(i.e. compare VDGE to the VIGEs)

# This can be done i.e. by
# 1. calculating evolvability. Two ways: 
#     - CVa, where CVa = (sqrt(Va)/Trait_mean))
#     - CViges, where CViges = (sqrt(V_ige)/Trait_mean))

# 2. mean standardized variances:
#     - Ia, where Ia = Va/(Trait_mean^2)
#     - Iiges, where Iiges = V_ige/(Trait_mean^2)

# 3. in addition to simply comparing the ratios  Va/Total_v_phen2 and 
# V_ige/Total_v_phen2, that is, H2.calculated vs. Social_h2.calculated


# What is the number of effect sizes with NO Va?
sum(is.na(dataset.IGE$Va))

# What is the number of effect sizes with NO H2?
sum(is.na(dataset.IGE$H2)) 


# we can proceed to estimate the ratio of variance explained by DGE 
# (i.e. H2.calculated) directly
dataset.IGE$H2_calc <- dataset.IGE$Va/dataset.IGE$Total_v_phen3 

with(dataset.IGE, plot(H2, H2_calc))
abline(0,1) # most squarely on the line, but some variation around that

# Use  H2 if the authors reported it, or calculate it ourselves if they did not
dataset.IGE$H2_2 <- ifelse(is.na(dataset.IGE$H2), 
                           dataset.IGE$H2_calc,
                           dataset.IGE$H2)

hist(dataset.IGE$H2_2) # None over 1. 


################################################################################
# ANALYSES
################################################################################

################################################################################
# Q1: What is the proportion of variance in traits explained by IGEs?       ####

# OTHER  META-ANALYSES ANALYZING VARIANCE/VP did:
# Moore et al. take m2 or they do the vM/VP ratio, did not mention lack of VP 
# Dochtermann gets directly h2, didn't have to do ratios to VP
# Holtmann et al. also gets directly repeatability (r or ICC) 


# transforming Social_h2.calculated using Fisher's Z transformation, which, 
# contrary to the (genetic) correlations, is unbounded and normally distributed
# (though keep in mind that Social_h2.calculated can only be positive [0,1])

# quick check: No effect sizes larger than 1
dataset.IGE[dataset.IGE$Social_h2_2>1 & 
              !(is.na(dataset.IGE$Social_h2_2)), ]


# transforming social_h2 (which is a type of ICC) to its Fisher's Z equivalent
dataset.IGE$Social_h2_2_Zr <- r.to.Zr(dataset.IGE$Social_h2_2)
hist(dataset.IGE$Social_h2_2_Zr,breaks=100)


# calculating the sampling variance in Zr (i.e. VZr) using the number of 
# individuals with records as our unit of replication. VZr corresponds to the
# sampling variance for all the analyses 
dataset.IGE$VZr <- VZr(dataset.IGE$N_id_w_records)
hist(dataset.IGE$VZr,breaks=100)
hist(dataset.IGE$N_id_w_records,breaks=100)


# first subset the dataset so that only full entries are available: by full
# we mean those with both Zr and VZr. This subset is used for both Q1 and Q2.
# Note that we did not have access the N_id_w_records of two effect sizes, which 
# means that we had to exclude one Social_h2.calculated effect size due to this.

dataset.IGE.subset1A <- as.data.frame(dataset.IGE 
                                      %>% filter(!(is.na(Social_h2_2_Zr)) & 
                                                   !(is.na(VZr))))

# some numbers and exploration
nrow(dataset.IGE.subset1A)
length(unique(dataset.IGE.subset1A$Paper_id))
length(unique(dataset.IGE.subset1A$Species_name.2))
summary(dataset.IGE.subset1A)
# 146 effect sizes, 40 papers, 19 species


# The next step if we want to do robust analyses that account for the 
# non-independence of sampling variances from the same study is to fit
# var-covar for the sampling variances, which can be calculated as follows

################################################################################
# Sampling variance-covariance matrix

# We specified sampling variance as a variance-covariance matrix, with the 
# sampling variance for each effect size on the diagonal, and the covariance 
# between these measures as off-diagonal elements. The model assumes a 0.5 
# correlation between the effect size sample variances with the same Paper_id. 
# For a similar approach, see O'Dea et al., 2019
# (https://doi.org/10.1111/faf.12394)or Kim et al. 2022 
# (https://doi.org/10.1111/1365-2656.13554), from which this code was taken from

# Creating a var-covar matrix assuming a 0.5 correlation between effect sizes 
# from the same Paper_id covariance = (0.5 * sqrt(vi.1) * sqrt(vi.2))

# Creates a matrix (called 'subset1A_VCV_ESVar') with the dimensions =  
# n(effect_sizes) x n(effect_sizes)
subset1A_VCV_ESVar <- matrix(0, 
                             nrow=nrow(dataset.IGE.subset1A), 
                             ncol=nrow(dataset.IGE.subset1A))

# Names rows and columns for each Record_id
rownames(subset1A_VCV_ESVar) <- dataset.IGE.subset1A[,"Record_id"]
colnames(subset1A_VCV_ESVar) <- dataset.IGE.subset1A[,"Record_id"]

# Finds effect sizes that come from the same study
shared_coord.subset1A <- which(dataset.IGE.subset1A[,"Paper_id"] %in% 
                                 dataset.IGE.subset1A[
                                   duplicated(dataset.IGE.subset1A[,"Paper_id"]), 
                                   "Paper_id"]==TRUE)

combinations.subset1A <- do.call("rbind", tapply(shared_coord.subset1A, 
                                                 dataset.IGE.subset1A[
                                                   shared_coord.subset1A, 
                                                   "Paper_id"], 
                                                 function(x) t(utils::combn(x, 
                                                                            2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1:dim(combinations.subset1A)[1]) {
  p1 <- combinations.subset1A[i, 1]
  p2 <- combinations.subset1A[i, 2]
  p1_p2_cov <- 0.5*
    sqrt(dataset.IGE.subset1A[p1, "VZr"])*
    sqrt(dataset.IGE.subset1A[p2, "VZr"])
  subset1A_VCV_ESVar[p1, p2] <- p1_p2_cov
  subset1A_VCV_ESVar[p2, p1] <- p1_p2_cov
} 

# Enters previously calculated effect size sampling variances into diagonals 
diag(subset1A_VCV_ESVar) <- dataset.IGE.subset1A[,"VZr"]

# # In case you want to visually double check the matrix outside of R
# write.csv(subset1A_VCV_ESVar, 'subset1A_VCV_ESVar.csv')


# Raw data funnel plot 

# # No more outliers after removing non-Gaussian traits and using parameters 
# # if they are reported and calculating them if not 
# par(mfrow=c(1, 2))
# funnel(dataset.IGE.subset1A$Social_h2_2_Zr, 
#        dataset.IGE.subset1A$VZr,
#        yaxis="sei",
#        xlab="Effect size (Zr)", 
#        digits=2, las=1) 
# funnel(dataset.IGE.subset1A$Social_h2_2_Zr, 
#        dataset.IGE.subset1A$VZr, 
#        yaxis="seinv",
#        xlab="Effect size (Zr)",  
#        digits=2, las=1) 


# creating a copy of Species_name.2 for phylogenetic effect
dataset.IGE.subset1A$Species_name.2.phylo <- dataset.IGE.subset1A$Species_name.2

# saving the subset for script 006_figures.R
write.csv(dataset.IGE.subset1A, file = "data/subsets/dataset_IGE_subset1A.csv",
          row.names = F)

# subsetting phylo_cor for this subset of the data
phylo_cor.subset1A <- phylo_cor[rownames(phylo_cor) %in% 
                                  unique(as.character(dataset.IGE.subset1A$Species_name.2)),
                                colnames(phylo_cor) %in% 
                                  unique(as.character(dataset.IGE.subset1A$Species_name.2))]


# MAIN MODEL = multilevel meta-analysis (i.e. the intercept-only model)
meta.model.IGE.subset1A <- rma.mv(Social_h2_2_Zr,
                                  subset1A_VCV_ESVar, 
                                  mods = ~ 1,
                                  random = list(~ 1 | Paper_id,
                                                ~ 1 | Group_id,
                                                ~ 1 | Population2,
                                                ~ 1 | Species_name.2,
                                                ~ 1 | Species_name.2.phylo,
                                                ~ 1 | Record_id),
                                  R = list(Species_name.2.phylo = phylo_cor.subset1A),
                                  method = "REML", 
                                  test = "t", 
                                  data = dataset.IGE.subset1A, 
                                  control=list(rel.tol=1e-9))

# saving the model for script 006_figures.R
save(meta.model.IGE.subset1A, file = "data/models/meta_model_IGE_subset1A.Rdata")

# model can be loaded instead of run using the following
# load("data/models/meta_model_IGE_subset1A.Rdata")


# Printing the summary results of the model
print(meta.model.IGE.subset1A, digits=3)

# Printing the results again, but adding the 95% prediction interval, 
# which uses the heterogeneity to generate an interval that should contain 95%
# of the effect sizes of any future or unknown studies with similar features
# as those included in the current database
predict(meta.model.IGE.subset1A, digits=3)

# # Model funnel plots
# par(mfrow = c(1, 1))
# funnel(meta.model.IGE.subset1A)


# # We can look for potential outliers based on the residuals of the model
# # any points below -2 or above 2 could be potential outliers
# resid <- rstandard(meta.model.IGE.subset1A)
# plot(1:nrow(dataset.IGE.subset1A), resid$z, type="b")
# abline(h = 2)
# 
# dataset.IGE.subset1A[resid$z > 2,]
# #A chicken paper provides most of these estimates and can't see anything 
# obviously suspicious. We used the VP they reported rather than calculating it 
# ourselves (although did use the social H2 we calculated) probably just higher 
# values in that context i.e. feather pecking 
# other value is a gull laying date paper - nothing obviously wrong


# using pluralistic approach to explore heterogeneity (Yang et al. 2023)

# why is the following not the same as above? Because the following is the 
# unstandardized raw variance whereas the previous is the sampling variance
round(sum(meta.model.IGE.subset1A$sigma2),4)
round(i2_ml(meta.model.IGE.subset1A),4)
round(cv_ml(meta.model.IGE.subset1A),4)
round(m1_ml(meta.model.IGE.subset1A),4)


# Output as table
ige_1a_res <- hux(predict(meta.model.IGE.subset1A, digits=3)) %>%
  add_rows(c("Prediction", "Standard error", "CI lower", "CI upper", "PI lower",
             "PI upper"), after = 0)
ige_1a_res


# SENSITIVITY ANALYSIS
# What happens when we repeat this only using values supplied by the authors 
# (so not calculating any ourselves)? The variable for that is "Social_h2"
# Note that the dataset for this analysis is much smaller (~1/4 of the original)
dataset.IGE$Social_h2_Zr <- r.to.Zr(dataset.IGE$Social_h2)

dataset.IGE.subset1Aii <- as.data.frame(dataset.IGE 
                                        %>% filter(!(is.na(Social_h2_Zr)) & 
                                                     !(is.na(VZr))))

# some numbers and exploration
nrow(dataset.IGE.subset1Aii)
length(unique(dataset.IGE.subset1Aii$Paper_id))
length(unique(dataset.IGE.subset1Aii$Species_name.2))
summary(dataset.IGE.subset1Aii)
# 40 effect sizes, 14 papers, 8 species
subset1Aii_VCV_ESVar <- matrix(0,
                               nrow=nrow(dataset.IGE.subset1Aii), 
                               ncol=nrow(dataset.IGE.subset1Aii))

# Names rows and columns for each Record_id
rownames(subset1Aii_VCV_ESVar) <- dataset.IGE.subset1Aii[,"Record_id"]
colnames(subset1Aii_VCV_ESVar) <- dataset.IGE.subset1Aii[,"Record_id"]

# Finds effect sizes that come from the same study
shared_coord.subset1Aii <- which(dataset.IGE.subset1Aii[,"Paper_id"] %in% 
                                   dataset.IGE.subset1Aii[
                                     duplicated(dataset.IGE.subset1Aii[,"Paper_id"]), 
                                     "Paper_id"]==TRUE)

combinations.subset1Aii <- do.call("rbind", tapply(shared_coord.subset1Aii, 
                                                   dataset.IGE.subset1Aii[
                                                     shared_coord.subset1Aii, 
                                                     "Paper_id"], 
                                                   function(x) t(utils::combn(x, 
                                                                              2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1:dim(combinations.subset1Aii)[1]) {
  p1 <- combinations.subset1Aii[i, 1]
  p2 <- combinations.subset1Aii[i, 2]
  p1_p2_cov <- 0.5*
    sqrt(dataset.IGE.subset1Aii[p1, "VZr"])*
    sqrt(dataset.IGE.subset1Aii[p2, "VZr"])
  subset1Aii_VCV_ESVar[p1, p2] <- p1_p2_cov
  subset1Aii_VCV_ESVar[p2, p1] <- p1_p2_cov
} 

# Enters previously calculated effect size sampling variances into diagonals 
diag(subset1Aii_VCV_ESVar) <- dataset.IGE.subset1Aii[,"VZr"]

# # Raw data funnel plot 
# par(mfrow=c(1, 2))
# funnel(dataset.IGE.subset1Aii$Social_h2_2_Zr, 
#        dataset.IGE.subset1Aii$VZr,
#        yaxis="sei",
#        xlab="Effect size (Zr)", 
#        digits=2, las=1) 
# funnel(dataset.IGE.subset1Aii$Social_h2_2_Zr, 
#        dataset.IGE.subset1Aii$VZr, 
#        yaxis="seinv",
#        xlab="Effect size (Zr)",  
#        digits=2, las=1) 


# creating a copy of Species_name.2 for phylogenetic effect
dataset.IGE.subset1Aii$Species_name.2.phylo <- dataset.IGE.subset1Aii$Species_name.2

# subsetting phylo_cor for this subset of the data
phylo_cor.subset1Aii <- phylo_cor[rownames(phylo_cor) %in% 
                                    unique(as.character(dataset.IGE.subset1Aii$Species_name.2)),
                                  colnames(phylo_cor) %in% 
                                    unique(as.character(dataset.IGE.subset1Aii$Species_name.2))]

# Fitting main model for this sensitivity analysis
meta.model.IGE.subset1Aii <- rma.mv(Social_h2_Zr,
                                    subset1Aii_VCV_ESVar, 
                                    mods = ~ 1,
                                    random = list(~ 1 | Paper_id,
                                                  ~ 1 | Group_id,
                                                  ~ 1 | Population2,
                                                  ~ 1 | Species_name.2,
                                                  ~ 1 | Species_name.2.phylo,
                                                  ~ 1 | Record_id),
                                    R = list(Species_name.2.phylo = phylo_cor.subset1Aii),
                                    method = "REML", 
                                    test = "t", 
                                    data = dataset.IGE.subset1Aii)

# Printing the summary results of the model
print(meta.model.IGE.subset1Aii, digits=3)


# Printing the results again, but adding the 95% prediction interval, 
# which uses the heterogeneity to generate an interval that should contain 95%
# of the effect sizes of any future or unknown studies with similar features
# as those included in the current database
predict(meta.model.IGE.subset1Aii, digits=3)

# Despite the much smaller dataset, the results are rather consistent both in
# terms of the meta-analytic mean and the uncertainty around it


################################################################################
# Q2: Does the magnitude of IGEs vary among with factors? ####
# 1.	B) Are there differences between: (i.e. moderators for the meta-analysis 
# on point 1A)
################################################################################

################################################################################
# Fixed effect of partner
summary(dataset.IGE.subset1A$Fixed_eff_of_partner_trait)
# 97 no, 49 yes
sum(summary(dataset.IGE.subset1A$Fixed_eff_of_partner_trait))

IGEmeta.regression.FEPartner <- rma.mv(Social_h2_2_Zr,
                                       subset1A_VCV_ESVar, 
                                       mods = ~ Fixed_eff_of_partner_trait,
                                       random = list(~ 1 | Paper_id,
                                                     ~ 1 | Group_id,
                                                     ~ 1 | Population2,
                                                     ~ 1 | Species_name.2,
                                                     ~ 1 | Species_name.2.phylo,
                                                     ~ 1 | Record_id),
                                       R = list(Species_name.2.phylo = phylo_cor.subset1A),
                                       method = "REML", 
                                       test = "t", 
                                       data = dataset.IGE.subset1A)

# saving the model for script 006_figures.R
save(IGEmeta.regression.FEPartner, 
     file = "data/models/IGEmeta_regression_FEPartner.Rdata")

# model can be loaded instead of run using the following
# load("data/models/IGEmeta_regression_FEPartner.Rdata")

# Printing the summary results of the model
print(IGEmeta.regression.FEPartner, digits=3)

# Calculate marginal R2 with r2_ml
R2.mr.FEPartner <- r2_ml(IGEmeta.regression.FEPartner) 
round(R2.mr.FEPartner*100, 1)

# Result: Although the estimate is negative as predicted, there is no
# statistically significant effect, with the estimate being very small and the
# uncertainty rather large. Also, very little variance explained.


# removing the intercept for easy visualization
IGEmeta.regression.FEPartner.NI <- rma.mv(Social_h2_2_Zr,
                                          subset1A_VCV_ESVar,
                                          mods = ~ Fixed_eff_of_partner_trait-1,
                                          random = list(~ 1 | Paper_id,
                                                        ~ 1 | Group_id,
                                                        ~ 1 | Population2,
                                                        ~ 1 | Species_name.2,
                                                        ~ 1 | Species_name.2.phylo,
                                                        ~ 1 | Record_id),
                                          R = list(Species_name.2.phylo = phylo_cor.subset1A),
                                          method = "REML",
                                          test = "t",
                                          data = dataset.IGE.subset1A, 
                                          control=list(rel.tol=1e-9))

print(IGEmeta.regression.FEPartner.NI, digits=3)


################################################################################
# Trait category
summary(dataset.IGE.subset1A$Trait_category)
table(dataset.IGE.subset1A$Trait_category)
# 60 behaviour, 37 morphology, 24 development, 12 survival, 10 reproduction, 
# 3 metabolism and physiology
sum(summary(dataset.IGE.subset1A$Trait_category))

# Although there are only 3 data points for "metabolism & physiology", we have
# decided to keep this level but keep this low sample size in mind when 
# interpreting the results for this level
IGEmeta.regression.TraitCat <- rma.mv(Social_h2_2_Zr,
                                      subset1A_VCV_ESVar, 
                                      mods = ~ Trait_category,
                                      random = list(~ 1 | Paper_id,
                                                    ~ 1 | Group_id,
                                                    ~ 1 | Population2,
                                                    ~ 1 | Species_name.2,
                                                    ~ 1 | Species_name.2.phylo,
                                                    ~ 1 | Record_id),
                                      R = list(Species_name.2.phylo = phylo_cor.subset1A),
                                      method = "REML", 
                                      test = "t", 
                                      data = dataset.IGE.subset1A, 
                                      control=list(rel.tol=1e-9)) # to allow true convergence, details explained in https://www.metafor-project.org/doku.php/tips:convergence_problems_rma_mv

# saving the model for script 006_figures.R
save(IGEmeta.regression.TraitCat, 
     file = "data/models/IGEmeta_regression_TraitCat.Rdata")

# model can be loaded instead of run using the following
# load("data/models/IGEmeta_regression_TraitCat.Rdata")

# Printing the summary results of the model
print(IGEmeta.regression.TraitCat, digits=3)


# Post-hoc Wald tests to test for statistically significant differences between 
# all levels
car::linearHypothesis(IGEmeta.regression.TraitCat, rbind(c(0,1,-1,0,0,0))) # development vs. metabolism & physiology
car::linearHypothesis(IGEmeta.regression.TraitCat, rbind(c(0,1,0,-1,0,0))) # development vs. morphology 
car::linearHypothesis(IGEmeta.regression.TraitCat, rbind(c(0,1,0,0,-1,0))) # development vs. reproduction
car::linearHypothesis(IGEmeta.regression.TraitCat, rbind(c(0,1,0,0,0,-1))) # development vs. survival

car::linearHypothesis(IGEmeta.regression.TraitCat, rbind(c(0,0,1,-1,0,0))) # metabolism & physiology vs. morphology
car::linearHypothesis(IGEmeta.regression.TraitCat, rbind(c(0,0,1,0,-1,0))) # metabolism & physiology vs. reproduction
car::linearHypothesis(IGEmeta.regression.TraitCat, rbind(c(0,0,1,0,0,-1))) # metabolism & physiology vs. survival

car::linearHypothesis(IGEmeta.regression.TraitCat, rbind(c(0,0,0,1,-1,0))) # morphology vs. reproduction
car::linearHypothesis(IGEmeta.regression.TraitCat, rbind(c(0,0,0,1,0,-1))) # morphology vs. survival

car::linearHypothesis(IGEmeta.regression.TraitCat, rbind(c(0,0,0,0,1,-1))) # reproduction vs. survival


# Calculate marginal R2 with r2_ml
R2.mr.TraitCat <- r2_ml(IGEmeta.regression.TraitCat) 
round(R2.mr.TraitCat*100, 1)


# removing the intercept for easy plotting and visualization
IGEmeta.regression.TraitCat.NI <- rma.mv(Social_h2_2_Zr,
                                         subset1A_VCV_ESVar, 
                                         mods = ~ Trait_category-1,
                                         random = list(~ 1 | Paper_id,
                                                       ~ 1 | Group_id,
                                                       ~ 1 | Population2,
                                                       ~ 1 | Species_name.2,
                                                       ~ 1 | Species_name.2.phylo,
                                                       ~ 1 | Record_id),
                                         R = list(Species_name.2.phylo = phylo_cor.subset1A),
                                         method = "REML", 
                                         test = "t", 
                                         data = dataset.IGE.subset1A, 
                                         control=list(rel.tol=1e-9)) 

print(IGEmeta.regression.TraitCat.NI, digits=3)


################################################################################
# Age
summary(dataset.IGE.subset1A$Age)
# 107 Adult, 35 juvenile, 4 both
sum(summary(dataset.IGE.subset1A$Age))

IGEmeta.regression.Age <- rma.mv(Social_h2_2_Zr,
                                 subset1A_VCV_ESVar, 
                                 mods = ~ Age,
                                 random = list(~ 1 | Paper_id,
                                               ~ 1 | Group_id,
                                               ~ 1 | Population2,
                                               ~ 1 | Species_name.2,
                                               ~ 1 | Species_name.2.phylo,
                                               ~ 1 | Record_id),
                                 R = list(Species_name.2.phylo = phylo_cor.subset1A),
                                 method = "REML", 
                                 test = "t", 
                                 data = dataset.IGE.subset1A)

# saving the model for script 006_figures.R
save(IGEmeta.regression.Age, 
     file = "data/models/IGEmeta_regression_Age.Rdata")

# model can be loaded instead of run using the following
# load("data/models/IGEmeta_regression_Age.Rdata")

# Printing the summary results of the model
print(IGEmeta.regression.Age, digits=3)


# Post-hoc Wald tests to test for statistically significant differences between 
# all levels
car::linearHypothesis(IGEmeta.regression.Age, rbind(c(0,1,-1))) # both vs. juv


# Calculate marginal R2 with r2_ml
R2.mr.Age <- r2_ml(IGEmeta.regression.Age) 
round(R2.mr.Age*100, 1)


# removing the intercept for easy plotting and visualization
IGEmeta.regression.Age.NI <- rma.mv(Social_h2_2_Zr,
                                    subset1A_VCV_ESVar, 
                                    mods = ~ Age-1,
                                    random = list(~ 1 | Paper_id,
                                                  ~ 1 | Group_id,
                                                  ~ 1 | Population2,
                                                  ~ 1 | Species_name.2,
                                                  ~ 1 | Species_name.2.phylo,
                                                  ~ 1 | Record_id),
                                    R = list(Species_name.2.phylo = phylo_cor.subset1A),
                                    method = "REML", 
                                    test = "t", 
                                    data = dataset.IGE.subset1A)

print(IGEmeta.regression.Age.NI, digits=3)


################################################################################
# Sex
summary(dataset.IGE.subset1A$Sex)
# 40 females only, 7 males only, 99 both
sum(summary(dataset.IGE.subset1A$Sex))

IGEmeta.regression.Sex <- rma.mv(Social_h2_2_Zr,
                                 subset1A_VCV_ESVar, 
                                 mods = ~ Sex,
                                 random = list(~ 1 | Paper_id,
                                               ~ 1 | Group_id,
                                               ~ 1 | Population2,
                                               ~ 1 | Species_name.2,
                                               ~ 1 | Species_name.2.phylo,
                                               ~ 1 | Record_id),
                                 R = list(Species_name.2.phylo = phylo_cor.subset1A),
                                 method = "REML", 
                                 test = "t", 
                                 data = dataset.IGE.subset1A)

# saving the model for script 006_figures.R
save(IGEmeta.regression.Sex, 
     file = "data/models/IGEmeta_regression_Sex.Rdata")

# model can be loaded instead of run using the following
# load("data/models/IGEmeta_regression_Sex.Rdata")

# Printing the summary results of the model
print(IGEmeta.regression.Sex, digits=3)


# Post-hoc Wald tests to test for statistically significant differences between 
# all levels
car::linearHypothesis(IGEmeta.regression.Sex, rbind(c(0,1,-1))) # F vs. M


# Calculate marginal R2 with r2_ml
R2.mr.Sex <- r2_ml(IGEmeta.regression.Sex) 
round(R2.mr.Sex*100, 1)

# removing the intercept for easy plotting and visualization
IGEmeta.regression.Sex.NI <- rma.mv(Social_h2_2_Zr,
                                    subset1A_VCV_ESVar, 
                                    mods = ~ Sex-1,
                                    random = list(~ 1 | Paper_id,
                                                  ~ 1 | Group_id,
                                                  ~ 1 | Population2,
                                                  ~ 1 | Species_name.2,
                                                  ~ 1 | Species_name.2.phylo,
                                                  ~ 1 | Record_id),
                                    R = list(Species_name.2.phylo = phylo_cor.subset1A),
                                    method = "REML", 
                                    test = "t", 
                                    data = dataset.IGE.subset1A)

print(IGEmeta.regression.Sex.NI, digits=3)


################################################################################
# Study type
summary(dataset.IGE.subset1A$Study_type)
# 140 experimental 6 observational
sum(summary(dataset.IGE.subset1A$Study_type))

IGEmeta.regression.StudyType <- rma.mv(Social_h2_2_Zr,
                                       subset1A_VCV_ESVar, 
                                       mods = ~ Study_type,
                                       random = list(~ 1 | Paper_id,
                                                     ~ 1 | Group_id,
                                                     ~ 1 | Population2,
                                                     ~ 1 | Species_name.2,
                                                     ~ 1 | Species_name.2.phylo,
                                                     ~ 1 | Record_id),
                                       R = list(Species_name.2.phylo = phylo_cor.subset1A),
                                       method = "REML", 
                                       test = "t", 
                                       data = dataset.IGE.subset1A)

# saving the model for script 006_figures.R
save(IGEmeta.regression.StudyType, 
     file = "data/models/IGEmeta_regression_StudyType.Rdata")

# model can be loaded instead of run using the following
# load("data/models/IGEmeta_regression_StudyType.Rdata")

# Printing the summary results of the model
print(IGEmeta.regression.StudyType, digits=3)


# Calculate marginal R2 with r2_ml
R2.mr.StudyType <- r2_ml(IGEmeta.regression.StudyType) 
round(R2.mr.StudyType*100, 1)


# removing the intercept for easy plotting and visualization
IGEmeta.regression.StudyType.NI <- rma.mv(Social_h2_2_Zr,
                                          subset1A_VCV_ESVar, 
                                          mods = ~ Study_type-1,
                                          random = list(~ 1 | Paper_id,
                                                        ~ 1 | Group_id,
                                                        ~ 1 | Population2,
                                                        ~ 1 | Species_name.2,
                                                        ~ 1 | Species_name.2.phylo,
                                                        ~ 1 | Record_id),
                                          R = list(Species_name.2.phylo = phylo_cor.subset1A),
                                          method = "REML", 
                                          test = "t", 
                                          data = dataset.IGE.subset1A)

print(IGEmeta.regression.StudyType.NI, digits=3)


################################################################################
# Population type
summary(dataset.IGE.subset1A$Population_type)
# 111 captive, 29 semicaptive, 6 wild

IGEmeta.regression.PopType <- rma.mv(Social_h2_2_Zr,
                                     subset1A_VCV_ESVar, 
                                     mods = ~ Population_type,
                                     random = list(~ 1 | Paper_id,
                                                   ~ 1 | Group_id,
                                                   ~ 1 | Population2,
                                                   ~ 1 | Species_name.2,
                                                   ~ 1 | Species_name.2.phylo,
                                                   ~ 1 | Record_id),
                                     R = list(Species_name.2.phylo = phylo_cor.subset1A),
                                     method = "REML", 
                                     test = "t", 
                                     data = dataset.IGE.subset1A)

# saving the model for script 006_figures.R
save(IGEmeta.regression.PopType, 
     file = "data/models/IGEmeta_regression_PopType.Rdata")

# model can be loaded instead of run using the following
# load("data/models/IGEmeta_regression_PopType.Rdata")

# Printing the summary results of the model
print(IGEmeta.regression.PopType, digits=3)


# Post-hoc Wald tests to test for statistically significant differences between 
# all levels
car::linearHypothesis(IGEmeta.regression.PopType, rbind(c(0,1,-1))) # semicaptive vs. wild


# Calculate marginal R2 with r2_ml
R2.mr.PopType <- r2_ml(IGEmeta.regression.PopType) 
round(R2.mr.PopType*100, 1)


# removing the intercept for easy plotting and visualization
IGEmeta.regression.PopType.NI <- rma.mv(Social_h2_2_Zr,
                                        subset1A_VCV_ESVar, 
                                        mods = ~ Population_type-1,
                                        random = list(~ 1 | Paper_id,
                                                      ~ 1 | Group_id,
                                                      ~ 1 | Population2,
                                                      ~ 1 | Species_name.2,
                                                      ~ 1 | Species_name.2.phylo,
                                                      ~ 1 | Record_id),
                                        R = list(Species_name.2.phylo = phylo_cor.subset1A),
                                        method = "REML", 
                                        test = "t", 
                                        data = dataset.IGE.subset1A)

print(IGEmeta.regression.PopType.NI, digits=3)


################################################################################
# Livestock # exploratory analyses after seeing the results

# unique(dataset.IGE.subset1A$Species_name.2) #of these, 
# #Neovison vison (farmed mink), Sus scrofa (pig), Bos taurus (cattle), and Gallus gallus (chicken), are  livestock 
# Could consider  Oryctolagus cuniculus (hare) livestock in this context, but not including for now) 

dataset.IGE.subset1A$Livestock = factor(ifelse(dataset.IGE.subset1A$Species_name.2 %in% 
                                                 c("Neovison vison",
                                                   "Sus scrofa", 
                                                   "Bos taurus",
                                                   "Gallus gallus"), 
                                               "yes", 
                                               "no" ))

summary(dataset.IGE.subset1A$Livestock)
# 70 yes, 76 no

IGEmeta.regression.Livestock <- rma.mv(Social_h2_2_Zr,
                                       subset1A_VCV_ESVar, 
                                       mods = ~ Livestock,
                                       random = list(~ 1 | Paper_id,
                                                     ~ 1 | Group_id,
                                                     ~ 1 | Population2,
                                                     ~ 1 | Species_name.2,
                                                     ~ 1 | Species_name.2.phylo,
                                                     ~ 1 | Record_id),
                                       R = list(Species_name.2.phylo = phylo_cor.subset1A),
                                       method = "REML", 
                                       test = "t", 
                                       data = dataset.IGE.subset1A)

# saving the model for script 006_figures.R
save(IGEmeta.regression.Livestock, 
     file = "data/models/IGEmeta_regression_Livestock.Rdata")

# model can be loaded instead of run using the following
# load("data/models/IGEmeta_regression_Livestock.Rdata")

# Printing the summary results of the model
print(IGEmeta.regression.Livestock, digits=3)


# Calculate marginal R2 with r2_ml
R2.mr.Livestock <- r2_ml(IGEmeta.regression.Livestock) 
round(R2.mr.Livestock*100, 1)


# removing the intercept for easy plotting and visualization
IGEmeta.regression.Livestock.NI <- rma.mv(Social_h2_2_Zr,
                                          subset1A_VCV_ESVar, 
                                          mods = ~ Livestock-1,
                                          random = list(~ 1 | Paper_id,
                                                        ~ 1 | Group_id,
                                                        ~ 1 | Population2,
                                                        ~ 1 | Species_name.2,
                                                        ~ 1 | Species_name.2.phylo,
                                                        ~ 1 | Record_id),
                                          R = list(Species_name.2.phylo = phylo_cor.subset1A),
                                          method = "REML", 
                                          test = "t", 
                                          data = dataset.IGE.subset1A)

print(IGEmeta.regression.Livestock.NI, digits=3)



################################################################################
# Q3:. What's the relative and absolute importance of direct and indirect  
# genetic effects?
# a: compare VDGE to the VIGEs
# b: compare H2 and social H2
# c: compare direct evolvability vs indirect evolvability (for choices on this 
# see below)
################################################################################

# For each we need a subset where no rows have NA
# We then re-arrange so that the values of the 2 variables are in a single 
# column, with another column defining whether direct or indirect.
# This column is then used as a moderator


################################################################################
# 3A: Va vs Vi####
# First we compare Va and Vi 

dataset.IGE.subset3A <- as.data.frame(dataset.IGE 
                                      %>% filter(!(is.na(Va)) &
                                                   !(is.na(V_ige)) &
                                                   !(is.na(VZr))))


# some numbers and exploration
nrow(dataset.IGE.subset3A)
length(unique(dataset.IGE.subset3A$Paper_id))
length(unique(dataset.IGE.subset3A$Species_name.2))
summary(dataset.IGE.subset3A)
# 160 effect sizes, 43 papers, 20 species

# there are 9 V_ige values reported as 0 in the original articles. 
table(dataset.IGE.subset3A$V_ige==0)

# In order to # be able to use these values for the analyses involving Va and
# V_ige, we have decided to convert those 0 to the mininum possible rounded 
# value based on how the 0 was reported. For example, if the value was reported 
# as 0.00 we would transform it into 0.001, or if the value was reported as 0, 
# we would transform it into 0.1. Below are the 9 values and how they were
# reported

dataset.IGE.subset3A[dataset.IGE.subset3A$V_ige==0,]
# Record_id 3; IGE0501: Reported as 0.000
# Record_id 22; IGE0857: Reported as 0.00
# Record_id 75; IGE0468: Reported as 0.000
# Record_id 79; IGE0468: Reported as 0.000
# Record_id 87; IGE0468: Reported as 0.000
# Record_id 88; IGE0468: Reported as 0.000
# Record_id 128; IGE0406: Reported as 0.000
# Record_id 132; IGE0406: Reported as 0.000
# Record_id 181; IGE0203: Reported as 0

# substituting all those values for their corresponding minimum value
dataset.IGE.subset3A$V_ige_2 <- ifelse((dataset.IGE.subset3A$Paper_id %in% 
                                          c("IGE0501","IGE0468","IGE0406")) & dataset.IGE.subset3A$V_ige==0,
                                       0.0001,
                                       ifelse((dataset.IGE.subset3A$Paper_id %in% 
                                                 c("IGE0857")) & dataset.IGE.subset3A$V_ige==0,
                                              0.001,ifelse((dataset.IGE.subset3A$Paper_id %in% 
                                                              c("IGE0203")) & dataset.IGE.subset3A$V_ige==0,
                                                           0.1,
                                                           dataset.IGE.subset3A$V_ige)
                                       ))

dataset.IGE.subset3A_long_no_zeros <- dataset.IGE.subset3A %>%
  pivot_longer(cols = c("Va", "V_ige_2"), 
               values_to = "dat3A_var", names_to = "direct_social") %>%
  mutate(dat3A_var_l = log(dat3A_var), # using logged value as it helps the model fit a very skewed data distribution (see histograms below)
         dat_weights = log(N_id_w_records)) %>% # for this analysis, we considered using the log of sample size as weights for the effect sizes but decided that using VZr, which is 1/(N-3) makes the results of all analyses more comparable to each other given that log(sample size) and 1/(N-3) are non-linearly associated
  filter(is.finite(dat3A_var_l)) %>% as.data.frame()


# given that the dataset is now "double" the size despite being the same data
# we need to account for this new level of nonindependence. The easiest way at
# this point is to continue using all random effects we used before, but add
# an additional one to account for the within-study/residual variance, given 
# that "Record_id" is no longer a unit-level random effect. Although we kept
# the same names, we need to interpret the random effects slightly different
# to what we did for the "non-long" models
dataset.IGE.subset3A_long_no_zeros$Record_id_long <- 1:nrow(dataset.IGE.subset3A_long_no_zeros)


hist(dataset.IGE.subset3A_long_no_zeros$dat3A_var, breaks = 100) # very skewed
hist(dataset.IGE.subset3A_long_no_zeros$dat3A_var_l, breaks = 100) # log helps making it look gaussian, but it also means we lose 9 effect sizes due to V_ige being 0 (from IGE0501 IGE0857 IGE0468 IGE0468 IGE0468 IGE0468 IGE0406 IGE0406 IGE0406 IGE0203)

nrow(dataset.IGE.subset3A_long_no_zeros)
length(unique(dataset.IGE.subset3A_long_no_zeros$Paper_id))
length(unique(dataset.IGE.subset3A_long_no_zeros$Species_name.2))
# 320 effect sizes, 43 papers, 20 species

# we decided to use the log transformed variance estimates despite having to 
# exclude 9 effect sizes because otherwise the models below have a really hard
# time converging. Nonetheless, models with non-log transformed values lead
# to the same qualitative result: that is, Va is considerably higher than V_ige

subset3A_long_no_zeros_VCV_ESVar <- matrix(0, 
                                           nrow=nrow(dataset.IGE.subset3A_long_no_zeros), 
                                           ncol=nrow(dataset.IGE.subset3A_long_no_zeros))

# Names rows and columns for each Record_id_long
rownames(subset3A_long_no_zeros_VCV_ESVar) <- dataset.IGE.subset3A_long_no_zeros[,"Record_id_long"]
colnames(subset3A_long_no_zeros_VCV_ESVar) <- dataset.IGE.subset3A_long_no_zeros[,"Record_id_long"]

# Finds effect sizes that come from the same study
shared_coord.subset3A_long_no_zeros <- which(dataset.IGE.subset3A_long_no_zeros[,"Paper_id"] %in% 
                                               dataset.IGE.subset3A_long_no_zeros[
                                                 duplicated(dataset.IGE.subset3A_long_no_zeros[,"Paper_id"]), 
                                                 "Paper_id"]==TRUE)

combinations.subset3A_long_no_zeros <- do.call("rbind", tapply(shared_coord.subset3A_long_no_zeros, 
                                                               dataset.IGE.subset3A_long_no_zeros[
                                                                 shared_coord.subset3A_long_no_zeros, 
                                                                 "Paper_id"], 
                                                               function(x) t(utils::combn(x, 
                                                                                          2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1:dim(combinations.subset3A_long_no_zeros)[1]) {
  p1 <- combinations.subset3A_long_no_zeros[i, 1]
  p2 <- combinations.subset3A_long_no_zeros[i, 2]
  p1_p2_cov <- 0.5*
    sqrt(dataset.IGE.subset3A_long_no_zeros[p1, "VZr"])*
    sqrt(dataset.IGE.subset3A_long_no_zeros[p2, "VZr"])
  # sqrt(dataset.IGE.subset3A_long[p1, "dat_weights"])* # testing whether using log(N_id_w_records) as our estimate of sampling variance for the varcovar matrix works. It does not.
  # sqrt(dataset.IGE.subset3A_long[p2, "dat_weights"])
  subset3A_long_no_zeros_VCV_ESVar[p1, p2] <- p1_p2_cov
  subset3A_long_no_zeros_VCV_ESVar[p2, p1] <- p1_p2_cov
} 

# Enters previously calculated effect size sampling variances into diagonals 
diag(subset3A_long_no_zeros_VCV_ESVar) <- dataset.IGE.subset3A_long_no_zeros[,"VZr"]


# creating a copy of Species_name.2 for phylogenetic effect
dataset.IGE.subset3A_long_no_zeros$Species_name.2.phylo <- dataset.IGE.subset3A_long_no_zeros$Species_name.2

phylo_cor.subset3A_no_zeros <- phylo_cor[rownames(phylo_cor) %in% 
                                           unique(as.character(dataset.IGE.subset3A_long_no_zeros$Species_name.2)),
                                         colnames(phylo_cor) %in% 
                                           unique(as.character(dataset.IGE.subset3A_long_no_zeros$Species_name.2))]



#library(dplyr)
dataset.IGE.subset3A_long_no_zeros%>%group_by(direct_social)%>%summarise(Average=log(mean(dat3A_var, rm.na=TRUE)))
dataset.IGE.subset3A_long_no_zeros$direct_social<- recode_factor(dataset.IGE.subset3A_long_no_zeros$direct_social, Va  = "V_a", V_ige = "V_ige")

# saving the subset for script 006_figures.R
write.csv(dataset.IGE.subset3A_long_no_zeros, file = "data/subsets/dataset_IGE_subset3A_long_no_zeros.csv",
          row.names = F)

meta.model.IGE.subset3A_no_zeros <- rma.mv(dat3A_var_l,
                                           subset3A_long_no_zeros_VCV_ESVar,
                                           mods = ~ direct_social,
                                           random = list(~ 1 | Paper_id,
                                                         ~ 1 | Group_id,
                                                         ~ 1 | Population2,
                                                         ~ 1 | Species_name.2,
                                                         ~ 1 | Species_name.2.phylo,
                                                         ~ 1 | Record_id_long,
                                                         ~ 1 | Record_id),
                                           R = list(Species_name.2.phylo = phylo_cor.subset3A_no_zeros),
                                           method = "REML",
                                           test = "t",
                                           data = dataset.IGE.subset3A_long_no_zeros)

# saving the model for script 006_figures.R
save(meta.model.IGE.subset3A_no_zeros, 
     file = "data/models/IGEmeta_regression_Va_vs_Vige.Rdata")

# model can be loaded instead of run using the following
# load("data/models/IGEmeta_regression_Va_vs_Vige.Rdata")

# Printing the summary results of the model
print(meta.model.IGE.subset3A_no_zeros, digits=3)

# Calculate marginal R2 with r2_ml
R2.IGE.subset3A_no_zeros <- r2_ml(meta.model.IGE.subset3A_no_zeros) 
round(R2.IGE.subset3A_no_zeros*100, 1)


# Run without intercept
meta.model.IGE.subset3A_no_zeros.NI <- rma.mv(dat3A_var_l,
                                              subset3A_long_no_zeros_VCV_ESVar,
                                              mods = ~ direct_social-1,
                                              random = list(~ 1 | Paper_id,
                                                            ~ 1 | Group_id,
                                                            ~ 1 | Population2,
                                                            ~ 1 | Species_name.2,
                                                            ~ 1 | Species_name.2.phylo,
                                                            ~ 1 | Record_id_long,
                                                            ~ 1 | Record_id),
                                              R = list(Species_name.2.phylo = phylo_cor.subset3A_no_zeros),
                                              method = "REML", 
                                              test = "t", 
                                              data = dataset.IGE.subset3A_long_no_zeros)

# Printing the summary results of the model
print(meta.model.IGE.subset3A_no_zeros.NI, digits=3)
predict(meta.model.IGE.subset3A_no_zeros.NI,transf=exp,digits=3)[1] # direct effect
predict(meta.model.IGE.subset3A_no_zeros.NI,transf=exp,digits=3)[2] # social effect
#social effect much lower, exp(-3.775) = 0.02293709 vs direct of exp(-0.019) = 0.9811794

################
# same analysis as above but without transforming 0 into their minimum value
# dataset.IGE.subset3A_long <- dataset.IGE.subset3A %>%
#   pivot_longer(cols = c("Va", "V_ige"), 
#                values_to = "dat3A_var", names_to = "direct_social") %>%
#   mutate(dat3A_var_l = log(dat3A_var), # using logged value as it helps the model fit a very skewed data distribution (see histograms below)
#          dat_weights = log(N_id_w_records)) %>% # for this analysis, we considered using the log of sample size as weights for the effect sizes but decided that using VZr, which is 1/(N-3) makes the results of all analyses more comparable to each other given that log(sample size) and 1/(N-3) are non-linearly associated
#   filter(is.finite(dat3A_var_l)) %>% as.data.frame()
# 
# 
# # given that the dataset is now "double" the size despite being the same data
# # we need to account for this new level of nonindependence. The easiest way at
# # this point is to continue using all random effects we used before, but add
# # an additional one to account for the within-study/residual variance, given 
# # that "Record_id" is no longer a unit-level random effect. Although we kept
# # the same names, we need to interpret the random effects slightly different
# # to what we did for the "non-long" models
# dataset.IGE.subset3A_long$Record_id_long <- 1:nrow(dataset.IGE.subset3A_long)
# 
# 
# hist(dataset.IGE.subset3A_long$dat3A_var, breaks = 100) # very skewed
# hist(dataset.IGE.subset3A_long$dat3A_var_l, breaks = 100) # log helps making it look gaussian, but it also means we lose 9 effect sizes due to V_ige being 0 (from IGE0501 IGE0857 IGE0468 IGE0468 IGE0468 IGE0468 IGE0406 IGE0406 IGE0406 IGE0203)
# 
# nrow(dataset.IGE.subset3A_long)
# length(unique(dataset.IGE.subset3A_long$Paper_id))
# length(unique(dataset.IGE.subset3A_long$Species_name.2))
# # 311 effect sizes, 43 papers, 20 species
# 
# # we decided to use the log transformed variance estimates despite having to 
# # exclude 9 effect sizes because otherwise the models below have a really hard
# # time converging. Nonetheless, models with non-log transformed values lead
# # to the same qualitative result: that is, Va is considerably higher than V_ige
# 
# subset3A_long_VCV_ESVar <- matrix(0, 
#                                   nrow=nrow(dataset.IGE.subset3A_long), 
#                                   ncol=nrow(dataset.IGE.subset3A_long))
# 
# # Names rows and columns for each Record_id_long
# rownames(subset3A_long_VCV_ESVar) <- dataset.IGE.subset3A_long[,"Record_id_long"]
# colnames(subset3A_long_VCV_ESVar) <- dataset.IGE.subset3A_long[,"Record_id_long"]
# 
# # Finds effect sizes that come from the same study
# shared_coord.subset3A_long <- which(dataset.IGE.subset3A_long[,"Paper_id"] %in% 
#                                       dataset.IGE.subset3A_long[
#                                         duplicated(dataset.IGE.subset3A_long[,"Paper_id"]), 
#                                         "Paper_id"]==TRUE)
# 
# combinations.subset3A_long <- do.call("rbind", tapply(shared_coord.subset3A_long, 
#                                                       dataset.IGE.subset3A_long[
#                                                         shared_coord.subset3A_long, 
#                                                         "Paper_id"], 
#                                                       function(x) t(utils::combn(x, 
#                                                                                  2))))
# 
# # Calculates the covariance between effect sizes and enters them in each 
# # combination of coordinates
# for (i in 1:dim(combinations.subset3A_long)[1]) {
#   p1 <- combinations.subset3A_long[i, 1]
#   p2 <- combinations.subset3A_long[i, 2]
#   p1_p2_cov <- 0.5*
#     sqrt(dataset.IGE.subset3A_long[p1, "VZr"])*
#     sqrt(dataset.IGE.subset3A_long[p2, "VZr"])
#   # sqrt(dataset.IGE.subset3A_long[p1, "dat_weights"])* # testing whether using log(N_id_w_records) as our estimate of sampling variance for the varcovar matrix works. It does not.
#   # sqrt(dataset.IGE.subset3A_long[p2, "dat_weights"])
#   subset3A_long_VCV_ESVar[p1, p2] <- p1_p2_cov
#   subset3A_long_VCV_ESVar[p2, p1] <- p1_p2_cov
# } 
# 
# # Enters previously calculated effect size sampling variances into diagonals 
# diag(subset3A_long_VCV_ESVar) <- dataset.IGE.subset3A_long[,"VZr"]
# 
# 
# # creating a copy of Species_name.2 for phylogenetic effect
# dataset.IGE.subset3A_long$Species_name.2.phylo <- dataset.IGE.subset3A_long$Species_name.2
# 
# phylo_cor.subset3A <- phylo_cor[rownames(phylo_cor) %in% 
#                                   unique(as.character(dataset.IGE.subset3A_long$Species_name.2)),
#                                 colnames(phylo_cor) %in% 
#                                   unique(as.character(dataset.IGE.subset3A_long$Species_name.2))]
# 
# 
# meta.model.IGE.subset3A <- rma.mv(dat3A_var_l,
#                                   subset3A_long_VCV_ESVar,
#                                   mods = ~ direct_social,
#                                   random = list(~ 1 | Paper_id,
#                                                 ~ 1 | Group_id,
#                                                 ~ 1 | Population2,
#                                                 ~ 1 | Species_name.2,
#                                                 ~ 1 | Species_name.2.phylo,
#                                                 ~ 1 | Record_id_long,
#                                                 ~ 1 | Record_id),
#                                   R = list(Species_name.2.phylo = phylo_cor.subset3A),
#                                   method = "REML",
#                                   test = "t",
#                                   data = dataset.IGE.subset3A_long)
# 
# 
# # Printing the summary results of the model
# print(meta.model.IGE.subset3A, digits=3)
# 
# 
# # Calculate marginal R2 with r2_ml
# R2.IGE.subset3A <- r2_ml(meta.model.IGE.subset3A) 
# round(R2.IGE.subset3A*100, 1)
# 
# 
# # Run without intercept
# meta.model.IGE.subset3A.NI <- rma.mv(dat3A_var_l,
#                                      subset3A_long_VCV_ESVar,
#                                      mods = ~ direct_social-1,
#                                      random = list(~ 1 | Paper_id,
#                                                    ~ 1 | Group_id,
#                                                    ~ 1 | Population2,
#                                                    ~ 1 | Species_name.2,
#                                                    ~ 1 | Species_name.2.phylo,
#                                                    ~ 1 | Record_id_long,
#                                                    ~ 1 | Record_id),
#                                      R = list(Species_name.2.phylo = phylo_cor.subset3A),
#                                      method = "REML", 
#                                      test = "t", 
#                                      data = dataset.IGE.subset3A_long)
# 
# # Printing the summary results of the model
# print(meta.model.IGE.subset3A.NI, digits=3)
# predict(meta.model.IGE.subset3A.NI,transf=exp) #to easily back transform estimates for visualization
# 
# 
# # # sensitivity analysis using non-log transform Va and V_ige estimates, which
# # # which takes some play around with optimizer to make it work. The commented
# # # optimizers are those that did not converge, the remaining ran, results are
# # # similar but it is unclear how much we can trust the results given the little
# # # discrepancies between the different optimizers
# # # more in: https://www.metafor-project.org/doku.php/tips:convergence_problems_rma_mv
# # res <- list(
# #   # rma.mv(dat3A_var,subset3A_long_VCV_ESVar,mods = ~ direct_social,random = list(~ 1 | Paper_id,~ 1 | Group_id,~ 1 | Population2,~ 1 | Species_name.2,~ 1 | Species_name.2.phylo,~ 1 | Record_id),
# #   #        R = list(Species_name.2.phylo = phylo_cor.subset3A),method = "REML",test = "t",data = dataset.IGE.subset3A_long
# #   #        ,control=list(optimizer="nlminb", rel.tol=1e-10)),
# #   # rma.mv(dat3A_var,subset3A_long_VCV_ESVar,mods = ~ direct_social,random = list(~ 1 | Paper_id,~ 1 | Group_id,~ 1 | Population2,~ 1 | Species_name.2,~ 1 | Species_name.2.phylo,~ 1 | Record_id),
# #   #        R = list(Species_name.2.phylo = phylo_cor.subset3A),method = "REML",test = "t",data = dataset.IGE.subset3A_long
# #   #        ,control=list(optimizer="Nelder-Mead")),
# #   # rma.mv(dat3A_var,subset3A_long_VCV_ESVar,mods = ~ direct_social,random = list(~ 1 | Paper_id,~ 1 | Group_id,~ 1 | Population2,~ 1 | Species_name.2,~ 1 | Species_name.2.phylo,~ 1 | Record_id),
# #   #        R = list(Species_name.2.phylo = phylo_cor.subset3A),method = "REML",test = "t",data = dataset.IGE.subset3A_long
# #   #        ,control=list(optimizer="BFGS")),
# #   rma.mv(dat3A_var,subset3A_long_VCV_ESVar,mods = ~ direct_social,random = list(~ 1 | Paper_id,~ 1 | Group_id,~ 1 | Population2,~ 1 | Species_name.2,~ 1 | Species_name.2.phylo,~ 1 | Record_id),
# #          R = list(Species_name.2.phylo = phylo_cor.subset3A),method = "REML",test = "t",data = dataset.IGE.subset3A_long
# #          ,control=list(optimizer="bobyqa")),
# #   rma.mv(dat3A_var,subset3A_long_VCV_ESVar,mods = ~ direct_social,random = list(~ 1 | Paper_id,~ 1 | Group_id,~ 1 | Population2,~ 1 | Species_name.2,~ 1 | Species_name.2.phylo,~ 1 | Record_id),
# #          R = list(Species_name.2.phylo = phylo_cor.subset3A),method = "REML",test = "t",data = dataset.IGE.subset3A_long
# #          ,control=list(optimizer="nloptr", maxeval=1000)),
# #   rma.mv(dat3A_var,subset3A_long_VCV_ESVar,mods = ~ direct_social,random = list(~ 1 | Paper_id,~ 1 | Group_id,~ 1 | Population2,~ 1 | Species_name.2,~ 1 | Species_name.2.phylo,~ 1 | Record_id),
# #          R = list(Species_name.2.phylo = phylo_cor.subset3A),method = "REML",test = "t",data = dataset.IGE.subset3A_long
# #          ,control=list(optimizer="nlm")),
# #   rma.mv(dat3A_var,subset3A_long_VCV_ESVar,mods = ~ direct_social,random = list(~ 1 | Paper_id,~ 1 | Group_id,~ 1 | Population2,~ 1 | Species_name.2,~ 1 | Species_name.2.phylo,~ 1 | Record_id),
# #          R = list(Species_name.2.phylo = phylo_cor.subset3A),method = "REML",test = "t",data = dataset.IGE.subset3A_long
# #          ,control=list(optimizer="hjk")),
# #   # rma.mv(dat3A_var,subset3A_long_VCV_ESVar,mods = ~ direct_social,random = list(~ 1 | Paper_id,~ 1 | Group_id,~ 1 | Population2,~ 1 | Species_name.2,~ 1 | Species_name.2.phylo,~ 1 | Record_id),
# #   #        R = list(Species_name.2.phylo = phylo_cor.subset3A),method = "REML",test = "t",data = dataset.IGE.subset3A_long
# #   #        ,control=list(optimizer="nmk")),
# #   rma.mv(dat3A_var,subset3A_long_VCV_ESVar,mods = ~ direct_social,random = list(~ 1 | Paper_id,~ 1 | Group_id,~ 1 | Population2,~ 1 | Species_name.2,~ 1 | Species_name.2.phylo,~ 1 | Record_id),
# #          R = list(Species_name.2.phylo = phylo_cor.subset3A),method = "REML",test = "t",data = dataset.IGE.subset3A_long
# #          ,control=list(optimizer="mads")),
# #   # rma.mv(dat3A_var,subset3A_long_VCV_ESVar,mods = ~ direct_social,random = list(~ 1 | Paper_id,~ 1 | Group_id,~ 1 | Population2,~ 1 | Species_name.2,~ 1 | Species_name.2.phylo,~ 1 | Record_id),
# #   #        R = list(Species_name.2.phylo = phylo_cor.subset3A),method = "REML",test = "t",data = dataset.IGE.subset3A_long
# #   #        ,control=list(optimizer="ucminf")),
# #   # rma.mv(dat3A_var,subset3A_long_VCV_ESVar,mods = ~ direct_social,random = list(~ 1 | Paper_id,~ 1 | Group_id,~ 1 | Population2,~ 1 | Species_name.2,~ 1 | Species_name.2.phylo,~ 1 | Record_id),
# #   #        R = list(Species_name.2.phylo = phylo_cor.subset3A),method = "REML",test = "t",data = dataset.IGE.subset3A_long
# #   #        ,control=list(optimizer="lbfgsb3c")),
# #   rma.mv(dat3A_var,subset3A_long_VCV_ESVar,mods = ~ direct_social,random = list(~ 1 | Paper_id,~ 1 | Group_id,~ 1 | Population2,~ 1 | Species_name.2,~ 1 | Species_name.2.phylo,~ 1 | Record_id),
# #          R = list(Species_name.2.phylo = phylo_cor.subset3A),method = "REML",test = "t",data = dataset.IGE.subset3A_long
# #          ,control=list(optimizer="subplex")),
# #   rma.mv(dat3A_var,subset3A_long_VCV_ESVar,mods = ~ direct_social,random = list(~ 1 | Paper_id,~ 1 | Group_id,~ 1 | Population2,~ 1 | Species_name.2,~ 1 | Species_name.2.phylo,~ 1 | Record_id),
# #          R = list(Species_name.2.phylo = phylo_cor.subset3A),method = "REML",test = "t",data = dataset.IGE.subset3A_long
# #          ,control=list(optimizer="BBoptim")))
# # 
# # 
# # tab <- lapply(res, function(x) data.frame(optimizer=x$control$optimizer, ll=logLik(x),
# #                                           sigma2.1=x$sigma2[1],
# #                                           sigma2.2=x$sigma2[2],
# #                                           sigma2.2=x$sigma2[3],
# #                                           sigma2.2=x$sigma2[4],
# #                                           sigma2.2=x$sigma2[5],
# #                                           sigma2.2=x$sigma2[6]))
# # tab <- do.call(rbind, tab)
# # dfround(tab, digits=5)



################################################################################
# 3B: h2 vs social h2


#### Before we compare the two values in the same model, we can run an
# intercept-only model for h2 as we did for social h2 in section 1

# MAIN MODEL = multilevel meta-analysis (i.e. the intercept-only model)
meta.model.DGE.subset1A <- rma.mv(r.to.Zr(dataset.IGE.subset1A$H2_2),
                                  subset1A_VCV_ESVar, 
                                  mods = ~ 1,
                                  random = list(~ 1 | Paper_id,
                                                ~ 1 | Group_id,
                                                ~ 1 | Population2,
                                                ~ 1 | Species_name.2,
                                                ~ 1 | Species_name.2.phylo,
                                                ~ 1 | Record_id),
                                  R = list(Species_name.2.phylo = phylo_cor.subset1A),
                                  method = "REML", 
                                  test = "t", 
                                  data = dataset.IGE.subset1A)


# Printing the summary results of the model
print(meta.model.DGE.subset1A, digits=3)

# Printing the results again, but adding the 95% prediction interval, 
# which uses the heterogeneity to generate an interval that should contain 95%
# of the effect sizes of any future or unknown studies with similar features
# as those included in the current database
predict(meta.model.DGE.subset1A, digits=2)


# Estimating heterogeneity as I2 (Nakagawa and Santos 2012)
I2.model.1A.DGE <- i2_ml(meta.model.DGE.subset1A)
round(I2.model.1A.DGE, 1)


######

dataset.IGE$H2_2_Zr <- r.to.Zr(dataset.IGE$H2_2) 
par(mfrow=c(1,2))
hist(dataset.IGE$H2_2_Zr, breaks=100) 
hist(dataset.IGE$Social_h2_2_Zr, breaks=100)

# more comparable
par(mfrow=c(1,1))

# We need H2, social h2, and sampling variance
dataset.IGE.subset3B <- as.data.frame(dataset.IGE 
                                      %>% filter(!(is.na(Social_h2_2_Zr)) &
                                                   !(is.na(H2_2_Zr)) &
                                                   !(is.na(VZr))))


# some numbers and exploration
nrow(dataset.IGE.subset3B)
length(unique(dataset.IGE.subset3B$Paper_id))
length(unique(dataset.IGE.subset3B$Species_name.2))
summary(dataset.IGE.subset3B)
# 146 effect sizes, 40 papers, 19 species, so the same as subset 1A? YES
setdiff(dataset.IGE.subset1A$H2_2,dataset.IGE.subset3B$H2_2)
setdiff(dataset.IGE.subset3B$H2_2,dataset.IGE.subset1A$H2_2)
setdiff(dataset.IGE.subset1A$Social_h2_2,dataset.IGE.subset3B$Social_h2_2)
setdiff(dataset.IGE.subset3B$Social_h2_2,dataset.IGE.subset1A$Social_h2_2)


# convert to long format to compare direct and social heritability using a 
# moderator
dataset.IGE.subset3B_long <- data.frame(dataset.IGE.subset3B %>%
                                          pivot_longer(cols = c("H2_2_Zr", 
                                                                "Social_h2_2_Zr"), 
                                                       values_to = "dat3B_Zr", 
                                                       names_to = "direct_social"))

# given that the dataset is now "double" the size despite being the same data
# we need to account for this new level of nonindependence. The easiest way at
# this point is to continue using all random effects we used before, but add
# an additional one to account for the within-study/residual variance, given 
# that "Record_id" is no longer a unit-level random effect. Although we kept
# the same names, we need to interpret the random effects slightly different
# to what we did for the "non-long" models
dataset.IGE.subset3B_long$Record_id_long <- 1:nrow(dataset.IGE.subset3B_long)


subset3B_long_VCV_ESVar <- matrix(0, 
                                  nrow=nrow(dataset.IGE.subset3B_long), 
                                  ncol=nrow(dataset.IGE.subset3B_long))

# Names rows and columns for each Record_id_long
rownames(subset3B_long_VCV_ESVar) <- dataset.IGE.subset3B_long[,"Record_id_long"]
colnames(subset3B_long_VCV_ESVar) <- dataset.IGE.subset3B_long[,"Record_id_long"]

# Finds effect sizes that come from the same study
shared_coord.subset3B_long <- which(dataset.IGE.subset3B_long[,"Paper_id"] %in% 
                                      dataset.IGE.subset3B_long[
                                        duplicated(dataset.IGE.subset3B_long[,"Paper_id"]), 
                                        "Paper_id"]==TRUE)

combinations.subset3B_long <- do.call("rbind", tapply(shared_coord.subset3B_long, 
                                                      dataset.IGE.subset3B_long[
                                                        shared_coord.subset3B_long, 
                                                        "Paper_id"], 
                                                      function(x) t(utils::combn(x, 
                                                                                 2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1:dim(combinations.subset3B_long)[1]) {
  p1 <- combinations.subset3B_long[i, 1]
  p2 <- combinations.subset3B_long[i, 2]
  p1_p2_cov <- 0.5*
    sqrt(dataset.IGE.subset3B_long[p1, "VZr"])*
    sqrt(dataset.IGE.subset3B_long[p2, "VZr"])
  subset3B_long_VCV_ESVar[p1, p2] <- p1_p2_cov
  subset3B_long_VCV_ESVar[p2, p1] <- p1_p2_cov
} 


# Enters previously calculated effect size sampling variances into diagonals 
diag(subset3B_long_VCV_ESVar) <- dataset.IGE.subset3B_long[,"VZr"]


# As the species are identical I assume we don't need to alter the phylogeny?
# Still need to create a copy of Species_name.2 for phylogenetic effect
dataset.IGE.subset3B_long$Species_name.2.phylo <- dataset.IGE.subset3B_long$Species_name.2

# saving the subset for script 006_figures.R
write.csv(dataset.IGE.subset3B_long, file = "data/subsets/dataset_IGE_subset3B_long.csv",
          row.names = F)

meta.model.IGE.subset3B <- rma.mv(dat3B_Zr,
                                  subset3B_long_VCV_ESVar,
                                  mods = ~ direct_social,
                                  random = list(~ 1 | Paper_id,
                                                ~ 1 | Group_id,
                                                ~ 1 | Population2,
                                                ~ 1 | Species_name.2,
                                                ~ 1 | Species_name.2.phylo,
                                                ~ 1 | Record_id_long,
                                                ~ 1 | Record_id),
                                  R = list(Species_name.2.phylo = phylo_cor.subset1A),
                                  method = "REML", 
                                  test = "t", 
                                  data = dataset.IGE.subset3B_long)

# saving the model for script 006_figures.R
save(meta.model.IGE.subset3B, 
     file = "data/models/IGEmeta_regression_h2_vs_socialh2.Rdata")

# model can be loaded instead of run using the following
# load("data/models/IGEmeta_regression_h2_vs_socialh2.Rdata")

# Printing the summary results of the model
print(meta.model.IGE.subset3B, digits=3)


# Calculate marginal R2 with r2_ml
R2.IGE.subset3B <- r2_ml(meta.model.IGE.subset3B) 
round(R2.IGE.subset3B*100, 1)


# Run without intercept
meta.model.IGE.subset3B.NI <- rma.mv(dat3B_Zr,
                                     subset3B_long_VCV_ESVar, 
                                     mods = ~ direct_social-1,
                                     random = list(~ 1 | Paper_id,
                                                   ~ 1 | Group_id,
                                                   ~ 1 | Population2,
                                                   ~ 1 | Species_name.2,
                                                   ~ 1 | Species_name.2.phylo,
                                                   ~ 1 | Record_id_long,
                                                   ~ 1 | Record_id),
                                     R = list(Species_name.2.phylo = phylo_cor.subset1A),
                                     method = "REML", 
                                     test = "t", 
                                     data = dataset.IGE.subset3B_long)

# Printing the summary results of the model
print(meta.model.IGE.subset3B.NI, digits=3)


################################################################################
# 3C: Ia vs Ii

# For evolvability #

#     - Ia, where Ia = Va/(Trait_mean^2)
#     - Iiges, where Iiges = V_ige/(Trait_mean^2)

# some comments to keep in mind, to calculate evolvability
# raw variance divided by the mean^2 whenever Mean_standardized==no, and if yes
# then the variance as reported

nrow(dataset.IGE[dataset.IGE$Variance_standardized == "no" & 
                   !(is.na(dataset.IGE$Trait_mean2)),])

#105 rows where trait is reported and variance was not standardised 

dataset.IGE$I_a <- ifelse(dataset.IGE$Variance_standardized == "no" & 
                            dataset.IGE$Trait_mean2 > 0.001 & # this excludes a very small value (0.00074) and 3 negative means that lead to a very skewed distribution otherwise. Whether to include these four or not lead to similar qualitative results, but the estimates seem more reliable when these four values are not included
                            !(is.na(dataset.IGE$Trait_mean2)),
                          dataset.IGE$Va / (dataset.IGE$Trait_mean2^2 ),
                          NA )


dataset.IGE$I_ige = ifelse(dataset.IGE$Variance_standardized == "no" & 
                             dataset.IGE$Trait_mean2 > 0.001 & # this excludes a very small value (0.00074) and 3 negative means that lead to a very skewed distribution otherwise. Whether to include these four or not lead to similar qualitative results, but the estimates seem more reliable when these four values are not included
                            !(is.na(dataset.IGE$Trait_mean2)),
                           dataset.IGE$V_ige / (dataset.IGE$Trait_mean2^2 ),
                           NA )

par(mfrow=c(1,2))
hist(dataset.IGE$I_a,breaks=40)
# Large values are a fish paper but they don't seem wacky
# just small means, so not removing them
hist(dataset.IGE$I_ige,breaks=40)
#range from 0-0.2

dataset.IGE.subset3C <- as.data.frame(dataset.IGE 
                                      %>% filter(!(is.na(I_a)) &
                                                   !(is.na(I_ige)) &
                                                   !(is.na(VZr))))

# some numbers and exploration
nrow(dataset.IGE.subset3C)
length(unique(dataset.IGE.subset3C$Paper_id))
length(unique(dataset.IGE.subset3C$Species_name.2))
summary(dataset.IGE.subset3C)
# 97 effect sizes, 30 papers, 16 species

dataset.IGE.subset3C_long = dataset.IGE.subset3C %>%
  pivot_longer(cols = c("I_a", "I_ige"), 
               values_to = "dat3C_I", 
               names_to = "direct_social") %>%
  mutate(dat_weights = log(N_id_w_records))%>% # for this analysis, we considered using the log of sample size as weights for the effect sizes but decided that using VZr, which is 1/(N-3) makes the results of all analyses more comparable to each other given that log(sample size) and 1/(N-3) are non-linearly associated
  filter(is.finite(dat3C_I), 
         is.finite(dat_weights)) %>% as.data.frame()


# given that the dataset is now "double" the size despite being the same data
# we need to account for this new level of nonindependence. The easiest way at
# this point is to continue using all random effects we used before, but add
# an additional one to account for the within-study/residual variance, given 
# that "Record_id" is no longer a unit-level random effect. Although we kept
# the same names, we need to interpret the random effects slightly different
# to what we did for the "non-long" models
dataset.IGE.subset3C_long$Record_id_long <- 1:nrow(dataset.IGE.subset3C_long)

nrow(dataset.IGE.subset3C_long) # lost no values, double 97 = 194!
length(unique(dataset.IGE.subset3C_long$Paper_id))
length(unique(dataset.IGE.subset3C_long$Species_name.2))

subset3C_long_VCV_ESVar <- matrix(0, 
                                  nrow=nrow(dataset.IGE.subset3C_long), 
                                  ncol=nrow(dataset.IGE.subset3C_long))

# Names rows and columns for each Record_id_long
rownames(subset3C_long_VCV_ESVar) <- dataset.IGE.subset3C_long[,"Record_id_long"]
colnames(subset3C_long_VCV_ESVar) <- dataset.IGE.subset3C_long[,"Record_id_long"]

# Finds effect sizes that come from the same study
shared_coord.subset3C_long <- which(dataset.IGE.subset3C_long[,"Paper_id"] %in% 
                                      dataset.IGE.subset3C_long[
                                        duplicated(dataset.IGE.subset3C_long[,"Paper_id"]), 
                                        "Paper_id"]==TRUE)

combinations.subset3C_long <- do.call("rbind", tapply(shared_coord.subset3C_long, 
                                                      dataset.IGE.subset3C_long[
                                                        shared_coord.subset3C_long, 
                                                        "Paper_id"], 
                                                      function(x) t(utils::combn(x, 
                                                                                 2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1:dim(combinations.subset3C_long)[1]) {
  p1 <- combinations.subset3C_long[i, 1]
  p2 <- combinations.subset3C_long[i, 2]
  p1_p2_cov <- 0.5*
    sqrt(dataset.IGE.subset3C_long[p1, "VZr"])*
    sqrt(dataset.IGE.subset3C_long[p2, "VZr"])
  subset3C_long_VCV_ESVar[p1, p2] <- p1_p2_cov
  subset3C_long_VCV_ESVar[p2, p1] <- p1_p2_cov
} 

# Enters previously calculated effect size sampling variances into diagonals 
diag(subset3C_long_VCV_ESVar) <- dataset.IGE.subset3C_long[,"VZr"]

# creating a copy of Species_name.2 for phylogenetic effect
dataset.IGE.subset3C_long$Species_name.2.phylo <- dataset.IGE.subset3C_long$Species_name.2

# saving the subset for script 006_figures.R
write.csv(dataset.IGE.subset3C_long, file = "data/subsets/dataset_IGE_subset3C_long.csv",
          row.names = F)

phylo_cor.subset3C <- phylo_cor[rownames(phylo_cor) %in% 
                                  unique(as.character(dataset.IGE.subset3C_long$Species_name.2)),
                                colnames(phylo_cor) %in% 
                                  unique(as.character(dataset.IGE.subset3C_long$Species_name.2))]

meta.model.IGE.subset3C <- rma.mv(dat3C_I,
                                  subset3C_long_VCV_ESVar,
                                  mods = ~ direct_social,
                                  random = list(~ 1 | Paper_id,
                                                ~ 1 | Group_id,
                                                ~ 1 | Population2,
                                                ~ 1 | Species_name.2,
                                                ~ 1 | Species_name.2.phylo,
                                                ~ 1 | Record_id_long,
                                                ~ 1 | Record_id),
                                  R = list(Species_name.2.phylo = phylo_cor.subset3C),
                                  method = "REML", 
                                  test = "t", 
                                  data = dataset.IGE.subset3C_long,
                                  control=list(optimizer="nlminb", rel.tol=1e-8)) # to allow true convergence, details explained in https://www.metafor-project.org/doku.php/tips:convergence_problems_rma_mv

# saving the model for script 006_figures.R
save(meta.model.IGE.subset3C, 
     file = "data/models/IGEmeta_regression_Ia_vs_Iige.Rdata")

# model can be loaded instead of run using the following
# load("data/models/IGEmeta_regression_Ia_vs_Iige.Rdata")

# Printing the summary results of the model
print(meta.model.IGE.subset3C, digits=3)


# Calculate marginal R2 with r2_ml
R2.IGE.subset3C <- r2_ml(meta.model.IGE.subset3C) 
round(R2.IGE.subset3C*100, 1)


# Run without intercept
meta.model.IGE.subset3C.NI <- rma.mv(dat3C_I,
                                     subset3C_long_VCV_ESVar,
                                     mods = ~ direct_social-1,
                                     random = list(~ 1 | Paper_id,
                                                   ~ 1 | Group_id,
                                                   ~ 1 | Population2,
                                                   ~ 1 | Species_name.2,
                                                   ~ 1 | Species_name.2.phylo,
                                                   ~ 1 | Record_id_long,
                                                   ~ 1 | Record_id),
                                     R = list(Species_name.2.phylo = phylo_cor.subset3C),
                                     method = "REML", 
                                     test = "t", 
                                     data = dataset.IGE.subset3C_long,
                                     control=list(optimizer="nlminb", rel.tol=1e-8)) 

# Printing the summary results of the model
print(meta.model.IGE.subset3C.NI, digits=3)



################################################################################
# Q4: Do IGEs typically alter evolutionary trajectories?

# We  get at this by calculating T2 which includes IGEs and their covariance 
# with DGEs (and group size) and comparing to H2, see if it is larger or smaller
# We also look at the DGE-IGE correlation and see the overall direction

################################################################################
# 4A: What is the magnitude of the DGE-IGE correlation?

hist(dataset.IGE$R_a_ige) #quite even either side of 0, with additional peaks at -1 and 1
sum(!(is.na(dataset.IGE$R_a_ige))) # 115 effect sizes in total

dataset.IGE$R_a_ige_calc <- dataset.IGE$Cov_a_ige / sqrt(dataset.IGE$Va*dataset.IGE$V_ige_2)
# there were 10 V_ige values = 0 that ultimately lead to 9 Inf estimates that
# would have been excluded from the analyses. As explained above, we decided
# to transform those 0's to minimal values, and we are using them throughout.

with(dataset.IGE, plot(R_a_ige,R_a_ige_calc))
abline(0,1) 

dataset.IGE$R_a_ige_2 <- ifelse(is.na(dataset.IGE$R_a_ige), 
                                dataset.IGE$R_a_ige_calc,
                                dataset.IGE$R_a_ige)

hist(dataset.IGE$R_a_ige_2)
# slightly flatter around 0 but otherwise similar dist. Crucially none beyond 
# -1 and 1 

sum(!(is.na(dataset.IGE$R_a_ige_2)))
# 148 effect sizes, so we managed to add 33, mostly just around zero

dataset.IGE$R_a_ige_2_Zr <- r.to.Zr(dataset.IGE$R_a_ige_2) 

dataset.IGE.subset4B <- as.data.frame(dataset.IGE 
                                      %>% filter(!(is.na(R_a_ige_2_Zr)) & 
                                                   !(is.na(VZr))))

# some numbers and exploration
nrow(dataset.IGE.subset4B)
length(unique(dataset.IGE.subset4B$Paper_id))
length(unique(dataset.IGE.subset4B$Species_name.2))
summary(dataset.IGE.subset4B)
# 146 effect sizes, 40 papers, 17 species

subset4B_VCV_ESVar <- matrix(0, 
                             nrow=nrow(dataset.IGE.subset4B), 
                             ncol=nrow(dataset.IGE.subset4B))

# Names rows and columns for each Record_id
rownames(subset4B_VCV_ESVar) <- dataset.IGE.subset4B[,"Record_id"]
colnames(subset4B_VCV_ESVar) <- dataset.IGE.subset4B[,"Record_id"]

# Finds effect sizes that come from the same study
shared_coord.subset4B <- which(dataset.IGE.subset4B[,"Paper_id"] %in% 
                                 dataset.IGE.subset4B[
                                   duplicated(dataset.IGE.subset4B[,"Paper_id"]), 
                                   "Paper_id"]==TRUE)

combinations.subset4B <- do.call("rbind", tapply(shared_coord.subset4B, 
                                                 dataset.IGE.subset4B[
                                                   shared_coord.subset4B, 
                                                   "Paper_id"], 
                                                 function(x) t(utils::combn(x, 
                                                                            2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1:dim(combinations.subset4B)[1]) {
  p1 <- combinations.subset4B[i, 1]
  p2 <- combinations.subset4B[i, 2]
  p1_p2_cov <- 0.5*
    sqrt(dataset.IGE.subset4B[p1, "VZr"])*
    sqrt(dataset.IGE.subset4B[p2, "VZr"])
  subset4B_VCV_ESVar[p1, p2] <- p1_p2_cov
  subset4B_VCV_ESVar[p2, p1] <- p1_p2_cov
} 

# Enters previously calculated effect size sampling variances into diagonals 
diag(subset4B_VCV_ESVar) <- dataset.IGE.subset4B[,"VZr"]

# # Raw data funnel plot 
# par(mfrow=c(1, 2))
# funnel(dataset.IGE.subset4B$Social_h2_2_Zr, 
#        dataset.IGE.subset4B$VZr,
#        yaxis="sei",
#        xlab="Effect size (Zr)", 
#        digits=2, las=1) 
# funnel(dataset.IGE.subset4B$Social_h2_2_Zr, 
#        dataset.IGE.subset4B$VZr, 
#        yaxis="seinv",
#        xlab="Effect size (Zr)",  
#        digits=2, las=1) 
# # Looks OK

# creating a copy of Species_name.2 for phylogenetic effect
dataset.IGE.subset4B$Species_name.2.phylo <- dataset.IGE.subset4B$Species_name.2

# saving the subset for script 006_figures.R
write.csv(dataset.IGE.subset4B, file = "data/subsets/dataset_IGE_subset4B.csv",
          row.names = F)

phylo_cor.subset4B <- phylo_cor[rownames(phylo_cor) %in% 
                                  unique(as.character(dataset.IGE.subset4B$Species_name.2)),
                                colnames(phylo_cor) %in% 
                                  unique(as.character(dataset.IGE.subset4B$Species_name.2))]

meta.model.IGE.subset4B <- rma.mv(R_a_ige_2_Zr,
                                  subset4B_VCV_ESVar, 
                                  mods = ~ 1,
                                  random = list(~ 1 | Paper_id,
                                                ~ 1 | Group_id,
                                                ~ 1 | Population2,
                                                ~ 1 | Species_name.2,
                                                ~ 1 | Species_name.2.phylo,
                                                ~ 1 | Record_id),
                                  R = list(Species_name.2.phylo = phylo_cor.subset4B),
                                  method = "REML", 
                                  test = "t", 
                                  data = dataset.IGE.subset4B)

# saving the model for script 006_figures.R
save(meta.model.IGE.subset4B, file = "data/models/meta_model_IGE_DGE-IGE_correlation.Rdata")

# model can be loaded instead of run using the following
# load("data/models/meta_model_IGE_DGE-IGE_correlation.Rdata")

# Printing the summary results of the model
print(meta.model.IGE.subset4B, digits=3)


# Printing the results again, but adding the credibility/prediction interval, 
# which uses the heterogeneity to generate an interval that should contain 95%
# of the effect sizes of any future or unknown studies with similar features
# as those included in the current database
predict(meta.model.IGE.subset4B, digits=3)

# # Model funnel plots
# par(mfrow = c(1, 1))
# funnel(meta.model.IGE.subset4B) # lovely spread


# # We can look for potential outliers based on the residuals of the model
# # any points below -2 or above 2 could be potential outliers
# resid<-rstandard(meta.model.IGE.subset4B)
# plot(1:nrow(dataset.IGE.subset4B), resid$z, type="b")
# abline(h = 2) #limited issues
# 
# dataset.IGE.subset4B[resid$z > 2,]
# these are all very strong positive corrs (0.98+), 3 reported, 2 calculated


# using pluralistic approach to explore heterogeneity (Yang et al. 2023)

# why is the following not the same as above? Because the following is the 
# unstandardized raw variance whereas the previous is the sampling variance
round(sum(meta.model.IGE.subset4B$sigma2),4)
round(i2_ml(meta.model.IGE.subset4B),4)
round(cv_ml(meta.model.IGE.subset4B),4)
round(m1_ml(meta.model.IGE.subset4B),4)

# Output as table
ige_4b_res <- hux(predict(meta.model.IGE.subset4B, digits=3)) %>%
  add_rows(c("Prediction", "Standard error", 
             "CI lower", "CI upper", "
             PI lower", "PI upper"), 
           after = 0)
ige_4b_res


################################################################################
# 4B: Total heritability (T2) vs narrow-sense heritability (h2)

sum(!(is.na(dataset.IGE$V_tbv))) #97 effect sizes presented

hist(dataset.IGE$V_tbv) # wild range
hist(dataset.IGE$V_tbv/dataset.IGE$Total_v_phen3) # this is T2, its better, ranges from 0 to 5

# Calculate Vtbv:
dataset.IGE$V_tbv_calc <- dataset.IGE$Va + 
  dataset.IGE$Cov_a_ige*2*(dataset.IGE$Mean_group_size-1) + 
  dataset.IGE$V_ige*((dataset.IGE$Mean_group_size-1)^2)

# Calculate Vtbv but using using the transform V_ige 0 values (see above):
# substituting all those values for their corresponding minimum value
dataset.IGE$V_ige_2 <- ifelse((dataset.IGE$Paper_id %in% 
                                 c("IGE0501","IGE0468","IGE0406")) & dataset.IGE$V_ige==0,
                              0.0001,
                              ifelse((dataset.IGE$Paper_id %in% 
                                        c("IGE0857")) & dataset.IGE$V_ige==0,
                                     0.001,ifelse((dataset.IGE$Paper_id %in% 
                                                     c("IGE0203")) & dataset.IGE$V_ige==0,
                                                  0.1,
                                                  dataset.IGE$V_ige)
                              ))


dataset.IGE$V_tbv_calc_2 <- dataset.IGE$Va + 
  dataset.IGE$Cov_a_ige*2*(dataset.IGE$Mean_group_size-1) + 
  dataset.IGE$V_ige_2*((dataset.IGE$Mean_group_size-1)^2)

par(mfrow=c(1,2))

with(dataset.IGE, plot(V_tbv, V_tbv_calc, xlim=c(0,8000),ylim=c(0,8000)))
abline(a = 0, b = 1) # pretty good correlation -> reliable

with(dataset.IGE, plot(V_tbv, V_tbv_calc_2, xlim=c(0,8000),ylim=c(0,8000)))
abline(a = 0, b = 1) # pretty good correlation -> reliable


# from now on, we proceed using V_tbv_calc_2, which includes the extra values
# reported as 0 and transformed to a minimum value by us
dataset.IGE$V_tbv_2 <- ifelse(is.na(dataset.IGE$V_tbv), 
                              dataset.IGE$V_tbv_calc_2,
                              dataset.IGE$V_tbv)

sum(!(is.na(dataset.IGE$V_tbv_2))) # 116 effect sizes now, so we manage to add 19 additional effect sizes

hist(dataset.IGE$V_tbv_2)
hist(dataset.IGE$V_tbv_2/dataset.IGE$Total_v_phen3) # better, ranges from 0 to 5

max(dataset.IGE$T2, na.rm = T) #2.49

sum(!(is.na(dataset.IGE$T2))) # 103 effect sizes

nrow(dataset.IGE[complete.cases(dataset.IGE$T2,
                                dataset.IGE$V_tbv),]) #85 rows with both

dataset.IGE$T2_calc <- dataset.IGE$V_tbv_2 / dataset.IGE$Total_v_phen3 

with(dataset.IGE, plot(T2, T2_calc))
abline(0,1) # rather good correlation, only some more var between 0.1 and 0.5 -> reliable

dataset.IGE$T2_2 <- ifelse(is.na(dataset.IGE$T2), 
                           dataset.IGE$T2_calc,
                           dataset.IGE$T2)

sum(!(is.na(dataset.IGE$T2_2))) # 116 effect sizes now, so we added 13 additional effect sizes
max(dataset.IGE$T2_2, na.rm = T) # some large values: 2.49, 1.35, 2.22, 2.31, 1.29, 1.37

with(dataset.IGE, plot(T2_2, H2_2, xlim = c(0,2.5), ylim = c(0,2.5)))
abline(0,1)
# T2 typically higher than h2, a lot of variation around line 

dataset.IGE.subset4A <- as.data.frame(dataset.IGE 
                                      %>% filter(!(is.na(T2_2)) &
                                                   !(is.na(H2_2)) &
                                                   !(is.na(VZr))))


# some numbers and exploration
nrow(dataset.IGE.subset4A)
length(unique(dataset.IGE.subset4A$Paper_id))
length(unique(dataset.IGE.subset4A$Species_name.2))
summary(dataset.IGE.subset4A)
# 110 effect sizes, 34 papers, 15 species

hist(dataset.IGE.subset4A$Mean_group_size,breaks=50)
summary(dataset.IGE.subset4A$Mean_group_size)


dataset.IGE.subset4A_long <- dataset.IGE.subset4A %>%
  pivot_longer(cols = c("H2_2", "T2_2"), 
               values_to = "dat4A_herit", 
               names_to = "direct_social") %>%
  mutate(dat_weights = log(N_id_w_records)) %>% # for this analysis, we considered using the log of sample size as weights for the effect sizes but decided that using VZr, which is 1/(N-3) makes the results of all analyses more comparable to each other given that log(sample size) and 1/(N-3) are non-linearly associated
  filter(is.finite(dat_weights)) %>% 
  as.data.frame()

hist(dataset.IGE.subset4A_long$dat4A_herit, breaks = 100)

# given that the dataset is now "double" the size despite being the same data
# we need to account for this new level of nonindependence. The easiest way at
# this point is to continue using all random effects we used before, but add
# an additional one to account for the within-study/residual variance, given 
# that "Record_id" is no longer a unit-level random effect. Although we kept
# the same names, we need to interpret the random effects slightly different
# to what we did for the "non-long" models
dataset.IGE.subset4A_long$Record_id_long <- 1:nrow(dataset.IGE.subset4A_long)


nrow(dataset.IGE.subset4A_long) #double 110
length(unique(dataset.IGE.subset4A_long$Paper_id))
length(unique(dataset.IGE.subset4A_long$Species_name.2))
# 220 effect sizes, 34 papers, 15 species

subset4A_long_VCV_ESVar <- matrix(0, nrow=nrow(dataset.IGE.subset4A_long), 
                                  ncol=nrow(dataset.IGE.subset4A_long))

# Names rows and columns for each Record_id_long
rownames(subset4A_long_VCV_ESVar) <- dataset.IGE.subset4A_long[,"Record_id_long"]
colnames(subset4A_long_VCV_ESVar) <- dataset.IGE.subset4A_long[,"Record_id_long"]

# Finds effect sizes that come from the same study
shared_coord.subset4A_long <- which(dataset.IGE.subset4A_long[,"Paper_id"] %in% 
                                      dataset.IGE.subset4A_long[
                                        duplicated(dataset.IGE.subset4A_long[,"Paper_id"]), 
                                        "Paper_id"]==TRUE)

combinations.subset4A_long <- do.call("rbind", tapply(shared_coord.subset4A_long, 
                                                      dataset.IGE.subset4A_long[
                                                        shared_coord.subset4A_long, 
                                                        "Paper_id"], 
                                                      function(x) t(utils::combn(x, 
                                                                                 2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1:dim(combinations.subset4A_long)[1]) {
  p1 <- combinations.subset4A_long[i, 1]
  p2 <- combinations.subset4A_long[i, 2]
  p1_p2_cov <- 0.5*
    sqrt(dataset.IGE.subset4A_long[p1, "VZr"])*
    sqrt(dataset.IGE.subset4A_long[p2, "VZr"])
  subset4A_long_VCV_ESVar[p1, p2] <- p1_p2_cov
  subset4A_long_VCV_ESVar[p2, p1] <- p1_p2_cov
} 

# Enters previously calculated effect size sampling variances into diagonals 
diag(subset4A_long_VCV_ESVar) <- dataset.IGE.subset4A_long[,"VZr"]

# creating a copy of Species_name.2 for phylogenetic effect
dataset.IGE.subset4A_long$Species_name.2.phylo <- dataset.IGE.subset4A_long$Species_name.2

# saving the subset for script 006_figures.R
write.csv(dataset.IGE.subset4A_long, file = "data/subsets/dataset_IGE_subset4A_long.csv",
          row.names = F)

phylo_cor.subset4A <- phylo_cor[rownames(phylo_cor) %in% 
                                  unique(as.character(dataset.IGE.subset4A_long$Species_name.2)),
                                colnames(phylo_cor) %in% 
                                  unique(as.character(dataset.IGE.subset4A_long$Species_name.2))]

meta.model.IGE.subset4A <- rma.mv(dat4A_herit,
                                  subset4A_long_VCV_ESVar,
                                  mods = ~ direct_social,
                                  random = list(~ 1 | Paper_id,
                                                ~ 1 | Group_id,
                                                ~ 1 | Population2,
                                                ~ 1 | Species_name.2,
                                                ~ 1 | Species_name.2.phylo,
                                                ~ 1 | Record_id_long,
                                                ~ 1 | Record_id),
                                  R = list(Species_name.2.phylo = phylo_cor.subset4A),
                                  method = "REML", 
                                  test = "t", 
                                  data = dataset.IGE.subset4A_long)

# saving the model for script 006_figures.R
save(meta.model.IGE.subset4A, 
     file = "data/models/IGEmeta_regression_h2_vs_Totalh2.Rdata")

# model can be loaded instead of run using the following
# load("data/models/IGEmeta_regression_h2_vs_Totalh2.Rdata")

print(meta.model.IGE.subset4A, digits=3)


# Calculate marginal R2 with r2_ml
R2.IGE.subset4A <- r2_ml(meta.model.IGE.subset4A) 
round(R2.IGE.subset4A*100, 1)


# Run without intercept
meta.model.IGE.subset4A.NI <- rma.mv(dat4A_herit,
                                     subset4A_long_VCV_ESVar,
                                     mods = ~ direct_social-1,
                                     random = list(~ 1 | Paper_id,
                                                   ~ 1 | Group_id,
                                                   ~ 1 | Population2,
                                                   ~ 1 | Species_name.2,
                                                   ~ 1 | Species_name.2.phylo,
                                                   ~ 1 | Record_id_long,
                                                   ~ 1 | Record_id),
                                     R = list(Species_name.2.phylo = phylo_cor.subset4A),
                                     method = "REML", 
                                     test = "t", 
                                     data = dataset.IGE.subset4A_long)

# Printing the summary results of the model
print(meta.model.IGE.subset4A.NI, digits=3) 



################################################################################
# 4Bii - what contributes more to the large increase of T2, covariance or IGE variance? 
# T2 is Va + 2(n-1)Cov + (n-1)^2*Vige

with(dataset.IGE.subset4A, plot(Cov_a_ige*2*(Mean_group_size-1),
                                V_ige_2*((Mean_group_size-1)^2)))
abline(a = 0, b = 1) #above the line means the 2nd term is larger
# pretty much all above the line

# Limit to see most of data better
with(dataset.IGE.subset4A, plot(Cov_a_ige*2*(Mean_group_size-1),
                                V_ige_2*((Mean_group_size-1)^2), 
                                xlim = c(-10, 1500),
                                ylim=c(0, 4000)))
abline(a = 0, b = 1)
# few below the line

# Limit further to see most of data better
with(dataset.IGE.subset4A, plot(Cov_a_ige*2*(Mean_group_size-1),
                                V_ige_2*((Mean_group_size-1)^2), 
                                xlim = c(-10, 15),
                                ylim=c(0, 10)))
abline(a = 0, b = 1)
# few below the line

summary(dataset.IGE.subset4A$Cov_a_ige)
summary((dataset.IGE.subset4A$V_ige_2*((dataset.IGE.subset4A$Mean_group_size-1)^2)))
correlation=(dataset.IGE.subset4A$Cov_a_ige*2*(dataset.IGE.subset4A$Mean_group_size-1))
summary(correlation)
group.size=(dataset.IGE.subset4A$V_ige_2*((dataset.IGE.subset4A$Mean_group_size-1)^2))
summary(group.size)
aa= group.size>correlation
median(group.size-correlation, na.rm=TRUE)
median(group.size, na.rm=TRUE)-mean(correlation, na.rm=TRUE)

table(aa)
(57*100)/78


################################################################################
# 4Biii: "response to direct selection and unrelated individuals" vs narrow-sense heritability (h2)

#we want to calculate Va + [n-1]*Cov_DI
dataset.IGE$RDSUR=dataset.IGE$Va + ((dataset.IGE$Mean_group_size-1)*dataset.IGE$Cov_a_ige)
summary(dataset.IGE$RDSUR)

sum(!(is.na(dataset.IGE$RDSUR))) #92 effect sizes presented

hist(dataset.IGE$RDSUR) # wide range
hist(dataset.IGE$RDSUR/dataset.IGE$Total_v_phen3) # this is T2, its better, ranges from 0 to 5

with(dataset.IGE, plot(V_tbv, RDSUR, xlim=c(0,8000),ylim=c(0,8000)))
abline(a = 0, b = 1) # pretty good correlation

dataset.IGE$T2_new <- dataset.IGE$RDSUR / dataset.IGE$Total_v_phen3 

with(dataset.IGE, plot(T2, T2_new))
abline(0,1) # rather good correlation, only some more var between 0.1 and 0.5 -> reliable

sum(!(is.na(dataset.IGE$T2_new))) # 79 effect sizes now, so we added 13 additional effect sizes
max(dataset.IGE$T2_new, na.rm = T) # some large values: 2.49, 1.35, 2.22, 2.31, 1.29, 1.37

with(dataset.IGE, plot(T2_new, H2_2, xlim = c(0,2.5), ylim = c(0,2.5)))
abline(0,1)
# T2 typically higher than h2, a lot of variation around line 

dataset.IGE.subset4ARD <- as.data.frame(dataset.IGE 
                                        %>% filter(!(is.na(T2_new)) &
                                                     !(is.na(H2_2)) &
                                                     !(is.na(VZr))))

# some numbers and exploration
nrow(dataset.IGE.subset4ARD)
length(unique(dataset.IGE.subset4ARD$Paper_id))
length(unique(dataset.IGE.subset4ARD$Species_name.2))
summary(dataset.IGE.subset4ARD)
# 78 effect sizes, 24 papers, 12 species

hist(dataset.IGE.subset4ARD$Mean_group_size,breaks=50)
summary(dataset.IGE.subset4ARD$Mean_group_size)

dataset.IGE.subset4ARD_long <- dataset.IGE.subset4ARD %>%
  pivot_longer(cols = c("H2_2", "T2_new"), 
               values_to = "dat4A_herit", 
               names_to = "direct_social") %>%
  mutate(dat_weights = log(N_id_w_records)) %>% # for this analysis, we considered using the log of sample size as weights for the effect sizes but decided that using VZr, which is 1/(N-3) makes the results of all analyses more comparable to each other given that log(sample size) and 1/(N-3) are non-linearly associated
  filter(is.finite(dat_weights)) %>% 
  as.data.frame()

hist(dataset.IGE.subset4ARD_long$dat4A_herit, breaks = 100)

# given that the dataset is now "double" the size despite being the same data
# we need to account for this new level of nonindependence. The easiest way at
# this point is to continue using all random effects we used before, but add
# an additional one to account for the within-study/residual variance, given 
# that "Record_id" is no longer a unit-level random effect. Although we kept
# the same names, we need to interpret the random effects slightly different
# to what we did for the "non-long" models
dataset.IGE.subset4ARD_long$Record_id_long <- 1:nrow(dataset.IGE.subset4ARD_long)

nrow(dataset.IGE.subset4ARD_long) #double 156
length(unique(dataset.IGE.subset4ARD_long$Paper_id))
length(unique(dataset.IGE.subset4ARD_long$Species_name.2))
# 156 effect sizes, 24 papers, 12 species

subset4ARD_long_VCV_ESVar <- matrix(0, nrow=nrow(dataset.IGE.subset4ARD_long), 
                                    ncol=nrow(dataset.IGE.subset4ARD_long))

# Names rows and columns for each Record_id_long
rownames(subset4ARD_long_VCV_ESVar) <- dataset.IGE.subset4ARD_long[,"Record_id_long"]
colnames(subset4ARD_long_VCV_ESVar) <- dataset.IGE.subset4ARD_long[,"Record_id_long"]

# Finds effect sizes that come from the same study
shared_coord.subset4ARD_long <- which(dataset.IGE.subset4ARD_long[,"Paper_id"] %in% 
                                        dataset.IGE.subset4ARD_long[
                                          duplicated(dataset.IGE.subset4ARD_long[,"Paper_id"]), 
                                          "Paper_id"]==TRUE)

combinations.subset4ARD_long <- do.call("rbind", tapply(shared_coord.subset4ARD_long, 
                                                        dataset.IGE.subset4ARD_long[
                                                          shared_coord.subset4ARD_long, 
                                                          "Paper_id"], 
                                                        function(x) t(utils::combn(x, 
                                                                                   2))))

# Calculates the covariance between effect sizes and enters them in each 
# combination of coordinates
for (i in 1:dim(combinations.subset4ARD_long)[1]) {
  p1 <- combinations.subset4ARD_long[i, 1]
  p2 <- combinations.subset4ARD_long[i, 2]
  p1_p2_cov <- 0.5*
    sqrt(dataset.IGE.subset4ARD_long[p1, "VZr"])*
    sqrt(dataset.IGE.subset4ARD_long[p2, "VZr"])
  subset4ARD_long_VCV_ESVar[p1, p2] <- p1_p2_cov
  subset4ARD_long_VCV_ESVar[p2, p1] <- p1_p2_cov
} 

# Enters previously calculated effect size sampling variances into diagonals 
diag(subset4ARD_long_VCV_ESVar) <- dataset.IGE.subset4ARD_long[,"VZr"]

# creating a copy of Species_name.2 for phylogenetic effect
dataset.IGE.subset4ARD_long$Species_name.2.phylo <- dataset.IGE.subset4ARD_long$Species_name.2

phylo_cor.subset4ARD <- phylo_cor[rownames(phylo_cor) %in% 
                                    unique(as.character(dataset.IGE.subset4ARD_long$Species_name.2)),
                                  colnames(phylo_cor) %in% 
                                    unique(as.character(dataset.IGE.subset4ARD_long$Species_name.2))]

meta.model.IGE.subset4ARD <- rma.mv(dat4A_herit,
                                    subset4ARD_long_VCV_ESVar,
                                    mods = ~ direct_social,
                                    random = list(~ 1 | Paper_id,
                                                  ~ 1 | Group_id,
                                                  ~ 1 | Population2,
                                                  ~ 1 | Species_name.2,
                                                  ~ 1 | Species_name.2.phylo,
                                                  ~ 1 | Record_id_long,
                                                  ~ 1 | Record_id),
                                    R = list(Species_name.2.phylo = phylo_cor.subset4ARD),
                                    method = "REML", 
                                    test = "t", 
                                    data = dataset.IGE.subset4ARD_long)

print(meta.model.IGE.subset4ARD, digits=3)

# Calculate marginal R2 with r2_ml
R2.IGE.subset4AARD <- r2_ml(meta.model.IGE.subset4ARD) 
round(R2.IGE.subset4AARD*100, 1)


# Run without intercept
meta.model.IGE.subset4ARD.NI <- rma.mv(dat4A_herit,
                                       subset4ARD_long_VCV_ESVar,
                                       mods = ~ direct_social-1,
                                       random = list(~ 1 | Paper_id,
                                                     ~ 1 | Group_id,
                                                     ~ 1 | Population2,
                                                     ~ 1 | Species_name.2,
                                                     ~ 1 | Species_name.2.phylo,
                                                     ~ 1 | Record_id_long,
                                                     ~ 1 | Record_id),
                                       R = list(Species_name.2.phylo = phylo_cor.subset4ARD),
                                       method = "REML", 
                                       test = "t", 
                                       data = dataset.IGE.subset4ARD_long)

# Printing the summary results of the model
print(meta.model.IGE.subset4ARD.NI, digits=3) 



################################################################################
# PUBLICATION BIAS ANALYSES
################################################################################

# Following recommendations in Nakagawa et al. 2022, MEE

################################################################################
# Is there evidence of small-study effects for Social_h2_2_Zr?

# calculating SE for publication bias
dataset.IGE.subset1A$sei <- sqrt(dataset.IGE.subset1A$VZr)

# fitting uni-moderator metaregression: SE
meta.model.IGE.subset1A.SSE <- rma.mv(Social_h2_2_Zr,
                                      subset1A_VCV_ESVar, 
                                      mods = ~ sei,
                                      random = list(~ 1 | Paper_id,
                                                    ~ 1 | Group_id,
                                                    ~ 1 | Population2,
                                                    ~ 1 | Species_name.2,
                                                    ~ 1 | Species_name.2.phylo,
                                                    ~ 1 | Record_id),
                                      R = list(Species_name.2.phylo = phylo_cor.subset1A),
                                      method = "REML", 
                                      test = "t", 
                                      data = dataset.IGE.subset1A)

# saving the model for script 006_figures.R
save(meta.model.IGE.subset1A.SSE, file = "data/models/meta_model_IGE_subset1A_SSE.Rdata")

# model can be loaded instead of run using the following
# load("data/models/meta_model_IGE_subset1A_SSE.Rdata")

summary(meta.model.IGE.subset1A.SSE,digits = 3) 
# since the slope is not different from 0, there is no evidence of asymmetry in
# the funnel plot, and thus, no evidence of small-study effects in this dataset

# Calculate marginal R2 with r2_ml
R2.IGE.subset1A.SSE <- r2_ml(meta.model.IGE.subset1A.SSE) 
round(R2.IGE.subset1A.SSE*100, 1)

# in addition, adding the SE as moderator explains close to 0% of the
# heterogeneity

# since the intercept is statistically significant from 0, it is recommended
# to fit the same model but using sampling variance instead of SE to obtain
# a less potentially biased intercept. In this case, however, given that the
# p-value of the slope is rather high, it is unlikely that the intercept of the 
# following model would be much different than that of the previous. Run but
# not shown

# fitting uni-moderator metaregression: SV
meta.model.IGE.subset1A.SSE.SV <- rma.mv(Social_h2_2_Zr,
                                         subset1A_VCV_ESVar,
                                         mods = ~ 1 + VZr,
                                         random = list(~ 1 | Paper_id,
                                                       ~ 1 | Group_id,
                                                       ~ 1 | Population2,
                                                       ~ 1 | Species_name.2,
                                                       ~ 1 | Species_name.2.phylo,
                                                       ~ 1 | Record_id),
                                         R = list(Species_name.2.phylo = phylo_cor.subset1A),
                                         method = "REML",
                                         test = "t",
                                         data = dataset.IGE.subset1A)

summary(meta.model.IGE.subset1A.SSE.SV, digits = 3)


################################################################################
# Is there evidence of decline effects for Social_h2_2_Zr?

# mean-centring Year_of_Publication
dataset.IGE.subset1A$Year_of_Publication.c <- scale(dataset.IGE.subset1A$Year_of_Publication,
                                                    scale=F)

# fitting uni-moderator metaregression: SE
meta.model.IGE.subset1A.DE <- rma.mv(Social_h2_2_Zr,
                                     subset1A_VCV_ESVar, 
                                     mods = ~ Year_of_Publication.c,
                                     random = list(~ 1 | Paper_id,
                                                   ~ 1 | Group_id,
                                                   ~ 1 | Population2,
                                                   ~ 1 | Species_name.2,
                                                   ~ 1 | Species_name.2.phylo,
                                                   ~ 1 | Record_id),
                                     R = list(Species_name.2.phylo = phylo_cor.subset1A),
                                     method = "REML", 
                                     test = "t", 
                                     data = dataset.IGE.subset1A)

# saving the model for script 006_figures.R
save(meta.model.IGE.subset1A.DE, file = "data/models/meta_model_IGE_subset1A_DE.Rdata")

# model can be loaded instead of run using the following
# load("data/models/meta_model_IGE_subset1A_DE.Rdata")

summary(meta.model.IGE.subset1A.DE,3) 
# since the slope is not different from 0, there is no evidence of decline
# effects in this dataset

# Calculate marginal R2 with r2_ml
R2.IGE.subset1A.DE <- r2_ml(meta.model.IGE.subset1A.DE) 
round(R2.IGE.subset1A.DE*100, 1)


