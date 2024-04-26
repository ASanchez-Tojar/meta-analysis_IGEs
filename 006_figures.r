################################################################################
# Authors: 
#
# Alfredo Sanchez-Tojar (@ASanchez_Tojar)
# Profile: https://goo.gl/PmpPEB
# Bielefeld University
# Email: alfredo.tojar@gmail.com

################################################################################
# Description of script and Instructions
################################################################################

# Script first created on the 19th of January 2024

# This script is to create the figures of a meta-analysis on Indirect Genetic
# Effects across species.

################################################################################
# Packages needed
################################################################################

# load packages
pacman::p_load(metafor,orchaRd, ggpubr,tidyr, extrafont, ggplot2, DT, RColorBrewer)

# cleaning up
rm(list=ls())

################################################################################
# Functions needed
################################################################################

# none

################################################################################
# Importing data and models
################################################################################

# importing main dataset
dataset.IGE.subset1A <- read.csv("data/subsets/dataset_IGE_subset1A.csv",
                                 header=T)

# importing subsets
dataset.IGE.subset3A_long_no_zeros <- read.csv("data/subsets/dataset_IGE_subset3A_long_no_zeros.csv",
                                               header=T)
dataset.IGE.subset3B_long <- read.csv("data/subsets/dataset_IGE_subset3B_long.csv",
                                      header=T)
dataset.IGE.subset3C_long <- read.csv("data/subsets/dataset_IGE_subset3C_long.csv",
                                      header=T)
dataset.IGE.subset4A_long <- read.csv("data/subsets/dataset_IGE_subset4A_long.csv",
                                      header=T)
dataset.IGE.subset4B <- read.csv("data/subsets/dataset_IGE_subset4B.csv",
                                 header=T)

# models imported from script 005_data_analysis.R
load("data/models/meta_model_IGE_subset1A.Rdata") # meta.model.IGE.subset1A
load("data/models/IGEmeta_regression_FEPartner.Rdata") # IGEmeta.regression.FEPartner
load("data/models/IGEmeta_regression_TraitCat.Rdata") # IGEmeta.regression.TraitCat
load("data/models/IGEmeta_regression_Age.Rdata") # IGEmeta.regression.Age
load("data/models/IGEmeta_regression_Sex.Rdata") # IGEmeta.regression.Sex
load("data/models/IGEmeta_regression_StudyType.Rdata") # IGEmeta.regression.StudyType
load("data/models/IGEmeta_regression_PopType.Rdata") # IGEmeta.regression.PopType
load("data/models/IGEmeta_regression_Livestock.Rdata") # IGEmeta.regression.Livestock
load("data/models/IGEmeta_regression_Va_vs_Vige.Rdata") # meta.model.IGE.subset3A_no_zeros
load("data/models/IGEmeta_regression_h2_vs_socialh2.Rdata") # meta.model.IGE.subset3B
load("data/models/IGEmeta_regression_Ia_vs_Iige.Rdata") # meta.model.IGE.subset3C
load("data/models/IGEmeta_regression_h2_vs_Totalh2.Rdata") # meta.model.IGE.subset4A
load("data/models/meta_model_IGE_DGE-IGE_correlation.Rdata") # meta.model.IGE.subset4B
load("data/models/meta_model_IGE_subset1A_SSE.Rdata") # meta.model.IGE.subset1A.SSE
load("data/models/meta_model_IGE_subset1A_DE.Rdata") # meta.model.IGE.subset1A.DE

################################################################################
# Q1: What is the proportion of variance in traits explained by IGEs?

# Printing the summary results of the model
print(meta.model.IGE.subset1A, digits=3)
predict(meta.model.IGE.subset1A, digits=3)

meta.model.IGE.subset1A_plot <- orchard_plot(mod_results(meta.model.IGE.subset1A, 
                                                         mod = "1",
                                                         group = "Paper_id", 
                                                         data = dataset.IGE.subset1A), 
                                             xlab = "Effect size (r)", 
                                             trunk.size = 2, 
                                             branch.size = 2,
                                             alpha = 0.3,
                                             transfm = "tanh",
                                             fill = T)+
  scale_fill_manual(values="grey") +
  scale_colour_manual(values="grey")+
  theme(legend.direction="horizontal", legend.title = element_text(size =8),
        legend.text = element_text(size = 10), 
        axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_blank()) +
  labs( x=expression(Social~h^2))

#meta.model.IGE.subset1A_plot

ggsave("figures/Fig2.png", plot = meta.model.IGE.subset1A_plot, 
       height = 100, width = 190, units = "mm", dpi = 600)


################################################################################
# Q2: Does the magnitude of IGEs vary among with factors? ####

# 1.	B) Are there differences between: (i.e. moderators for the meta-analysis 
# on point 1A)


################################################################################
# Fixed effect of partner

# Printing the summary results of the model
print(IGEmeta.regression.FEPartner, digits=3)

IGEmeta.regression.FEPartner_plot <- orchard_plot(mod_results(IGEmeta.regression.FEPartner,
                                                             mod = "Fixed_eff_of_partner_trait",
                                                             group = "Paper_id",
                                                             data = dataset.IGE.subset1A),
                                                 xlab = "Effect size (r)",
                                                 angle = 0,
                                                 trunk.size = 2,
                                                 branch.size = 1.5,
                                                 alpha = 0.3,
                                                 transfm = "tanh",
                                                 fill = T)

IGEmeta.regression.FEPartner_plot = IGEmeta.regression.FEPartner_plot +
 scale_x_discrete(labels = c('Partner effect: no','Partner effect: yes'))

# IGEmeta.regression.FEPartner_plot

ggsave("figures/FigS4.png", plot = IGEmeta.regression.FEPartner_plot, 
       height = 100, width = 190, units = "mm", dpi = 600)


################################################################################
# Trait category

# Printing the summary results of the model
print(IGEmeta.regression.TraitCat, digits=3)

IGEmeta.regression.TraitCat_plot <- orchard_plot(mod_results(IGEmeta.regression.TraitCat, 
                                                             mod = "Trait_category", 
                                                             group = "Paper_id", 
                                                             data = dataset.IGE.subset1A), 
                                                 xlab = "Effect size (r)",
                                                 angle = 0, 
                                                 trunk.size = 2,
                                                 branch.size = 1.5,
                                                 alpha = 0.3,
                                                 transfm = "tanh",
                                                 fill = T)

IGEmeta.regression.TraitCat_plot = IGEmeta.regression.TraitCat_plot + 
  theme(plot.title = element_text(face="bold", size = 12,hjust = 0.5),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 10))+
  scale_x_discrete(labels = c('Behaviour','Development',
                              'Metabolism & Physiology', "Morphology", 
                              "Reproduction", "Survival"))

#IGEmeta.regression.TraitCat_plot

################################################################################
# Age

# Printing the summary results of the model
print(IGEmeta.regression.Age, digits=3)

IGEmeta.regression.Age_plot <- orchard_plot(mod_results(IGEmeta.regression.Age, 
                                                        mod = "Age", 
                                                        group = "Paper_id", 
                                                        data = dataset.IGE.subset1A), 
                                            xlab = "Effect size (r)",
                                            angle = 0, 
                                            trunk.size = 2,
                                            branch.size = 1.5,
                                            alpha = 0.3,
                                            transfm = "tanh",
                                            fill = T)

IGEmeta.regression.Age_plot = IGEmeta.regression.Age_plot+
  scale_x_discrete(labels=c("Juv" = "      Juvenile", "Ad" = "Adult"))+
  theme(plot.title = element_text(face="bold", size = 12,hjust = 0.5),
                                     axis.title.y = element_text(size = 15),
                                     axis.title.x = element_blank(),
                                     axis.text.y = element_text(size = 15), 
                                     axis.text.x = element_text(size = 10))

#IGEmeta.regression.Age_plot

################################################################################
# Sex

# Printing the summary results of the model
print(IGEmeta.regression.Sex, digits=3)

IGEmeta.regression.Sex_plot <- orchard_plot(mod_results(IGEmeta.regression.Sex, 
                                                        mod = "Sex", 
                                                        group = "Paper_id", 
                                                        data = dataset.IGE.subset1A), 
                                            xlab = "Effect size (r)",
                                            angle = 0, 
                                            trunk.size = 2,
                                            branch.size = 1.5,
                                            alpha = 0.3,
                                            transfm = "tanh",
                                            fill = T)

IGEmeta.regression.Sex_plot = IGEmeta.regression.Sex_plot + 
  scale_x_discrete(labels=c("M" = "Male", "F" = "Female"))+
  theme(plot.title = element_text(face="bold", size = 12,hjust = 0.5),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 10))
  # theme(plot.title = element_text(face="bold", size = 12,hjust = 0.5),
  #                                 axis.title.y = element_text(size = 15),
  #                                 axis.title.x = element_blank(size = 15),
  #                                 axis.text.y = element_text(size = 15), 
  #                                 axis.text.x = element_text(size = 10))

IGEmeta.regression.Sex_plot

# ################################################################################
# # Study type
# 
# # Printing the summary results of the model
# print(IGEmeta.regression.StudyType, digits=3)
# 
# IGEmeta.regression.StudyType_plot <- orchard_plot(mod_results(IGEmeta.regression.StudyType, 
#                                                               mod = "Study_type", 
#                                                               group = "Paper_id", 
#                                                               data = dataset.IGE.subset1A), 
#                                                   xlab = "Effect size (r)",
#                                                   angle = 0, 
#                                                   trunk.size = 2,
#                                                   branch.size = 1.5,
#                                                   alpha = 0.3,
#                                                   transfm = "tanh",
#                                                   fill = T)
# 
# IGEmeta.regression.StudyType_plot = IGEmeta.regression.StudyType_plot+
# theme(plot.title = element_text(face="bold", size = 12,hjust = 0.5),
#       axis.title.y = element_text(size = 15),
#       axis.title.x = element_text(size = 15),
#       axis.text.y = element_text(size = 15),
#       axis.text.x = element_text(size = 10))
# IGEmeta.regression.StudyType_plot

################################################################################
# Population type

# Printing the summary results of the model
print(IGEmeta.regression.PopType, digits=3)

IGEmeta.regression.PopType_plot <- orchard_plot(mod_results(IGEmeta.regression.PopType, 
                                                            mod = "Population_type", 
                                                            group = "Paper_id", 
                                                            data = dataset.IGE.subset1A), 
                                                xlab = "Effect size (r)",
                                                angle = 0, 
                                                trunk.size = 2,
                                                branch.size = 1.5,
                                                alpha = 0.3,
                                                transfm = "tanh",
                                                fill = T)

IGEmeta.regression.PopType_plot = IGEmeta.regression.PopType_plot+
  theme(plot.title = element_text(face="bold", size = 12,hjust = 0.5),
                                     axis.title.y = element_text(size = 15),
                                     axis.title.x = element_text(size = 15),
                                     axis.text.y = element_text(size = 15), 
                                     axis.text.x = element_text(size = 10))

#IGEmeta.regression.PopType_plot

################################################################################
# Livestock # exploratory analyses after seeing the results

# Printing the summary results of the model
print(IGEmeta.regression.Livestock, digits=3)

IGEmeta.regression.Livestock_plot <- orchard_plot(mod_results(IGEmeta.regression.Livestock,
                                                              mod = "Livestock",
                                                              group = "Paper_id",
                                                              data = dataset.IGE.subset1A),
                                                  xlab = "Effect size (r)",
                                                  angle = 0, 
                                                  trunk.size = 2,
                                                  branch.size = 1.5,
                                                  alpha = 0.3,
                                                  transfm = "tanh",
                                                  fill = T)

IGEmeta.regression.Livestock_plot <- IGEmeta.regression.Livestock_plot + 
scale_x_discrete(labels=c("Yes" = "Livestock", "No" = "Non-livestock"))

ggsave("figures/FigS5.png", plot = IGEmeta.regression.Livestock_plot, height = 100, width = 190, units = "mm", dpi = 600)

###################

fig3<- ggarrange(#IGEmeta.regression.FEPartner_plot,
                  IGEmeta.regression.TraitCat_plot,
                  IGEmeta.regression.Age_plot,
                  IGEmeta.regression.Sex_plot,
                  #IGEmeta.regression.StudyType_plot,
                  IGEmeta.regression.PopType_plot,
                  # IGEmeta.regression.Livestock_plot,
                  nrow = 2, ncol = 2,  
                  labels = c("(a)", "(b)","(c)","(d)"),common.legend = T)
                  # nrow = 3, ncol = 2,  
                  # labels = c("(a)", "(b)","(c)","(d)", "(e)"),common.legend = T)

fig3T=annotate_figure(fig3, top = text_grob(expression(bold(Social~h^2)), 
                                            color = "black", 
                                            face = "bold", size = 20))
#fig3T

ggsave("figures/Fig3A-E.png", plot = last_plot(), 
       height = 250, width = 400, units = "mm", dpi = 600,bg ="white")

# pdf("Fig.3.pdf", width=14, height=8, onefile=F)
# dev.off()


################################################################################
# Q3:. What's the relative and absolute importance of direct and indirect  
# genetic effects?
# a: compare VDGE to the VIGEs
# b: compare H2 and social H2
# c: compare direct evolvability vs indirect evolvability (for choices on this 
# see below)
################################################################################

################################################################################
# 3A: Va vs Vi####
# First we compare Va and Vi 

# Printing the summary results of the model
print(meta.model.IGE.subset3A_no_zeros, digits=3)

IGEmeta.regression.Va.vs.Vige_plot <- orchard_plot(mod_results(meta.model.IGE.subset3A_no_zeros, 
                                                               mod = "direct_social", 
                                                               group = "Paper_id", 
                                                               data = dataset.IGE.subset3A_long_no_zeros), 
                                                   xlab =  "Effect size (log(var))",
                                                   angle = 0, 
                                                   trunk.size = 1,
                                                   branch.size = 1.5,
                                                   alpha = 0.3,
                                                   fill = T)+
  ggtitle("Unstandardized") +
  scale_fill_manual(values=c("#f2d116","aquamarine4")) +
  scale_colour_manual(values=c("#f2d116","aquamarine4"))+
  scale_x_discrete(labels=c("V_a" = "DGE", "V_ige_2" = "IGE"))+
  theme(plot.title = element_text(face="bold", size = 12,hjust = 0.5),
              legend.direction="horizontal", 
              legend.title = element_text(size =10),
              legend.text = element_text(size = 10), 
              axis.text.y = element_text(size = 15, angle = 90, vjust = 0.5, hjust=0.5),
              axis.title.x = element_text(size = 15),
              axis.text.x = element_text(size = 10))

IGEmeta.regression.Va.vs.Vige_plot

################################################################################
# 3B: h2 vs social h2

# Printing the summary results of the model
print(meta.model.IGE.subset3B, digits=3)

IGEmeta.regression.h2.vs.socialh2_plot <- orchard_plot(mod_results(meta.model.IGE.subset3B, 
                                                                   mod = "direct_social", 
                                                                   group = "Paper_id", 
                                                                   data = dataset.IGE.subset3B_long), 
                                                       xlab =  "Effect size (r)",
                                                       angle = 0, 
                                                       trunk.size = 1,
                                                       branch.size = 1.5,
                                                       alpha = 0.3,
                                                       transfm = "tanh",
                                                       fill = T)+
  ggtitle("Variance-standardized") +
  scale_fill_manual(values=c("#f2d116","aquamarine4")) +
  scale_colour_manual(values=c("#f2d116","aquamarine4"))+
  scale_x_discrete(labels=c("Social_h2_2_Zr" = "IGE", "H2_2_Zr" = "DGE"))+
  theme(plot.title = element_text(face="bold", size = 12,hjust = 0.5),
        axis.text.y = element_blank(), 
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 10))

IGEmeta.regression.h2.vs.socialh2_plot

################################################################################
# 3C: Ia vs Ii

# Printing the summary results of the model
print(meta.model.IGE.subset3C, digits=3)

IGEmeta.regression.Ia.vs.Iige_plot <- orchard_plot(mod_results(meta.model.IGE.subset3C, 
                                                               mod = "direct_social", 
                                                               group = "Paper_id", 
                                                               data = dataset.IGE.subset3C_long), 
                                                   xlab =  "Effect size (I)",
                                                   angle = 0, 
                                                   trunk.size = 1,
                                                   branch.size = 1.5,
                                                   alpha = 0.3,
                                                   fill = T)+
  ggtitle("Mean-standardized") +
  scale_fill_manual(values=c("#f2d116","aquamarine4")) +
  scale_colour_manual(values=c("#f2d116","aquamarine4"))+
  scale_x_discrete(labels=c("I_a" = "DGE","I_ige" = "IGE"))+
  theme(plot.title = element_text(face="bold", size = 12,hjust = 0.5),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 10))

IGEmeta.regression.Ia.vs.Iige_plot


fig4 <- ggarrange(IGEmeta.regression.Va.vs.Vige_plot,IGEmeta.regression.h2.vs.socialh2_plot,
                  IGEmeta.regression.Ia.vs.Iige_plot, nrow = 1, ncol = 3, 
                  labels = c("(a)", "(b)","(c)"),common.legend = T,
                  font.label=list(color="black",size=10)) + bgcolor("white") 

fig4

ggsave("figures/Fig4A-C.png", plot = fig4, 
      height = 150, width = 250, units = "mm", dpi = 600)


################################################################################
# Q4: Do IGEs typically alter evolutionary trajectories?


################################################################################
# 4A: Total heritability (T2) vs narrow-sense heritability (h2)

# Printing the summary results of the model
print(meta.model.IGE.subset4A, digits=3)

IGEmeta.regression.h2.vs.Totalh2_plot <- orchard_plot(mod_results(meta.model.IGE.subset4A, 
                                                                  mod = "direct_social", 
                                                                  group = "Paper_id", 
                                                                  data = dataset.IGE.subset4A_long), 
                                                      xlab = "Effect size (r)",
                                                      angle = 0,
                                                      trunk.size = 2,
                                                      branch.size = 1.5,
                                                      alpha = 0.3,
                                                      transfm = "tanh",
                                                      fill = T)+
  scale_fill_manual(values=c("#f2d116","aquamarine4")) +
  scale_colour_manual(values=c("#f2d116","aquamarine4"))+
  scale_x_discrete(labels=c("T2_2" = expression(~tau^2), "H2_2" = expression(~h^2)))+
  theme(plot.title = element_text(face="bold", size = 15,hjust = 0.5),
        legend.direction="horizontal", legend.title = element_text(size =10),
        legend.text = element_text(size = 10), 
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 15, angle = 90, vjust = 0.5, hjust=0.5),
        axis.text.x = element_text(size = 15))

IGEmeta.regression.h2.vs.Totalh2_plot


################################################################################
# 4B: What is the magnitude of the DGE-IGE correlation?

# Printing the summary results of the model
print(meta.model.IGE.subset4B, digits=3)
predict(meta.model.IGE.subset4B, digits=3)

meta.model.IGE.DGE.IGE.correlation_plot <- orchard_plot(mod_results(meta.model.IGE.subset4B, 
                                                                    mod = "1", 
                                                                    group = "Paper_id", 
                                                                    data = dataset.IGE.subset4B), 
                                                        xlab = "Effect size (r)",
                                                        trunk.size = 2, 
                                                        branch.size = 1.5,
                                                        alpha = 0.5,
                                                        transfm = "tanh",
                                                        fill = T)+
  #ggtitle("DGE-IGE correlation") +
  scale_fill_manual(values="grey") +
  scale_colour_manual(values="grey")+
  theme(plot.title = element_text(face="bold", size = 15,hjust = 0.5),
        legend.direction="horizontal", legend.title = element_text(size =10),
      legend.text = element_text(size = 10), 
      axis.title = element_text(size = 15),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_blank()) +
  labs( x="DGE-IGE correlation")

meta.model.IGE.DGE.IGE.correlation_plot

fig5 <- ggarrange(IGEmeta.regression.h2.vs.Totalh2_plot,meta.model.IGE.DGE.IGE.correlation_plot, ncol = 2,  nrow = 1, 
                  labels = c("(a)", "(b)"), common.legend = T,
                  font.label=list(color="black",size=10)) 

fig5

ggsave("figures/Fig5A-B.png", 
       plot = fig5, height = 100, width = 200, units = "mm", dpi = 600, bg="white")



################################################################################
# PUBLICATION BIAS ANALYSES
################################################################################

# Following recommendations in Nakagawa et al. 2022, MEE

################################################################################
# Is there evidence of small-study effects for Social_h2_2_Zr?

# calculating SE for publication bias
dataset.IGE.subset1A$sei <- sqrt(dataset.IGE.subset1A$VZr)

# Printing the summary results of the model
print(meta.model.IGE.subset1A.SSE, digits=3)

IGEmeta.regression.small.study.effects_plot <- orchaRd::bubble_plot(meta.model.IGE.subset1A.SSE,
                                                                    mod = "sei", 
                                                                    group = "Paper_id", 
                                                                    xlab = "SE", 
                                                                    ylab = "Indirect Genetic Effects: Social h2 (Zr)", 
                                                                    legend.pos = "top.left")

IGEmeta.regression.small.study.effects_plot


################################################################################
# Is there evidence of decline effects for Social_h2_2_Zr?

# mean-centring Year_of_Publication
dataset.IGE.subset1A$Year_of_Publication.c <- scale(dataset.IGE.subset1A$Year_of_Publication,
                                                    scale=F)


# Printing the summary results of the model
print(meta.model.IGE.subset1A.DE, digits=3)

IGEmeta.regression.decline.effects_plot <- orchaRd::bubble_plot(meta.model.IGE.subset1A.DE,
                                                                    mod = "Year_of_Publication.c", 
                                                                    group = "Paper_id", 
                                                                    xlab = "Year of Publication (mean-centred; mean = 2014)", 
                                                                    ylab = "Indirect Genetic Effects: Social h2 (Zr)", 
                                                                    legend.pos = "top.left")

IGEmeta.regression.decline.effects_plot
