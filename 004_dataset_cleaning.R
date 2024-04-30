
##############################################################
# Authors: 
#
# Francesca Santostefano
# Profile: https://www.researchgate.net/profile/Francesca-Santostefano
# University du Quebec Montreal
# Email: F.Santostefano2@exeter.ac.uk
#
# David Fisher
# Profile: https://evoetholab.com/
# Univeristy of Aberdeen
# Email: david.fisher@abdn.ac.uk
#
# Alfredo Sanchez-Tojar (@ASanchez_Tojar)
# Profile: https://goo.gl/PmpPEB
# Bielefeld University
# Email: alfredo.tojar@gmail.com

##############################################################
# Description of script and Instructions
##############################################################

# This script is to import and clean the data collected for a systematic review
# on indirect genetic effects conducted in Web of Science and
# Scopus.

#necessary files:
#"Data from papers-Grid view_FS.csv"
#"Data from papers AST-Grid view_FS.csv"
#"Data from papers DF-Grid view_FS.csv"
#"Data from papers FS-Grid view_FS.csv"
#"Data from papers MM-Grid view_FS.csv"


##############################################################

# Packages needed 

##############################################################

pacman::p_load(data.table,stringr,ape,plyr,rotl,treebase,ggplot2)

# cleaning up
rm(list=ls())

##############################################################

#Reading the data

##############################################################

dataTR <- fread("data/Data from papers-Grid view_FS.csv")
dataAST <- fread("data/Data from papers AST-Grid view_FS.csv")
dataDF <- fread("data/Data from papers DF-Grid view_FS.csv")
dataFS <- fread("data/Data from papers FS-Grid view_FS.csv")
dataMM <- fread("data/Data from papers MM-Grid view_FS.csv")

#checking structure of each dataset
#str(dataTR) #15 obs. of  62 variables
#str(dataAST) #38 obs. of  62 variables
#str(dataDF) #71 obs. of  62 variables
#str(dataFS) #46 obs. of  62 variables
#str(dataMM) #43 obs. of  62 variables

#converting to data frame
dataTR<-as.data.frame(dataTR)
dataAST<-as.data.frame(dataAST)
dataDF<-as.data.frame(dataDF)
dataFS<-as.data.frame(dataFS)
dataMM<-as.data.frame(dataMM)

#making one df
data<-rbind(dataTR, dataAST, dataDF, dataFS, dataMM)

#renaming columns
names(data)<-str_replace_all(names(data), c(" " = "_" ))
names(data)<-str_replace_all(names(data), c("/" = "_" ))
names(data)<-str_replace_all(names(data), c("#" = "N" ))
names(data)<-str_replace_all(names(data), c("-" = "_" ))
names(data)[names(data) == "VPE_(foc)" ] <- "VPE_foc"
names(data)[names(data) == "VPE_(soc)" ] <- "VPE_soc"
names(data)[names(data) == "r_PE_PE(S)" ] <- "r_PE_PE_Soc"
names(data)[names(data) == "cov__PE_PE(S)" ] <- "cov_PE_PE_Soc"

#capitalize first letter of column names
names(data)<-str_to_title(names(data)) 
names(data)

#replace blank cells with NAs
data[data==""]<-NA

#####going through each column for formatting

#removing papers excluded during data extraction
data[is.na(data$Excluded_during_data_extraction),"Excluded_during_data_extraction"] <- "no" 
#####length(unique(data[data$Excluded_during_data_extraction=="yes","Paper_id"]))
##### exclude them
data <- data[data$Excluded_during_data_extraction=="no",]

#Record_id
#giving a unique record id to each row
data$Record_id<-1:nrow(data)
data$Record_id<-as.factor(data$Record_id) #####convert to factor

#Paper_id
#convert as factor
data$Paper_id<-as.factor(data$Paper_id)

#Year
data$Year <- str_replace_all(data$Year, " - ", "-") ##### removing spaces from ranges
data$Year <- str_replace_all(data$Year, ", ", " and ") ##### transforming , to and to standardize the variable
is.na(data$Year) <- data$Year == "unclear" ##### unclear should be NA to standardize the variable

#Species_name
##### adjust names to match their own species, but first create a second Species_name variable to keep the original data as it was
data$Species_name.2 <- data$Species_name
data[data$Species_name.2=="Bos taurus (Aosta Chestnut and Aosta Black Pied breeds)" & !(is.na(data$Species_name.2)),"Species_name.2"] <- "Bos taurus"

# ##### checking that all species names are valid and no synonyms exist (i.e. same species named differently in the dataset)
# ##### obtaining data frame listing the Open Tree identifiers potentially matching our list of species (be aware that this will take a few minutes, and you can load the data below)
taxa <- tnrs_match_names(names = unique(data$Species_name.2))
 
# ##### saving the taxonomic data created on the 3rd of Dec 2022 to speed the process in the future
# save(taxa,file = "data/taxa_Open_Tree_of_Life.RData")

##### loading the taxonomic data created on the 3rd of Dec 2022
#load("data/taxa_Open_Tree_of_Life.RData") #taxa

##### check approximate matches
taxa[taxa$approximate_match==TRUE,]

##### check approximate matches
taxa[taxa$is_synonym==TRUE,]

##### using the same names for Species_name.2 as in the Tree of Life
data$Species_name.2 <- ifelse(data$Species_name.2=="Larus novaehollandiae scopulinus",
                              "Chroicocephalus scopulinus",
                              data$Species_name.2)

data$Species_name.2 <- ifelse(data$Species_name.2=="Litopenaeus vannamei",
                              "Penaeus vannamei",
                              data$Species_name.2)

##### check number of matches
taxa[taxa$number_matches>1,]#####Gadus morhua is the correct one because in our dataset it was included as a fish


##### further checks
ott_id_tocheck <- taxa[taxa$number_matches != 1,"ott_id"]

for(i in 1:length(ott_id_tocheck)){
  print(inspect(taxa, ott_id = ott_id_tocheck[i])) #####Gadus morhua is the correct one because in our dataset it was included as a fish
}

#Error in if (!check_numeric(ott_id)) { : 
#missing value where TRUE/FALSE needed 
#check this error


##### otherwise everything seems in order regarding the species names

#convert as factor
data$Species_name.2<-as.factor(data$Species_name.2)
levels(data$Species_name.2)

#check plots that can be deleted

#frequencies of species (by lines of data, not by study)
#sp_fq<-plyr::count(data, 'Species_name.2')
#sp_fq

#quick plot
#sp_fqp<-ggplot(sp_fq, aes(x=Species_name.2, y=freq, fill=Species_name.2)) +
#  geom_bar(stat="identity")+theme_minimal() +
#  ylab("Frequence") +
#  xlab("Species") +
#  geom_text(aes(label=freq))
#sp_fqp<-sp_fqp + coord_flip() 
#sp_fqp

#Taxon
#convert as factor
data$Taxon<-as.factor(data$Taxon)
levels(data$Taxon)

#standardize name
data$Taxon [data$Taxon == 'insects' ]<- 'arthropods'

##### reset factor levels so that insects dissapears
data$Taxon <- factor(data$Taxon)


#check plots that can be deleted

#frequencies of taxon (by lines of data, not by study)
#ta_fq<-plyr::count(data, 'Taxon')
#ta_fq

#quick plot
#ta_fqp<-ggplot(ta_fq, aes(x=Taxon, y=freq, fill=Taxon)) +
#  geom_bar(stat="identity")+theme_minimal() +
#  ylab("Frequence") +
#  xlab("Taxon") +
#  geom_text(aes(label=freq))
#ta_fqp<-ta_fqp + coord_flip() 
#ta_fqp

#population 

#we create a second column where we shorten and standardize the names by hand
data$Population2 <-NA

#We give the location + the country code, or when location not available, name of institute or breed, and country of first author affiliation
data$Population2 [data$Population == 'Capalbio, Italy' ]<- 'Capalbio_IT'
data$Population2 [data$Population == 'White Mountain Research Station, California, US' ]<- 'WhiteMountain_US'
data$Population2 [data$Population == 'Kluane Red Squirrel Project, Yukon, Canada' ]<- 'Kluane_CAN'
data$Population2 [data$Population == 'Kaikoura Peninsula, New Zealand' ]<- 'Kaikoura_NZ'
data$Population2 [data$Population == 'Research Centre Foulum, Denmark' ]<- 'Foulum_DK'
data$Population2 [data$Population == 'Manjimup, Western Australia' ]<- 'Manjimup_AUS'
data$Population2 [data$Population == 'Mandarte Island, BC, Canada' ]<- 'MandarteIsland_CAN'
data$Population2 [data$Population == 'Rivelin Valley, Sheffield, UK' ]<- 'RivelinValley_UK'
data$Population2 [data$Population == 'Jitra Aquaculture Extension Center, Kedah, Malaysia' ]<- 'Jitra_MY'
data$Population2 [data$Population == 'Hendrix Genetics' ]<- 'HendrixGenetics_NL'
data$Population2 [data$Population == 'Valle Aosta, Italy' ]<- 'ValleAosta_IT'
data$Population2 [data$Population == 'Aquaculture Extension Center, Jitra, Kedah State,Malaysia' ]<- 'Jitra_MY'
data$Population2 [data$Population == 'Yukon Territory, Canada' ]<- 'Kluane_CAN'
data$Population2 [data$Population == 'National Institute of Animal Science, Cheonan, South Korea)' ]<- 'Cheonan_KOR'
data$Population2 [data$Population == 'White Leghorns, Denmark' ]<- 'WhiteLeghorns_DK'
data$Population2 [data$Population == 'Institute of Pig Genetics, Netherlands' ]<- 'PigGenetics_NL'
data$Population2 [data$Population == 'Dahomey stock' ]<- 'Dahomey_UK'
data$Population2 [data$Population == 'Hendrix genetics' ]<- 'HendrixGenetics_NL'
data$Population2 [data$Population == 'IRTA, Spain' ]<- 'Irta_ESP'
data$Population2 [data$Population == 'White Leghorns, Netherlands' ]<- 'HendrixGenetics_NL'
data$Population2 [data$Population == 'Purdue Poultry Research Centre' ]<- 'Purdue_US'
data$Population2 [data$Population == 'Wheel running mice' ]<- 'WheelRunning_CAN'
data$Population2 [data$Population == 'Goland sire line, Gorzagri, Italy' ]<- 'Gorzagri_IT'
data$Population2 [data$Population == 'Landrace pigs, Landrace Sunin, Republic of Korea' ]<- 'LandraceSunin_KOR'
data$Population2 [data$Population == 'Aosta Region, Italy' ]<- 'ValleAosta_IT'
data$Population2 [data$Population == 'Wageningen, the Netherlands.' ]<- 'Wageningen_NL'
data$Population2 [data$Population == 'Pig Improvement Company, Franklin, KY, US' ]<- 'PigImprovementCompany_US'
data$Population2 [data$Population == 'University of Nebraska Agricultural Research Development Center Swine Research Farm' ]<- 'AgriculturalResearch_US'
data$Population2 [data$Population == 'Gotland, Sweden' ]<- 'Gotland_SWE'
data$Population2 [data$Population == 'Hainan Higene Aquaculture Technology, Wenchang City, China' ]<- 'WenchangCity_CHN'
data$Population2 [data$Population == 'Forster stock, University of Queensland, AUS' ]<- 'Forster_AUS'
data$Population2 [data$Population == 'University of Queensland, AUS' ]<- 'UniversityQueensland_AUS'
data$Population2 [data$Population == 'Institute for Pig Genetics, Beilen, The Netherlands' ]<- 'PigGenetics_NL'
data$Population2 [data$Population == 'Institut de Selection Animale B.V., Hendrix Genetics' ]<- 'HendrixGenetics_NL'
data$Population2 [data$Population == 'Linthal and Elm, Glarus, Switzerland' ]<- 'Glarus_SUI'
data$Population2 [data$Population == 'Matsalu National Park, Estonia' ]<- 'Matsalu_EST'
data$Population2 [data$Population == 'Mandarte island, Canada' ]<- 'MandarteIsland_CAN'
data$Population2 [data$Population == 'laboratory at the Ludwig Maximilians University of Munich' ]<- 'Capalbio_IT'
data$Population2 [data$Population == 'Isle of Rum, Scotland' ]<- 'IsleRum_UK'
data$Population2 [data$Population == 'Tromso, Norway' ]<- 'Tromso_NOR'
data$Population2 [data$Population == 'Smithfield Premium Genetics Group (Roanoke Rapids, NC)' ]<- 'Smithfield_US'
data$Population2 [data$Population == 'nucleus breeding farm in Korea' ]<- 'Farm_KOR'
data$Population2 [data$Population == 'IRTA experimental facilities (Spain)' ]<- 'Irta_ESP'
data$Population2 [data$Population == 'Institut de Sélection Animale B.V., the layer breeding division of Hendrix Genetics' ]<- 'HendrixGenetics_NL'
data$Population2 [data$Population == 'Danish pig breeding program, DanAvl, Denmark' ]<- 'DanAvl_DK'
data$Population2 [data$Population == 'Yorkshire Sunjin, Danyang, Korea,' ]<- 'YorkshireSunjin_KOR'
data$Population2 [data$Population == 'experimental farm of Institute for Pig Genetics (IPG) and experimental farm of the Animal Science Group (ASG), Netherlands' ]<- 'PigGenetics_NL'


#check  all the HendrixGenetics and WhiteLeghorn papers, if they are the same
#check<- data [which(data$Population2 == 'HendrixGenetics_NL' |data$Population2 == 'WhiteLeghorns_DK'),]
#check
#looks ok, white leghorn is the name of the breed, one paper is in DK all the others are in hendrix genetics or its other division (institute de selection animale) 
#so I assume they are the same pop and we name them Hendrix except for the DK pop

#check the Korean pigs papers if they are the same
#Farm_KOR Cheonan_KOR YorkshireSunjin_KOR
#check2<- data [which(data$Population2 == 'Farm_KOR' |data$Population2 == 'Cheonan_KOR'|data$Population2 == 'YorkshireSunjin_KOR'),]
#check2
#IGE 0846 has the same values as IGE 0786, We remove it as it is a duplicate study
#the paper only used social breeding values of individuals estimated with an IGE model to create different treatment groups to test the effect of high and low SBVs on behaviours
#we realized this is a duplicate at this point, but the removal could be included early at the beginning of the script

data$Excluded_during_data_extraction[data$Paper_id == "IGE0846"] = "yes" 

##### let's then remove it from the dataset straight away
data<- data[data$Excluded_during_data_extraction=="no",]


#check for NAs
sum(is.na (data$Population2))
##31 NA after having excluded studies with no data (see above). 
#rename NAs based on location or institute of first author
#check3<- data [is.na(data$Population2),]
#check3

data$Population2[data$Paper_id == 'IGE0199']<- 'Bielefeld_DE'
data$Population2[data$Paper_id == 'IGE1038']<- 'Sichuan_CHN'
data$Population2[data$Paper_id == 'IGE0242']<- 'AnimalResearch_US'
data$Population2[data$Paper_id == 'IGE0195']<- 'PigGenetics_NL' #should be same as other papers we have
data$Population2[data$Paper_id == 'IGE0144']<- 'Gorzagri_IT'  #should be same as other papers we have
data$Population2[data$Paper_id == 'IGE0920']<- 'PigGenetics_NL' #should be same as other papers we have
data$Population2[data$Paper_id == 'IGE0927']<- 'NordicGenetics_SWE'
data$Population2[data$Paper_id == 'IGE0507']<- 'HendrixGenetics_NL'

#study_coordinates
#there are very few records with this variable so we won't use it.
###### let's then delete the variable to keep our working dataset cleaner
data <- subset(data, select = -(Study_coordinates))

#Population_type
#convert as factor
data$Population_type<-as.factor(data$Population_type)
#levels(data$Population_type)

#frequencies of population type (by data lines, not study)
#pp_fq<-plyr::count(data, 'Population_type')
#pp_fq

#quick plot
#pp_fqp<-ggplot(pp_fq, aes(x=Population_type, y=freq, fill=Population_type)) +
#  geom_bar(stat="identity")+theme_minimal() +
#  ylab("Frequence") +
#  xlab("Population type") +
#  geom_text(aes(label=freq))
#pp_fqp<-pp_fqp + coord_flip() 
#pp_fqp


#Study_type
#convert as factor
data$Study_type<-as.factor(data$Study_type)
#levels(data$Study_type)

#frequencies of Study type (by data lines, not study)
#st_fq<-plyr::count(data, 'Study_type')
#st_fq

#quick plot
#st_fqp<-ggplot(st_fq, aes(x=Study_type, y=freq, fill=Study_type)) +
#  geom_bar(stat="identity")+theme_minimal() +
#  ylab("Frequence") +
#  xlab("Study type") +
#  geom_text(aes(label=freq))
#st_fqp<-st_fqp + coord_flip() 
#st_fqp

#Sex
##### IGE0494 had NA. We checked the original paper and they indeed do not mention anything about sex, but we then think it is safe to assume that both sexes were analyzed together.
data[data$Paper_id=="IGE0494","Sex"] <- "both"

#convert as factor
data$Sex<-as.factor(data$Sex)
#levels(data$Sex)

#frequencies of sex (by data lines, not study)
#sex_fq<-plyr::count(data, 'Sex') 
#sex_fq

#quick plot
#sex_fqp<-ggplot(sex_fq, aes(x=Sex, y=freq, fill=Sex)) +
#  geom_bar(stat="identity")+theme_minimal() +
#  ylab("Frequence") +
#  xlab("Sex") +
#  geom_text(aes(label=freq))
#sex_fqp<-sex_fqp + coord_flip() 
#sex_fqp

#Age
#convert as factor
data$Age<-as.factor(data$Age)
#levels(data$Age)

#frequencies of age (by data lines, not study)
#age_fq<-plyr::count(data, 'Age')
#age_fq

#quick plot
#age_fqp<-ggplot(age_fq, aes(x=Age, y=freq, fill=Age)) +
#  geom_bar(stat="identity")+theme_minimal() +
#  ylab("Frequence") +
#  xlab("Age") +
#  geom_text(aes(label=freq))
#age_fqp<-age_fqp + coord_flip() 
#age_fqp

#N_Id_Pedigree
#data$N_id_pedigree
#ok
#hist(data$N_id_pedigree,  breaks = 50)

#N_Records
#data$N_records
#ok
#sum(is.na (data$N_records))
##### 9 after having excluded studies with no data (see above)

#hist(data$N_records,  breaks = 50)

#N_id_w_records
#data$N_id_w_records
#ok

#hist(data$N_id_w_records,  breaks = 50)
#sum(is.na (data$N_id_w_records))
##### 10 after having excluded studies with no data (see above)

#N_Sires
#data$N_sires
#ok

#hist(data$N_sires,  breaks = 50)

#N_Dams
#data$N_dams
#ok

#hist(data$N_dams,  breaks = 50)

#N_Families
#data$N_families
#record 80 to 193 10 per generation" needs to be replaced; 78 depth of ped so 780
data$N_families[data$N_families == "10 per generation" ] <- "780"
#convert as numeric
data$N_families<-as.numeric(data$N_families)
#hist(data$N_families,  breaks = 50)

#Depth_of_ped
#data$Depth_of_ped
#ok
#hist(data$Depth_of_ped,  breaks = 50)
#sum(is.na (data$Depth_of_ped))
##### 63 after having excluded studies with no data (see above)

#Trait_name 
#I don't think we'll do anything with this, too heterogeneous/many levels
#convert as factor
data$Trait_name <-as.factor(data$Trait_name )
#levels(data$Trait_name )

#Trait_category
##### we checked the two entries where the trait type was categorised as other, and both were hoarding behaviour in North American red squirrels, which was now classify as a behaviour as otherwise the category is not large enough to be useful
data$Trait_category <- ifelse(data$Trait_category=="other",
                              "behaviour",
                              data$Trait_category)
  
#convert as factor
data$Trait_category<-as.factor(data$Trait_category)
#levels(data$Trait_category)

#frequencies of Trait_category
#tr_fq<-plyr::count(data, 'Trait_category')
#tr_fq

#quick plot
#tr_fqp<-ggplot(tr_fq, aes(x=Trait_category, y=freq, fill=Trait_category)) +
#  geom_bar(stat="identity")+theme_minimal() +
#  ylab("Frequence") +
#  xlab("Trait category") +
#  geom_text(aes(label=freq))
#tr_fqp<-tr_fqp + coord_flip() 
#tr_fqp

#Mean_standardized
#ok
##### transforming to lower case just because the other yes/no variables are in lowercase (standardization)
data$Mean_standardized <- tolower(data$Mean_standardized)
data$Mean_standardized<-as.factor(data$Mean_standardized)
#levels(data$Mean_standardized)

#frequencies of Mean standardized
#ms_fq<-plyr::count(data, 'Mean_standardized')
#ms_fq

#Trait_mean
#data$Trait_mean
#ok

#hist(data$Trait_mean,  breaks = 50)
#sum(is.na (data$Trait_mean))
#75 after having excluded studies with no data (see above)

#Variance_standardized
#ok
##### transforming to lower case just because the other yes/no variables are in lowercase (standardization)
data$Variance_standardized <- tolower(data$Variance_standardized)
data$Variance_standardized<-as.factor(data$Variance_standardized)
#levels(data$Variance_standardized)

#frequencies of Variance standardized
#vs_fq<-plyr::count(data, 'Variance_standardized')
#vs_fq

#Trait_sd
#data$Trait_sd
#ok

#hist(data$Trait_sd,  breaks = 50)
#sum(is.na (data$Trait_sd))
##### 97 after having excluded studies with no data (see above), and therefore, there are 22 Trait_mean without an estimate of SD

#quickly plotting ln(mean)~ln(SD) to spot any potential outliers; but keep in mine that values =<0 are not plotted because of the log transformation
#plot(log(data$Trait_sd),log(data$Trait_mean))

#frequencies of Treatment_Group
#I don't think we'll do anything with this, too heterogeneous/many levels and mostly NAs
##### we need to transform this variable into an ID for us to have the possibility to include it as a random effect
data$Group_id <- as.numeric(as.factor(paste(data$Paper_id,data$Treatment_group,sep="_")))
# data$Treatment_group<-as.factor(data$Treatment_group)
# sum(is.na(data$Treatment_group))
# # 152
# 
# #frequencies of Treatment_Group
# tg_fq<-plyr::count(data, 'Treatment_group')
# tg_fq


#Fixed_eff_of_partner_trait
data$Fixed_eff_of_partner_trait<-as.factor(data$Fixed_eff_of_partner_trait)
#levels(data$Fixed_eff_of_partner_trait)

#frequencies of Fixed_eff_of_partner_trait
#pt_fq<-plyr::count(data, 'Fixed_eff_of_partner_trait')
#pt_fq

#Other_fixed_eff
data$Other_fixed_eff<-as.factor(data$Other_fixed_eff)
#levels(data$Other_fixed_eff)

#frequencies of Other_fixed_eff
#of_fq<-plyr::count(data, 'Other_fixed_eff')
#of_fq

#Mean_group_size 
data$Mean_group_size
#record 28 to 30 "9-11" needs to be replaced (with mean, so 10) 
data$Mean_group_size[data$Paper_id == 'IGE0243'] <- "10"
data$Mean_group_size<-as.numeric(data$Mean_group_size)
#hist(data$Mean_group_size ,  breaks = 50)

#Va 
#data$Va
#replace , with nothing because those "," were meant to represent the thousand limit in the original paper
data$Va<-gsub(',','',data$Va)
#ok (huge range variation depending on whether it's std or not)
#data$Va <-as.numeric(data$Va )
#hist(data$Va,  breaks = 50)

##### for those standardized, checking whether values are as expected (i.e. between 0 and 1?)
#summary(data[data$Variance_standardized=="yes","Va"]) ##### all as expected

#V_ige
#data$V_ige
#replace then only <0.001 value with 0.0001 to make it usable for the analyses
data$V_ige[data$V_ige == "<0.001" ] <- "0.0001"
#ok (huge range variation depending on whether it's std or not)
#data$V_ige <-as.numeric(data$V_ige )
#hist(data$V_ige,  breaks = 50)

##### for those standardized, checking whether values are as expected (i.e. between 0 and 1?)
#summary(data[data$Variance_standardized=="yes","V_ige"]) ##### all as expected

#Vpe_foc 
#data$Vpe_foc
#replace then only <0.001 value with 0.0001 to make it usable for the analyses
data$Vpe_foc[data$Vpe_foc == "<0.001" ] <- "0.0001"
#ok (huge range variation depending on whether it's std or not, lots of NAs)
data$Vpe_foc<-as.numeric(data$Vpe_foc)
#hist(data$Vpe_foc,  breaks = 50)

##### for those standardized, checking whether values are as expected (i.e. between 0 and 1?)
summary(data[data$Variance_standardized=="yes","Vpe_foc"]) ##### all as expected

#Vpe_soc 
#data$Vpe_soc
#replace then only <0.001 value with 0.0001 to make it usable for the analyses
data$Vpe_soc[data$Vpe_soc == "<0.001" ] <- "0.0001"
#replace , with nothing because those "," were meant to represent the thousand limit in the original paper
data$Vpe_soc<-gsub(',','',data$Vpe_soc)
#ok (huge range variation depending on whether it's std or not, lots of NAs)
data$Vpe_soc<-as.numeric(data$Vpe_soc)
#hist(data$Vpe_soc,  breaks = 50)

##### for those standardized, checking whether values are as expected (i.e. between 0 and 1?)
#summary(data[data$Variance_standardized=="yes","Vpe_soc"]) ##### all as expected

#other variances will only really be used to add up to VP  
#data$V_other_1
#replace , with nothing because those "," were meant to represent the thousand limit in the original paper
data$V_other_1<-gsub(',','',data$V_other_1)
#ok
data$V_other_1<-as.numeric(data$V_other_1)

##### for those standardized, checking whether values are as expected (i.e. between 0 and 1?)
summary(data[data$Variance_standardized=="yes","V_other_1"]) ##### all as expected

#data$V_other_2
data$V_other_2<-as.numeric(data$V_other_2)

##### for those standardized, checking whether values are as expected (i.e. between 0 and 1?)
summary(data[data$Variance_standardized=="yes","V_other_2"]) ##### all as expected

#data$V_other_3
data$V_other_3<-as.numeric(data$V_other_3)

##### for those standardized, checking whether values are as expected (i.e. between 0 and 1?)
summary(data[data$Variance_standardized=="yes","V_other_3"]) ##### all as expected, there is only one left with var standardized

#we have very few values left for these 
#data$V_other_4
data$V_other_4<-as.numeric(data$V_other_4)

#data$V_other_5
data$V_other_5<-as.numeric(data$V_other_5)

#data$V_other_6
data$V_other_6<-as.numeric(data$V_other_6)

#V_residual 
#data$V_residual
#replace (pi^2)/3 with 3.289868
data$V_residual[data$V_residual == "(pi^2)/3" ] <- "3.289868"
#replace , with nothing
data$V_residual<-gsub(',','',data$V_residual)
#ok (huge range variation depending on whether it's std or not)
data$V_residual<-as.numeric(data$V_residual)
#hist(data$V_residual,  breaks = 50)

##### for those standardized, checking whether values are as expected (i.e. between 0 and 1?)
summary(data[data$Variance_standardized=="yes","V_residual"]) ##### all as expected except value 2286, but that is correct because the authors of that paper provided all values but the residual variance as standardized

#Total_v_phen
#data$Total_v_phen
#replace , with nothing because those "," were meant to represent the thousand limit in the original paper
data$Total_v_phen<-gsub(',','',data$Total_v_phen)
#ok (huge range variation depending on whether it's std or not)
data$Total_v_phen<-as.numeric(data$Total_v_phen)
#hist(data$Total_v_phen,  breaks = 50)
sum(is.na (data$Total_v_phen))
#96 records where it's not already provided 

#Total_v_phen got swapped (during data entry I assume): change it 
data[data$Paper_id=="IGE1038" & data$Trait_name =="residual feed intake (RFI)","Total_v_phen"] = 33449.05
data[data$Paper_id=="IGE1038" & data$Trait_name=="feed conversion ratio (FCR)","Total_v_phen"] = 307.61


#H2
#data$H2
#replace then only <0.001 value with 0.0001 to make it usable for the analyses
data$H2[data$H2 == "<0.001" ] <- "0.0001"
data$H2<-as.numeric(data$H2)
#hist(data$H2,  breaks = 50)
#check whether record 171 is a typo, it should be between 0 and 1 
data[data$Paper_id=="IGE1050","H2"] <- 0.145

#Social_h2
#data$Social_h2
#hist(data$Social_h2,  breaks = 50)
#check whether record 171 is a typo,it should be between 0 and 1 
#yes, replace 4.8 with 0.048
#data[data$Record_id==171, "Social_h2"] <- 0.048
data[data$Paper_id=="IGE1050","Social_h2"] <- 0.048


#R_a_ige
#data$R_a_ige
#some reformatting of values
data$R_a_ige<-as.numeric(data$R_a_ige)
#hist(data$R_a_ige,  breaks = 50)
#check whether records 202 and 192 are a typo, it should be between -1 and 1
#record 192 is not a typo, however this value gets dropped anyway with the code assigning NA to those estimates for which V_ige == 0 (just below)
#202 is a typo
data[data$R_a_ige==-37.0000 & !(is.na(data$R_a_ige)), "R_a_ige"] <- -0.37 
#also when variances are reported as 0, then the correlation should be set to undefined, so we replace those instances with NA
sum(data$V_ige == 0, na.rm=T)
data$R_a_ige [data$V_ige == 0]<- "NA"
data$R_a_ige<-as.numeric(data$R_a_ige)
#hist(data$R_a_ige,  breaks = 50)

#R_pe_pe_soc
data$R_pe_pe_soc
#ok 
#very few values, so we will not do anything with this  
#hist(data$R_pe_pe_soc,  breaks = 50)


#Cov_a_ige
#data$Cov_a_ige
#some reformatting of values
data$Cov_a_ige<-gsub(' ','',data$Cov_a_ige)
#replace then only <0.001 value with 0.0001 to make it usable for the analyses
data$Cov_a_ige[data$Cov_a_ige == "<0.001" ] <- "0.0001"
#replace , with nothing
data$Cov_a_ige<-gsub(',','',data$Cov_a_ige)

data$Cov_a_ige<-as.numeric(data$Cov_a_ige)
#hist(data$Cov_a_ige,  breaks = 50)

#Cov_pe_pe_soc
#data$Cov_pe_pe_soc
#very few values
#hist(data$Cov_pe_pe_soc,  breaks = 50)

#V_tbv
#data$V_tbv
#replace , with nothing
data$V_tbv<-gsub(',','',data$V_tbv)
data$V_tbv<-as.numeric(data$V_tbv)
#hist(data$V_tbv,  breaks = 50)

#T2 can be outside 0 and 1
#data$T2
#ok
#hist(data$T2,  breaks = 50)

#"Data_location"
#data$Data_location
##### just adding some small changes to make it look a bit cleaner and more standardized
data$Data_location <- str_replace_all(data$Data_location, c("text" = "Text" ))
data$Data_location <- str_replace_all(data$Data_location, c("page" = "pg" ))
data$Data_location <- str_replace_all(data$Data_location, c("and" = "&" ))
data$Data_location <- str_replace_all(data$Data_location, c("p " = "pg " ))
data$Data_location <- str_replace_all(data$Data_location, c("p2" = "pg 2" ))
data$Data_location <- str_replace_all(data$Data_location, c(" -" = "-" ))

#Screener_id
#just for us but no need to manipulate
#data$Screener_id

#Notes
#for us to help data  post processing, no need to manipulate, will use them one by one 
data$Notes

#Screener_id
#just for us but no need to manipulate
data$Second_screener_id

#Second_screener_notes
#for us to help data post processing, no need to manipulate, will use them one by one 
data$Second_screener_notes

##### I have excluded these studies at the beginning of the script to (1) one make coding easier for some particular case, (2) because I don't think we need to keep them for anything else at this point. 
# #Excluded_during_data_extraction
# #convert as factor
# #where yes, remove from dataset? or keep for summary analyses of literature? for now we keep them
# data$Excluded_during_data_extraction<-as.factor(data$Excluded_during_data_extraction)
# levels(data$Excluded_during_data_extraction)

#Decision
#this was for revision in extraction phase, no needed here 
##### We rename this variable to avoid confusion later on
data$Second_screener_original_decision <- data$Decision
data <- subset(data,select =  -c(Decision))
data$Second_screener_original_decision<-as.factor(data$Second_screener_original_decision)
#levels(data$Second_screener_original_decision)

#Changes_after_second_screening
#this was for us to keep track of changes in extraction, not needed here 
data$Changes_after_second_screening

#Final_comments
#to help with data post processing no need to manipulate, will use them one by one; only two screeners had it, 
# use along with Notes and Second_screener_notes from previous rounds to finalize and approve the extraction
data$Final_comments




######################################################################################################
#data processing

###we decided that we would rather contact all the authors who didn't provide VP than calculate ourselves, but in some cases calculating VP seems OK.

#QUESTIONS OF THE PROJECT AND VARIABLES WE NEED TO ANSWER
# 1.      A) What is the proportion of variance in traits explained by IGEs?

#we need the IGE/VP ratio 

#papers with NO v_ige
sum(is.na (data$V_ige)) 
#12 

#papers with no VP
sum(is.na (data$Total_v_phen)) 
#96

#some papers provide the ratio already 
sum(!is.na (data$Social_h2)) 
#46

#so for those who provide the ratio, or both variance components, we don't need to ask (we should safely calculate the ratio if we have the two components)
#check how many VIGE/VP we still don't have with 

missingSocial_h2<-complete.cases(data[, "Social_h2"])|(complete.cases(data[, "Total_v_phen"]) & complete.cases(data[, "V_ige"]))
data$missingSocial_h2<-(!missingSocial_h2)
sum(!missingSocial_h2)
#58 Vige/VP still missing 
#where true, contact authors

#alternatively, we ask directly for VIGE and VP since they may be useful for other points too (e.g. CVs, mean std variances, etc)
data$missingVP<-is.na (data$Total_v_phen)
data$missingVIGE<-is.na (data$V_ige)


#1.	B) Are there differences between: (i.e. moderators for the meta analysis on point 1A)
#this is all already done

#2.      What's the relative and absolute importance of direct and indirect genetic effects?
#(i.e. compare VDGE to the VIGEs)
 
#This can be done i.e. by
#1. calculating evolvability (CVa, where CVa = (sqrt(va)/m)), and the respective for IGEs (CViges) 
#2. or mean standardized variances, i.e. Va/mean^2 and Vige/mean^2
#3. in addition to comparing the ratios  va/vp and vige/vp

#parameters we need:
#Va, V_ige, VP 
#Trait mean 



sum(is.na (data$Trait_mean)) #38% of traits don't have a mean provided (75) 
data$missingMean<-is.na (data$Trait_mean)
#where true, contact authors

#knowing also that we have missing values in vige and va, again some could be backtransformed from h2 and social_h2 (see above)?
#we would have to be sure at least for how h2 is calculated
sum(is.na (data$V_ige)) #12 this we should have covered above 
sum(is.na (data$Va))  #19

#it seems at least for VA we should ask directly
data$missingVA<-is.na (data$Va)


#3.    A) Do IGEs change evolutionary trajectories?
#This is done by comparing the t2 (VTBV/Vp) estimated with IGEs and the DGE-IGE correlation vs h2 without IGEs. So to extract VTBV  (and Vp ),  t2 and h2. **TBV = total breeding value
#parameters we need:
#use T2, when available, or estimate if VTBV and VP are presented both (if reliable)
sum(is.na(data$T2))
#84 missing T2 
sum(is.na(data$V_tbv))
#93 missing VTBV 

#compare with h2
sum(is.na(data$H2))
#116 NA, again we could ask for the data since back-transformation didn't seem reliable from va


#3. 	B)  More specifically, do IGEs constrain or speed up the evolution of traits (through their correlation with DGEs)? 
#parameters we need: 
#here just use R_a_ige when provided and get an overall estimate, don't ask
sum(is.na (data$R_a_ige))
#74


#missing id with records (we need to ask for sample size)
sum(is.na (data$N_id_w_records)) #10
data$missingIdRecords<-(is.na (data$N_id_w_records))


#writing this file to check
write.csv(data, file = 'data/fulldataset.csv')


###################################################

#building dataset for contacting the authors


#we had planned to ask social h2, mean trait and # of records when we thought we could extract VP by back-transforming
setDT(data)
authors <- data[, .(Paper_id, Record_id, Trait_name, Excluded_during_data_extraction, missingSocial_h2, missingMean, missingIdRecords, missingVA, missingVIGE, missingVP)]
names(authors)

paperslist<-fread("data/animal_model_papers_to_be_assigned_full.csv")
names(paperslist)[names(paperslist) == "studyID" ] <- "Paper_id"

totlist<-merge(x=authors,y=paperslist, by="Paper_id", all.x=T)
str(totlist)
#View(totlist)

data2<-totlist
#combine column of missing data and remove those with all records
data2$combo<-paste(data2$missingSocial_h2, data2$missingMean, data2$missingIdRecords, sep = "_")
data2<-data2[!data2$combo == "FALSE_FALSE_FALSE"]

str(data2)
#100 obs
data2$Paper_id<-as.factor(data2$Paper_id)
levels(data2$Paper_id)
#32

#write.csv(data2, file = 'authorscontact.csv')


#if we would ask vp, va, vige directly: 

data3<-totlist
#combine column of missing data and remove those with all records
data3$combo<-paste(data3$missingVP, data3$missingVA, data3$missingVIGE, data3$missingMean, data3$missingIdRecords, sep = "_")
data3<-data3[!data3$combo == "FALSE_FALSE_FALSE_FALSE_FALSE"]

str(data3)
#129
data3$Paper_id<-as.factor(data3$Paper_id)
levels(data3$Paper_id)
#41


names(data3)
data3<-subset(data3, select =  -c(combo, Excluded_during_data_extraction, missingSocial_h2))

#write.csv(data3, file = 'authorscontact2.csv')

#we ended up using this but only asking for VP, # of records, trait mean



#################################
#################################
#################################
#Code for filling in data from the authors contacted


dataset = read.csv("data/fulldataset.csv", header=T)
summary(dataset)

#IGE0455
dataset$Trait_mean[dataset$Paper_id == "IGE0455" & dataset$Trait_name == "age at first birth (AFB)"] = 25.57
dataset$Total_v_phen[dataset$Paper_id == "IGE0455" & dataset$Trait_name == "age at first birth (AFB)"] = 23.81

#IGE0857
dataset$Trait_mean[dataset$Paper_id == "IGE0857" & dataset$Trait_name == "feeding rates"] = 2.6

#IGE0253
dataset$Trait_mean[dataset$Paper_id == "IGE0253" & dataset$Trait_name == "dyadic fights (fighting ability)"] = 23.0091
dataset$Total_v_phen[dataset$Paper_id == "IGE0253" & dataset$Trait_name == "dyadic fights (fighting ability)"] = 11.2414

#IGE0502
dataset$Trait_mean[dataset$Paper_id == "IGE0502" & dataset$Trait_name == "social dominance"] = 0.4365 ##### although we had a value already assigned to the Trait_mean of this study (0.5), the correct value is indeed 0.4365
dataset$Total_v_phen[dataset$Paper_id == "IGE0502" & dataset$Trait_name == "social dominance"] = 0.2460

#IGE0572 (note author supplied updated values for all variance components, which we are using to maintain relevance to the total VP)
dataset$Total_v_phen[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 1"] = 47.3456
dataset$Total_v_phen[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 2"] = 54.2079
dataset$Total_v_phen[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 3"] = 61.1060
dataset$Total_v_phen[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 4"] = 68.4242
dataset$Total_v_phen[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 5"] = 79.3523

dataset$Va[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 1"] = 7.8753
dataset$Va[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 2"] = 12.9368
dataset$Va[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 3"] = 13.3603
dataset$Va[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 4"] = 10.0450
dataset$Va[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 5"] = 5.7880

dataset$V_ige[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 1"] = 0.3939
dataset$V_ige[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 2"] = 0.2360
dataset$V_ige[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 3"] = 0.2034
dataset$V_ige[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 4"] = 0.1899
dataset$V_ige[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 5"] = 0.2290

dataset$Cov_a_ige[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 1"] = -1.4437
dataset$Cov_a_ige[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 2"] = -1.5163
dataset$Cov_a_ige[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 3"] = -1.5430
dataset$Cov_a_ige[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 4"] = -1.2356
dataset$Cov_a_ige[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 5"] = -0.8612

#V other 1 is litter, V other 2 is ""pseudo-environmental variance"", V other 3 is group

dataset$V_other_1[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 1"] = 8.6388
dataset$V_other_1[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 2"] = 5.8548
dataset$V_other_1[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 3"] = 4.3509
dataset$V_other_1[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 4"] = 3.7080
dataset$V_other_1[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 5"] = 3.7664

dataset$V_other_2[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 1"] = 23.0578
dataset$V_other_2[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 2"] = 28.1216
dataset$V_other_2[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 3"] = 34.3635
dataset$V_other_2[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 4"] = 42.4814
dataset$V_other_2[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 5"] = 52.6743

dataset$V_other_3[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 1"] = 5.6029
dataset$V_other_3[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 2"] = 7.4531
dataset$V_other_3[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 3"] = 9.5971
dataset$V_other_3[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 4"] = 12.3518
dataset$V_other_3[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 5"] = 15.9109

#variance ratios, based on straight forward division but adding any way

dataset$H2[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 1"] = 0.1663
dataset$H2[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 2"] = 0.2387
dataset$H2[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 3"] = 0.2186
dataset$H2[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 4"] = 0.1468
dataset$H2[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 5"] = 0.0729

dataset$Social_h2[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 1"] = 0.0083
dataset$Social_h2[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 2"] = 0.0004
dataset$Social_h2[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 3"] = 0.0033
dataset$Social_h2[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 4"] = 0.0028
dataset$Social_h2[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 5"] = 0.0029

dataset$T2[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 1"] = 0.1471
dataset$T2[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 2"] = 0.0604
dataset$T2[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 3"] = 0.0282
dataset$T2[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 4"] = 0.0300
dataset$T2[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 5"] = 0.0624

#IGE0348
dataset$Trait_mean[dataset$Paper_id == "IGE0348" & dataset$Trait_name == "female liability for extra-pair reproduction"] = 0.285

#IGE0199
dataset$Total_v_phen[dataset$Paper_id == "IGE0199" & dataset$Trait_name == "egg length"] = 0.224
dataset$Total_v_phen[dataset$Paper_id == "IGE0199" & dataset$Trait_name == "egg number"] = 6.757
dataset$Total_v_phen[dataset$Paper_id == "IGE0199" & dataset$Trait_name == "egg pod length"] = 2.643
dataset$Total_v_phen[dataset$Paper_id == "IGE0199" & dataset$Trait_name == "egg pod number"] = 0.04
dataset$Total_v_phen[dataset$Paper_id == "IGE0199" & dataset$Trait_name == "latency to first egg pod"] = 0.599

#IGE0389
dataset$N_id_w_records[dataset$Paper_id == "IGE0389" & dataset$Trait_name == "Diameter at breast height (2 years)"] = 9107
dataset$N_id_w_records[dataset$Paper_id == "IGE0389" & dataset$Trait_name == "Diameter at breast height (4 years)"] = 8913
dataset$N_id_w_records[dataset$Paper_id == "IGE0389" & dataset$Trait_name == "Mycosphaerella leaf disease"] = 9107

#updating this as the trait was analysed after an arcsine transformation, which is then 48.03 
dataset$Trait_mean[dataset$Paper_id == "IGE0389" & dataset$Trait_name == "Mycosphaerella leaf disease"] = 48.03

#IGE0418 (note same dataset as above)
dataset$N_id_w_records[dataset$Paper_id == "IGE0418" & dataset$Trait_name == "Over-bark diameter at breast height at 2 years"] = 9107
dataset$N_id_w_records[dataset$Paper_id == "IGE0418" & dataset$Trait_name == "Over-bark diameter at breast height at 4 years"] = 8913
dataset$N_id_w_records[dataset$Paper_id == "IGE0418" & dataset$Trait_name == "Over-bark diameter at breast height at 8 years"] = 8767
dataset$N_id_w_records[dataset$Paper_id == "IGE0418" & dataset$Trait_name == "percentage of juvenile foliage affected by Mycosphaerella leaf disease"] = 9107

dataset$Trait_mean[dataset$Paper_id == "IGE0418" & dataset$Trait_name == "Over-bark diameter at breast height at 4 years"] = 135.1
dataset$Trait_mean[dataset$Paper_id == "IGE0418" & dataset$Trait_name == "percentage of juvenile foliage affected by Mycosphaerella leaf disease"] = 48.03 

dataset$Total_v_phen[dataset$Paper_id == "IGE0418" & dataset$Trait_name == "Over-bark diameter at breast height at 2 years"] = 245.6
dataset$Total_v_phen[dataset$Paper_id == "IGE0418" & dataset$Trait_name == "Over-bark diameter at breast height at 4 years"] = 524.0
dataset$Total_v_phen[dataset$Paper_id == "IGE0418" & dataset$Trait_name == "Over-bark diameter at breast height at 8 years"] = 1464.2
dataset$Total_v_phen[dataset$Paper_id == "IGE0418" & dataset$Trait_name == "percentage of juvenile foliage affected by Mycosphaerella leaf disease"] = 188.6

#IGE0537 (note traits were variance standardised, so VP is effectively 1) ##### we are adding these values for completeness, but note that VP was standardized so, we will actually not use these values but rather assume VP = 1 by using the variable "Total_v_phen2" created below
dataset$Total_v_phen[dataset$Paper_id == "IGE0537" & 
                       dataset$Trait_name == "average daily gain (ADG)" & 
                       dataset$Treatment_group == "full feeding (F)"] = 123.8874

dataset$Total_v_phen[dataset$Paper_id == "IGE0537" & 
                       dataset$Trait_name == "average daily gain (ADG)" & 
                       dataset$Treatment_group == "restricted feeding (R)"] = 170.6174

#IGE0509
dataset$Trait_mean[dataset$Paper_id == "IGE0509" & dataset$Trait_name == "laying date"] = 306.2

#IGE0483
dataset$Trait_mean[dataset$Paper_id == "IGE0483" & 
                     dataset$Trait_name == "Aggression" & 
                     dataset$Treatment_group == "nymph low-density"] = 2.5065

dataset$Trait_mean[dataset$Paper_id == "IGE0483" & 
                     dataset$Trait_name == "Aggression" & 
                     dataset$Treatment_group == "nymph high-density"] = 2.4942

dataset$Trait_mean[dataset$Paper_id == "IGE0483" & 
                     dataset$Trait_name == "Aggression" & 
                     dataset$Treatment_group == "Young adult low-density"] = 2.9130

dataset$Trait_mean[dataset$Paper_id == "IGE0483" & 
                     dataset$Trait_name == "Aggression" & 
                     dataset$Treatment_group == "Young adult high-density"] = 2.9269

dataset$Trait_mean[dataset$Paper_id == "IGE0483" & 
                     dataset$Trait_name == "Aggression" & 
                     dataset$Treatment_group == "old adult low-density"] = 2.6466

dataset$Trait_mean[dataset$Paper_id == "IGE0483" & 
                     dataset$Trait_name == "Aggression" & 
                     dataset$Treatment_group == "old adult high-density"] = 2.5545

dataset$Total_v_phen[dataset$Paper_id == "IGE0483" & 
                       dataset$Trait_name == "Aggression" & 
                       dataset$Treatment_group == "nymph low-density"] = 0.6042

dataset$Total_v_phen[dataset$Paper_id == "IGE0483" & 
                       dataset$Trait_name == "Aggression" & 
                       dataset$Treatment_group == "nymph high-density"] = 0.5376

dataset$Total_v_phen[dataset$Paper_id == "IGE0483" & 
                       dataset$Trait_name == "Aggression" & 
                       dataset$Treatment_group == "Young adult low-density"] = 0.5210

dataset$Total_v_phen[dataset$Paper_id == "IGE0483" & 
                       dataset$Trait_name == "Aggression" & 
                       dataset$Treatment_group == "Young adult high-density"] = 0.4705

dataset$Total_v_phen[dataset$Paper_id == "IGE0483" & 
                       dataset$Trait_name == "Aggression" & 
                       dataset$Treatment_group == "old adult low-density"] = 0.7853

dataset$Total_v_phen[dataset$Paper_id == "IGE0483" & 
                       dataset$Trait_name == "Aggression" & 
                       dataset$Treatment_group == "old adult high-density"] = 0.8452

#IGE0292
dataset$Total_v_phen[dataset$Paper_id == "IGE0292" & dataset$Trait_name == "neck bite marks" ] = 4.4754
dataset$Total_v_phen[dataset$Paper_id == "IGE0292" & dataset$Trait_name == "body bite marks" ] = 5.9418
dataset$Total_v_phen[dataset$Paper_id == "IGE0292" & dataset$Trait_name == "tail bite marks" ] = 6.7858

#IGE0501
dataset$Total_v_phen[dataset$Paper_id == "IGE0501" & dataset$Trait_name == "approach rate" ] = 0.0775
dataset$Total_v_phen[dataset$Paper_id == "IGE0501" & dataset$Trait_name == "naso-anal contact rate" ] = 0.0471
dataset$Total_v_phen[dataset$Paper_id == "IGE0501" & dataset$Trait_name == "mounting rate" ] = 0.0793
dataset$Total_v_phen[dataset$Paper_id == "IGE0501" & dataset$Trait_name == "rearing rate" ] =  0.0707
dataset$Total_v_phen[dataset$Paper_id == "IGE0501" & dataset$Trait_name == "reciprocal latency to fight" ] = 0.1968

dataset$Trait_mean[dataset$Paper_id == "IGE0501" & dataset$Trait_name == "approach rate" ] = 0.4925
dataset$Trait_mean[dataset$Paper_id == "IGE0501" & dataset$Trait_name == "naso-anal contact rate" ] = 0.1820
dataset$Trait_mean[dataset$Paper_id == "IGE0501" & dataset$Trait_name == "mounting rate" ] = 0.1660
dataset$Trait_mean[dataset$Paper_id == "IGE0501" & dataset$Trait_name == "rearing rate" ] =  0.2461
dataset$Trait_mean[dataset$Paper_id == "IGE0501" & dataset$Trait_name == "reciprocal latency to fight" ] = -2.2224

#IGE0516
dataset$Total_v_phen[dataset$Paper_id == "IGE0516" & dataset$Trait_name == "Social dominance" ] = 0.25 
dataset$Trait_mean[dataset$Paper_id == "IGE0516" & dataset$Trait_name == "Social dominance" ] = 0.5  

#IGE0196
dataset$Total_v_phen[dataset$Paper_id == "IGE0196" & dataset$Trait_name == "breeding date" ] = 165.2161

#IGE0517
dataset$Total_v_phen[dataset$Paper_id == "IGE0517" & dataset$Trait_name == "divorce" ] = 0.0857
dataset$Trait_mean[dataset$Paper_id == "IGE0517" & dataset$Trait_name == "divorce" ] = 0.0947

#IGE494
dataset$Total_v_phen[dataset$Paper_id == "IGE0494" & dataset$Trait_name == "average daily gain" ] = 0.75

#IGE0746
dataset$Trait_mean[dataset$Paper_id == "IGE0746" & 
                     dataset$Trait_name == "neck feather condition score" & 
                     dataset$Treatment_group == "line W1"] = 1.38

dataset$Trait_mean[dataset$Paper_id == "IGE0746" & 
                     dataset$Trait_name == "back feather condition score" & 
                     dataset$Treatment_group == "line W1"] = 1.15

dataset$Trait_mean[dataset$Paper_id == "IGE0746" & 
                     dataset$Trait_name == "rump feather condition score" & 
                     dataset$Treatment_group == "line W1"] = 1.14

dataset$Trait_mean[dataset$Paper_id == "IGE0746" & 
                     dataset$Trait_name == "belly feather condition score" & 
                     dataset$Treatment_group == "line W1"] = 1.31

dataset$Trait_mean[dataset$Paper_id == "IGE0746" & 
                     dataset$Trait_name == "neck feather condition score" & 
                     dataset$Treatment_group == "line WB"] = 2.30

dataset$Trait_mean[dataset$Paper_id == "IGE0746" & 
                     dataset$Trait_name == "back feather condition score" & 
                     dataset$Treatment_group == "line WB"] = 1.63

dataset$Trait_mean[dataset$Paper_id == "IGE0746" & 
                     dataset$Trait_name == "rump feather condition score" & 
                     dataset$Treatment_group ==  "line WB"] = 2.11

dataset$Trait_mean[dataset$Paper_id == "IGE0746" & 
                     dataset$Trait_name == "belly feather condition score" & 
                     dataset$Treatment_group == "line WB"] = 1.60

#IGE0493
dataset$Total_v_phen[dataset$Paper_id == "IGE0493" & dataset$Trait_name == "old cone hoard size" ] = 191371304
dataset$Total_v_phen[dataset$Paper_id == "IGE0493" & dataset$Trait_name == "new cone hoard size" ] =  314210694
dataset$Total_v_phen[dataset$Paper_id == "IGE0493" & dataset$Trait_name == "lifetime reproductive success (LRS)" ] = 3.982

dataset$Trait_mean[dataset$Paper_id == "IGE0493" & dataset$Trait_name == "old cone hoard size" ] = 6969.252
dataset$Trait_mean[dataset$Paper_id == "IGE0493" & dataset$Trait_name == "new cone hoard size" ] = 9937.322
dataset$Trait_mean[dataset$Paper_id == "IGE0493" & dataset$Trait_name == "lifetime reproductive success (LRS)" ] = 1.384

#could also add distribution specific variance (as V other 2) for these 3 traits as they are Poisson dist:
dataset$V_other_2[dataset$Paper_id == "IGE0493" & dataset$Trait_name == "old cone hoard size" ] = 0.000997
dataset$V_other_2[dataset$Paper_id == "IGE0493" & dataset$Trait_name == "new cone hoard size" ] = 0.150
dataset$V_other_2[dataset$Paper_id == "IGE0493" & dataset$Trait_name == "lifetime reproductive success (LRS)" ] = 4.947

#IGE0852
dataset$Trait_mean[dataset$Paper_id == "IGE0852" & dataset$Trait_name == "parturition date" ] = 5.13 


#IGE0144 
dataset$Va[dataset$Paper_id == "IGE0144" & dataset$Trait_name == "carcass weight (CW)"] = 0.34

dataset$V_ige[dataset$Paper_id == "IGE0144" & dataset$Trait_name == "carcass weight (CW)"] = 0.007

dataset$Variance_standardized[dataset$Paper_id == "IGE0144" & dataset$Trait_name == "carcass weight (CW)"] = "yes" ##### "Yes" substituted by "yes"

#Note for traits that are variance standardised VP is effectively 1. 

#IGE0824
#Trait "oviposition pre mating" is allegedly missing V_IGE (and social h2), 
#but this trait was expressed prior to the social interaction, so should have no IGE anyway.
#Suggest removing trait from analysis:
dataset$Excluded_during_data_extraction[dataset$Paper_id == "IGE0824" & dataset$Trait_name == "oviposition pre mating"] = "yes" ##### "Yes" substituted by "yes"

##### let's then remove it from the dataset straight away
dataset <- dataset[dataset$Excluded_during_data_extraction=="no",]


#IGE0458
dataset$Trait_mean[dataset$Paper_id == "IGE0458" & dataset$Trait_name == "nesting site preference" ] = 0.497
dataset$Total_v_phen[dataset$Paper_id == "IGE0458" & dataset$Trait_name == "nesting site preference" ] = 0.25
dataset$N_id_w_records[dataset$Paper_id == "IGE0458" & dataset$Trait_name == "nesting site preference" ] = 1537

#JM [original author] said: "1005 pairs for which we identified both the male and the female. Among those, there were 791 unique females, and 746 unique males." so have summed n males and n females to get n unique indivs

##################################

#correcting one record after comments inspection
#h2 set to NA in IGE 0197 following comments 
dataset$H2[dataset$Paper_id == "IGE0197"] = NA 


# adding VP = 1 if trait is variance standardized, else copy value of reported VP
dataset$Total_v_phen2<-ifelse(dataset$Variance_standardized == "yes", 1, dataset$Total_v_phen) ##### "Yes" substituted by "yes"

# adding mean = 0 if trait is mean standardized, else copy value of reported mean
dataset$Trait_mean2<-ifelse(dataset$Mean_standardized == "yes", 0, dataset$Trait_mean)

#we should use these Total_v_phen2 and Trait_mean2 as VP and mean, or eventually replace the original ones
head(dataset)

#removing first column
dataset$X<-NULL

#here we can decide to remove other columns that are not important like the columns with "missing" that we used to figure out the authors contacs
dataset$Excluded_during_data_extraction<-NULL #this is all no now

#Final_comments can be removed l (it was just a summary of first and second screener comments, but not all screeners made comments so it is heterogeneous)
dataset$Final_comments<-NULL


### last checks after comparing extracted h2 and social h2 with  V_ige/Total_v_phen2 and Va/Total_v_phen2
#there are discrepancies in 14 papers

#from here, we decide to use what they report instead of calculating the ratios ourselves in some cases, going through them one by one
#I add this variable to flag studies where we should keep their values instead

dataset$Reported_values_flag<-NA

#a few of the cases could not be resolved so we decide to discard the values for that study
#it is not the same as excluded during data extraction 

dataset$Excluded_after_control<-NA

#in the following 9 studies it's where we contacted the authors
#"IGE0418" "IGE0199" "IGE0253" "IGE0493" "IGE0144" "IGE0572" "IGE0502" "IGE0458" "IGE0517"

#IGE0418
#the disagreement in h2 is 2% so I wouldn’t worry about it.

#IGE0199
#disagreement in h2 is 5% and for social h2 0% so fine

#IGE0253
#Summing their variance components gives 1.613, which used as the VP gives the h2 of 0.009 they report. 
#Authors told us VP was 11.24 which is way out
#we decide to use the VP summing up the variance components as it gives the right h2 reported in the paper and override what the authors communicated

dataset$Total_v_phen[dataset$Paper_id == "IGE0253" ] = 1.613

#"IGE0493" 
#used a Poisson error structure with log-link function for old cones
#so VP is on original scale but variance estimates are not
#use the h2 reported in the paper for all traits, we need to flag this
#we can’t divide the VIGE by trait mean or variance to get social h2 or CV IGE due to the differing scales, how to flag this?
#Set to VP2 to NA

dataset$Reported_values_flag[dataset$Paper_id == "IGE0493" ] = "YES"
dataset$Total_v_phen2[dataset$Paper_id == "IGE0493" ] = NA


#"IGE0144"
#typo  in R code, VIGE should be 0.007 not 0.07
dataset$V_ige[dataset$Paper_id == "IGE0144" & dataset$Trait_name == "carcass weight (CW)" ] = 0.007

#"IGE0572"
#IGE0572 I did not spot that they multiplied the VIGE by 49 (square of number of group mates) to calculate their social h2. 
#Since we haven’t done this with any other papers, we should not here. So update (matches our calculated values): 
  
dataset$social_h2[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 1"] = 0.0083
dataset$social_h2[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 2"] = 0.0004
dataset$social_h2[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 3"] = 0.0033
dataset$social_h2[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 4"] = 0.0028
dataset$social_h2[dataset$Paper_id == "IGE0572" & dataset$Trait_name == "Average daily gain week 5"] = 0.0029

#"IGE0502" 
#seems similar issue to 493, where maybe we have discrepancies between scales?
#also in this case we should use what they report and not our ratio
#again if we add variable, flag this study
dataset$Reported_values_flag[dataset$Paper_id == "IGE0502" ] = "YES"

#"IGE0458" 
#not clear what to use, discrepancies between pdf and online table for values
# Binary model may be causing issues. Told by authors VP was 0.25 but they use a VR of 10 in the models, which goes some way to explaining why they get a much smaller h2 and social h2 than us. 
#Values in dataset do not match  Table 3 in paper 
#May be issue with non gaussian traits

#I didn't do anything here, we can try to extract the data from the pdf instead, or just decide to use the values they report
#and so this as a study where we don't recalculate things ourselves

#"IGE0517"
#binary model may be causing issues, told us VP was 0.086 but report VA of 0.07. Their h2 of 0.04 must be correct, or more correct than ours.
#again if we add variable, flag this study

dataset$Reported_values_flag[dataset$Paper_id == "IGE0517" ] = "YES"

### 5 remaining studies
#0243, 1050, 0920, 0746, 1038 

#"IGE0243" 
#diff in h2 is 4% and perhaps down to rounding so not worth worrying about.

#"IGE1050" –
#the lower est in our calculated values could be due to fixed effects as of focal and partner being included – 
#these mean VP is higher than summing the rest of the variance components.
#I summed the other variance components and using that as VP gives the same h2 as the authors report, 
#so I think that is what is going on.


#"IGE0920" –
#was a mistake in our data extraction, we extracted h2 from a model without IGE in the paper
#we should set to NA and use one we calculate for all traits in that paper

dataset$H2[dataset$Paper_id == "IGE0920"] = NA 

# IGE0746
#No idea how they get such low values. For neck feathers in the W1 line for example they report VA clearly as 0.135, and VP as 0.306 (Table S2),
#and they say h2 is Va / VP in methods (eq. 3), so our calc of 0.44 should be correct. 
#But they report 0.059 in Table 3, which is way off and not even something simple like a change from % to proportion and so x10 or x100 …
#so not sure we understand where values come from, either we accept what authors report or we discard study

#Decided to accept what authors report

# "IGE 1038"
# they report both “phenotypic variance” [VP] and “total phenotypic variance” [TVP; always much larger than VP] – they use VP for the denominator in their h2, but TVP as the denominator for their T2. 
#We extracted TVP as our VP, hence why we calculated a lower h2 than they did. 
#As noted in the comments by Alfredo during extraction, it’s not clear what TVP is – they don’t ever describe it in the paper.
#I would replace the VP values we have with the VPs they give (Table 1)
#and then ignore their T2, we can calculate it ourselves from the TBV and VP they give.

dataset$Total_v_phen[dataset$Paper_id == "IGE1038" & dataset$Trait_name == "average daily gain (ADG)"] = 5335.67
dataset$Total_v_phen[dataset$Paper_id == "IGE1038" & dataset$Trait_name == "days to 100 kg (D100)"] = 631.40
dataset$Total_v_phen[dataset$Paper_id == "IGE1038" & dataset$Trait_name == "backfat thickness to 100kg (B100)"] = 292.49
dataset$Total_v_phen[dataset$Paper_id == "IGE1038" & dataset$Trait_name == "average daily feed intake (ADFI)"] = 38278.14
dataset$Total_v_phen[dataset$Paper_id == "IGE1038" & dataset$Trait_name == "feed conversion ratio (FCR)"] = 0.0029
dataset$Total_v_phen[dataset$Paper_id == "IGE1038" & dataset$Trait_name == "residual feed intake (RFI)"] = 25993.18




#################

write.csv(dataset, file = 'data/dataset_final_after_cleaning_and_adding_author_contact_FS_MM.csv') 


