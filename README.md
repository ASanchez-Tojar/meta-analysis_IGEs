DOI: to be added via the Zenodo integration upon acceptance of this work in a scientific journal.

# Meta-analysis of indirect genetic effects

This repository contains the data, code and other materials used in the following study:

---

Francesca Santostefano*, Maria Moiron*, Alfredo Sánchez-Tójar*, David N. Fisher*. 2024. *Indirect genetic effects increase the heritable variation available to selection and are largest for behaviours: a meta-analysis*. Submitted. Preprint available [here](LINK)

*All authors contributed equally to this work

---

The repository consists of an Rproject ([meta-analysis_IGEs.Rproj](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/meta-analysis_IGEs.Rproj)) with 6 scripts and 3 folders containing the data either needed to run these scripts or created by the scripts, including the figures. For a detailed description of all the files, please see below. For any further information about this repository, please contact: [Alfredo Sánchez-Tójar](https://scholar.google.co.uk/citations?hl=en&user=Sh-Rjq8AAAAJ&view_op=list_works&sortby=pubdate), email: alfredo.tojar@gmail.com. 

## Information about scripts, folders and files within:

Scripts:
-	[001_systematic_review.R ](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/001_systematic_review.R): used to import the results of the systematic literature review conducted in Web of Science and Scopus.
  
    * Input: Web of Science and Scopus keyword search, and the 8 .bib files corresponding to the snowballing search
      - Web of Science keyword search results: [WoS_IGEs_search.bib](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/literature_review/systematic_search/WoS_IGEs_search.bib)
      - Scopus keyword search results: [Scopus_IGEs_search.bib](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/literature_review/systematic_search/Scopus_IGEs_search.bib)
      - Snowballing search results: 8 .bib files found [here](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/literature_review/snowballing)
      - [Deduplicated list of unique references](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/literature_review/search_unique_references_extracted.csv) to be formatted for importing into [Rayyan](https://rayyan.qcri.org/)
    
    * Output:
      - [Preliminary list of unique references](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/literature_review/search_unique_references_extracted.csv) for further deduplication outside R.
      - [Deduplicated full reference list](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/literature_review/IGE_search_unique_references_rayyan.csv) to be imported into [Rayyan](https://rayyan.qcri.org/) for the title-and-abstract screening
      - [R session information](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/literature_review/deduplicating_Rpackages_session.txt) for reproducibility purposes
  
-	[002_fulltext_screening.R ](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/002_fulltext_screening.R): used to import the results of the title-and-abstract screening conducted using [Rayyan](https://rayyan.qcri.org/) and to generate the database needed to proceed with the full-text screening phase.
  
    * Input:
      - [List of unique reference list with title-and-abstract decisions](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/literature_review/title_and_abstract_screening/title-and-abstract_decisions_rayyan_studyID.csv)
    
    * Output:
      - Four full-text screening subsets (as .xlsx files) including the list of references assigned to each of the 4 observers. Files available [here](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/literature_review/fulltext_screening)
      - [R session information](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/literature_review/fulltext_screening/fulltext_templates_Rpackages_session.txt) for reproducibility purposes
      
-	[003_data_extraction.R ](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/003_data_extraction.R): used to import the results of the full-text screening and to generate the databases needed to proceed to the data extraction phase.
  
    * Input:
      - [List of unique reference list with title-and-abstract decisions](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/literature_review/title_and_abstract_screening/title-and-abstract_decisions_rayyan_studyID.csv)
      - [Final list of studies included after full-text screening](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/literature_review/fulltext_screening/Final_fulltext_screening_responses_including_conflict_resolution_google_form_data.xlsx)
      - [Final list of studies included after full-text screening and using animal models](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/literature_review/data_extraction/animal_model_papers_to_be_assigned_20200703.csv)
    
    * Output:
      - Four data-extraction subsets (as .xlsx files) including the list of references assigned to each of the 4 observers. Files available [here](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/literature_review/data_extraction)
      - Four data double-extraction subsets (as .xlsx files) including the list of references assigned to each of the 4 observers. Files available [here](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/literature_review/data_extraction/double-checking)
        
-	[004_dataset_cleaning.R ](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/004_dataset_cleaning.R): used to clean the meta-analytic dataset and generate the version for the analyses that is imported by [005_data_analysis.r](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/005_data_analysis.r).
  
    * Input:
      - Four datasets with the full data extraction of each of the 4 observers. The four files are .csv and their name starts with '*Data from papers*'. All available [here]([https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/literature_review/data_extraction](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/data))
      - [Final list of studies included after full-text screening and using animal models](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/data/animal_model_papers_to_be_assigned_full.csv)
      - [List of references created by this script before contacting authors, edited outside R, and imported back to be filled after having contacted the authors](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/data/fulldataset.csv)
    
    * Output:
      - [List of references before contacting authors](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/data/fulldataset.csv)
      - [Full list of references with data, including data obtained by author correspondence](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/data/dataset_final_after_cleaning_and_adding_author_contact_FS_MM.csv). This list will be imported by [005_data_analysis.r](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/005_data_analysis.r) for the analyses
        
-	[005_data_analysis.r ](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/005_data_analysis.r): used to perform all the analyses.
  
    * Input: 
      - [Final dataset for analyses](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/data/dataset_final_after_cleaning_and_adding_author_contact_FS_MM.csv): the metadata providing a description for each of the variables is available [here](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/data/dataset_final_after_cleaning_and_adding_author_contact_FS_MM_METADATA.txt)
      - [Taxonomic data extracted from the Open Tree of Life](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/data/taxa_Open_Tree_of_Life.RData) created by this script and then imported for reproducibility purposes
      - [Tree file](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/data/tree.Rdata) created by this script and then imported for reproducibility purposes
      - [Phylogenetic correlation matrix](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/data/phylo_cor.Rdata) created by this script and then imported for reproducibility purposes
    
    * Output:
      - [Taxonomic data extracted from the Open Tree of Life](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/data/taxa_Open_Tree_of_Life.RData) using the function 'tnrs_match_names()'
      - [Tree file](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/data/tree.Rdata)
      - [Phylogenetic correlation matrix](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/data/phylo_cor.Rdata)
      - Meta-analytic and meta-regression models. All available as .RData [here](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/data/models) and will be imported by [006_figures.r](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/006_figures.r) for generating the figures.
      - Model-specific datasets. All available as .csv [here](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/data/subsets) and will be imported by [006_figures.r](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/006_figures.r) for generating the figures.
        
-	[006_figures.r ](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/006_figures.r): used to generate the all figures included in the main text and the supplementary materials (except Figure 1 from the main text).
  
    * Input:
      - Meta-analytic and meta-regression models from [here](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/data/models)
      - Model-specific datasets from [here](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/data/subsets)
    
    * Output:
      - [All figures](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/figures) included in the main text and the supplementary materials (except Figure 1 from the main text)

Folders:
-	[literature_review](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/literature_review):
    * [systematic_search](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/literature_review/systematic_search): contains 2 .bib files imported by [001_systematic_review.R](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/001_systematic_review.R)
    * [snowballing](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/literature_review/snowballing): contains 8 .bib files imported by [001_systematic_review.R](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/001_systematic_review.R)
    * [title_and_abstract_screening](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/literature_review/title_and_abstract_screening): containing the [list of unique reference list with title-and-abstract decisions](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/literature_review/title_and_abstract_screening/title-and-abstract_decisions_rayyan_studyID.csv)
    * [fulltext_screening](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/literature_review/fulltext_screening): containing all files generated by [003_data_extraction.R](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/003_data_extraction.R) as well as the [final list of references with all the full-text screening decisions](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/literature_review/fulltext_screening/Final_fulltext_screening_responses_including_conflict_resolution_google_form_data.xlsx)
    * [data_extraction](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/literature_review/data_extraction): containing all files generated by [003_data_extraction.R](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/003_data_extraction.R)

- [data](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/data): contains all data files necessary for the following three scripts: [004_dataset_cleaning.R](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/004_dataset_cleaning.R), [005_data_analysis.r](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/005_data_analysis.r) and [006_figures.r](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/006_figures.r)

- [figures](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/figures): contains all figures generated by [006_figures.r](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/006_figures.r) and included in the main text and the supplementary materials (except Figure 1 from the main text)
