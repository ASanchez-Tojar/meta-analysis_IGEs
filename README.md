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
-	[001_systematic_review.R](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/001_systematic_review.R): used to import the results of the systematic literature review conducted in Web of Science and Scopus.
  
    * Input: Web of Science and Scopus keyword search, and the 8 .bib files corresponding to the snowballing search
      - Web of Science keyword search results: [WoS_IGEs_search.bib](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/literature_review/systematic_search/WoS_IGEs_search.bib)
      - Scopus keyword search results: [Scopus_IGEs_search.bib](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/literature_review/systematic_search/Scopus_IGEs_search.bib)
      - Snowballing search results: 8 .bib files found [here](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/literature_review/snowballing)
      - [Deduplicated list of unique references](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/literature_review/search_unique_references_extracted.csv) to be formatted for importing into [Rayyan](https://rayyan.qcri.org/)
    
    * Output:
      - [Preliminary list of unique references](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/literature_review/search_unique_references_extracted.csv) for further deduplication outside R.
      - [Deduplicated full reference list](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/literature_review/IGE_search_unique_references_rayyan.csv) to be imported into [Rayyan](https://rayyan.qcri.org/) for the title-and-abstract screening
      - [R session information](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/literature_review/deduplicating_Rpackages_session.txt) for reproducibility purposes
  
-	[002_fulltext_screening.R](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/002_fulltext_screening.R): used to import the results of the title-and-abstract screening conducted using [Rayyan](https://rayyan.qcri.org/) and to generate the database needed to proceed with the full-text screening phase.
-	[003_data_extraction.R](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/003_data_extraction.R): used to import the results of the full-text screening and to generate the databases needed to proceed to the data extraction phase.
-	[004_dataset_cleaning.R](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/004_dataset_cleaning.R): used to clean the meta-analytic dataset and generate the version for the analyses that is imported by [005_data_analysis.r](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/005_data_analysis.r).
-	[005_data_analysis.r](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/005_data_analysis.r): used to perform all the analyses.
-	[006_figures.r](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/006_figures.r): used to generate the all figures included in the main text and the supplementary materials (except Figure 1 from the main text).

Folders:
-	[literature_review](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/literature_review):
    * [systematic_search](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/literature_review/systematic_search): contains 2 .bib files imported by [001_systematic_review.R](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/001_systematic_review.R)
    * [snowballing](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/literature_review/snowballing): contains 8 .bib files imported by [001_systematic_review.R](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/001_systematic_review.R)
    * [title_and_abstract_screening](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/literature_review/title_and_abstract_screening):
    * [fulltext_screening](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/literature_review/fulltext_screening):
    * [data_extraction](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/literature_review/data_extraction):
- [data](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/data):
- [figures](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/tree/main/figures): contains all figures generated by [006_figures.r](https://github.com/ASanchez-Tojar/meta-analysis_IGEs/blob/main/006_figures.r) and included in the main text and the supplementary materials (except Figure 1 from the main text).
