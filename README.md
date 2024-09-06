# Proteomics-Reliability-Analysis
Pipeline for a novel way for the imputation of missing values in LC-MS/MS

Example study: 4 biological replicates, 6 timepoints
Each protein was classified for each timepoint as reliably or unreliably detected or undetected depending on its number of not-assigned values (NAs), and the number of NAs of its immediately adjacent days (neighbors). 
A day for a specific protein is considered as a supporting neighbor if it has 0, 1 or 2 NAs (out of 4 replicates). 
A day for a specific protein is considered as an unsopporting neighbor if it has 3 or 4 NAs (out of 4 replicates).

In our study: 
Days 0 and 5 were considered as Reliably Undetected when all replicates were NAs, and days 1 – 4, besides that, must had at least one unsupporting neighbor (with 3-4 NAs). NAs were replaced by the minimum of detection of the dataset (Deterministic Minimum Imputation method (Meleth, Deshane, and Kim 2005)). 
Days with one or no NAs were defined as Reliably Detected and their abundance values were kept. 
Finally, days with two or more NAs were classified as Unreliably Detected when they had at least one supporting neighbor (with 0, 1 or 2 NAs), keeping their quantification values; otherwise, they were classified as Unreliably Undetected, and its quantification values were replaced by NAs in all replicates. 

All those proteins which were Reliably or Unreliably Undetected in every timepoint were discarded. The remaining NA values were estimated by k-Nearest Neighbors (kNN) imputation (k = 10) (Troyanskaya et al. 2001) 
