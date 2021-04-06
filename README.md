# MVPAemoDys_Analyses
 
R scripts and preprocessed data to reproduce the mega-analysis in "Affective neural signatures do not distinguish women with emotion dysregulation from healthy controls: A mega-analysis across three task-based fMRI studies" by Sicorello et al. (2021). For a preprint, visit: https://www.medrxiv.org/content/10.1101/2021.02.03.21251077v1

-Main analysis scripts are:
1) within-person analyses: 'AffDysPattern_MainAnalysisWithin_AnalysisScript.R'
2) between-person analyses: 'AffDysPattern_MainAnalysisBetween_AnalysisScript.R'

-'AffDysPattern_functions' contains some data processing helper functions

-'AffDysPattern_reliability' was used to calculate internal consistency, as reported in the methods section

MCMC results of stan models were too large and led to error messages on upload. If you want to reproduce our exact
Bayes factors, please send a message to mauriziosicorello@gmail.com and/or via twitter/researchgate.

Paths are set relative to script location, but this might only work on windows machines.
If it doesn't work, check the "setwd()" commands, especially the path separator ("/" or "\")


