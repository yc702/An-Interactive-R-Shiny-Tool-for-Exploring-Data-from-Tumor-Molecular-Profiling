# An-Interactive-R-Shiny-Tool-for-Exploring-Data-from-Tumor-Molecular-Profiling
An Interactive R/Shiny Tool for Exploring Data from Tumor Molecular Profiling

**This Shiny App allows for the following displays:**

1. **Heatmaps/tables of biomarkers for _disease_ and _disease category_, by _technology_ and _actionability_:** 
    * Heatmap
        + Sort the top 50 biomarkers that have been most mutated among patients
        + Color can be Red/Green or Blue/Yellow
    * Table
        + Table has first 12 columns from raw data
        + Table allows for sorting and searching
  
2. **Biomarker correlation circle plot by _technology_:**
    * Only consider predetermined set of biomarkers (may expand in the future) filtered by technology
    
3. **Steps to run R Shiny App:**
    * The dataset has not been uploaded so the app does not work at this time;
    * The app should be open using the _Shiny_Markdown.Rmd_ file and the rest of the files are the helper functions;
    * The datasets of _Actionability Assessment.csv_ and _Functional Category.csv_ contain information of the gene actionability levels and functional categories that the genes belong to.

