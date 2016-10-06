# 2016-paper-ACMG-penetrance

========================
ACMG-Penetrance Pipeline
========================

Date: October 5, 2016

Here are some brief notes on how to run the pipeline for our 2016 paper on the penetrance of the ACMG-56.  
Using these commands, you should be able to completely recreate all figures in the paper.


-------------------------------------------------------------

Before we get started, we need to install all the necessary software, including:
 - R
 - RStudio (contains RMarkdown and knitr): https://github.com/rstudio/rstudio
 - SAMtools (the code uses tabix to download various ACMG-specific regions from 1000 Genomes): https://github.com/samtools/samtools
 
 
-------------------------------------------------------------
Check out the source repository, which contains our Markdown scripts and data:: 
 - git clone https://github.com/jamesdiao/2016-paper-ACMG-penetrance/


-------------------------------------------------------------

There are 4 versions of the Markdown script: 
 - **ACMG_scratch_10_5.Rmd**: compiles everything from scratch (runtime: 30-50 minutes)
 - **ACMG_load_VCF_10_5.Rmd**: uses pre-downloaded VCFs from 1000 Genomes (runtime: 20-40 minutes). 
This skips the step where VCFs are downloaded directly from 1000 Genomes. SAMtools and tabix are no longer required. 
 - **ACMG_load_Rdata_10_5.Rmd**: uses a saved data objects: (runtime: 5-10 minutes)
This skips the step where the VCFs are downloaded directly from 1000 Genomes, and the step where the VCFs are unnested, processed, and imported. SAMtools and tabix are no longer required. 
 - **ACMG_fast_10_5.Rmd**: uses cached files and uses cached data objects in the repo to recreate figures. 
This skips almost all data processing tasks and allows the figures to be generated immediately. SAMtools and tabix are no longer required. 

Open any of these .Rmd files in Rstudio. Click 'Knit', and R will compile the code as written. 
At this point, '**ACMG_*_10_5.html**' will open, containing the paper with the figures you just created.
All 4 Markdown scripts should give identical results. 
