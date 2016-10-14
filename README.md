# 2016-paper-ACMG-penetrance

========================
ACMG-Penetrance Analysis
========================

#### Date: October 14, 2016

Here are some brief notes on how to run the pipeline for our 2016 paper on the penetrance of the ACMG-56. <br />
Using these commands, you should be able to completely recreate all figures in the paper.


-------------------------------------------------------------

### SOFTWARE

Before we get started, we need to install all the necessary software, including:

 - R <br />
 - RStudio (contains RMarkdown and knitr): https://github.com/rstudio/rstudio <br />
 - SAMtools/tabix (tabix is used to download various ACMG-specific regions from 1000 Genomes): https://github.com/samtools/samtools
 
Note: SAMtools and tabix are only required for the downloading step. Using zipped VCFs or data objects from the source repo allows you to generate the figures without tabix. 
 
-------------------------------------------------------------

### CLONE REPO

Navigate to an appropriate directory and clone the source repository, which contains our Markdown scripts and data: <br />
- `$ git clone https://github.com/jamesdiao/2016-paper-ACMG-penetrance/`


-------------------------------------------------------------
### RECREATING FIGURES

To run the code, open `ACMG_Penetrance.Rmd` and click "Knit" to compile everything into `ACMG_Penetrance.html`. 

There are four different ways to run the code: <br />
 - **Running from scratch**: This is the default method. Everything will be compiled on the spot (runtime: 30-50 minutes). <br />
 - **Running with pre-downloaded 1000 Genomes VCFs**: Unzip `1000G.zip` before knitting the .Rmd, and the program will skip the step where the VCFs are downloaded directly from 1000 Genomes. SAMtools and tabix are no longer required. <br />
 - **Running with processed 1000 Genomes data frame**: Unzip `ACMG_1000G.RData.zip` before knitting the .Rmd, and the program will skip the step where the VCFs are downloaded from 1000 Genomes and then unnested, processed, and imported. SAMtools and tabix are no longer required. <br />
 - **Running with cached data objects**: Unzip `ACMG_Penetrance_cache.zip` and `ACMG_Penetrance_files.zip` before knitting the .Rmd, and the program will use cached data objects to recreate the figures. This skips all data processing tasks and allows the figures to be generated immediately (identical to those in the included .html). SAMtools and tabix are no longer required. <br />
<br />
All 4 procedures give the same results, which should all be identical to `ACMG_Penetrance.html`


<br />
<br />
<br />

