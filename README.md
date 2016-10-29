# 2016-paper-ACMG-penetrance

#### Date: October 16, 2016

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
 - **Running with pre-downloaded 1000 Genomes VCFs**: Unzip `1000G.zip` before knitting the .Rmd, and the program will skip the step where the VCFs are downloaded directly from 1000 Genomes. SAMtools and tabix are no longer required. (20-40 minutes) <br />
 - **Running with processed 1000 Genomes data frame**: Unzip `ACMG_1000G.RData.zip` before knitting the .Rmd, and the program will skip the step where the VCFs are downloaded from 1000 Genomes and then unnested, processed, and imported. SAMtools and tabix are no longer required. (5-10 minutes) <br />
 - **Running with cached data objects**: Unzip `ACMG_Penetrance_cache.zip` and `ACMG_Penetrance_files.zip` before knitting the .Rmd, and the program will use cached data objects to recreate the figures. This skips all data processing tasks and allows the figures to be generated immediately (identical to those in the included .html). SAMtools and tabix are no longer required. (<1 minute) <br />
<br />

All 4 procedures give the same results, which should all be identical to `ACMG_Penetrance.html`


-------------------------------------------------------------

### SUPPLEMENTARY FILES 
`Environ_2016-10-16` contains all relevant data objects retained in the environment after running all code chunks. Loading this into memory allows you to modify and execute individual code chunks in the .Rmd file. <br />
`ACMG_Lit_Small.csv` contains all prevalence estimates for ACMG-56 related diseases from the medical literature. <br />
**Supplementary_Files** contains intermediary files generated while executing the script. <br />
- `ACMG-56_Panel.txt` contains the names of the ACMG-56 genes. <br />
- `clinvar_[date].vcf` contains the raw ClinVar VCF from the ClinVar FTP from that particular date. <br />
- `phase3map.txt` contains the ancestral population of each of the 2,504 individuals in 1000 Genomes. <br />
- `clinvar_query.txt` contains the table of results that are matched by the search query "(APC[GENE] OR MYH11[GENE]... OR WT1[GENE]) AND (clinsig_pathogenic[prop] OR clinsig_likely_pathogenic[prop])' from the ClinVar website. 
**Processed_Files** contains the key data frames used in the analysis as plain text (csv or tab-delimited). Key modifications from raw data are described in the .Rmd file. <br />
- `ACMG_1000G_Processed.csv` contains info for all variants found in ACMG-56 genes in the 1000 Genomes cohort. <br />
- `ACMG_ExAC_Processed.csv` contains info for all variants found in ACMG-56 genes in the ExAC cohort. <br />
- `ClinVar_Processed.csv` contains info for all variants in the ClinVar VCF. <br />
-----------------------------------------------------------------

### CONTACT  

Please contact James Diao (james [dot] diao [at] yale [dot] edu) with any questions. 

<br />
<br />
<br />

