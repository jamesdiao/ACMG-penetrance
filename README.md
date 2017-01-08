# 2017-ACMG-penetrance

#### Date: 07 January 2017

Here are some brief notes on how to run the pipeline for our 2017 project on the penetrance of the ACMG-59. <br />
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

Clone this source repository, containing our Markdown scripts and data: <br />
- `$ git clone https://github.com/jamesdiao/2016-paper-ACMG-penetrance/`


-------------------------------------------------------------

### RECREATING FIGURES FOR ACMG-PENETRANCE

1. Navigate to the directory `/ACMG_Penetrance/`.  
2. Open `/ACMG_Penetrance/ACMG_Penetrance.Rmd`.  
3. Click "Knit" to compile everything into `/ACMG_Penetrance/ACMG_Penetrance.html`. 

There are four different ways to run the code: <br />
 - **Running from scratch**: This is the default method. Everything will be compiled on the spot (runtime: 30-50 minutes). <br />
 - **Running with pre-downloaded 1000 Genomes VCFs**: Unzip `1000G.zip` before knitting the .Rmd, and the program will skip the step where the VCFs are downloaded directly from 1000 Genomes. SAMtools and tabix are no longer required. (20-40 minutes) <br />
 - **Running with processed 1000 Genomes data frame**: Unzip `/ACMG_Penetrance/ACMG_1000G.rds.zip` before knitting the .Rmd, and the program will skip the step where the VCFs are downloaded from 1000 Genomes and then unnested, processed, and imported. SAMtools and tabix are no longer required. (5-10 minutes) <br />
 - **Running with cached data objects**: Unzip `/ACMG_Penetrance/ACMG_Penetrance_cache.zip` and `/ACMG_Penetrance/ACMG_Penetrance_files.zip` before knitting the .Rmd, and the program will use cached data objects to recreate the figures. This skips all data processing tasks and allows the figures to be generated immediately (identical to those in the included .html). SAMtools and tabix are no longer required. (<1 minute) <br />
<br />

All 4 procedures give the same results, which should all be identical to `/ACMG_Penetrance/ACMG_Penetrance.html`

-----------------------------------------------------------------

### RECREATING FIGURES FOR CLINVAR_REPORTS

1. Navigate to the directory `/ClinVar_Reports/`. 
2. Edit input parameters loaded from `ClinVar_Parameters.txt`. <br /> 
More readme details in the file header. 
3. Open `/ClinVar_Reports/ClinVar_Script.Rmd` and click "Knit".  
 - Individual reports for each date will be generated and saved in `\ClinVar_Reports\Reports`.  
 - All reports will be combined into a final report in `/ClinVar_Reports/ClinVar_Penetrance.html`.  

-------------------------------------------------------------

### DESCRIPTION OF FILES/FOLDERS
1. **1000G.zip** contains all VCF files collected from 1000 Genomes, saved as comma-separated values. Unzipping this allows `ACMG_Penetrance.Rmd` to skip the step where the VCFs are downloaded, processed, and imported. 
2. **ACMG_Penetrance/** contains all files needed to generate `ACMG_Penetrance.html`. <br />
3. **ClinVar_Reports/** contains all files needed to generate `ACMG_Penetrance.html`. <br />
`Environ_[date]` contains all relevant data objects retained in the environment from after running all code chunks (from given [date]). Loading this into memory allows you to modify and execute individual code chunks in the .Rmd file.  <br />
4. **ExAC/** contains all VCF files collected from ExAC, saved as comma-separated values.  
5. **gnomAD/** contains all VCF files collected from gnomAD, saved as comma-separated values.  
6. **Key_Figures/** contains the following:  
- `Figure_0.pdf`: Barplot of Min/Point/Max Penetrance (using gnomAD)
- `Figure_1.pdf`: Heatmap of Max Penetrance by Ancestry (using gnomAD)
- `Table_1.csv`: CSV table documenting population parameters for each disease. Columns include: Phenotype, GeneReviews, Prevalence, Case Allele Frequency (CAF), Population Allele Frequency (PAF), Max_Penetrance. 
- `Table_1.pdf`: Same as `Table_1.csv`, but pander-formatted and saved as pdf. 
- `Table_1.Rmd`: Code used to generate `Table_1.pdf` from `Table_1.csv`.
7. **Shiny_App** contains the app (and dependencies) that allows users to generate `Figure_0.pdf` and `Figure_1.pdf` for custom inputs. 
8. **Supplementary_Files/** contains the following input files and intermediary files. <br />
- `ACMG_SF_v2.0.txt` contains the table from pages 4-5 of the Kalia et al. paper establishing ACMG SD v2.0. <br />
- `clinvar_result_[date].txt` contains the table downloaded from the ClinVar website matched by the search query: <br />"(APC[GENE] OR MYH11[GENE]... OR WT1[GENE]) AND (clinsig_pathogenic[prop] OR clinsig_likely_pathogenic[prop])".  <br />
- `cvquery_hg19_[date].bed` contains chromosomal locations of all variants in hg19 (as given by ClinVar). Ex: "chr1	17349110	17349111."  
- `cvquery_hg38_[date].bed` contains chromosomal locations of all variants in hg38 (converted by liftOver).
- `download_output.txt` contains the most recent output from downloading of VCFs from 1000 Genomes. It lists the selected transcription region and associated chromosomal locations for all genes. 
- `ind_test.rds` contains the output of 1000 trials, each sampling 1000 variants from 1000 Genomes. This was used to evaluate independence assumptions across genes. 
- `Literature_Prevalence_Estimates.csv` contains all prevalence estimates for ACMG-59 related diseases from the medical literature. Sources, quoted sections, phenotype, abbreviations, inheritance mode, and reporting guidelines are listed in additional columns. <br />
- `phase3map.txt` contains the ancestral population of each of the 2,504 individuals in 1000 Genomes. <br />
- `query_input.txt` contains the search query used to collected `clinvar_result_[date].txt`.  Line 1 is the query itself; line 2 is the generated URL. Only line 2 will give results (due to the size of the query). 

9. **Supplementary_Files/Processed_Files/** contains the key data frames used in the analysis as plain text (csv or tab-delimited). Key modifications from raw data are described in the .Rmd file. <br />
- `ACMG_1000G.csv.zip` contains info for all variants found in ACMG-59 genes in the 1000 Genomes cohort. <br />
- `ACMG_ExAC.csv.zip` contains info for all variants found in ACMG-59 genes in the ExAC cohort. <br />
- `ACMG_gnomAD.csv.zip` contains info for all variants found in ACMG-59 genes in the gnomAD cohort. <br />
- `ClinVar_Processed.txt.zip` contains info for all variants in the ClinVar VCF. <br />

-----------------------------------------------------------------

### CONTACT  

Please contact James Diao (james [dot] diao [at] yale [dot] edu) with any questions. 

<br />
<br />
<br />

