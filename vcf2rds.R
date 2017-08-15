get_clinvar <- function(clinvar_file) {
  
  extract_element <- function(phrase) {
    str_match_all(input$INFO, sprintf('%s=([^;]*);', phrase)) %>% 
      lapply('[[', 2) %>% unlist
  }
  file.by.line <- readLines(clinvar_file)
  #file_date <- as.Date(strsplit(file.by.line[2],"=")[[1]][2], "%Y%m%d")
  #system(sprintf("mv %s ClinVar_Reports/clinvar_%s.vcf", clinvar_file, file_date))
  clean.lines <- file.by.line[!grepl("##.*", file.by.line)] #Remove ## comments
  clean.lines[1] <- sub('.', '', clean.lines[1]) #Remove # from header
  input <- read.table(text = paste(clean.lines, collapse = "\n"), header = T, stringsAsFactors = F, 
                      comment.char = "", quote = "", sep = "\t")
  input <- input[nchar(input$REF)==1,] #deletions
  alt_num <- sapply(strsplit(input$ALT,","),length) #number of alts
  acceptable_nchar <- 2*alt_num-1 #adds in the length from commas, if each alt is 1 nt.
  input <- input[nchar(input$ALT)==acceptable_nchar,] #insertions
  input$ALT <- strsplit(input$ALT,",")
  input$CLNALLE <- extract_element('CLNALLE') %>% strsplit(',', fixed = T) %>% lapply(as.integer)
  input$CLNSIG <- extract_element('CLNSIG') %>% strsplit(',', fixed = T)
  input$CLNDBN <- extract_element('CLNDBN') %>% strsplit(',', fixed = T)
  
  #CLNALLE has 0,-1,3,4 --> CLNSIG has 1,2,3,4 --> ALT has 1. 
  taking <- sapply(input$CLNALLE, function(x) x[x>0] ) #Actual elements > 0. Keep these in CLNSIG and ALT 
  taking_loc <- sapply(input$CLNALLE, function(x) which(x>0) )#Tracks locations for keeping in CLNALLE
  keep <- sapply(taking, length)>0 #reduce everything to get rid of 0 and -1
  # Reduce, reduce, reduce. 
  taking <- taking[keep]
  taking_loc <- taking_loc[keep]
  input <- input[keep,]
  
  #Make this more readable
  input$ALT <- sapply(1:nrow(input), function(row) {
    input$ALT[[row]][taking[[row]]]
  })
  
  col_subset <- function(name) {
    sapply(1:nrow(input), function(row) {
      unlist(input[row,name])[taking_loc[[row]]]
    })
  }
  input$CLNSIG <- col_subset("CLNSIG")
  input$CLNALLE <- col_subset("CLNALLE")
  input$CLNDBN <- col_subset("CLNDBN")
  filter_condition <- input[,unlist(lapply(input, typeof))=="list"] %>% 
    apply(1,function(row) !any(is.na(row)))
  input <- input %>% filter(filter_condition) %>%
    unnest %>% unite(VAR_ID, CHROM, POS, REF, ALT, sep = "_", remove = F) %>%
    select(VAR_ID, CHROM, POS, REF, ALT, ID, CLNSIG, CLNDBN) %>% 
    mutate(CLNSIG = strsplit(CLNSIG,"|",fixed = T)) %>% 
    mutate(CLNDBN = strsplit(CLNDBN,"|",fixed = T)) %>% 
    mutate(POS = as.integer(POS))
  input$CLNSIG <- sapply(input$CLNSIG, function(x) as.integer(x))
  input$pathogenic_incl_conflicts <- sapply(input$CLNSIG, function(x) any(x %in% c(4,5)))
  input$pathogenic_no_conflicts <- sapply(input$CLNSIG, function(x) any(x %in% c(4,5)) & !(any(x %in% c(2,3)))) 
  #input$INTERP <- input$pathogenic_no_conflicts
  #input$LMM <- grepl("Laboratory_for_Molecular_Medicine",input$INFO)
  return(input)
}

require(stringr)
setwd('/Users/jamesdiao/Documents/Kohane_Lab/clinvaR/inst/extdata')
get_date_list <- function() {
  clinvar_reports <- system("ls Archive_Tables", intern = T)
  clinvar_reports <- clinvar_reports[grep(".tsv",clinvar_reports)]
  do.call("rbind", lapply(clinvar_reports, function(tsv) {
    read.table(file = sprintf("Archive_Tables/%s", tsv), sep = "\t", 
               header = T, stringsAsFactors = F)
  })) -> archive
  regmatches(archive$Name,regexpr("^File:clinvar_(20.{6})\\.vcf\\.gz$",archive$Name)) %>% 
    str_extract("20.{6}") %>% as.Date(format = "%Y%m%d") %>% unique()
}
all_dates <- get_date_list()

get_dates <- function(date_list, range) {
  date_list <- date_list %>% as.Date()
  if (!exists("all_dates"))
    all_dates <- get_date_list()
  if (range) {
    keep <- between(all_dates, min(date_list), max(date_list))
  } else {
    keep <- all_dates %in% date_list
  }
  all_dates[keep] %>% as.character()
}
date_list <- c("2012-06-16","2013-01-14","2014-02-11","2015-01-06","2016-01-04",as.character(Sys.Date()))
range = TRUE
clinvar_dates <- get_dates(date_list, range)

Sys.setlocale('LC_ALL','C')
for (cv_date in clinvar_dates) {
    # Find file (this is input to ClinVar_Report.Rmd)
    clinvar_file <- sprintf("clinvar_%s.vcf", cv_date)
    # Copy into the main folder
    system(sprintf("cp VCF/clinvar_%s.vcf.gz .", cv_date))
    # Gunzip the VCF
    system(sprintf("gunzip clinvar_%s.vcf.gz", cv_date))
    # Download ClinVar 
    orig_clinvar <- download_clinvar(clinvar_file)
    orig_clinvar <- orig_clinvar[!duplicated(orig_clinvar$VAR_ID),]
    saveRDS(orig_clinvar, sprintf('clinvar_%s.rds', cv_date))
}

if (any(grepl("clinvar_.*\\.vcf", system("ls", intern = T))))
  system("rm clinvar_*.vcf*")

temp <- lapply(clinvar_dates[1:3], function(cv_date) {
  readRDS(sprintf('clinvar_%s.rds', cv_date))
})

temp[[1]]$VAR_ID -> x
mean(x %in% temp[[2]]$VAR_ID)

for (cv_date in clinvar_dates) {
  readRDS(sprintf('clinvar_%s.rds', cv_date)) %>% 
    select()
  
}
  

