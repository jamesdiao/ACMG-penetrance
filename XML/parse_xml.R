
pkg_list <- c("scrapeR","RMySQL", "knitr","pander","ggplot2","ggrepel", "tibble","curl","tidyr","dplyr","XML")
installed <- pkg_list %in% installed.packages()[,"Package"]
if (!all(installed))
  install.packages(pkg_list[!installed])
sapply(pkg_list, require, character.only = T)
setwd("/Users/jamesdiao/Documents/Kohane_Lab/2016-paper-ACMG-penetrance")
data <- xmlParse(file = "RCV000077146.xml")
xml_data <- xmlToList(data)
xml_data$Title
xml_data$ReferenceClinVarAssertion$ClinicalSignificance
xml_data$ReferenceClinVarAssertion$MeasureSet$Measure$SequenceLocation
