# 2016-paper-ACMG-penetrance

====================================
Running the ACMG-penetrance pipeline
====================================

:Date: October 4, 2016

Here are some brief notes on how to run the pipeline for our 2016 paper on the penetrance of the ACMG-56.  
Using these commands you should be able to completely recreate the figures in the paper.

The instructions below will reproduce all of the figures in the paper, and will then compile the paper from scratch using the new figures.

1. Collect necessary data
-------------------------------------------------------------
Check out the source repository and grab the initial data sets::
 git clone https://github.com/jamesdiao/2016-paper-ACMG-penetrance/
 
2. Installing necessary software
-------------------------------------------------------------

Before we get started, we need to install all the necessary software, including:
 - R
 - RStudio
 - SAMtools


Recreating the figures
-------------------------------------------------------------

Now go into the pipeline directory and open Master_9_30.Rmd in Rstudio. Click 'Knit', and wait for it to fully load. 
At this point, 'Master_9_30.html' will open, containing the paper with the figures you just created.
