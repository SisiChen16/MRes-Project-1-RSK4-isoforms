# Aim: To correlate the RSK4 isoform expression with overall survival of cancer patients

# Load packages 
library(survminer)
library(survival) # survival analysis 
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)

# Import data 
setwd("/Users/chensisi/Documents/RNAseq/")
load("merge_combine.RData")

# Univariate survival analysis by TPM

isoform_tpm <- merge_combine %>%
  select("ENST00000620340.4","ENST00000262752.4")


