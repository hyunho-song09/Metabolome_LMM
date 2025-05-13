# Load required libraries
library(dplyr)
library(ggplot2)
library(xlsx)
library(mixOmics)
library(vegan)

# Load metabolome and microbiome data
raw.metabolome <- read.csv("example_data01.csv")
raw.microbiome <- read.csv("example_data02.csv")