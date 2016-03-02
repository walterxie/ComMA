library(vegan)
library(vegetarian)
library(grid)
library(gridExtra)
library(data.table)
library(xtable)
library(tools)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(colorspace)
#library(Matrix)
library(reshape2)
library(ggdendro)
library(plyr)
library(cluster)
library(ape)

#if(!exists("sourcePath")) stop("source path to initiate modules is missing !")
#if(!exists("workingPath")) stop("working path containing data is missing !")

#if(!exists("verbose")) verbose <- TRUE

# utils and fundmental functions, run them before analysis
source("R/IO.R")
source("R/Utils.R")
source("R/UtilsCM.R")
source("R/Diversities.R")

# most abundant 150 OTUs
threshold = 150
