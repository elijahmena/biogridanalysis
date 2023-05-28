This repository includes the script for the BioGRID comparison that was performed in Timms, Mena, et al., 2023, Nature Cell Biology.

# Required software
Analysis was performed using R 4.1.0 and uses the packages dplyr, ggplot2, stringr.

# Required data
Raw data from the paper (same data an in Supplementary files) are in the raw_data folder. BioGRID data can be downloaded from https://downloads.thebiogrid.org/BioGRID. The paper used the Homo sapiens 4.4.220 dataset, which is included in the raw_data folder and is trimmed to only include the Gene Symbol columns.

# Output
**allhits.csv** This file contains all hits above a specified threshold as well as whether they are found in the BioGRID database or not.
**allhits.eps** This plot shows a histogram of the number of hits shared between random simulated datasets and the BioGRID interaction database.

