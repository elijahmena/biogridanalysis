# Script written by Elijah Mena for Timms, Mena et al., Nature Cell Biology, 2023.
# Script finds hits from screen data that are also present in the BioGRID database of physical interactions
#   and it compares that to the number of hits expected from simulated screens of random hits

# Need screen data as csv file containing statistical scores, the substrate column, and the E3 (usually labeled "id")
# Also need BioGRID datasets (used 4.4.220 dataset for analysis)
# BioGRID data is accessible here:
# https://downloads.thebiogrid.org/BioGRID

# Timms, Mena et al. used the following cutoffs for classifying hits: p<0.01 for 6-bin screens and Mageck>2.5 for 1-bin screens.

# dplyr, ggplot2, and stringr are required
library(dplyr)
library(ggplot2)
library(stringr)

# Read in Biogrid data and remove interactions
biogrid <- read.csv("raw_data/BIOGRID-Homo_sapiens_4.4.220_onlyGeneSymbols.csv", header = TRUE)
biogrid$interaction <- paste(pmin(toupper(biogrid$Official.Symbol.Interactor.A), toupper(biogrid$Official.Symbol.Interactor.B)),
                             pmax(toupper(biogrid$Official.Symbol.Interactor.A), toupper(biogrid$Official.Symbol.Interactor.B)),
                             sep = "-")

# Read CSV files and clean up columns (For C-term screens, need to try to fill blank gene names)

## CUL2 C-terminal peptide screen
cterm_cul2 <- read.csv("raw_data/cul2_cterm.csv", header = TRUE)
cterm_cul2['Screen'] <- 'C-term CUL2'
pairings <- read.csv("raw_data/extra_gene_names.csv", header = TRUE)
cterm_cul2 <- merge(cterm_cul2, pairings, by = "Substrate_ID", all.x = TRUE) # Merge the dataframes based on Substrate_ID
cterm_cul2$Substrate.x[cterm_cul2$Substrate.x == ""] <- 
  cterm_cul2$Substrate.y[cterm_cul2$Substrate.x == ""] # Fill in the blank values in the Substrate column
cterm_cul2 <- subset(cterm_cul2, select = -c(Uniprot_ID,AA_seq.y,Substrate.y)) # Remove unnecessary columns
cterm_cul2 <- cterm_cul2 %>% # Rename columns
  rename("Substrate" = "Substrate.x", "AA_seq" = "AA_seq.x")
complete_rows <- !(cterm_cul2$Substrate == "" | is.na(cterm_cul2$Substrate))  # Identify rows with non-blank values
cterm_cul2 <- cterm_cul2[complete_rows, ] # Subset the data frame
cterm_cul2 <- cterm_cul2  %>% select(c('Substrate','id','Score','Screen'))
cterm_cul2_hits <- cterm_cul2[cterm_cul2$Score > 2.5, ]

## CUL4 C-terminal peptide screen
cterm_cul4 <- read.csv("raw_data/cul4_cterm.csv", header = TRUE)
cterm_cul4['Screen'] <- 'C-term CUL4'
pairings <- read.csv("raw_data/extra_gene_names.csv", header = TRUE)
cterm_cul4 <- merge(cterm_cul4, pairings, by = "Substrate_ID", all.x = TRUE) # Merge the dataframes based on Substrate_ID
cterm_cul4$Substrate.x[cterm_cul4$Substrate.x == ""] <- 
  cterm_cul4$Substrate.y[cterm_cul4$Substrate.x == ""] # Fill in the blank values in the Substrate column
cterm_cul4 <- subset(cterm_cul4, select = -c(Uniprot_ID,AA_seq.y,Substrate.y)) # Remove unnecessary columns
cterm_cul4 <- cterm_cul4 %>% # Rename columns
  rename("Substrate" = "Substrate.x", "AA_seq" = "AA_seq.x")
complete_rows <- !(cterm_cul4$Substrate == "" | is.na(cterm_cul4$Substrate))  # Identify rows with non-blank values
cterm_cul4 <- cterm_cul4[complete_rows, ] # Subset the data frame
cterm_cul4 <- cterm_cul4  %>% select(c('Substrate','id','Score','Screen'))
cterm_cul4_hits <- cterm_cul4[cterm_cul4$Score > 2.5, ]

## CUL3 ORF 1-bin screen
orf_1bin <- read.csv("raw_data/cul3_orfs_1bin.csv", header = TRUE)
orf_1bin['Screen'] <- 'ORF CUL3 1-bin'
orf_1bin <- orf_1bin  %>% select(c('Substrate','id','Score','Screen'))
orf_1bin_hits <- orf_1bin[orf_1bin$Score > 2.5, ]

## CUL3 ORF 6-bin screen
orf_6bin <- read.csv("raw_data/cul3_orfs_6bin.csv", header = TRUE)
orf_6bin['Screen'] <- 'ORF CUL3 6-bin'
orf_6bin <- orf_6bin %>%
  rename(Score = P.value..Mann.Whitney.U.,
         id = Target) #Change column names to be consistent
orf_6bin <- orf_6bin  %>% select(c('Substrate','id','Score','Screen'))
orf_6bin_hits <- orf_6bin[orf_6bin$Score < 0.01, ]
orf_6bin_hits <- orf_6bin_hits[!duplicated(orf_6bin_hits),] #Remove duplicates

## Internal degron screen (group1, sort1)
internal_g1s1 <- read.csv("raw_data/internal_g1s1.csv", header = TRUE)
internal_g1s1['Screen'] <- 'Internal (Group1, Sort1)'
internal_g1s1 <- internal_g1s1  %>% select(c('Substrate','id','Score','Screen'))
internal_g1s1$Substrate <- sub("_.*", "", internal_g1s1$Substrate)
internal_g1s1_hits <- internal_g1s1[internal_g1s1$Score > 2.5, ]

## Internal degron screen (group1, sort2)
internal_g1s2 <- read.csv("raw_data/internal_g1s2.csv", header = TRUE)
internal_g1s2['Screen'] <- 'Internal (Group1, Sort2)'
internal_g1s2 <- internal_g1s2  %>% select(c('Substrate','id','Score','Screen'))
internal_g1s2$Substrate <- sub("_.*", "", internal_g1s2$Substrate)
internal_g1s2_hits <- internal_g1s2[internal_g1s2$Score > 2.5, ]

## Internal degron screen (group2, sort1)
internal_g2s1 <- read.csv("raw_data/internal_g2s1.csv", header = TRUE)
internal_g2s1['Screen'] <- 'Internal (Group2, Sort1)'
internal_g2s1 <- internal_g2s1  %>% select(c('Substrate','id','Score','Screen'))
internal_g2s1$Substrate <- sub("_.*", "", internal_g2s1$Substrate)
internal_g2s1_hits <- internal_g2s1[internal_g2s1$Score > 2.5, ]

## Internal degron screen (group2, sort2)
internal_g2s2 <- read.csv("raw_data/internal_g2s2.csv", header = TRUE)
internal_g2s2['Screen'] <- 'Internal (Group2, Sort2)'
internal_g2s2 <- internal_g2s2  %>% select(c('Substrate','id','Score','Screen'))
internal_g2s2$Substrate <- sub("_.*", "", internal_g2s2$Substrate)
internal_g2s2_hits <- internal_g2s2[internal_g2s2$Score > 2.5, ]

# Combine into screen data and hit data into respective pooled datasets
allscreened <- rbind(cterm_cul2, cterm_cul4, orf_1bin, orf_6bin, 
                     internal_g1s1, internal_g1s2, internal_g2s1, internal_g2s2)
allhits <- rbind(cterm_cul2_hits, cterm_cul4_hits, orf_1bin_hits, orf_6bin_hits, 
                 internal_g1s1_hits, internal_g1s2_hits, internal_g2s1_hits, internal_g2s2_hits)

# Make columns of interactions (pmin/pmax makes sure the interaction order is consistent w/ Biogrid)
allscreened$interaction <- paste(pmin(toupper(allscreened$Substrate), toupper(allscreened$id)),
                           pmax(toupper(allscreened$Substrate), toupper(allscreened$id)),
                           sep = "-")
allhits$interaction <- paste(pmin(toupper(allhits$Substrate), toupper(allhits$id)),
                           pmax(toupper(allhits$Substrate), toupper(allhits$id)),
                           sep = "-")

# Create a new column in the allhits dataframe to indicate presence in the Biogrid dataframe
allhits$InBiogrid <- allhits$interaction %in% biogrid$interaction

#Rename and re-order columns before export
allhits <- allhits %>% # Rename columns
  rename("E3" = "id")
allhits <- allhits[, c('Substrate','E3','interaction','Screen','Score','InBiogrid')]
allscreened <- allscreened %>% # Rename columns
  rename("E3" = "id")

###### Optionally Remove Cullin hits from data ######
#allhits <- allhits[!grepl('CUL', allhits$E3),]
#allscreened <- allscreened[!grepl('CUL', allscreened$E3),]

###### EXPORT CSV ######
write.csv(allhits, "analyzed_data/allhits.csv", row.names = FALSE)

# Remove all duplicate interactions from data before comparison with Biogrid
allscreened_unique <- allscreened[!duplicated(allscreened$interaction), ]
allhits_unique <- allhits[!duplicated(allhits$interaction), ]

# Calculate number of hits (should yield 1013 hits, with 31 in common with BioGRID)
cat("The total number of unique hits found is:", nrow(allhits_unique), "\n")
num_yes <- sum(allhits_unique$InBiogrid == "TRUE")
cat("Number of interactions present on biogrid is:", num_yes, "\n")


#Function determines number of protein pairs in common between two dataframes
# Note: Interaction headings must be labeled 'interaction'
find_common_protein_pairs <- function(df1, df2) {
  InteractionPresent <- df1$interaction %in% df2$interaction
  return(sum(InteractionPresent == "TRUE"))
}


# Function that produces a new random dataframe with n number of hits
random_sample <- function(data, n) {
  # check if n is greater than number of rows in data
  if(nrow(data) < n){
    stop("n is greater than the number of rows in the data")
  }
  
  # sample n rows at random from the data
  sample_rows <- sample(nrow(data), n, replace = FALSE)
  
  # return a new data frame containing the sampled rows
  return(data[sample_rows, ])
}


# Function that creates m mock random datasets (each with n number of hits), 
# and for each random dataset calculates how many protein pairs are in common with Biogrid data
# Returns a dataframe containing the number of protein pairs in common for each random screen
simulate <- function(m,n) {
  i <- 0
  all_simulations <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(all_simulations) <- c('i', 'num_in_common')
  
  while (i < m) {
    # For each of m datasets, find number of pairs in common between random sampling of n hits and biogrid data
    num_hits <- find_common_protein_pairs(random_sample(allscreened_unique,n),biogrid)
    
    # Add data to a dataframe
    all_simulations[nrow(all_simulations) + 1,] = c(i,num_hits)
    i <- i + 1
  }
  return(all_simulations)
}

# Perform simulations
# First argument is the number of simulated datasets (increase number if desired)
# Second argument is the number of hits within each screen
result <- simulate(100,1013)

# Graph data in histogram
# Need to change geom_vline to correspond to actual screen data
plot <- ggplot(result, aes(x=num_in_common)) +
  geom_histogram(binwidth=1,color="black", fill="grey40") +
  theme_classic() +
  xlim(-1,35) +
  ggtitle('Simulations of random screen data') +
  xlab('Number of hits that are also present in Biogrid repository')+
  ylab('Number of simulations')+
  geom_vline(xintercept = 31, linetype="dashed", #Red line is the number of actual hits in common observed
             color = "red", size=1.5)

##### Save Plot ######
ggsave(filename = 'analyzed_data/allhits.eps', plot = plot, device = "eps", width = 5, height = 3)

