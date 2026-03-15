# The Core Microbiome

library(phyloseq)
library(microbiome)
library(tidyverse)

# loading data

taxonomy = read.delim('Datasets/taxonomy.tsv')
tree = read_tree('Datasets/tree.nwk')

counts = read.delim('Datasets/anemia-feature-table.txt')
metadata = read.delim("Datasets/anemia_metadata.txt",
                      header = TRUE,
                      sep = "",
                      row.names = 1)