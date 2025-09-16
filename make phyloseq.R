library(tidyverse)
library(phyloseq)
library(writexl) 
library(tidyr)

# read shared file
shared <- read.table("/Users/caseymeili/Desktop/26454R/Fastq/cbd/cbd.opti_mcc.shared", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Drop metadata columns (keep just OTU counts)
otu_table <- shared %>%
  select(-label, -numOtus) %>%
  pivot_longer(-Group, names_to = "OTU", values_to = "Count") %>%
  pivot_wider(names_from = Group, values_from = Count)

write.csv(otu_table, "/Users/caseymeili/Downloads/cbd_otu_table.csv", row.names = FALSE)


