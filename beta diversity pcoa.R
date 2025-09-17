# beta diversity for the comparison of two treatment groups
# PCoA plots colored by variable (treatment and if the animal was infected)
# check dispersion and run PERMANOVA

library(readxl)
library(phyloseq)
library(ape)
library(dplyr)
library(vegan)
library(ggplot2)
library(ggrepel)

# create phyloseq object for all samples
otu_mat <-read_excel("/Users/caseymeili/Desktop/26454R/Fastq/cbd/cbd_phyloseq.xlsx", sheet="OTU")
tax_mat<- read_excel("/Users/caseymeili/Desktop/26454R/Fastq/cbd/cbd_phyloseq.xlsx", sheet="taxon")
Meta <-read_excel("/Users/caseymeili/Desktop/26454R/Fastq/cbd/cbd_phyloseq.xlsx", sheet="Samples")
Meta <- Meta %>%
  tibble::column_to_rownames("Sample")
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("#OTU ID")
tax_mat <- tax_mat %>%
  tibble::column_to_rownames("#OTU ID")
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(Meta)
samples
Physeq <-phyloseq(OTU, TAX, samples)

# calculate bray curtis
pcoa_bc = ordinate(Physeq, "PCoA", "bray") 

# rename levels
# check current levels
levels(sample_data(Physeq)$Group)
levels(samples$Group)

# rename levels
sample_data(Physeq)$Day <- factor(sample_data(Physeq)$Day,
                                    levels = c("B", "A"),
                                    labels = c("-1", "7 DPI"))

# plot PCoA
plot_ordination(Physeq, pcoa_bc, color = "Treatment", shape = "Infection") +
  geom_point(size = 3) +
  stat_ellipse(aes(group = Treatment, color = Treatment), 
               type = "t", linetype = 1, linewidth = .5) +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("#06b5ce", "#ee00a8")) +
  scale_shape_manual(values = c("TMEV" = 16, "PBS" = 17)) +
  labs(title = "7 dpi (All Samples)", color = "Treatment")

# extract distance matrix
bc_dist <- phyloseq::distance(Physeq, method = "bray")

# run PERMANOVA
adonis2(bc_dist ~ Treatment, data = as(sample_data(Physeq), "data.frame"))

# check dispersion
disp <- betadisper(bc_dist, sample_data(Physeq)$Treatment)

# test for homogeneity of dispersions
anova(disp)  

# permutation test
permutest(disp)  

# export figure
ggsave("/Users/caseymeili/Downloads/all-samples.jpg",
       plot = last_plot(),
       width = 4, dpi = 300)

