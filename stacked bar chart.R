# for creating a stacked bar chart of the 10 most abundant taxa 
# one bar per sample, seperated by group

library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)

# load data and create phyloseq object
otu_mat <- read_excel("/Users/caseymeili/Desktop/26454R/Fastq/cbd/cbd_phyloseq_infected.xlsx", sheet="OTU")
tax_mat <- read_excel("/Users/caseymeili/Desktop/26454R/Fastq/cbd/cbd_phyloseq_infected.xlsx", sheet="taxon")
Meta <- read_excel("/Users/caseymeili/Desktop/26454R/Fastq/cbd/cbd_phyloseq_infected.xlsx", sheet="Samples")
Meta <- Meta %>% tibble::column_to_rownames("Sample")
otu_mat <- otu_mat %>% tibble::column_to_rownames("#OTU ID")
tax_mat <- tax_mat %>% tibble::column_to_rownames("#OTU ID")
OTU <- otu_table(as.matrix(otu_mat), taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(tax_mat))
SAMPLES <- sample_data(Meta)
Physeq <- phyloseq(OTU, TAX, SAMPLES)

# replace NA genus with "Unknown"
tax_mat <- as(tax_table(Physeq), "matrix")
tax_mat[, "Genus"] <- ifelse(is.na(tax_mat[, "Genus"]) | tax_mat[, "Genus"] == "", 
                             "Unknown", tax_mat[, "Genus"])
tax_table(Physeq) <- tax_table(tax_mat)

# melt phyloseq to long format
df_long <- psmelt(Physeq)
df_long$Abundance <- as.numeric(df_long$Abundance)

# identify top N taxa overall
top_n <- 10
top_taxa <- df_long %>%
  mutate(Genus = as.character(Genus)) %>%
  group_by(Genus) %>%
  summarize(TotalAbundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = top_n) %>%
  pull(Genus)

df_long$Genus_plot <- ifelse(df_long$Genus %in% top_taxa, df_long$Genus, "Other")

# create combined grouping variable for ordering
df_long$Group <- paste(df_long$Treatment, df_long$Day, sep = "_")
df_long$Group <- factor(df_long$Group, levels = c("Vehicle_B", "Vehicle_A", "CBD_B", "CBD_A"))

# rename group levels for clarity on plot
df_long$Group <- factor(df_long$Group,
                        levels = c("Vehicle_B", "Vehicle_A", "CBD_B", "CBD_A"),
                        labels = c("Vehicle Day -1", "Vehicle 7 DPI",
                                   "CBD Day -1", "CBD 7 DPI"))

# reorder Genus_plot so "Other" is last
df_long$Genus_plot <- factor(df_long$Genus_plot,
                             levels = c(setdiff(unique(df_long$Genus_plot), "Other"), "Other"))

# define custom color palette
# you can also easily use the default colors, I just find them incredibly ugly 

# Assign colors to taxa
taxa_levels <- levels(df_long$Genus_plot)

custom_colors <- c(
  "#3B82F6", # blue
  "#EC4899", # pink
  "#FFD700", # yellow
  "#10B981", # green
  "#F97316", # orange
  "#8B5CF6", # purple
  "#F43F5E", # red
  "#00CED1", # cyan 
  "#F59E0B", # gold-orange
  "#4F63D7", # indigo
  "Other" = "#A3A3A3" 
)

taxa_levels <- levels(df_long$Genus_plot)
names(custom_colors) <- c(setdiff(taxa_levels, "Other"), "Other")

# plot bar chart
ggplot(df_long, aes(x = Sample, y = Abundance, fill = Genus_plot)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ Group, scales = "free_x", nrow = 1, strip.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(x = "", y = "Relative Abundance (%)", fill = "Genus",
       title = "Stacked bar chart of top taxa per sample (percent abundance)")

ggsave("/Users/caseymeili/Downloads/cbd-d.jpg",
       plot = last_plot(),
       width = 11, dpi = 300)
