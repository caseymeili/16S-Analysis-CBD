# differential abundance analysis using DESeq2

library(phyloseq)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)

# create phyloseq object for infected animals
otu_mat <-read_excel("/Users/caseymeili/Desktop/26454R/Fastq/cbd/cbd_phyloseq_infected.xlsx", sheet="OTU")
tax_mat<- read_excel("/Users/caseymeili/Desktop/26454R/Fastq/cbd/cbd_phyloseq_infected.xlsx", sheet="taxon")
Meta <-read_excel("/Users/caseymeili/Desktop/26454R/Fastq/cbd/cbd_phyloseq_infected.xlsx", sheet="Samples")
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

# filter low-abundance taxa
Physeq_filt <- filter_taxa(Physeq, function(x) sum(x) > 5, TRUE)

# convert to DESeq2 object
dds <- phyloseq_to_deseq2(Physeq_filt, ~ Treatment)

# run DESeq2
dds <- DESeq(dds)

# extract results
res <- results(dds, contrast = c("Treatment", "Vehicle", "CBD"))
res <- res[order(res$padj, na.last = NA), ]

# add taxonomy
res_tax <- as.data.frame(res)
tax_table_df <- as(tax_table(Physeq_filt), "matrix") %>% as.data.frame()
res_tax <- cbind(res_tax, tax_table_df[rownames(res_tax), ])

# use Genus or Family as labels
res_tax$label <- ifelse(!is.na(res_tax$Genus), as.character(res_tax$Genus),
                        ifelse(!is.na(res_tax$Family), as.character(res_tax$Family), rownames(res_tax)))

# add significance & enrichment columns
res_tax$significant <- ifelse(res_tax$padj < 0.05, "Significant", "Not Significant")
res_tax$enriched <- ifelse(res_tax$log2FoldChange > 0, "Vehicle", "CBD")

# create a color column: gray for non-significant, colors for significant
res_tax$color <- ifelse(res_tax$significant == "Significant",
                        ifelse(res_tax$enriched == "Vehicle", "#ee00a8", "#06b5ce"),
                        "gray80")

# select top taxa for labeling (significant only)
top_taxa <- res_tax %>%
  filter(significant == "Significant") %>%
  slice_max(order_by = abs(log2FoldChange), n = 10)

# plot volcano plot
ggplot(res_tax, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = color), alpha = 0.6) +
  scale_color_identity() +  # use exact colors in 'color' column
  geom_text_repel(data = top_taxa,
                  aes(label = label),
                  size = 3,
                  max.overlaps = 20) +
  theme_minimal() +
  xlab("log2 Fold Change (CBD vs. Vehicle)") +
  ylab("-log10 adjusted p-value") +
  ggtitle("Differentially Abundant Taxa - Infected Only") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  theme(plot.title = element_text(hjust = 0.5))

# save plot
ggsave("/Users/caseymeili/Downloads/volcano_plo_cbd_d7t_infected.jpg", plot = last_plot(), width = 6, height = 5, dpi = 300)          

# filter all significant taxa and sort by significance
sig_taxa <- res_tax %>%
  filter(padj < 0.05) %>%
  arrange(padj)  

# export to CSV
write.csv(sig_taxa, "/Users/caseymeili/Downloads/significant_taxa_cbd_7dpi.csv",
          row.names = TRUE)

