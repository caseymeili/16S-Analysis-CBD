# TMEV 16S Analysis

### Raw Data Processing and Alpha Diversity
Raw 16S data processed using mothur v.1.48.3. Pipeline adapted from mothur MiSeq SOP (accessed August 2025), https://mothur.org/wiki/miseq_sop/. Code from recent analysis of the V3-V4 region (titled cbd-pipeline.sh) is included and can be easily modified for other regions and experiments. 

**Required File for Sequence Processings**
- SILVA database (release 138.2) obtained from: https://mothur.org/wiki/silva_reference_files/
- PDS reference files (version 19) obtained from: https://mothur.org/wiki/rdp_reference_files/


### Beta Diversity Analysis
Shared file from mothur was converted into a Phyloseq object (using the "make phyloseq.R" script). Beta diversity and differential taxa analysis were conducted with R v and RStudio. Primary packages include Phyloseq, vegan, and ggplot2. All packages are listed in the corresponding R scripts.



### Questions 
Questions can be directed to Casey Meili (casey.meili@pharm.utah.edu)



**References**

Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41

Kozich JJ, Westcott SL, Baxter NT, Highlander SK, Schloss PD. (2013): Development of a dual-index sequencing strategy and curation pipeline for analyzing amplicon sequence data on the MiSeq Illumina sequencing platform. Applied and Environmental Microbiology. 79(17):5112-20.


