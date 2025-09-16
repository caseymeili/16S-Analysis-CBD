#!/bin/bash

# 16S Data Processing (V3-V4)
# mothur v.1.48.2
# adapted from MiSeq SOP (accessed August 2025), https://mothur.org/wiki/miseq_sop/

####################################################################################

# launch mothur, set working directory, and set processors
ml mothur/1.48.2
set.dir(input=/Users/caseymeili/Desktop/26454R/Fastq/cbd) 
set.current(processors=12)

# make stability file (file containing gz files in directory)
make.file(inputdir=/Users/caseymeili/Desktop/26454R/Fastq/cbd, type=gz, prefix=stability) 

# combine two sets of reads for each sample, and then combine the data from all of the samples
make.contigs(file=stability.files)
summary.seqs(fasta=stability.trim.contigs.fasta, count=stability.contigs.count_table) 

# screen to remove sequences longer than 460 bp, shorter than 400 bp, with ambiguous bases, or more than 8 homopolymers
screen.seqs(fasta=stability.trim.contigs.fasta, count=stability.contigs.count_table, maxambig=0, minlength=400, maxlength=460, maxhomop=8) 
summary.seqs(fasta=stability.trim.contigs.good.fasta, count=stability.contigs.good.count_table)

# merge duplicated sequences
unique.seqs(fasta=stability.trim.contigs.good.fasta, count=stability.contigs.good.count_table)
summary.seqs(fasta=stability.trim.contigs.good.unique.fasta, count=stability.trim.contigs.good.count_table)

# create database customized to the region of interest (V3-V4)
# SILVA database (silva.nr_v138_2.align & silva.nr_v138_2.tax) obtained from mothur silva reference files, release 138.2

# determine where primers sit on the alignment using oligos file
pcr.seqs(fasta=silva.nr_v138_2.align, taxonomy=silva.nr_v138_2.tax, oligos=zymo_primers.oligos, keepdots=F)

# rename output files to friendlier names
rename.file(input=silva.nr_v138_2.pcr.align, new=silva.V3V4.pcr.align)
rename.file(input=silva.nr_v138_2.pcr.tax,   new=silva.V3V4.pcr.tax)

# align primer-trimmed sequences to new V3–V4 SILVA reference
align.seqs(fasta=stability.trim.contigs.good.unique.fasta, reference=silva.V3V4.pcr.align)
summary.seqs(fasta=stability.trim.contigs.good.unique.align)

# remove sequences that align outside of expected positions (start=4965, end=21977), usually due to poor alignment or non-specific amplification
screen.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table, start=4965, end=21977)                      

# remove overhangs (probably not necessary since sequencing was paired-end)
filter.seqs(fasta=stability.trim.contigs.good.unique.good.align, vertical=T, trump=.)

# remove redundancy (potentially created by trimming the ends)
unique.seqs(fasta=stability.trim.contigs.good.unique.good.filter.fasta, count=stability.trim.contigs.good.good.count_table)
summary.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.fasta, count=stability.trim.contigs.good.unique.good.filter.count_table)

# precluster allowing for up to 2 differences between sequences
pre.cluster(fasta=stability.trim.contigs.good.unique.good.filter.unique.fasta, count=stability.trim.contigs.good.unique.good.filter.count_table, diffs=2)

# remove chimeras using VSEARCH algorithm
chimera.vsearch(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)

# classify sequences using reference files
# PDS reference files (trainset19_072023.pds.fasta & trainset19_072023.pds.tax) were obtained from mothur RDP reference files version 19 (https://mothur.org/wiki/rdp_reference_files/)
classify.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.count_table, reference=trainset19_072023.pds.fasta, taxonomy=trainset19_072023.pds.tax)

# remove unwanted sequences (18S, archaea, chloroplasts, mitochondria, unknown)
remove.lineage(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pds.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
summary.tax(taxonomy=current, count=current)

# rename files before analysis (so the file names aren't so yucky)
rename.file(fasta=current, count=current, taxonomy=current, prefix=cbd)

# clustering into OTUs
# default cutoff used for clustering (0.03)
dist.seqs(fasta=cbd.fasta, cutoff=0.03)
cluster(column=cbd.dist, count=cbd.count_table)

# determine how many sequences are in each OTU
make.shared(list=cbd.opti_mcc.list, count=cbd.count_table, label=0.03)

# consensus taxonomy for each OTU
# taxonomy file is the output of the classify.seqs command
classify.otu(list=cbd.opti_mcc.list, count=cbd.count_table, taxonomy=cbd.taxonomy, label=0.03)

# make shared file for ASV analysis
make.shared(count=cbd.count_table)
# output is cbd.asv.list and cbd.asv.shared

# generate consensus taxonomy for ASVs
classify.otu(list=cbd.asv.list, count=cbd.count_table, taxonomy=cbd.taxonomy, label=ASV)

# get number of sequences in each sample
count.groups(shared=cbd.opti_mcc.shared)
count.groups(shared=cbd.asv.shared)

# generate subsampled file for analysis
# size of smallest group from count.groups was 3317 (for both OTUs and ASVs)
sub.sample(shared=cbd.opti_mcc.list, size=3371)
sub.sample(shared=cbd.asv.list, size=3371)

# alpha diversity
# generate rarefaction curves for plotting in R
rarefaction.single(shared=cbd.opti_mcc.shared, calc=sobs)

# calculate coverage, number of seqs, and diversity indices (Simpson, inverse Simpson, Shannon)
# be sure to subsample to the size of the smallest group for all alpha diversity analysis
# if the smallest group is too small, discard sample
summary.single(shared=cbd.opti_mcc.shared, calc=nseqs-coverage-sobs-simpson-invsimpson-shannon, subsample=T)
summary.single(shared=cbd.asv.shared, calc=nseqs-coverage-sobs-simpson-invsimpson-shannon, subsample=T)

