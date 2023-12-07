#This repository contains scripts for analyses related to the study "The Role of the Gut Microbiome in Anorexia Nervosa: Insights from a Conventional, Antibiotic-treated, and Germ-free Mouse Model".

01.demultiplexing.txt - Bash commands for demultiplexing fastq files provided by the sequencing facility.

02.dada_script.R - Quality filtering and denoising of demultiplexed fastq files with data2.

03.prep_phyloseq.R - Assigns the taxonomy and creates a bacterial phylogeny. Identifies and eliminates contaminating ASVs (chimeras).

04.Duplicate merger.R - Merges the dada2 outputs of the two sequencing runs into a single phyloseq object. 

05.Statistics.R - Script used for statistical analysis.


