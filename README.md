# Klebsiella-spp
This folder contains the scripts that were used for the Klebsiella spp. study.


# calculate_pangenome_repetetive.py

This script performs repetetive the size of the Pan and Core genome of bacterial strains. For every repetition the order of the added strains is randomize.

As input gene_presence_absence.Rtab from roary (https://sanger-pathogens.github.io/Roary/) is needed.

# AMR_heatmap.R
This script generates a heatmap displaying the occurrence of AMR genes screened using Kleborate
Input files are 
(i) the phylogenetic tree of dereplicated Klebsiella assemblies (n=548) (core genome alignment from roary, tree by Fasttee) and 
(ii) the output from Kleborate screening for AMR genes (Kleborate was run  as two genomes were missing in the first screen)
(iii) a file 'species.csv' displaying whcih strain belongs to which Klebsiella species (needed for coloring)

# heatmap_masses.R
This script generates a heatmap, displaying the frequency of ribosomal mass alleles, predicted from Klebsiella species genomic data.
Inputs are 
(i) The predicted masses in binary form, displaying which Klebsiella assembly (rows) encodes which mass alleles (columns) and 
(ii) The file 'Table_S1_Assemblies' which lists all assemblies and also includes the information which strains were included in the genomic analysis (dereplicated, excluding strains which shared > 99.9% ANIm with another strain in the selection)

