# Klebsiella-spp
This folder contains the scripts that were used for the Klebsiella spp. study.


# calculate_pangenome_repetetive.py

This script performs repetetive the size of the Pan and Core genome of bacterial strains. For every repetition the order of the added strains is randomize.

As input gene_presence_absence.Rtab from roary (https://sanger-pathogens.github.io/Roary/) is needed.

# AMR_heatmap.R
This script generates a heatmap displaying the occurrence of AMR genes screened using Kleborate
Input files are 
(i) the phylogenetic tree of dereplicated Klebsiella assemblies (n=548) (core genome alignment from roary, tree by Fasttee) and 
(ii) the output from Kleborate screening for AMR genes 
(iii) a file 'species.csv' displaying whcih strain belongs to which Klebsiella species (needed for coloring)

# heatmap_masses.R
This script generates a heatmap, displaying the frequency of ribosomal mass alleles, predicted from Klebsiella species genomic data.
Inputs are 
(i) The predicted masses in binary form, displaying which Klebsiella assembly (rows) encodes which mass alleles (columns) and 
(ii) The file 'Table_S1_Assemblies' which lists all assemblies and also includes the information which strains were included in the genomic analysis (dereplicated, excluding strains which shared > 99.9% ANIm with another strain in the selection)

# KandO_visualize.R
This script visualises which K and O loci have been identified per Klebsiella species. It uses the file KandO_final_table_fused.tsv as input containg the species ID and the K and O loci as identified by Kleborate (v0.3.0) (https://github.com/katholt/Kleborate) 

# plasmid_count_visualize.R
This script visualises how many plasmid replicons have been identified per Klebsiella species. It uses plasmid_count_final_table.tsv as input which lists the species and the number of plasmid replicons identified by Plasmidfinder (https://cge.cbs.dtu.dk/services/PlasmidFinder/)

# plasmid_presence_visualize.R
This script visualises which plasmid replicons have been identified per Klebsiella species.  It uses plasmid_presence_final_table.tsv as input table which lists which plasmids have been identified per strain using Plasmidfinder (https://cge.cbs.dtu.dk/services/PlasmidFinder/)
