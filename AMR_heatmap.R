# Aline Cuénod
# This script generates a heatmap displaying the occurrence of AMR genes screened using Kleborate
# Input files are 
# (i) the phylogenetic tree of dereplicated Klebsiella assemblies (n=548) (core genome alignment from roary, tree by Fasttee) and 
# (ii) the output from Kleborate screening for AMR genes (Kleborate was run  as two genomes were missing in the first screen)
# (iii) a file 'species.csv' displaying whcih strain belongs to which Klebsiella species (needed for coloring)

# load pakages
library('ape')
library('stringr')
library('purrr')
library('ggtree')
library('dplyr')

# read in cg - phylogenetic tree (by roary)
nwk<-("./tree.newick")
tree <- read.tree(nwk)

# read Kleborate output
kleborate<-read.csv('./Kleborate_output_all.tab', sep = '\t', header = F)
colnames(kleborate)<-c("strain", "species", "species_match", "contig_count", "N50", "largest_contig", "ST", "virulence_score", "resistance_score", "num_resistance_classes", "num_resistance_genes", "Yersiniabactin", "YbST", "Colibactin", "CbST", "Aerobactin", "AbST", "Salmochelin", "SmST", "rmpA", "rmpA2", "wzi", "K_locus", "K_locus_problems", "K_locus_confidence", "K_locus_identity", "K_locus_missing_genes", "O_locus", "O_locus_problems", "O_locus_confidence", "O_locus_identity", "O_locus_missing_genes", "Chr_ST", "gapA", "infB", "mdh", "pgi", "phoE", "rpoB", "tonB", "ybtS", "ybtX", "ybtQ", "ybtP", "ybtA", "irp2", "irp1", "ybtU", "ybtT", "ybtE", "fyuA", "clbA", "clbB", "clbC", "clbD", "clbE", "clbF", "clbG", "clbH", "clbI", "clbL", "clbM", "clbN", "clbO", "clbP", "clbQ", "AGly", "Col", "Fcyn", "Flq", "Gly", "MLS", "Ntmdz", "Phe", "Rif", "Sul", "Tet", "Tmt", "Bla", "Bla_Carb", "Bla_ESBL", "Bla_ESBL_inhR", "Bla_broad", "Bla_broad_inhR")

# remove empty rows
kleborate<-kleborate[!kleborate$strain == 'strain',]
# the two missing genomes were repeated, add these outputs speperately
kleborate_rep<-read.csv('./Kleborate_output_repeat_2021.tab', sep = '\t', header = T)
# add
kleborate<-merge(kleborate, kleborate_rep, by = colnames(kleborate), all = T)

# list columns associated to antibiotic resistance
res_columns<-c("AGly", "Bla", "Bla_broad", "Bla_broad_inhR", "Bla_Carb", "Bla_ESBL", "Bla_ESBL_inhR", "Fcyn", "Flq", "Gly", "MLS", "Ntmdz", "Phe", "Rif", "Sul", "Tet", "Tmt")
# select these from kleborate output
kleborate_res<-kleborate[,c("strain", "species",res_columns)]

# define function to extract all unique genes per group
uniquegenes<-function(column){
  genes<-paste0(kleborate_res[,column], collapse = ';')
  # * stand for not extact nucleotide hit, but exact aa hit, ? stands for inclomplete hit (https://github.com/katholt/Kleborate/wiki/Antimicrobial-resistance). remove these for figure. remove '-'
  genes<-gsub('\\*|\\?', '',genes)
  genes<-unique(strsplit(genes, ';')[[1]])
  # remove '-'
  genes<-genes[!genes == '-']
  return(genes)
}

# extract all unique genes per group of antibiotics
genes_per_group<-list()
for (i in 1:length(res_columns)){
  genes_per_group[[res_columns[i]]]<-list(sort(uniquegenes(res_columns[i])))
#  print(length(genes_per_group[[res_columns[i]]][[1]]))
}

# For two groups of antibiotics, no genes were found in the selection of assemblies, remove these.
genes_per_group$Gly<-NULL
genes_per_group$Ntmdz<-NULL

# define function to convert the kleborate output into a binary occurence table
to_binary_df<-function(res){
  binary_list<-list()
  # query each gene in df
  for (i in 1:length(genes_per_group[[res]][[1]])){
    binary_list[[genes_per_group[[res]][[1]][i]]]<-grepl(paste0(genes_per_group[[res]][[1]][i],'[^[:alnum:]]*'), kleborate_res[,res])
  }
  df<-as.data.frame(do.call(cbind, binary_list))
  df['strain']<-kleborate_res$strain
  return(df)
}

# remove the antibiotic groups for which no genes were found from the list of columns, associateed to AMR 
no_empty_res<-res_columns[!res_columns %in% c('Gly', 'Ntmdz')]
binary_df_list<-list()

# For all other antibiotics, convert finding too binary output using the function defines above
for (i in 1:length(no_empty_res)){
  binary_df_list[[i]]<-to_binary_df(no_empty_res[i])
}

# merge all binary files using the 'strain' tag as unique identifier
binary_df_all<-  Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "strain", all.x = TRUE),
                        binary_df_list)

# add the information which strain belongs to which species
species<-read.csv('./species.csv')
species['species']<-gsub('\\-.*', '', species$species_used)
species$species_used<-NULL
species$genus<-NULL
# add species information to binary output
binary_df_all<-merge(binary_df_all, species, by.x = 'strain', by.y = 'accession')
# in order to link to tree, use strainnames as rownames
rownames(binary_df_all)<-binary_df_all$strain
binary_df_all$strain<-NULL

# check if grepl worked for all
for (i in 1:length(colnames(binary_df_all))){
  print(paste(i, table(binary_df_all[, colnames(binary_df_all)[i]])))
}

# remove the columns which are only found in 5 or less assemblies
freq<-as.data.frame(colSums(binary_df_all != FALSE))
freq['res']<-rownames(freq)

rare<-freq[freq$`colSums(binary_df_all != FALSE)` <= 5, 'res']


# replace all 'TRUE' enrties with the species name, in order to colour them accordingly in the plot
binary_df_all<-binary_df_all %>% dplyr::mutate_at(vars(!contains('species')), ~ ifelse(.== TRUE, species, 0))
# change the order of columns, bring species to front
binary_df_all<-binary_df_all[,c(ncol(binary_df_all), 2:(ncol(binary_df_all)-1))]

# for the figure, remove all columns with 5 or less occurrences
binary_df_all_sub<-binary_df_all[,!(colnames(binary_df_all) %in% rare)]

#plot the phylogenetic tree 
p <- ggtree(tree) 

# check the nodelabels in order to find which node to reroot on 
# p + geom_text(aes(label=node))
# root on node 779
tree1<-phytools::reroot(tree, node = 779)
# convert into ggtree
p1 <- ggtree(tree1) 
# in order to harmonise occurence of species with the text, flip branches
p2 <- p1 %>% flip(770, 146) %>% flip(858, 857)
p3 <- p2 %>% flip(852, 856) %>% flip(1006, 966)

# define colouring of the species
heatmap.colours <- c("white","darkmagenta","green","turquoise",
                     "darkblue","deeppink","red4","orange","yellow")
names(heatmap.colours) <- c('0', "grimontii", "huaxiensis", "michiganensis", "oxytoca", "pneumoniae", "quasipneumoniae", "quasivariicola", "variicola")

# Plot the heatmap to the phylogenetic tree
heatmap<-gheatmap(p3, binary_df_all_sub, color=NULL, 
         colnames_position="top", 
         colnames_angle=90, offset=0.001, width=1.5,
         hjust=0, font.size=2) +
  scale_fill_manual(values=heatmap.colours, breaks=names(heatmap.colours)) + 
  theme(legend.position = 'none')+ ylim(NA, 565) + 
  geom_treescale()

# export the figure
pdf('./heatmap.pdf', width = 12, height = 7)
heatmap
dev.off()

# count genes for per species
kleborate_res_count<-kleborate[,c("strain", res_columns)]
kleborate_res_count<-kleborate_res_count %>% dplyr::mutate_at(vars(!contains('strain')), ~ ifelse(.!= '-', 1, 0))
kleborate_res_count<-merge(kleborate_res_count, species, by.x = 'strain', by.y = 'accession')

# sum up to how many AB classes genes that confer res were detected
kleborate_res_count['sum']<-rowSums(kleborate_res_count[,!colnames(kleborate_res_count) %in% c('species', 'strain')])
# add column 'group' to summarise
kleborate_res_count['group']<-ifelse(kleborate_res_count$species %in% c('pneumoniae', 'quasipneumoniae', 'quasivariicola', 'variicola'), 'pneumoniae',
                                     ifelse(kleborate_res_count$species == 'huaxiensis', NA, 'oxytoca')) 

# find median and IQR of occurence of resistance per Klebsiella 'groups'
kleborate_res_count_group<-kleborate_res_count %>% 
  group_by(group) %>%
  summarise(median_AB_res = median(sum),  Q1=quantile(sum, probs = 0.25), Q3=quantile(sum, probs = 0.75))

# find median and IQR of occurence of resistance per Klebsiella 'species'
kleborate_res_count_species<-kleborate_res_count %>% 
  group_by(species) %>%
  summarise(median_AB_res = median(sum),  Q1=quantile(sum, probs = 0.25), Q3=quantile(sum, probs = 0.75))

# count the of resistance per Klebsiella 'groups'
check_group_all<-kleborate_res_count %>%
  group_by(group) %>% 
  summarise(total = n(), 
            across(!starts_with('s'), ~ mean(. == 1)))

# count the of resistance per Klebsiella species
check_species_all<-kleborate_res_count %>%
  group_by(species) %>% 
  summarise(total = n(), 
            across(!starts_with('s'), ~ mean(. == 1)))
