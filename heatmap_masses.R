# Aline Cuénod
# This script generates a heatmap, displaying the frequency of ribosomal mass alleles, predicted from Klebsiella species genomic data.
# Inputs are 
# (i) The predicted masses in binary form, displaying which Klebsiella assembly (rows) encodes which mass alleles (columns) and 
# (ii) The file 'Table_S1_Assemblies' which lists all assemblies and also includes the information which strains were included in the genomic analysis (dereplicated, excluding strains which shared > 99.9% ANIm with another strain in the selection)

#Load packages
library('dplyr')
library('ggplot2')
library('tidyr')

# import predicted ribosomal mass alleles
binary_full<-read.csv('./2021-05-Table_S4_predicted_masses.csv', sep = ',')
# change columnname
colnames(binary_full)[1] <- "Strain"
# harmonise strainnames
binary_full$Strain<-gsub('Klebsiella\\_huaxensis_WCHKl090001\\(T\\)', 'Klebsiella_huaxensis_WCHKl090001', binary_full$Strain)

# import list of assemblies
assemblies<-read.csv('./Table_S1_Assemblies.csv', sep = ',')
colnames(assemblies)<-assemblies[1,]
assemblies<-assemblies[-1,]
# select dereplicated strains which have been included in the genomic analysis
sel <- assemblies[assemblies$`Included in comperative genomics analysis` == TRUE,]

# Check which are missing in binary table 
# These had an incomplete set of predicted ribosomal mass alleles and where therefor excluded from the binary table
assemblies['in_binary_table_ribos']<-assemblies$Strain %in% binary_full$Strain

# reduce binary table to strains also included in genomic analysis in order to avoid bias being introduced by closely related strains
binary_sel<-binary_full[binary_full$Strain %in% sel$Strain,]

# from the strainname, exctract which strain belongs to which species
binary_sel<-binary_full[binary_full$Strain %in% sel$Strain,]
binary_sel['species'] <- gsub('(^.+\\_.+\\_)(.*)', '\\1', binary_sel$Strain)
binary_sel['species'] <- gsub('-subsp\\-.*$|\\_$', '', binary_sel$species)
binary_sel['species'] <- gsub('\\_', ' ', binary_sel$species)

# change the class of all binary columns to numeric
binary_col <- colnames(binary_sel)[!colnames(binary_sel) %in% c("species", "Asssembly","Strain")]
binary_sel[binary_col] <- sapply(binary_sel[binary_col],as.numeric)

# summarise the occurences of a mass allele by species by calculating the mean value.
binary_sel_sum <- binary_sel %>% 
  group_by(species) %>% 
  summarise(across(all_of(binary_col), mean))

# remove the mass alleles which do not vary between the strains and are the same throughout the entire set of genomes
binary_sel_sum_ <- binary_sel_sum[, colSums(binary_sel_sum != 0) > 0]
binary_sel_sum <-binary_sel_sum_
# set the species as rownames 
rownames(binary_sel_sum)<-binary_sel_sum$species
binary_sel_sum$species<-NULL

# convert to 'longer' format in order to plot
binary_sel_sum_long<-pivot_longer(binary_sel_sum_, cols = colnames(binary_sel_sum_)[!colnames(binary_sel_sum_) %in% 'species'])

# convert 'mass' to numeric in order to display all mass alleles with ascending mass
binary_sel_sum_long['mass']<-as.numeric(as.character(sub('(^.*\\.)(\\d{4,}\\.*\\d*)', '\\2', binary_sel_sum_long$name)))
# order the mass alleles according to their mass
binary_sel_sum_long$name<-factor(binary_sel_sum_long$name, levels = unique(binary_sel_sum_long$name[order(binary_sel_sum_long$mass)]))
# order the species
binary_sel_sum_long$species<-factor(binary_sel_sum_long$species, levels = c('Klebsiella pneumoniae', 'Klebsiella quasipneumoniae', 'Klebsiella variicola', 'Klebsiella quasivariicola', 'Klebsiella oxytoca', 'Klebsiella michiganensis', 'Klebsiella grimontii', 'Klebsiella huaxensis'), 
                                    labels =  c('K. pneumoniae', 'K. quasipneumoniae', 'K. variicola', 'K. quasivariicola', 'K. oxytoca', 'K. michiganensis', 'K. grimontii', 'K. huaxensis'))
# plot the heatmap, displaying the occurence of the differential mass alleles per species
frq_ribos<-ggplot(binary_sel_sum_long, aes(name, species, fill= value)) + 
  geom_tile(colour = 'grey') +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1), axis.text.y = element_text(face = "italic"), legend.position = 'bottom') +
  scale_y_discrete(limits=rev) + scale_fill_gradient(low = 'white', high = 'black') + 
  xlab('') +
  ylab('') + labs(fill = "Fraction of assemblies\nencoding mass allele")

# export the plot
pdf('./freq_ribos.pdf', width = 16, height =3.5)
frq_ribos
dev.off()
