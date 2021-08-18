# Daniel Wüthrich

#source("https://bioconductor.org/biocLite.R")
#biocLite("ComplexHeatmap")
library(ggplot2)
library(gridExtra)
library(grid)


mat=read.table("./plasmid_count_final_table.tsv", sep="\t",header = TRUE)


pdf(paste("plasmid.pdf", sep = ""),12,3.5)

ggplot(mat, aes(x=species, y= Plasmids, color=species )) +geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.title.x=element_blank())+scale_color_manual(values=c("gene" = "black","pneumoniae" = "#FF00BF","quasipneumoniae" = "#8A0808","quasivariicola" = "#DF7401","variicola" = "#FACC2E","oxytoca" = "#29088A","michiganensis" = "#0489B1","grimontii" = "#01DFD7","huaxensis" = "#01DF3A"))+
  xlim("pneumoniae","quasipneumoniae","quasivariicola","variicola","oxytoca","michiganensis" ,"grimontii","huaxensis" )
dev.off()