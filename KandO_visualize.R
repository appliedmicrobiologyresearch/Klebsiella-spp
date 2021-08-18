#Daniel Wüthrich

#source("https://bioconductor.org/biocLite.R")
#biocLite("ComplexHeatmap")
library(ggplot2)
library(gridExtra)
library(grid)


mat=read.table("./KandO_final_table_fused.tsv", sep="\t",header = TRUE)

pdf(paste("O_K.pdf", sep = ""),17,3)

ggplot(mat, aes(x=locus, y=species, color = species))+geom_count()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5, hjust = 0)) +theme(legend.position="bottom",legend.title = element_blank()) +scale_color_manual(values=c("gene" = "black","pneumoniae" = "#FF00BF","quasipneumoniae" = "#8A0808","quasivariicola" = "#DF7401","variicola" = "#FACC2E","oxytoca" = "#29088A","michiganensis" = "#0489B1","grimontii" = "#01DFD7","huaxensis" = "#01DF3A"))+
   ylim("huaxensis","grimontii","michiganensis","oxytoca","variicola","quasivariicola","quasipneumoniae","pneumoniae" )

dev.off()

