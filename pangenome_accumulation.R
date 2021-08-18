# Daniel Wüthrich

require(grid)
require(ggplot2)

args = commandArgs(trailingOnly=TRUE)

pdf(paste(args[1],".pdf", sep = ""),15,10)
df=read.table(args[1],sep="\t",header=TRUE)


dodge <- position_dodge(width = 0.9)

g <- ggplot(df, aes(x=number_of_strains, y=number_of_orthologous_cluster, color=replicate, group=replicate, fill=replicate))
g <- g + geom_line(position=dodge, stat="identity")
g <- g + geom_errorbar(aes(ymax=number_of_orthologous_cluster+sd, ymin=number_of_orthologous_cluster-sd), position = dodge,width = 0.25)
#g <- g +  facet_grid(. ~ sample)
#g <- g + stat_function(fun=function(x)x^2, geom="line", aes(colour="blue")) # addaditionalyline
g

dev.off()
