# Daniel Wüthrich

library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
require(scales)



myfunction <- function(name){
 
file_name=name
dat=read.table(file_name,sep="\t",header=TRUE)
ggplot(dat, aes(order = (group),factor=group,color=group)) + 
    geom_segment(aes(x=start, xend=end, y=-strain, yend=-strain), size=0.5)  +
    xlab("Position [bp]") + ylab("") +theme(axis.text.y = element_text(size=0)) +theme(axis.text.x = element_text(angle = 45, hjust = 1))+
theme(legend.position="bottom", legend.direction="horizontal",legend.title=element_blank())+
   theme(panel.background = element_blank()) +scale_color_manual(values=c("gene" = "black","pneumoniae" = "#FF00BF","quasipneumoniae" = "#8A0808","quasivariicola" = "#DF7401","variicola" = "#FACC2E","oxytoca" = "#29088A","michiganensis" = "#0489B1","grimontii" = "#01DFD7","huaxensis" = "#01DF3A"))+ ggtitle(name)+ 
  theme(plot.title = element_text(hjust = 0.5))+ theme(legend.position="none",axis.ticks.y = element_blank())
}

myfunction2 <- function(name){
  
  file_name=paste ("",name,".tab",sep="")
  dat=read.table(file_name,sep="\t",header=TRUE)
  ggplot(dat, aes(order = (group),factor=group,color=group)) + 
    geom_segment(aes(x=start, xend=end, y=strain, yend=strain), size=1.65)  +
    xlab("Position [bp]") + ylab("") +theme(axis.text.y = element_text(size=0))+
    theme(legend.position="bottom", legend.direction="horizontal",legend.title=element_blank())+
    theme(panel.background = element_blank()) +scale_color_manual(values=c("red","green4","purple","darkorange2","blue"))+ ggtitle(name)+ 
    theme(plot.title = element_text(hjust = 0.5))+scale_x_continuous(limits = c(0, 700000))+ theme(legend.position="none",axis.ticks.y = element_blank(),plot.margin=unit(c(0.5,-1,0.5,-1), "cm"))
}

pdf("myOut.pdf",15,15)
#grid.arrange(myfunction2("N1115"),arrangeGrob(myfunction("BD_II")
#  ,myfunction("CAUH35")
#  ,myfunction("LC2W")
#  ,myfunction("LcA")
#  ,myfunction("LcY")
#  ,myfunction("W56"), nrow=1)
#  ,ncol=1, nrow=2)

myfunction("./virulence_data_reordered.txt")

dev.off()

#pdf("myOut_N1115.pdf",15,6.5)
#grid.arrange(myfunction2("N1115")
#             ,ncol=1, nrow=1)

#dev.off()
