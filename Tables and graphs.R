# This R file is used to generate Figures and table related to simulation studies.
# You need to run simulations first and write out the results with the last line of Rcode
savepdf <- function(file, width=16, height=10, pointsize = 10) {
  fname <- file
  pdf(fname, width=width/2.54, height=height/2.54,
      pointsize=pointsize)
  par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
}


#please set working directory to the parent file first.

library(RColorBrewer)
library(reshape2)
library(gridExtra)
library(latex2exp)
library(ggplot2)
library(cowplot)
library(stargazer)


###load results for no outliers
data <- read.table('p50n100no.txt',header = F)
colnames(data)=rep(c("Lasso","FL","AFL","RFL.H","RAFL.H","RFL.T","RAFL.T"),5)

data1 <- data[,1:7]
data2 <- data[,8:14]
data3 <- data[,15:21]
data4 <- data[,22:28]
data5 <- data[,29:35]

meltdata1 <- melt(data1)
meltdata2 <- melt(data2)
meltdata3 <- melt(data3)
meltdata4 <- melt(data4)
meltdata5 <- melt(data5[,-1])

p1 <- ggplot(meltdata1, aes(x=factor(variable), y=value)) + 
  geom_boxplot(fill=c("grey","#1B9E77","#1B9E77","pink","pink","#619CFF","#619CFF")) +
  ggtitle("") + 
  labs(y = "Euclidean Distance") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45,hjust=1), axis.text=element_text(size=12), axis.title.x=element_blank(), axis.title.y=element_text(size=16))

p2 <- ggplot(meltdata4, aes(x=factor(variable), y=value)) + 
  geom_boxplot(fill=c("grey","#1B9E77","#1B9E77","pink","pink","#619CFF","#619CFF")) +
  ggtitle("") + 
  labs(y = "MSPE") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45,hjust=1), axis.text=element_text(size=12), axis.title.x=element_blank(), axis.title.y=element_text(size=16))

p3<- ggplot(meltdata5, aes(x=factor(variable), y=value)) + 
  geom_boxplot(fill=c("#1B9E77","#1B9E77","pink","pink","#619CFF","#619CFF")) +
  ggtitle("") + 
  labs(y = "ARI") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45,hjust=1), axis.text=element_text(size=12), axis.title.x=element_blank(), axis.title.y=element_text(size=16))


noout <- plot_grid(p1, p2, p3, nrow = 1,align="v")
ggsave("p50n100no.pdf", plot = noout , width = 30, height = 10, units = "cm")

part1 <- format(round(apply(data,2,mean),3), nsmall=3)
table1<- matrix(part1,nrow=7,byrow=F)


###load results for y outliers
data <- read.table('p50n100y.txt',header = F)
colnames(data)=rep(c("Lasso","FL","AFL","RFL.H","RAFL.H","RFL.T","RAFL.T"),5)

data1 <- data[,1:7]
data2 <- data[,8:14]
data3 <- data[,15:21]
data4 <- data[,22:28]
data5 <- data[,29:35]

meltdata1 <- melt(data1)
meltdata2 <- melt(data2)
meltdata3 <- melt(data3)
meltdata4 <- melt(data4)
meltdata5 <- melt(data5[,-1])



p4 <- ggplot(meltdata1, aes(x=factor(variable), y=value)) + 
  geom_boxplot(fill=c("grey","#1B9E77","#1B9E77","pink","pink","#619CFF","#619CFF")) +
  ggtitle("") + 
  labs(y = "Euclidean Distance") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45,hjust=1), axis.text=element_text(size=12), axis.title.x=element_blank(), axis.title.y=element_text(size=16))

p5 <- ggplot(meltdata4, aes(x=factor(variable), y=value)) + 
  geom_boxplot(fill=c("grey","#1B9E77","#1B9E77","pink","pink","#619CFF","#619CFF")) +
  ggtitle("") + 
  labs(y = "MSPE") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45,hjust=1), axis.text=element_text(size=12), axis.title.x=element_blank(), axis.title.y=element_text(size=16))

p6<- ggplot(meltdata5, aes(x=factor(variable), y=value)) + 
  geom_boxplot(fill=c("#1B9E77","#1B9E77","pink","pink","#619CFF","#619CFF")) +
  ggtitle("") + 
  labs(y = "ARI") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45,hjust=1), axis.text=element_text(size=12), axis.title.x=element_blank(), axis.title.y=element_text(size=16))

yout <- plot_grid(p4, p5, p6, nrow = 1,align="v")
ggsave("p50n100y.pdf", plot = yout , width = 30, height = 10, units = "cm")


part1 <- format(round(apply(data,2,mean),3), nsmall=3)
table2<- matrix(part1,nrow=7,byrow=F)


###load results for xy outliers

data <- read.table('p50n100xy.txt',header = F)
colnames(data)=rep(c("Lasso","FL","AFL","RFL.H","RAFL.H","RFL.T","RAFL.T"),5)

data1 <- data[,1:7]
data2 <- data[,8:14]
data3 <- data[,15:21]
data4 <- data[,22:28]
data5 <- data[,29:35]

meltdata1 <- melt(data1)
meltdata2 <- melt(data2)
meltdata3 <- melt(data3)
meltdata4 <- melt(data4)
meltdata5 <- melt(data5[,-1])



p7 <- ggplot(meltdata1, aes(x=factor(variable), y=value)) + 
  geom_boxplot(fill=c("grey","#1B9E77","#1B9E77","pink","pink","#619CFF","#619CFF")) +
  ggtitle("") + 
  labs(y = "Euclidean Distance") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45,hjust=1), axis.text=element_text(size=12), axis.title.x=element_blank(), axis.title.y=element_text(size=16))

p8 <- ggplot(meltdata4, aes(x=factor(variable), y=value)) + 
  geom_boxplot(fill=c("grey","#1B9E77","#1B9E77","pink","pink","#619CFF","#619CFF")) +
  ggtitle("") + 
  labs(y = "MSPE") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45,hjust=1), axis.text=element_text(size=12), axis.title.x=element_blank(), axis.title.y=element_text(size=16))

p9<- ggplot(meltdata5, aes(x=factor(variable), y=value)) + 
  geom_boxplot(fill=c("#1B9E77","#1B9E77","pink","pink","#619CFF","#619CFF")) +
  ggtitle("") + 
  labs(y = "ARI") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45,hjust=1), axis.text=element_text(size=12), axis.title.x=element_blank(), axis.title.y=element_text(size=16))

xyout <- plot_grid(p7, p8, p9, nrow = 1,align="v")
ggsave("p50n100xy.pdf", plot = xyout , width = 30, height = 10, units = "cm")


allplot <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, align="v")
ggsave("p50n100all.pdf", plot = allplot , width = 35, height = 30, units = "cm")


part1 <- format(round(apply(data,2,mean),3), nsmall=3)
table3<- matrix(part1,nrow=7,byrow=F)



##################################################
#Table
###############################################
Name <- c("Scenario 1",rep(NA,6),"Scenario 2",rep(NA,6),"Scenario 3",rep(NA,6))
Method <- rep(c("Lasso","FL","AFL","RFL.H","RAFL.H","RFL.T","RAFL.T"),3)

tablefinal <- data.frame(Name,Method,rbind(table1,table2,table3))
colnames(tablefinal)=c("","Method","ED","TPR","TNR","MSPE","ARI")

stargazer(tablefinal, digits=3, summary = F, rownames=F, notes.align = 'r')

