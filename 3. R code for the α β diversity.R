# 20230101
# Ludi Liu
### ---------------------------
#--------diversity----------
#title: "Shannon Index"
# library
library(ggplot2)
library(ggpubr)
library(agricolae)
library(vegan)

# file read & reshape
pri <- read.csv()  # Genus file
pri <- as.data.frame(t(pri))
# shannon calculate
shannon <- data.frame(diversity(pri,"shannon"))
names(shannon) <- "shannon"
shannon$MetaID <- row.names(shannon)
# add group
group <- read.csv()  #Grouping file on HUA
mer <- merge(shannon, group, by = "MetaID")

# visible pic
library(gghalves)
#
p<-ggplot(mer,aes(x= group, y= shannon,fill=group,color=group))+
  scale_fill_manual(values = c("#44757a", "#d44c3c"))+
  scale_colour_manual(values = c("#44757a", "#d44c3c"))
p
#
p1<-p+geom_half_violin(position=position_nudge(x=0.1,y=0),
                       side='R',adjust=1.2,trim=F,color=NA,alpha=0.5)
p1
#
p2<-p1+geom_half_point(position=position_nudge(x=-0.35,y=0),size =2, 
                       shape =19,range_scale = 0.5,alpha=0.5)
p2
#
p3<-p2+geom_boxplot(outlier.shape = NA, #隐藏离群点；
                    width =0.1,
                    alpha=0.5)
p3
#
p4<-p3+coord_flip()
p4
#
p5 <- p4+theme_bw()+theme(panel.grid=element_blank())
p5

#----------------------------
#title: beta diversity
#llibrary
library(vegan)
library(ggplot2)
library(ggrepel)
library(magrittr)

#read
pri <- read.csv()  # Genus file
pri <- as.data.frame(t(pri))
group <- read.csv()  #Grouping file on HUA
pri$MetaID <- rownames(pri)
mer <- merge(group,pri, by = "MetaID")
pri <- mer[,c(5:(ncol(mer)))] 

group <- data.frame(mer[,c('group')])
row.names(group)=mer[,1]
colnames(group) <- ("group")

##------------NMDS-----------------
sol <- metaMDS(pri)
NMDS = data.frame(MDS1 = sol$points[,1], MDS2 = sol$points[,2])
site <- cbind.data.frame(NMDS, group)
group_average <- aggregate(cbind(MDS1, MDS2)~group, data = site, FUN = mean)

## visible pic
nmds <- site
Fig1a.taxa.pc1.density <-
  ggplot(nmds) +
  geom_density(aes(x=MDS1, group=group, fill=group),
               color="black", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#44757a","#d44c3c")) +
  theme_test() +
  labs(fill="")
Fig1a.taxa.pc1.density

Fig1a.taxa.pc2.density <-
  ggplot(nmds) +
  geom_density(aes(x=MDS2, group=group, fill=group),
               color="black", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#44757a","#d44c3c")) +
  theme_test() +
  labs(fill="") + 
  coord_flip()
Fig1a.taxa.pc2.density