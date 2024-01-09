# 20230101
# Ludi Liu

### ---------------------------
#read file
group <- read.csv()  #Grouping file on HUA/AI
pri <- read.csv()   #Species/pathway/KO/metabolites file

#prevalence筛选
pri$pre <- NA
for (i in 1:nrow(pri)) {
  pri[i,ncol(pri)] <- ncol(pri)-1-sum(pri[i,1:ncol(pri)-1]==0)
}
pri$pre <-  pri$pre/(ncol(pri)-1)*100
pri <- subset(pri,pri$pre>=10)
pri <- pri[,1:ncol(pri)-1]

#trans
pri <- as.data.frame(t(pri))
pri$MetaID<- rownames(pri)

#merge two files and wilcox test
wilcox.p <- NULL
n1 <-"CON"  #or AI_high
n2 <-"HUA"  #or AI_low

for (i in 1:(ncol(pri)-1)){
  mer <- merge(pri[,c(i,ncol(pri))], group, by = "MetaID")
  g1 <- mer[mer$group == n1, 2]
  g2 <- mer[mer$group == n2, 2]
  g1_not0 <- 82-sum(g1==0)  # 82 should be the number of the subjects in this group 
  g2_not0 <- 82-sum(g2==0)
  wilcox.p$metabo[i] <- colnames(mer)[2]
  wilcox.p$CON_not0[i] <- g1_not0 
  wilcox.p$CASE_not0[i] <- g2_not0
  wilcox.p$CON_max[i] <- max(g1)
  wilcox.p$CON_mean[i] <- sum(g1)/82
  wilcox.p$CASE_max[i] <- max(g2)
  wilcox.p$CASE_mean[i] <- sum(g2)/82
  
  ##engroup
  wil <- wilcox.test(x = g1,y = g2, paired = T)
  wilcox.p$P[i] <- wil$p.value
  wil <- wilcox.test(x = g1,y = g2, paired = T,alternative="greater")
  if(wil$p.value < 0.05)
    wilcox.p$engroup[i] <- n1
  else
    wilcox.p$engroup[i] <- "equal"
  wil <- wilcox.test(x = g1,y = g2, paired = T,alternative="less")
  if(wil$p.value < 0.05)
    wilcox.p$engroup[i] <- n2
}
#data.frame
p <- data.frame(wilcox.p)
wilcox.fdr<- p[order(p[,1]),]

#permutation 1000
data <- merge(group,pri, by = "MetaID")
data <- data[,4:ncol(data)]
permutated_res <- as.data.frame(array(dim = c((ncol(data)-1),k+1)))
colnames(permutated_res)[1] <- "MetaID"
permutated_res$MetaID <- colnames(data)[2:ncol(data)]

k = 1000
n <- c(1:(ncol(data)-1))
p <- c(1:k)
set.seed(615)

n1 <-"CON"  #or AI_high
n2 <-"HUA"  #or AI_low
for (i in n) {
  data_temp <- data
  for (j in p) {
    data_temp[,i+1] <- sample(data[,i+1])
    g1 <- data_temp[data_temp$group == n1, i+1]
    g2 <- data_temp[data_temp$group == n2, i+1]
    wil <- wilcox.test(x = g1,y = g2, paired = T)
    permutated_res[i,j+1] <- wil$p.value
  }
}
permutated_res<- permutated_res[order(permutated_res[,1]),]

wilcox.fdr$permutated_p <- NA
n <- c(1:nrow(wilcox.fdr))
for (i in n) {
  freq <- length(which(abs(permutated_res[i,2:ncol(permutated_res)])<abs(wilcox.fdr$P[i])))/k
  wilcox.fdr$permutated_p[i] <- freq
}

write.csv(wilcox.fdr,"permutated1000.csv")

#------- boxblot----------
library(ggplot2)
group <- read.csv()  #Grouping file on HUA/AI
data <- read.csv()  #Species/pathway/KO/metabolites file
pic <- read.csv()  #Features that need to be presented

# transformance
pri <- merge(pic, data, by = "MetaID")
row.names(pri) <- pri[,1]
pri <- pri[,-1]
pri <- as.data.frame(t(pri))
pri <- log10(pri)
pri$MetaID <- rownames(pri)
mer <- merge(pri, group,by = "MetaID")
mer <- mer[,-1]
mer.na <- rbind(colnames(mer),mer)

#add name
name <- NULL

for (i in 1:(ncol(mer)-1)){
  a <- NULL
  a$log10 <- mer[,i]
  a$group <- mer[,ncol(mer)]
  a$name <- mer[,ncol(mer)]
  a <- as.data.frame(a)
  for (j in 1:nrow(a))
    a[j,3] <- mer.na[1,i]
  name <- rbind.data.frame(name,a)
}

## draw pic(ggplot)
ggplot(name, aes(x=name, y=log10))+
  geom_boxplot(aes(fill = group),position = position_dodge(0.7),width=0.5,outlier.colour=NULL,outlier.shape=21,outlier.size=2) +
  labs(title="Species",x="", y="log10(abundance)")+
  scale_fill_manual(values=c("#44757a","#d44c3c"))+ 
  theme(plot.title=element_text(hjust=0.5), legend.title=element_blank(),axis.line = element_line(colour = "black"))+  
  theme_test()