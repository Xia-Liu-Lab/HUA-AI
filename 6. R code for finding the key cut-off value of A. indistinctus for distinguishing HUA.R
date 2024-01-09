# 20230101
# Ludi Liu
### ---------------------------
cli <- read.csv()  #clinical file of cohort
ai <- read.csv()  #0.0001 unit-wide window file ranged from 0-3.3478
a <- ai$AI
or <- as.data.frame(array(dim = c(length(a),5)))
colnames(or) <- c("cutoff", "OR", "LCI","UCI","P")

for (i in 1:length(a)) {
  cut <- a[i]
  or[i,1] <- cut
  temp <- cli
  temp$g_AI <- ifelse(temp$AI<cut,1,0)
  log <- glm(group~g_AI,family=binomial(link='logit'),data=temp) 
  x <- summary(log)
  or[i,2] <- exp(coef(log))[2]
  or[i,3] <- exp(confint(log))[2,1]
  or[i,4] <- exp(confint(log))[2,2]
  or[i,5] <- x[["coefficients"]][2,4]
}

write.csv(or,"cutoff.csv")

# figure
library(ggplot2)
library(ggalt)
data <- or[which(or$P<0.05),]
colnames(data)
p1<-ggplot(data, aes(x=cutoff, y=OR)) + 
  geom_xspline(size = 2, lineend = "round") + 
  geom_point(size=1)+
  theme(axis.title=element_text(size = 20),
        axis.text=element_text(size = 20),
        legend.text=element_text(size = 20),
        legend.title=element_text(size = 20),
        legend.position="top")+
  theme_test()
p1