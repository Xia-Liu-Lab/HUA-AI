# 20230101
# Ludi Liu
### ---------------------------
## ---- AI-related features ----
##---- UA_species -------
AI <- read.csv()  #AI abundance file
fea <- read.csv()  #pathway/metabolites file
data <- cbind(AI,fea)
#--glm
lm_result <- as.data.frame(array(dim = c(length(colnames(fea)),3)))
colnames(lm_result) <- c("NAME","coef","p")
lm_result$NAME <- colnames(fea)

lm_result <- as.data.frame(array(dim = c(length(colnames(fea)),3)))
colnames(lm_result) <- c("NAME","coef","p")
lm_result$NAME <- colnames(fea)

for (i in 1:length(colnames(fea))) {
  a<- glm(AI~ data[,i+1], data = data, family = gaussian) #"binomial"Age+Sex+BMI++Cr
  c <- summary(a)
  lm_result$coef[i] <-c[["coefficients"]][2,1]
  lm_result$p[i] <-c[["coefficients"]][2,4]
}

lm_result<- lm_result[order(lm_result$p),]

write.csv(lm_result,"AI_related.csv")


## ------ SHAP-random forest -------
data <- read.csv()  #pathway or metabolites file with group information
sub <- sample(1:164,123)
train <- data[sub,]
test <- data[-sub,]
# modeling
rf <- randomForest(group~., data=data)
# explanation
rf_exp <- DALEX::explain(rf,
                         data = data[,-1],
                         y=data$group==1,
                         label = "randomForest")data[1,]
# shap
shap_henry <- predict_parts(explainer = rf_exp, 
                            new_observation = data[1,], 
                            type = "shap",
                            B = 25)
plot(shap_henry, show_boxplots = FALSE) 

## shapviz
shp <- shapviz(shap_henry)

# pic
sv_importance(shp, kind = "bar", fill = "#187FB5")+
  theme_test()

s <- shp[["S"]]
write.csv(s,"shapvalue.csv")


## ------ multiblock partial least squares discriminant analysis -----
library(mixOmics)

group <- read.csv()  #group file
pri <- read.csv() #AI-related pathways and metabolites
pri <- as.matrix(pri)
colnames(pri)
path <- pri[3:26]
mebo <- pri[27:39]
X1 <- list(
  path = pri[3:26], 
  mebo = pri[27:39])
Y1 <- as.factor(group$group)
summary(Y1)

list.keepX1 <- list(species = c(1, 1), path = c(18,5), protein = c(5, 5))
MyResult.diablo1 <- block.splsda(X1, Y1)
MyResult.diablo2 <- block.plsda(X1, Y1)
plotIndiv(MyResult.diablo1) ## sample plot
plotVar(MyResult.diablo1,var.names = c(TRUE, TRUE),legend = TRUE, style = 'graphics') ## variable plot
circosPlot(MyResult.diablo1, cutoff=0.2,size.variables = 1) 
