# 20230101
# Ludi Liu

### ---------------------------
#--------table----------
#install.packages("tableone") 
library(tableone)
cli <- read.csv("1. clinical data of discovery cohort.csv",header=T,sep=",",row.names=1)
dput(names(cli))

myVars <- c("Sex", "Age", "BMI", "SBP", "DBP", "TG", "TC", "HDL", 
            "LDL", "GLU", "UA", "DD")
catVars <- c("Sex")
## Create a TableOne object(????̬)
tab1 <- CreateTableOne(vars = myVars, data = cli, factorVars = catVars)
print(tab1)
biomarkers <- c("Sex", "Age", "DBP", "TG", "TC", "HDL", 
                "LDL", "GLU", "UA")
table1 <- print(tab1, nonnormal = biomarkers)

## addOverall ????Overall??Ϣ
tab2 <-CreateTableOne(vars = myVars, strata = "group" , data = cli, factorVars = catVars,
                      addOverall = TRUE )
tab2

## exact????fisher??ȷ?????ı?��??wilcox???÷ǲμ????ı?��
table2 <- print(tab2, nonnormal = biomarkers, exact = catVars)
write.csv(table2,"table.csv")
### ---------------------------


### ---------------------------
# ----------- species ------------
# ------------ binary ------------
# ------ SHAP-random forest ------
data <- read.csv()  #Species file with group information
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

# ---------- continuous ------------
# -------- two-part model ----------
#### Packages loading
options("scipen"=999, "digits"=4)
library("dplyr")
library("ggplot2")
library("poolr")
library("tidyverse")
library("caret")
library("utils")
library("rsq")
library("ggthemes")
library("RColorBrewer")
library("ggsci")
library("patchwork")

#### Profiles loading
# Loading of group list and related metadata information
group_list <- read.csv("1. clinical data of discovery cohort.csv",row.names = 1,check.names = F)
# Loading of species data matrix
species_raw_table <- read.csv()  #Species file
# Transpose of the data matrix
species_raw_table <- as.data.frame(t(species_raw_table))

#### Construct a detect matrix (if 0 then "undetected", if 1 then detected)
# A detect matrix of species
species_detect <- species_raw_table
n <- c(1:ncol(species_raw_table))
p <- c(1:nrow(species_raw_table))
for (i in n) {
  for (j in p) {
    if(species_raw_table[j,i]>0){
      species_detect[j,i] <- 1
    }
  }
}

#### Select microbe with presence association with sleep score
dput(names(group_list))
UA_list <- subset(group_list,select = c("UA"))
species_detect <- merge(UA_list,species_detect,by="row.names")
row.names(species_detect) <- species_detect[,1]
species_detect <- species_detect[,-grep("Row.names",colnames(species_detect))]
species_assoscation_detect <- data.frame(ID=NA,beta_presence=NA,intercept=NA,t_detect=NA,p_species_raw_detect=NA)
n <- c(7:ncol(species_detect))
for (i in n) {
  a <- glm(UA ~ species_detect[,i]+Cr,data = species_detect)
  b <- summary(a)
  e1 <- b[["coefficients"]][1,1]
  c <- b[["coefficients"]][2,3]
  beta <- b[["coefficients"]][2,1]
  p_species <- b[["coefficients"]][2,4]
  ID <- colnames(species_detect)[i]
  species_assoscation_detect_temp <- data.frame(ID=ID,beta_presence=beta,intercept=e1,t_detect=c,p_species_raw_detect=p_species)
  species_assoscation_detect <- rbind(species_assoscation_detect,species_assoscation_detect_temp)
}
species_assoscation_detect <- species_assoscation_detect[-1,]

#### Select microbe with abundance association with sleep score
species_abundance <- species_raw_table
species_abundance <- merge(UA_list,species_abundance,by="row.names")
row.names(species_abundance) <- species_abundance[,1]
species_abundance <- species_abundance[,-grep("Row.names",colnames(species_abundance))]
species_assoscation_abundance <- data.frame(ID=NA,beta_abundance=NA,intercept=NA,t_abundance=NA,p_species_raw_abundance=NA)
n <- c(7:ncol(species_abundance))
for (i in n) {
  temp <- species_abundance
  temp <- temp[which(temp[,i]>0),]
  temp[,i] <- -(1/log10(temp[,i]))
  a <- glm(UA ~ temp[,i]+Cr,data = temp)
  b <- summary(a)
  e1 <- b[["coefficients"]][1,1]
  c <- b[["coefficients"]][2,3]
  beta <- b[["coefficients"]][2,1]
  p_species <- b[["coefficients"]][2,4]
  ID <- colnames(temp)[i]
  species_assoscation_abundance_temp <- data.frame(ID=ID,beta_abundance=beta,intercept=e1,t_abundance=c,p_species_raw_abundance=p_species)
  species_assoscation_abundance <- rbind(species_assoscation_abundance,species_assoscation_abundance_temp)
}
species_assoscation_abundance <- species_assoscation_abundance[-1,]

#### Calculate their combined p value
two_stage_p_species_temp1 <- subset(species_assoscation_detect,select = c("ID","p_species_raw_detect","beta_presence","t_detect"))
two_stage_p_species_temp2 <- subset(species_assoscation_abundance,select = c("ID","p_species_raw_abundance","beta_abundance","t_abundance"))
two_stage_p_species <- merge(two_stage_p_species_temp1,two_stage_p_species_temp2,by="ID")
two_stage_p_species$combined_p_species <- NA
n <- c(1:nrow(two_stage_p_species))
for (i in n) {
  z1 <- abs(qnorm(two_stage_p_species$p_species_raw_detect[i]/2))
  z2 <- abs(qnorm(two_stage_p_species$p_species_raw_abundance[i]/2))
  z <- (z1+z2)/sqrt(2)
  p <- (1-pnorm(abs(z)))*2
  temp <- p
  two_stage_p_species$combined_p_species[i] <- temp
}

#### Summary the associated microbe
length(which(two_stage_p_species$combined_p_species<0.05))
two_stage_p_species_diff <- two_stage_p_species[which(two_stage_p_species$combined_p_species<0.05),]
diff_list <- two_stage_p_species_diff$ID

#### Calculate the permutation importance of the final p
species_two_stage_diff_table <- species_raw_table[,colnames(species_raw_table)%in%diff_list]
permutated_results <- as.data.frame(array(dim=c(100,123)))
colnames(permutated_results) <- two_stage_p_species_diff$ID
n <- c(1:ncol(species_two_stage_diff_table))
p <- c(1:100)
for (i in n) {
  for (j in p) {
    permutated_data <- species_two_stage_diff_table
    permutated_data[,i] <- sample(species_two_stage_diff_table[,i])
    permutated_data_detect <- permutated_data
    permutated_data_abundance <- permutated_data
    q <- c(1:ncol(permutated_data_detect))
    r <- c(1:nrow(permutated_data_detect))
    for (l in q) {
      for (z in r) {
        if(permutated_data_detect[z,l]>0){
          permutated_data_detect[z,l] <- 1
        }
      }
    }
    #1.
    UA_list <- subset(group_list,select = c("UA"))
    permutated_data_detect <- merge(UA_list,permutated_data_detect,by="row.names")
    row.names(permutated_data_detect) <- permutated_data_detect[,1]
    permutated_data_detect <- permutated_data_detect[,-grep("Row.names",colnames(permutated_data_detect))]
    p1_temp <- summary(glm(UA ~ permutated_data_detect[,(i+1)],data = permutated_data_detect))[["coefficients"]][2,4]
    #2.
    UA_list <- subset(group_list,select = c("UA"))
    permutated_data_abundance <- merge(UA_list,permutated_data_abundance,by="row.names")
    row.names(permutated_data_abundance) <- permutated_data_abundance[,1]
    permutated_data_abundance <- permutated_data_abundance[,-grep("Row.names",colnames(permutated_data_abundance))]
    permutated_data_abundance <- permutated_data_abundance[-which(permutated_data_abundance[,i+1]==0),]
    permutated_data_abundance[,(i+1)] <- -(1/log10(permutated_data_abundance[,(i+1)]))
    p2_temp <- summary(glm(UA ~ permutated_data_abundance[,(i+1)],data = permutated_data_abundance))[["coefficients"]][2,4]
    #3.
    z1 <- abs(qnorm(p1_temp/2))
    z2 <- abs(qnorm(p2_temp/2))
    z <- (z1+z2)/sqrt(2)
    p3_temp <- (1-pnorm(abs(z)))*2
    #4.
    pfinal_temp <- min(p1_temp,p2_temp,p3_temp)
    permutated_results[j,i] <- pfinal_temp
  }
}

write.csv(two_stage_p_species,"two_stage_p_species.csv",row.names = F,na = "")
write.csv(two_stage_p_species_diff,"two_stage_p_species_sig.csv",row.names = F,na = "")

# ---------- validation ------------
#--------table----------
library(tableone)
cli <- read.csv("2. clinical data of Guangzhou validation cohort.csv",header=T,sep=",",row.names=1)

dput(names(cli))
myVars <- c("sex", "age", "BMI", "SBP", "DBP", "TG", "TC", "HDL", 
            "LDL", "Glu","UA")

catVars <- c("sex")
## Create a TableOne object
tab1 <- CreateTableOne(vars = myVars, data = cli, factorVars = catVars)
print(tab1)
biomarkers <- c("sex","age", "BMI", "SBP", "DBP","TG", "TC", "HDL","Glu", "UA") 

table1 <- print(tab1, nonnormal = biomarkers)

tab2 <-CreateTableOne(vars = myVars, strata = "group" , data = cli, factorVars = catVars,
                      addOverall = TRUE )
tab2

table2 <- print(tab2, nonnormal = biomarkers, exact = catVars)
write.csv(table2,"table_guangzhou.csv")

#------continuous-------
a <- glm(UA ~ AI , data = data, family = gaussian)
summary(a)
#figure
ggplot(data, aes(log10AI,UA))+
  geom_point(size=4,alpha=0.3,color="#778899")+
  geom_smooth(method = "lm", formula = y~x, color = "#1E90FF", fill = "#B0C4DE")+ #颜色选自https://colorbrewer2.org/
  theme_test()+
  theme(
    panel.grid.major = element_blank(),panel.grid.minor = element_blank()
  )