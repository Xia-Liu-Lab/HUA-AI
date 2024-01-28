# 20230101
# Ludi Liu
### -----------------------
## ------ mediation -------
library(mediation)
data <- read.csv()  #mediaton profile
a <- glm(HA~AI,data = data)
summary(a)
b <- glm(UA~HA+AI,data = data)
summary(b)

set.seed(123)
result = mediate(a,b,treat = "AI", mediator = "HA", boot = T)  #AI\HA can be changed according to the content
summary(result)