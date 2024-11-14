# load the package (only needed once per session)
library(lavaan)
library(MASS) 
source("Goodness_of_fit_test.R")

data1 <- read.table("Data1_HS1939.txt", header=T)
X <- as.matrix(data1)


# specify the model
HS.model <- ' visual  =~ x1 + x2 + x3      
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

# fit the model
fit1 <- cfa(HS.model, data = data1, likelihood = "wishart", estimator = "ML")
fit2 <- cfa(HS.model, data = data1, likelihood = "wishart", estimator = "MLM")
fit3 <- cfa(HS.model, data = data1, likelihood = "wishart", estimator = "MLMVS")

summary(fit1)
summary(fit2)
summary(fit3)

# Summary of goodness of fit tests and fit indices
GFT(fit1,fit2,fit3,X)

