# japan test
# using ridge and lasso regresssion
# load packages
library(olsrr)
library(glmnet)
library("dplyr")
library("ggpubr")
library(pls)
library(car)
library(leaps)

Fertility <- read.csv("/Users/karen/Desktop/664/project/japan_dataset.csv", header=T)
Fertilityx <- Fertility[,-4]

#regression attempt
lmfull <- lm(fertility~lifexpect+fem.enroll.+urb.popul.
             +inf.mortality+HHDI.CPI, data = Fertility)
summary(lmfull)
#life expectancy is found to be insignificant
ggqqplot(lmfull$residuals)
shapiro.test(lmfull$residuals)
#p=0.4817>0.05, the normality of residuals is significant

#DFFITS testing for influential points of the fullmode
dffitsfull<-as.data.frame(dffits(lmfull))
dffitsfull
thresholdfull<-2*sqrt(5/24) #p=5, n=24 for fullmode
plot(dffits(lmfull), type = 'h')
abline(h = thresholdfull, lty = 2)
abline(h = -thresholdfull, lty = 2)p
#outliers testing with full model
outlierTest(lmfull)
# No. 20 is an outlier
Fertilitynew <- Fertility[c(-20),]
Fertilityxnew <- Fertilitynew[,c(-1,-5)]

#regression attempt
lmfull1 <- lm(fertility~lifexpect+fem.enroll.+urb.popul.
             +inf.mortality+HHDI.CPI, data = Fertilitynew)
summary(lmfull1)
#life expectancy is found to be insignificant
ggqqplot(lmfull1$residuals)
shapiro.test(lmfull1$residuals)
#p=0.1919>0.05, the normality of residuals is significant

#multicollinearity test
car::vif(lmfull1)
# collinearity found

#backward selection
#regression attempt
lmback1 <- lm(fertility~lifexpect+fem.enroll.+urb.popul.
              +HHDI.CPI, data = Fertilitynew)
summary(lmback1)
#life expectancy is found to be insignificant
ggqqplot(lmback1$residuals)
shapiro.test(lmback1$residuals)
#p=0.1554>0.05, the normality of residuals is significant

lmback2 <- lm(fertility~fem.enroll.+urb.popul.
              +HHDI.CPI, data = Fertilitynew)
summary(lmback2)
#life expectancy is found to be insignificant
ggqqplot(lmback2$residuals)
shapiro.test(lmback2$residuals)
#p=0.01334<0.05, the normality of residuals is not significant


#lasso based on the full data
lambda_seq <- 10^seq(3, -3, by = -.1) # lambda sequence
lmlasso <- glmnet(Fertilityxnew, Fertilitynew$fertility, alpha = 1, lambda  = lambda_seq)
# Checking the model
summary(lmlasso)
# Using cross validation glmnet
Fertilityxnew <- data.matrix(Fertilityxnew)
lasso_cv <- cv.glmnet(Fertilityxnew, Fertilitynew$fertility, alpha = 1, lambda  = lambda_seq)
# Best lambda value
best_lambda <- lasso_cv$lambda.min
best_lambda
# the best lambda 0.001995262
lmlassobest <- glmnet(Fertilityxnew, Fertilitynew$fertility, alpha = 1, lambda  = best_lambda)
summary(lmlassobest)
coef(lmlassobest)
ols_vif_tol(lmlassobest)


e <- predict(lmlassobest, Fertilityxnew)-Fertilitynew$fertility
e <- as.numeric(e)
class(e) # list to vector transformation
shapiro.test(e) #0.1759
ggqqplot(e)
# normality is acceptable

#ridge based on the full data
lambda_seq <- 10^seq(4, -4, by = -.001) # lambda sequence
lmridge <- glmnet(Fertilityxnew, Fertilitynew$fertility, alpha = 0, lambda  = lambda_seq)
# Checking the model
summary(lmridge)
# Using cross validation glmnet
Fertilityxnew <- data.matrix(Fertilityxnew)
ridge_cv <- cv.glmnet(Fertilityxnew, Fertilitynew$fertility, alpha = 0, lambda  = lambda_seq)
# Best lambda value
best_lambda2 <- ridge_cv$lambda.min
best_lambda2
# the best lambda 0.0001592209
lmridgebest <- glmnet(Fertilityxnew, Fertilitynew$fertility, alpha = 0, lambda  = best_lambda2)
summary(lmridgebest)
coef(lmridgebest)


lmlassobest <- glmnet(fertilitytrainx, fertilitytrain$fertility, alpha = 1, lambda  = best_lambda2)
summary(lmlassobest)
coef(lmlassobest)


f <- predict(lmridgebest, Fertilityxnew)-Fertilitynew$fertility
f <- as.numeric(f)
class(f) # list to vector transformation
shapiro.test(f) #0.1979
ggqqplot(f)
# normality is acceptable

#PCR
set.seed (1000) # random
lmpcr <- pcr(fertility~lifexpect+fem.enroll.+urb.popul.
             +inf.mortality+HHDI.CPI, data = Fertilitynew,
             scale = TRUE, validation = "CV") 
summary(lmpcr)
# Plot the root mean squared error
validationplot(lmpcr)
# Plot the R2
validationplot(lmpcr, val.type = "R2")
# RMSEP test showed an minimum adjusted CV at ncomp=4, while the fifth component contributes little to %variance
# ncomp=5, prediction of PCR model
predict(lmpcr,Fertilityxnew,ncomp=5)
# best pcr model
lmpcrbest <- pcr(fertility~lifexpect+fem.enroll.+urb.popul.
                 +inf.mortality+HHDI.CPI, data = Fertilitynew,
                 scale = TRUE, validation = "CV", ncomp=5) 
a <- unlist(lmpcrbest$residuals)
class(a) # list to vector transformation
shapiro.test(a) #2.513e-05, not acceptable
ggqqplot(as.numeric(a))

#PLSR
set.seed (100)
lmplsr <- plsr(fertility~lifexpect+fem.enroll.+urb.popul.
               +inf.mortality+HHDI.CPI, data = Fertilitynew, 
               scale = TRUE, validation = "CV")
summary(lmplsr)
validationplot(lmplsr)
validationplot(lmplsr, val.type = "R2")
# Similar to pcr, ncomp=5 should be enough for the prediction
# best plsr model
lmplsrbest <- plsr(fertility~lifexpect+fem.enroll.+urb.popul.
                   +inf.mortality+HHDI.CPI, data = Fertilitynew,
                   scale = TRUE, validation = "CV", ncomp=5) 
# prediction of PLSR model
predict(lmplsrbest, Fertilityxnew, ncomp = 5)
b <- unlist(lmplsrbest$residuals)
class(b) # list to vector transformation
shapiro.test(as.numeric(b)) 
ggqqplot(as.numeric(b))
