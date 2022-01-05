Fertility <- read.csv('/Users/karen/Desktop/664/project/dataset.csv', header=T)

#regression without normality test of y variable
lmfull <- lm(fertility~lifexpect+univ.recruit.+urb.popul.+inf.mortality+income.CPI, data = Fertility)
summary(lmfull)
#life expectancy is found to be insignificant
library(ggpubr)
ggqqplot(lmfull$residuals)
shapiro.test(lmfull$residuals)
#p=0.5802>0.05, the normality of residuals is significant
#outliers testing with full model
install.packages('car')
library(car)
outlierTest(lmfull)
#
#DFFITS testing for influential points of the fullmode
dffitsfull<-as.data.frame(dffits(lmfull))
dffitsfull
thresholdfull<-2*sqrt(5/34) #p=5, n=34 for fullmode
plot(dffits(lmfull), type = 'h')
abline(h = thresholdfull, lty = 2)
abline(h = -thresholdfull, lty = 2)
#based on outlier and high-influential test, No. 33 and No. 34 should be excluded from the model

Fertilitynew <- Fertility[-c(33,34),]
lmnew <- lm(fertility~lifexpect+univ.recruit.+urb.popul.+inf.mortality+income.CPI, data = Fertilitynew)
summary(lmnew)
#based on backward selection, life expectancy and infant mortality are abandoned in the following model
lmnew1 <- update(lmnew,.~.-inf.mortality)
lmnew1 <- update(lmnew1,.~.-lifexpect)
summary(lmnew1)
# university recuitment, urban polulation, and income/CPI are significant
shapiro.test(lmnew1$residuals)
#residuals are normally distributed
shapiro.test(Fertilitynew$fertility)
#y variable is not normally distributed

#normality test
library(dplyr)
library(ggpubr)
ggqqplot(Fertility$fertility, main='fertility') # right skewed-long tail
ggqqplot(Fertility$lifexpect, main='lifexpect') # close to normal distribution
ggqqplot(Fertility$univ.recruit., main='univ.recruit.')
ggqqplot(Fertility$urb.popul.,main='urb.popul.')
ggqqplot(Fertility$inf.mortality,main='inf.mortality')
ggqqplot(Fertility$income.CPI,main='income/CPI') #right skewed
shapiro.test(Fertility$fertility)
shapiro.test(Fertility$lifexpect)
shapiro.test(Fertility$univ.recruit.)
shapiro.test(Fertility$urb.popul.)# normal
shapiro.test(Fertility$inf.mortality)
shapiro.test(Fertility$income.CPI)
#Obviously, the raw data cannot be regarded as normally distributed random variables
#(the normality tests of x variables are actually unnecessary)

res.lmfull <- residuals.lm(lmfull)

#data transformation/GLM regression

#BOX COX
library(MASS)
y <-log(Fertility$fertility)
x<-sqrt(Fertility$lifexpect);
shapiro.test(y)




bc <- boxcox(y ~ x);
lambda <- bc$x[which.max(bc$y)];
new_model <- lm(((y^lambda-1)/lambda) ~ x)




# university recruitment, urban population, and income/CPI are significant

#model selection
#training vs testing slpitting ratio = 5:1
ind <- seq(6, nrow(Fertilitynew), by=6) # an indicating vector from 6 to 30 by 6
fertilitytest <- Fertilitynew[c(ind),]
fertilitytrain <- Fertilitynew[-c(ind),]
fertilitytestx <- fertilitytest[,c(-2,-1)]
fertilitytrainx <- fertilitytrain[,c(-2,-1)]

#full model
lmf <- lm(fertility~lifexpect+univ.recruit.+urb.popul.
          +inf.mortality+income.CPI, data = fertilitytrain)
summary(lmf)
fertilitytestx <- data.frame(fertilitytestx) # data frame required
predict(lmf, fertilitytestx)
shapiro.test(lmf$residuals)
ggqqplot(lmf$residuals)

#backward selection
lmback1 <- update(lmf,.~.-inf.mortality) #remove infant mortality
summary(lmback1)
lmback2 <- update(lmback1,.~.-lifexpect) # remove life expectancy
summary(lmback2)
#normality
shapiro.test(lmback1$residuals)
ggqqplot(lmback1$residuals)
shapiro.test(lmback2$residuals)
ggqqplot(lmback2$residuals)
fertilitytestx <- data.frame(fertilitytestx) # data frame required
#multicollinearity
car::vif(lmback2)
# prediction
predict(lmback1, fertilitytestx)
predict(lmback2, fertilitytestx)

#PCR
set.seed (1000) # random
lmpcr <- pcr(fertility~lifexpect+univ.recruit.+urb.popul.
             +inf.mortality+income.CPI, data = fertilitytrain,
             scale = TRUE, validation = "CV") 
summary(lmpcr)
# Plot the root mean squared error
validationplot(lmpcr)
# Plot the R2
validationplot(lmpcr, val.type = "R2")
# RMSEP test showed an minimum adjusted CV at ncomp=4, while the fifth component contributes little to %variance
# ncomp=4, prediction of PCR model
predict(lmpcr,fertilitytestx,ncomp=4)
# best pcr model
lmpcrbest <- pcr(fertility~lifexpect+univ.recruit.+urb.popul.
                 +inf.mortality+income.CPI, data = fertilitytrain,
                 scale = TRUE, validation = "CV", ncomp=4) 
a <- unlist(lmpcrbest$residuals)
class(a) # list to vector transformation
shapiro.test(a)
ggqqplot(as.numeric(a))

#PLSR
set.seed (100)
lmplsr <- plsr(fertility~.lifexpect+univ.recruit.+urb.popul.
               +inf.mortality+income.CPI, data = fertilitytrain, 
               scale = TRUE, validation = "CV")
summary(lmplsr)
validationplot(lmplsr)
validationplot(lmplsr, val.type = "R2")
# Similar to pcr, ncomp=4 should be enough for the prediction
# best plsr model
lmplsrbest <- plsr(fertility~lifexpect+univ.recruit.+urb.popul.
                   +inf.mortality+income.CPI, data = fertilitytrain,
                   scale = TRUE, validation = "CV", ncomp=4) 
# prediction of PLSR model
predict(lmplsrbest, fertilitytestx, ncomp = 4)
b <- unlist(lmplsrbest$residuals)
b <- as.numeric(b)
class(b) # list to vector transformation
shapiro.test(b)
ggqqplot(b)

# Ridge regression
# Setting the range of lambda values
lambda_seq <- 10^seq(2, -2, by = -.1)
#data transformed to a matrix
fertilitytrainx <- data.matrix(fertilitytrainx) 
fertilitytestx <- data.matrix(fertilitytestx) 
# build up ridge regression
lmridge <- glmnet(fertilitytrainx, fertilitytrain$fertility, 
                  alpha = 0, lambda  = lambda_seq) 
# Checking the model
summary(lmridge)
# Using cross validation glmnet
ridge_cv <- cv.glmnet(fertilitytrainx, fertilitytrain$fertility, 
                      alpha = 0, lambda = lambda_seq)
# Best lambda value
best_lambda1 <- ridge_cv$lambda.min
best_lambda1
# the best lambda 0.01
lmridgebest <- glmnet(fertilitytrainx, fertilitytrain$fertility, 
                      alpha = 0, lambda = best_lambda1)
summary(lmridgebest)
coef(lmridgebest)
fertilitytestx <- data.matrix(fertilitytestx)
predict(lmridgebest, fertilitytestx)
# c is the residuals vector
c <- predict(lmridgebest, fertilitytrainx)-fertilitytrain$fertility
c <- as.numeric(c)
class(c) # list to vector transformation
shapiro.test(c)
ggqqplot(c)
# normality is acceptable

#LASSO regression
# build up LASSO regression
lmlasso <- glmnet(fertilitytrainx, fertilitytrain$fertility, alpha = 1, lambda  = lambda_seq)
# Checking the model
summary(lmlasso)
# Using cross validation glmnet
lasso_cv <- cv.glmnet(fertilitytrainx, fertilitytrain$fertility, alpha = 1, lambda  = lambda_seq)
# Best lambda value
best_lambda2 <- lasso_cv$lambda.min
best_lambda2
# the best lambda 0.01
lmlassobest <- glmnet(fertilitytrainx, fertilitytrain$fertility, alpha = 1, lambda  = best_lambda2)
summary(lmlassobest)
coef(lmlassobest)
fertilitytestx <- data.matrix(fertilitytestx)
predict(lmlassobest, fertilitytestx)
# c is the residuals vector
d <- predict(lmlassobest, fertilitytrainx)-fertilitytrain$fertility
d <- as.numeric(d)
class(d) # list to vector transformation
shapiro.test(d)
ggqqplot(d)
# normality is good
car::vif(lmlassobest)

# AIC
aicsub <- regsubsets(fertility~lifexpect+univ.recruit.+urb.popul.
                     +inf.mortality+income.CPI, data = fertilitytrain)
rs <- summary(aicsub)
rs$which
rs$rss
AIC <- 27*log(rs$rss/27) + (2:6)*2 #n=27, p from 2-6
plot(AIC ~ I(1:5), ylab="AIC", xlab="Number of Predictors")
AIC
#AIC completed
# all except for lifexpect should be included in this model
lmAIC <- lm(fertility~univ.recruit.+urb.popul.
            +inf.mortality+income.CPI, data = fertilitytrain)
summary(lmAIC)
#prediction of AIC model
fertilitytestx <- data.frame(fertilitytestx)
predaic <- predict.lm(lmAIC, fertilitytestx)
# collinearity test
car::vif(lmAIC)
# normality test
ggqqplot(lmAIC$residuals)
shapiro.test(lmAIC$residuals)

#lasso based on the full data
Fertilitynewx <- data.matrix(Fertilitynew[,c(-2,-1)])
class(Fertilitynewx)
lmlasso1 <- glmnet(Fertilitynewx, Fertilitynew$fertility, alpha = 1, lambda  = lambda_seq)
# Checking the model
summary(lmlasso1)
# Using cross validation glmnet
lasso1_cv <- cv.glmnet(Fertilitynewx, Fertilitynew$fertility, alpha = 1, lambda  = lambda_seq)
# Best lambda value
best_lambda3 <- lasso1_cv$lambda.min
best_lambda3
# the best lambda 0.01
lmlasso1best <- glmnet(Fertilitynewx, Fertilitynew$fertility, alpha = 1, lambda  = best_lambda3)
summary(lmlasso1best)
coef(lmlasso1best)
e <- predict(lmlasso1best, Fertilitynewx)-Fertilitynew$fertility
e <- as.numeric(e)
class(e) # list to vector transformation
shapiro.test(e)
ggqqplot(e)
# normality is acceptable

#ridge based on the full data
lmridge1 <- glmnet(Fertilitynewx, Fertilitynew$fertility, alpha = 0, lambda  = lambda_seq)
# Checking the model
summary(lmridge1)
# Using cross validation glmnet
ridge1_cv <- cv.glmnet(Fertilitynewx, Fertilitynew$fertility, alpha = 0, lambda  = lambda_seq)
# Best lambda value
best_lambda4 <- ridge1_cv$lambda.min
best_lambda4
# the best lambda 0.01
lmridge1best <- glmnet(Fertilitynewx, Fertilitynew$fertility, alpha = 0, lambda  = best_lambda4)
summary(lmridge1best)
coef(lmridge1best)
f <- predict(lmridge1best, Fertilitynewx)-Fertilitynew$fertility
f <- as.numeric(f)
class(f) # list to vector transformation
shapiro.test(f)
ggqqplot(f)
# normality is acceptable

# lmf, lmback1, lmback2 based on full data
lmf1 <- lm(fertility~lifexpect+univ.recruit.+urb.popul.+inf.mortality+income.CPI, data = Fertilitynew)
summary(lmf1)
shapiro.test(lmf1$residuals)
ggqqplot(lmf1$residuals)
#based on backward selection, life expectancy and infant mortality are abandoned in the following model
lmb1 <- update(lmf1,.~.-inf.mortality)
summary(lmb1)
shapiro.test(lmb1$residuals)
ggqqplot(lmb1$residuals)
lmb2 <- update(lmb1,.~.-lifexpect)
summary(lmb2)
shapiro.test(lmb2$residuals)
ggqqplot(lmb2$residuals)
#residuals are normally distributed
Anova(lmf1,lmb2)

car::vif(lmf1)
car::vif(lmb1)
car::vif(lmb2)
ols_coll_diag(lmf1)
ols_coll_diag(lmb1)
ols_coll_diag(lmb2)
# strong multicolinearity was found
#PCR based on lmback1
set.seed (1000)
lmpcrnew <- pcr(fertility~univ.recruit.+urb.popul.+income.CPI+lifexpect, data = Fertilitynew, scale = TRUE, validation = "CV")
summary(lmpcrnew)
# Plot the root mean squared error
validationplot(lmpcrnew)
# Plot the R2
validationplot(lmpcrnew, val.type = "R2")
# RMSEP test showed an minimum adjusted CV at ncomp=4, while the fifth component contributes little to %variance
# ncomp=4, prediction of PCR model
predict(lmpcrnew,fertilitytestx,ncomp=4)
anew <- unlist(lmpcrnew$residuals)
class(anew) # list to vector transformation
shapiro.test(anew)
ggqqplot(as.numeric(anew))

#PCR based on lmback2
set.seed (1000)
lmpcrnew2 <- pcr(fertility~univ.recruit.+urb.popul.+income.CPI, data = Fertilitynew, scale = TRUE, validation = "CV")
summary(lmpcrnew2)
# Plot the root mean squared error
validationplot(lmpcrnew2)
# Plot the R2
validationplot(lmpcrnew2, val.type = "R2")
# RMSEP test showed an minimum adjusted CV at ncomp=4, while the fifth component contributes little to %variance
# ncomp=4, prediction of PCR model
predict(lmpcrnew2,fertilitytestx,ncomp=4)
anew2 <- unlist(lmpcrnew2$residuals)
class(anew2) # list to vector transformation
shapiro.test(anew2)
ggqqplot(as.numeric(anew2))










