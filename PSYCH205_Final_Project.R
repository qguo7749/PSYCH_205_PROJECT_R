library(nlme)
library(lattice)
library(car)

data<-read.csv("ptcp_gnometrans.csv")


data_no_g11<-data[data$generation <11, ]
data_no_g11<-data_no_g11[data_no_g11$condition =="SSL", ]
drops <- c("network_id","replication","cloned","algorithm","algorithm_description")
data_no_g11<-data_no_g11[ , !(names(data_no_g11) %in% drops)]
data_no_g11<-na.omit(data_no_g11)

n <- nrow(data_no_g11)

summary(mod.lm0<-lm(p_trans ~ 1,data=data_no_g11))
summary(mod.lm1<-lm(p_trans ~ generation,data=data_no_g11))
summary(mod.lm2<-lm(p_trans ~ s_demo+generation,data=data_no_g11))
anova(mod.lm0,mod.lm1)
anova(mod.lm1,mod.lm2)
anova(mod.lm0,mod.lm2)
anova.1_2<-anova(mod.lm1,mod.lm2)


hist(data_no_g11$p_trans,main = "Histogram of Probability of Successful Transmission",xlab = 'Probability of Successful Transmission')
scatterplot(p_trans ~ generation, data=data_no_g11,xlab='generation',ylab = 'Probability of Successful Transmission')
title(main = "transmission probability vs generation ")
scatterplot(p_trans ~ s_demo, data=data_no_g11,xlab='demonstration quality',ylab = 'Probability of Successful Transmission')
title(main = "transmission probability vs demonstration score ")

(R2.1_2 <- 1 - anova.1_2$RSS[2]/anova.1_2$RSS[1])
(R2.adj.1_2 <- 1 - (anova.1_2$RSS[2]/mod.lm2$df.residual)/(anova.1_2$RSS[1]/mod.lm1$df.residual))


n.folds <- 10  
folds <- cut(seq(1,n),breaks=n.folds,labels=FALSE)
folds <- sample(folds, replace = FALSE)


MSE.0 <- array(data=0, dim = n.folds)
MSE.1 <- array(data=0, dim = n.folds)
MSE.2 <- array(data=0, dim = n.folds)


#Cross Validation
for(i in 1:n.folds){
  testIndexes <- which(folds==i,arr.ind=TRUE)
  testData <- data_no_g11[testIndexes, ]
  trainData <- data_no_g11[-testIndexes, ]
  
  model.0.cv <- lm(p_trans ~ 1, data = trainData)
  model.1.cv <- lm(p_trans ~ generation, data = trainData)
  model.2.cv <- lm(p_trans ~ s_demo+generation, data = trainData)
  
  pred.0 <- predict(model.0.cv, newdata = testData)
  pred.1 <- predict(model.1.cv, newdata = testData)
  pred.2 <- predict(model.2.cv, newdata = testData)
  
  MSE.0[i] <- sum((testData$p_trans - pred.0)^2)/nrow(testData)
  MSE.1[i] <- mean((testData$p_trans - pred.1)^2)
  MSE.2[i] <- mean((testData$p_trans - pred.2)^2)
  
}

# Now we can calculate all of the cross-validated R2 - notice that these are arrays that have the size of you cv folds.  You get one value of R2 per fold.  This will be useful to get a confidence interval on your R2CV

R2.cv.0_1 <- 1 - MSE.1/MSE.0
R2.cv.1_2 <- 1 - MSE.2/MSE.1
R2.cv.0_2 <- 1 - MSE.2/MSE.0

# Get means and SEs
R2.cv.m.0_1 <- mean(R2.cv.0_1)
R2.cv.m.1_2 <- mean(R2.cv.1_2)
R2.cv.m.0_2 <- mean(R2.cv.0_2)

R2.cv.se.0_1 <- sqrt(sum((R2.cv.0_1 - R2.cv.m.0_1)^2)/(n.folds-1))*sqrt(1/n.folds + 1/(n.folds-1))
R2.cv.se.1_2 <- sqrt(sum((R2.cv.1_2 - R2.cv.m.1_2)^2)/(n.folds-1))*sqrt(1/n.folds + 1/(n.folds-1))
R2.cv.se.0_2 <- sqrt(sum((R2.cv.0_2 - R2.cv.m.0_2)^2)/(n.folds-1))*sqrt(1/n.folds + 1/(n.folds-1))

hist(R2.cv.0_1, xlab="R2", ylab="Folds",
     main=paste("R2 scores across folds (mean = ",
                round(R2.cv.m.0_1, 3), ", se = ", round(R2.cv.se.0_1, 3), ")"))
abline(v=R2.cv.m.0_1, col="red")
abline(v=R2.cv.m.0_1 - 1.96*R2.cv.se.0_1, col="blue")
abline(v=R2.cv.m.0_1 + 1.96*R2.cv.se.0_1, col="blue")

hist(R2.cv.1_2, xlab="R2", ylab="Folds",
     main=paste("R2 scores across folds (mean = ",
                round(R2.cv.m.1_2, 3), ", se = ", round(R2.cv.se.1_2, 3), ")"))
abline(v=R2.cv.m.1_2, col="red")
abline(v=R2.cv.m.1_2 - 1.96*R2.cv.se.1_2, col="blue")
abline(v=R2.cv.m.1_2 + 1.96*R2.cv.se.1_2, col="blue")

hist(R2.cv.0_2, xlab="R2", ylab="Folds",
     main=paste("R2 scores across folds (mean = ",
                round(R2.cv.m.0_2, 3), ", se = ", round(R2.cv.se.0_2, 3), ")"))
abline(v=R2.cv.m.0_2, col="red")
abline(v=R2.cv.m.0_2 - 1.96*R2.cv.se.0_2, col="blue")
abline(v=R2.cv.m.0_2 + 1.96*R2.cv.se.0_2, col="blue")

# Print some results
print(sprintf('Model 1 R2=%.2f R2adj=%.2f R2cv=%.2f +- %.3f', summary(mod.lm1)$r.squared, summary(mod.lm1)$adj.r.squared, R2.cv.m.0_1, 1.96*R2.cv.se.0_1))
print(sprintf('Model 2 R2=%.2f R2adj=%.2f R2cv=%.2f +- %.3f', summary(mod.lm2)$r.squared, summary(mod.lm2)$adj.r.squared, R2.cv.m.0_2, 1.96*R2.cv.se.0_2))
print(sprintf('Model 1vs2 R2=%.2f R2adj=%.2f R2cv=%.2f +- %.3f', R2.1_2, R2.adj.1_2, R2.cv.m.1_2, 1.96*R2.cv.se.1_2))


mod.lme<-lme(p_trans ~ s_demo+generation, random = ~ 1 |participant_id , data=data_no_g11,method = 'ML')
summary(mod.lme)

compareCoefs(mod.lm2, mod.lme)

