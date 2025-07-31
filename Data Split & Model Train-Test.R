setwd("C:/Users/ASUS/Desktop/Snigi")
data<-read.csv(file="64 & 45 Top_ACC_Features_RF_Selected_with_Accession.csv")
V_p <- as.factor(c(rep("Antigen", 64), rep("Non-Antigen", 45)))
length(V_p)                   
length(V_p) 
set.seed(123)
####################
###Splitting data into training and testing sets
ids<-sample(1:nrow(data), round(0.8*nrow(data)),replace=FALSE)
train.x <- data[ids, ]
test.x <- data[-ids, ]
v_p.train<-V_p[ids]
v_p.test<- V_p[-ids]
###Rf model
library(randomForest)
model_rf<-randomForest(x=train.x, y=v_p.train,importance=TRUE,proximity=TRUE,ntreeTry=500,mtry=4,cv.fold=10)
####Predict on the test set
pred.test<-predict(model_rf, test.x)
table(pred.test)


####Confusion matrix
conf.rf <- table(v_p.test, pred.test)
summary(conf.rf)
plot(model_rf)
#####Calculate accuracy
accuracy <- sum(diag(conf.rf)/sum(conf.rf))
mean<-mean(accuracy)
err1 <- 1 -mean




