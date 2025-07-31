###Set working directory and load data
setwd("C:/Users/ASUS/Desktop/Snigi")
data <- read.csv("64 & 45 Top_ACC_Features_RF_Selected_with_Accession.csv")

###Extract class labels and feature matrix
V_p <- as.factor(data$Class)
features <- data[, !names(data) %in% "Class"]

###Train-test split (80-20)
set.seed(123)
ids <- sample(1:nrow(features), 0.8 * nrow(features))
train.x <- features[ids, ]
test.x  <- features[-ids, ]
train.y <- V_p[ids]
test.y  <- V_p[-ids]

###Load required package
library(randomForest)

###Train Random Forest model
model_rf <- randomForest(x = train.x, y = train.y,
                         importance = TRUE, proximity = TRUE,
                         ntree = 500, mtry = 4)

###Predict on test data
pred <- predict(model_rf, test.x)

###Confusion matrix
conf.mat <- table(Actual = test.y, Predicted = pred)
print(conf.mat)

###Accuracy
accuracy <- sum(diag(conf.mat)) / sum(conf.mat)
cat("Accuracy:", round(accuracy, 4), "\n")

###Error rate
error_rate <- 1 - accuracy
cat("Error Rate:", round(error_rate, 4), "\n")

###Plot model error rate across trees
plot(model_rf)

