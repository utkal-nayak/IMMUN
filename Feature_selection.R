####Load library
library(randomForest)

####Load full data including accession numbers
full_data <- read.csv("123&97 EZ_descriptor_ACC_features.csv")
accessions <- full_data[, 1]  # Save accession numbers
features <- full_data[, -1]  # Remove accession column

####Add class labels
features$Class <- as.factor(c(rep("Antigen", 123), rep("Non-Antigen", 97)))

####Separate features and labels
labels <- features$Class
feature_matrix <- features[, -ncol(features)]

###Run Random Forest 50 times to collect importance scores
set.seed(123)
importance_list <- replicate(50, {
  idx <- sample(1:nrow(features), 0.8 * nrow(features))
  rf <- randomForest(x = feature_matrix[idx, ], y = labels[idx], importance = TRUE, ntree = 500)
  importance(rf, type = 2)[, "MeanDecreaseGini"]
}, simplify = FALSE)

####Average importance and select top features
importance_mean <- rowMeans(do.call(cbind, importance_list))
sorted_features <- sort(importance_mean, decreasing = TRUE)
top_n <- 50
selected_features <- names(sorted_features)[1:min(top_n, length(sorted_features))]

####Subset top features and combine with accession + class
selected_data <- cbind(Accession = accessions, feature_matrix[, selected_features], Class = labels)

###Save to CSV
write.csv(selected_data, "123&97Top_ACC_Features_RF_Selected_with_Accession.csv", row.names = FALSE)

###Create full ranked feature table (optional)
ranked_features_df <- data.frame(
  Feature = names(sorted_features),
  MeanImportance = sorted_features,
  Rank = rank(-sorted_features, ties.method = "first")
)
