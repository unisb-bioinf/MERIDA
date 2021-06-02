library(caret)

args = commandArgs(trailingOnly=TRUE)

training_matrix_file= args[1]
test_matrix_file = args[2]
drug_name = args[3]
output_file = args[4]

training_matrix = read.table(file = training_matrix_file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
test_matrix = read.table(file = test_matrix_file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

filtered_features = which(colSums(training_matrix) > 0)
if(!(ncol(training_matrix) %in% filtered_features)){
  filtered_features = c(filtered_features, ncol(training_matrix))
}

training_matrix$weights = NULL
test_matrix$weights = NULL

training_matrix[, drug_name] = sapply(training_matrix[, drug_name], FUN = function(x){ if(x == 0){return('resistant')}else{return('sensitive')}})
test_matrix[, drug_name] = sapply(test_matrix[, drug_name], FUN = function(x){ if(x == 0){return('resistant')}else{return('sensitive')}})

training_matrix[, drug_name] = factor(training_matrix[, drug_name], levels =c ("sensitive", "resistant"))
test_matrix[, drug_name] = factor(test_matrix[, drug_name], levels = c("sensitive", "resistant"))


tr_response = training_matrix[, drug_name]
training_mtx_without_response = training_matrix[, !(colnames(training_matrix) %in% drug_name)]
k_to_be_tested = seq(from = 1, to = 50, by = 2)


#caret package
trctrl = trainControl(method = "cv", number = 5, summaryFunction = twoClassSummary, classProbs = TRUE)
set.seed(293746)
#fitting
knn_fit = train(x = training_mtx_without_response, y = tr_response, method = "knn", preProcess = c(), trControl = trctrl, tuneGrid = expand.grid(k = c(k_to_be_tested)),  metric = "ROC")
final_k = knn_fit$finalModel$k
cv_result_df = knn_fit$results
cv_sens_final_k = cv_result_df[cv_result_df$k == final_k, "Sens"]
cv_spec_final_k = cv_result_df[cv_result_df$k == final_k, "Spec"]

#prediction for test set
knn_predict = predict(knn_fit, newdata = test_matrix)
confusion_mtx = confusionMatrix(knn_predict, test_matrix[, drug_name])

test_sens = as.numeric(confusion_mtx$byClass["Sensitivity"])
test_spec = as.numeric(confusion_mtx$byClass["Specificity"])


output_vector = c("k" = final_k, "CV_sensitivity" = cv_sens_final_k, "CV_specificity" = cv_spec_final_k, "Test_sensitivity" = test_sens, "Test_specificity" = test_spec)

write.table(x = output_vector, file = output_file, quote = FALSE, row.names = TRUE, col.names = FALSE, sep = "\t")
