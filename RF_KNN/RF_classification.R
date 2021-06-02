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

training_matrix = training_matrix[,filtered_features]
test_matrix = test_matrix[,filtered_features]

training_matrix$weights = NULL
test_matrix$weights = NULL

training_matrix[, drug_name] = sapply(training_matrix[, drug_name], FUN = function(x){ if(x == 0){return('resistant')}else{return('sensitive')}})
test_matrix[, drug_name] = sapply(test_matrix[, drug_name], FUN = function(x){ if(x == 0){return('resistant')}else{return('sensitive')}})

training_matrix[, drug_name] = factor(training_matrix[, drug_name], levels =c ("sensitive", "resistant"))
test_matrix[, drug_name] = factor(test_matrix[, drug_name], levels = levels(training_matrix[,drug_name]))

#Need same factors in training and test set for prediction
same_factor_generation_matrix = rbind(training_matrix, test_matrix)

same_factor_generation_matrix[colnames(same_factor_generation_matrix)] = lapply(same_factor_generation_matrix[colnames(same_factor_generation_matrix)], factor)

training_matrix = same_factor_generation_matrix[rownames(training_matrix),]
test_matrix = same_factor_generation_matrix[rownames(test_matrix),]


tr_response = training_matrix[, drug_name]
training_mtx_without_response = training_matrix[, !(colnames(training_matrix) %in% drug_name)]


#caret package
trctrl = trainControl(method = "cv", number = 5, search = "grid",  summaryFunction = twoClassSummary, classProbs = TRUE)
mtry = as.integer(sqrt(ncol(training_mtx_without_response)))


grid = seq(from = max(mtry-50, 1), to = mtry + 50, by = 2)

set.seed(293746)
#fitting
rf_fit = train(x = training_mtx_without_response, y = tr_response, method = "rf", trControl = trctrl,  metric = "ROC", tuneGrid = expand.grid(.mtry = c(grid)))


final_m = rf_fit$bestTune$mtry
cv_result_df = rf_fit$results
cv_sens_final_m = cv_result_df[cv_result_df$mtry == final_m, "Sens"]
cv_spec_final_m = cv_result_df[cv_result_df$mtry == final_m, "Spec"]

#prediction for test set
rf_predict = predict(rf_fit, newdata = test_matrix)
confusion_mtx = confusionMatrix(rf_predict, test_matrix[, drug_name])

test_sens = as.numeric(confusion_mtx$byClass["Sensitivity"])
test_spec = as.numeric(confusion_mtx$byClass["Specificity"])


output_vector = c("mtry" = final_m, "CV_sensitivity" = cv_sens_final_m, "CV_specificity" = cv_spec_final_m, "Test_sensitivity" = test_sens, "Test_specificity" = test_spec)

write.table(x = output_vector, file = output_file, quote = FALSE, row.names = TRUE, col.names = FALSE, sep = "\t")
