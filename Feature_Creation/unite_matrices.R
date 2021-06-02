library(MASS)

args = commandArgs(trailingOnly=TRUE)
if(length(args) !=6){
  stop("Binarized gene expression matrix, pre-processed mutation and cnv matrix, drug name, gdsc version, output file for matrix and output file for features needed")
}

gene_expression_file = args[1]
mut_cnv_file = args[2]
drug = args[3]
gdsc_version = args[4]
output_file = args[5]
feature_output_file = args[6]

print(drug)
gene_expr_mat = read.table(file = gene_expression_file, header = TRUE, sep ="\t", stringsAsFactors = FALSE, check.names = FALSE)

feat_mat = read.table(mut_cnv_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

#Row and column names
colnames(feat_mat) = feat_mat[1,]
feat_mat = feat_mat[-1,]

rownames(feat_mat) = feat_mat[,1]
feat_mat = feat_mat[,-1]


#Unite all sensitive features
sensitive_matrix = feat_mat[-c(1,2), feat_mat["implication",] == "sensitive", drop = FALSE]
sensitive_features = colnames(sensitive_matrix)

#Unite all resistant features
resistant_matrix = feat_mat[-c(1,2), feat_mat["implication",] == "resistant", drop = FALSE]
resistant_features = colnames(resistant_matrix)


#Unite all "not sensitive" features
not_sensitive_matrix = feat_mat[-c(1,2), feat_mat["implication",] == "not sensitive", drop = FALSE]
not_sensitive_features = colnames(not_sensitive_matrix)


#Unite all "not resistant" features
not_resistant_matrix = feat_mat[-c(1,2), feat_mat["implication",] == "not resistant", drop= FALSE]
not_resistant_features = colnames(not_resistant_matrix)


feat_mat = feat_mat[-c(1,2),]


##########################################################################################
#Make a special matrix for the case that you are not interested in single already selected features
unite_columns = function(x_matrix, feat_mat, col_name){
  known_sensitivity_determinants = as.data.frame(apply(x_matrix, 1, function(x) {
    if(sum(as.numeric(x))>0){1} else{0}
  } ),drop = FALSE)
  
  colnames(known_sensitivity_determinants) = col_name
  remaining_features = setdiff(colnames(feat_mat),colnames(x_matrix))
  special_mat = feat_mat[, remaining_features]
  special_mat = cbind(special_mat, known_sensitivity_determinants)
  
  return(special_mat)
}

if(ncol(sensitive_matrix) != 0){
  feat_mat = unite_columns(sensitive_matrix, feat_mat, "known_sensitivity_determinants")
}

if(ncol(resistant_matrix) != 0){
  feat_mat = unite_columns(resistant_matrix, feat_mat, "known_resistance_determinants")
}

if(ncol(not_sensitive_matrix) != 0){
  feat_mat = unite_columns(not_sensitive_matrix, feat_mat, "known_not_sensitivity_determinants")
}

if(ncol(not_resistant_matrix) != 0){
  feat_mat = unite_columns(not_resistant_matrix, feat_mat, "known_not_resistance_determinants")
}


#Final matrix
final_matrix = apply(feat_mat, 2, function(x){as.numeric(x)})
rownames(final_matrix) = rownames(feat_mat)

final_matrix = cbind(final_matrix, gene_expr_mat[rownames(final_matrix),])


drug_response_discr = ""

if(gdsc_version == "GDSC2"){
  drug_response_discr="<path_to_discretized_drug_response_matrix_GDSC2>"
}else{
  drug_response_discr="<path_to_discretized_drug_response_matrix_GDSC1>"
}

drug_response_matrix = read.csv(drug_response_discr, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

drug_response_numbers = drug_response_matrix[rownames(final_matrix), drug, drop = FALSE]


final_matrix = as.data.frame(final_matrix)

final_matrix$weights = rep(1, nrow(final_matrix))

final_matrix = cbind(final_matrix, drug_response_numbers)

write.table(x = final_matrix, file = output_file, sep = "\t", quote = FALSE)



number_of_pre_selected_features = length(sensitive_features) + length(resistant_features) + length(not_sensitive_features) + length(not_resistant_features)

if(number_of_pre_selected_features != 0){
  pre_selected_features_matrix = data.frame(matrix(0, nrow = number_of_pre_selected_features, ncol = 4))
  colnames(pre_selected_features_matrix) = c("sensitive_selected",	"sensitive_not_selected",	"resistant_selected",	"resistant_not_selected")
  
  rownames_pre_selected_features_matrix = c(sensitive_features, resistant_features, not_sensitive_features, not_resistant_features)
  rownames(pre_selected_features_matrix) = rownames_pre_selected_features_matrix
  
  
  for(feature in sensitive_features){
    
    pre_selected_features_matrix[feature,1] = 1
  }
  for(feature in resistant_features){
    pre_selected_features_matrix[feature,3] = 1
  }
  for(feature in not_sensitive_features){
    pre_selected_features_matrix[feature,2] = 1
  }
  for(feature in not_resistant_features){
    pre_selected_features_matrix[feature, 4] = 1
  }
  
  write.table(x = pre_selected_features_matrix, file =feature_output_file, sep = "\t", quote = FALSE)
}else{
  print(drug)
  print("No features selected")
}


