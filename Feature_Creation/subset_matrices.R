args = commandArgs(trailingOnly = TRUE)

if(length(args) !=5){
  stop("5 arguments needed: full matrix, training set, test set, training matrix output file, test matrix output file")
}

full_matrix_file = args[1]
training_set_file= args[2]
test_set_file= args[3]
training_output = args[4]
test_output = args[5]

available_cl_tr = read.table(file = training_set_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
available_cl_te = read.table(file = test_set_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

full_matrix = read.table(file = full_matrix_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

training_mtx = full_matrix[as.character(available_cl_tr$V1),]
test_mtx = full_matrix[as.character(available_cl_te$V1),]


write.table(x = training_mtx, file = training_output, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(x = test_mtx, file = test_output, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)



