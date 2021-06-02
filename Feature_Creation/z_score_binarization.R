
args = commandArgs(trailingOnly=TRUE)
if(length(args) !=2){
  stop("Z score matrix needed and output file for binarized matrix")
}

matrixfile = args[1]
output_file = args[2]
print(matrixfile)
gene_expression_matrix = read.table(file = matrixfile, header = TRUE, sep ="\t", stringsAsFactors = FALSE, check.names = FALSE)

matrix = as.data.frame(t(as.matrix(gene_expression_matrix)))

#Make high/low expressed categories
threshold_down = -1.65
threshold_up = 1.65

down_matrix = data.frame(matrix(0, nrow = nrow(matrix), ncol = ncol(matrix)))
up_matrix = data.frame(matrix(0, nrow = nrow(matrix), ncol = ncol(matrix)))
rownames(down_matrix) = rownames(matrix)
colnames(down_matrix) = colnames(matrix)
rownames(up_matrix) = rownames(matrix)
colnames(up_matrix) = colnames(matrix)

for(i in 1:nrow(matrix)){
  
  for(j in 1:ncol(matrix)){
    
    if(matrix[i,j] < threshold_down){
      
      down_matrix[i,j] = 1
      
    }
    
    if(matrix[i,j]> threshold_up){
      up_matrix[i,j] = 1
    }
  }
}

down_to_be_kept = apply(down_matrix, 2, function(x){sum(x) != 0})
up_to_be_kept = apply(up_matrix, 2, function(x){sum(x) != 0})

down_matrix_final = down_matrix[,down_to_be_kept]
up_matrix_final = up_matrix[, up_to_be_kept]


colnames(up_matrix_final) = paste(colnames(up_matrix_final), "_high_expr", sep = "")
colnames(down_matrix_final) = paste(colnames(down_matrix_final), "_low_expr", sep = "")


final_matrix = cbind(up_matrix_final, down_matrix_final)

write.table(final_matrix, file = output_file, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE )


