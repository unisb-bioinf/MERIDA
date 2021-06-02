args = commandArgs(trailingOnly = TRUE)

if(length(args) !=6){
  stop("6 arguments needed: available cell lines file, drug_name, drug_file, training samples output file, test samples output file, percentage of test samples")
}

available_cl_file = args[1]
drug_name= args[2]
drug_file= args[3]
training_output = args[4]
test_output = args[5]
percentage_test = as.double(args[6])

available_cl = read.table(file = available_cl_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

#drug
drug_data = read.table(file = drug_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

if(!(drug_name %in% colnames(drug_data))){
  print("Drug not contained in drug response file")
  print(drug_name)
}

drug_data_specific_drug = drug_data[as.character(available_cl$V1), drug_name, drop = FALSE]


#cell line annotations
cell_line_annotation_file = '<path_to_Cell_Lines_Details.csv>'
cell_line_annotation_file2 = '<path_to_COSMIC_Tissue_Classification.csv>'


cell_line_annotations2 = read.csv(file = cell_line_annotation_file2, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(cell_line_annotations2) = cell_line_annotations2[1,]
cell_line_annotations2 = cell_line_annotations2[-1,]


needed_cell_lines_annotations = cell_line_annotations2[cell_line_annotations2$COSMIC_ID %in% as.character(available_cl$V1),]


library(splitstackshape)

set.seed(37360025)#usual seed

train_test = stratified(needed_cell_lines_annotations, group = c("Site"), size = percentage_test, bothSets = TRUE)

test_cls = as.character(train_test[[1]]$COSMIC_ID)
train_cls = as.character(train_test[[2]]$COSMIC_ID)

if(sum(drug_data_specific_drug[test_cls,1]) < 0.1 * length(test_cls)){#average number of sensitive cell lines corresponds to 10%
  print(drug_name)
  print("Fewer sensitive cell lines than to be expected, please check matrix")
  number_of_sensitive_cls = as.character(as.double(sum(drug_data_specific_drug[test_cls,1])))
  print("Number of sensitive cell lines:")
  print(number_of_sensitive_cls)
  number_of_expceted_sensitive_cls = as.character(as.double(0.1*length(test_cls)))
  print("Number of expected sensitive cell lines:")
  print(number_of_expceted_sensitive_cls)
}

write.table(x = train_cls, file = training_output, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(x = test_cls, file = test_output, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)



