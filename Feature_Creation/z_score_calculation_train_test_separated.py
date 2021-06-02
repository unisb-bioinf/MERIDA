#Author: Kerstin Lenhof
#Date: 12.03.2021
#Z-score calculation 
import sys
import numpy as np



def read_cell_line_set_file(set_file):
	
	cl_list = []
	with open(set_file, "r", encoding = "utf-8") as a_cls:
		
		for line in a_cls:
			
			sline = line.strip()
			
			if sline == "":
				continue
			else:
				cl_list.append(sline)
	return cl_list

def subset_matrix(input_matrix_file, train_set, test_set):
	
	#print(train_set)
	#print(test_set)
	
	train_set_indices = [0]#0 is the gene name
	test_set_indices = [0]
	
	with open(input_matrix_file, "r", encoding = "utf-8") as input_matrix:
		
		first_line = input_matrix.readline().strip('\n').split('\t')
		
		#print(first_line)
		for cell_line in train_set:
			if cell_line in test_set: 
				print("Cell line in train and test set: " + cell_line)
			else:
				
				train_set_indices.append(first_line.index(cell_line) +1)#first element is gene name
		for cell_line in test_set:
			if cell_line in train_set: 
				print("Cell line in train and test set: " + cell_line)
			else:
				
				test_set_indices.append(first_line.index(cell_line) +1)#first element is gene name
				
	
		training_matrix = [train_set]
		test_matrix = [test_set]
		
		for line in input_matrix:
			
			#print(line)
			sline = line.strip('\n').split('\t')
			
			train_sline = [sline[x] for x in train_set_indices]
			test_sline = [sline[x] for x in test_set_indices]
			
			training_matrix.append(train_sline)
			test_matrix.append(test_sline)
	
	#print(test_matrix)
	return([training_matrix, test_matrix])
			
			

def z_score_matrix_row(training_matrix, test_matrix, output_file):
	
	new_matrix = []
	
	
	number_of_genes = len(training_matrix)
	
	header = training_matrix[0]
	header.extend(test_matrix[0])
	
	for iterator in range(1,number_of_genes,1):
		
		line = training_matrix[iterator]
		
		gene_name = line[0]
		
		z_score_fields = [float(x) for x in line[1:]]
		
		#print(gene_name)
		if int(sum(z_score_fields)) == 0:
			print("zero row")
		else:
			z_score_np = np.array(z_score_fields)
			
			mean_res = np.mean(z_score_np)
			sample_sd = np.std(z_score_np, ddof=1)
			z_score_res = (z_score_np - mean_res)/sample_sd
			
			str_zscores = [str(x) for x in z_score_res]
			str_zscores.insert(0, gene_name)
			
			
			#Make transformation for test matrix with mean and sd from train matrix
			
			line_test = test_matrix[iterator]
			
			z_score_fields_test = [float(x) for x in line_test[1:]]
			z_score_np_test = np.array(z_score_fields_test)
			
			z_score_res_test = (z_score_np_test - mean_res)/sample_sd
			str_zscores_test = [str(x) for x in z_score_res_test]
			
			str_zscores.extend(str_zscores_test)
			
			new_matrix.append('\t'.join(str_zscores))
			
			
			
	
	with open(output_file, 'w', encoding='utf-8') as output_matrix:
		
		output_matrix.write("\t".join(header) + "\n")
		for entry in new_matrix:
			
			output_matrix.write(entry + "\n")
					

#def z_score_matrix_column():#TODO:implement
	

################Main function ##################################

def main(input_matrix_file, train_cl_file, test_cl_file, output_file, mode):
	
	if mode == "row":
		train_cls = read_cell_line_set_file(train_cl_file)
		test_cls = read_cell_line_set_file(test_cl_file)
		
		matrix_subsets = subset_matrix(input_matrix_file, train_cls, test_cls)
		z_score_matrix_row(matrix_subsets[0], matrix_subsets[1], output_file)
		return
	
	elif mode == "column":
		print("Mode not yet supported")
		return
	else:
		print("Unknown mode, supported modes are: 'row' and 'column'")
		return


#################################################################
if __name__ == "__main__":
	print(sys.argv)
	if len(sys.argv) < 6:
			sys.exit("This program needs the following arguments:\
					\n- input matrix file name (including path)\
					\n- training cell line file\
					\n- test cell line file\
					\n- output file name (including path)\
					\n- row or column (for row or column z-score calculation)") 

	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]) 	
	
