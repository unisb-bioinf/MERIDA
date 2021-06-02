#Author:Kerstin Lenhof
#Date: 17.01.2020
#Preparation script for GDSC data

import sys


copy_number_data = "<path_to_GDSC_copy_number_file>"
mutation_data ="<path_to_GDSC_mutation_file>"
gene_expr_zscore_data="<path_to_GDSC_gene_expression_file"
drug_response_gdsc1 = "<path_to_directory_for_GDSC1_drug_responses>"
drug_response_gdsc2 = "<path_to_directory_for_GDSC2_drug_responses>"


def available_cl_ge():

	ge_cl = []
	
	with open(gene_expr_zscore_data, 'r', encoding='utf-8') as ge_data:
		
		first_line = ge_data.readline()
		
		for line in ge_data:
            
			sline = line.strip('\n').split('\t')
            
			cl_name = sline[0].strip()
            
			if cl_name == "":
				continue
			else:
				if not cl_name in ge_cl:
					ge_cl.append(cl_name)
	
	return ge_cl
def available_cl_cnv():

	cnv_cl = []
	
	with open(copy_number_data, 'r', encoding='utf-8') as cnv_data:
		
		first_line = cnv_data.readline()
		
		for line in cnv_data:
			
			sline = line.strip('\n').split('\t')
			
			cosmic_id = sline[0].strip()
			
			if cosmic_id == "":
				continue
			else:
				if not cosmic_id in cnv_cl:
					cnv_cl.append(cosmic_id)
				
	return cnv_cl

def available_cl_mut():
	
	mut_cl = []
	
	with open(mutation_data, 'r', encoding='utf-8') as mut_data:
		
		first_line = mut_data.readline()
		
		for line in mut_data:
			sline = line.strip('\n').split('\t')
			
			cosmic_id =sline[1].strip()
			
			if cosmic_id == "":
				continue
			else:
				if not cosmic_id in mut_cl:
					mut_cl.append(cosmic_id)
			
	return mut_cl

def available_cl_drug(drug_name, gdsc1_or_2): #does not support synonyms yet
	
	drug_cl = []
	
	
	if(gdsc1_or_2 == "GDSC1"):
		filename = drug_response_gdsc1 + drug_name + ".txt" #it is necessary to use the filename of the drug in the directory as drug name
		
		with open(filename, 'r', encoding ='utf-8') as drug_data:
			
			first_line = drug_data.readline()
			
			for line in drug_data:
				
				sline = line.strip('\n').split('\t')
				
				if sline[0].strip() == "":
					continue
				else:
					drug_cl.append(sline[0].strip())
				
				
	else:
		
		filename = drug_response_gdsc2 + drug_name + ".txt"
		
		with open(filename, 'r', encoding = 'utf-8') as drug_data:
			
			first_line = drug_data.readline()
			
			for line in drug_data:
				
				sline = line.strip('\n').split('\t')
				
				
				if sline[0].strip() == "":
					continue
				else:
					drug_cl.append(sline[0].strip())
	
	return [drug_cl]		
			
			


def available_cell_lines_old_drug_data_cnv_mut_ge(drug_name, output_file, gdsc1_or_2):
	
	available_celllines = []
	available_celllines_cnv = available_cl_cnv()
	available_celllines_mutation = available_cl_mut()
	available_celllines_dr_response = available_cl_drug(drug_name, gdsc1_or_2)
	available_celllines_ge = available_cl_ge()
	
	#print(available_celllines_cnv)
	#print(available_celllines_dr_response)
	#print(available_celllines_ge)
	
	duplicate_measurements = len(available_celllines_dr_response)
	
	#print(str(len(available_celllines_dr_response)))
	
	output_string = output_file + drug_name + ".txt"
	available_celllines = []
	
	if not len(available_celllines_dr_response) ==1:
		print("Found drug twice")
	if not len(available_celllines_dr_response[0]) == 1:
		
		for element in available_celllines_cnv:
			if element in available_celllines_mutation:
				if element in available_celllines_dr_response[0]:
					if element in available_celllines_ge:
						available_celllines.append(element)
	else:		
		with open(output_string, 'w', encoding='utf-8') as output:
			output.write(available_celllines_dr_response[j][0])

	print(len(available_celllines))

	with open(output_string, 'w', encoding='utf-8') as output:
		for cl in available_celllines:
			
			output.write(cl + "\n")

def available_cell_lines_old_drug_data_mut_ge(drug_name, output_file, gdsc1_or_2):
	
	available_celllines = []
	available_celllines_mutation = available_cl_mut()
	available_celllines_dr_response = available_cl_drug(drug_name, gdsc1_or_2)
	available_celllines_ge = available_cl_ge()
	
	#print(available_celllines_cnv)
	#print(available_celllines_dr_response)
	#print(available_celllines_ge)
	
	duplicate_measurements = len(available_celllines_dr_response)
	
	for j in range(0, duplicate_measurements):
		if "/" in drug_name:
			drug_name = drug_name.replace("/", "--")
		output_string = output_file  + drug_name + "_" + str(j) + ".txt"
		available_celllines = []#empty cell lines for each duplicate measurement
		if not len(available_celllines_dr_response[j]) == 1:
			
			for element in available_celllines_mutations:
				if element in available_celllines_dr_response[j]:
					if element in available_celllines_ge:
						available_celllines.append(element)
		else:
			
			with open(output_string, 'w', encoding='utf-8') as output:
				output.write(available_celllines_dr_response[j][0])
		
		print(len(available_celllines))
		with open(output_string, 'w', encoding='utf-8') as output:
			for cl in available_celllines:
				
				output.write(cl + "\n")

def available_cell_lines_old_drug_data_cnv_mut(drug_name, output_file, gdsc1_or_2):
	
	available_celllines = []
	available_celllines_cnv = available_cl_cnv()
	available_celllines_mutation = available_cl_mut()
	available_celllines_dr_response = available_cl_drug(drug_name, gdsc1_or_2)

	
	#print(available_celllines_cnv)
	#print(available_celllines_dr_response)
	#print(available_celllines_ge)
	
	duplicate_measurements = len(available_celllines_dr_response)
	
	for j in range(0, duplicate_measurements):
		if "/" in drug_name:
			drug_name = drug_name.replace("/", "--")
		output_string = output_file + drug_name + "_" + str(j) + ".txt"
		available_celllines = []#empty cell lines for each duplicate measurement
		if not len(available_celllines_dr_response[j]) == 1:
			
			for element in available_celllines_cnv:
				if element in available_celllines_dr_response[j]:
					if element in available_celllines_mutation:
						available_celllines.append(element)
		else:
			
			with open(output_string, 'w', encoding='utf-8') as output:
				output.write(available_celllines_dr_response[j][0])
		
		print(len(available_celllines))
		with open(output_string, 'w', encoding='utf-8') as output:
			for cl in available_celllines:
				
				output.write(cl + "\n")

def available_cell_lines_old_drug_data_cnv_ge(drug_name, output_file, gdsc1_or_2):
	
	available_celllines = []
	available_celllines_cnv = available_cl_cnv()
	available_celllines_ge = available_cl_ge()
	available_celllines_dr_response = available_cl_drug(drug_name, gdsc1_or_2)

	
	#print(available_celllines_cnv)
	#print(available_celllines_dr_response)
	#print(available_celllines_ge)
	
	duplicate_measurements = len(available_celllines_dr_response)
	
	for j in range(0, duplicate_measurements):
		if "/" in drug_name:
			drug_name = drug_name.replace("/", "--")
		output_string = output_file  + drug_name + "_" + str(j) + ".txt"
		available_celllines = []#empty cell lines for each duplicate measurement
		if not len(available_celllines_dr_response[j]) == 1:
			
			for element in available_celllines_cnv:
				if element in available_celllines_dr_response[j]:
					if element in available_celllines_ge:
						available_celllines.append(element)
		else:
			
			with open(output_string, 'w', encoding='utf-8') as output:
				output.write(available_celllines_dr_response[j][0])
		
		print(len(available_celllines))
		with open(output_string, 'w', encoding='utf-8') as output:
			for cl in available_celllines:
				
				output.write(cl + "\n")


def available_cell_lines_old_drug_data_cnv(drug_name, output_file, gdsc1_or_2):
	
	available_celllines = []
	available_celllines_cnv = available_cl_cnv()
	available_celllines_dr_response = available_cl_drug(drug_name, gdsc1_or_2)

	
	#print(available_celllines_cnv)
	#print(available_celllines_dr_response)
	#print(available_celllines_ge)
	
	duplicate_measurements = len(available_celllines_dr_response)
	
	for j in range(0, duplicate_measurements):
		if "/" in drug_name:
			drug_name = drug_name.replace("/", "--")
		output_string = output_file  + drug_name + "_" + str(j) + ".txt"
		available_celllines = []#empty cell lines for each duplicate measurement
		if not len(available_celllines_dr_response[j]) == 1:
			
			for element in available_celllines_cnv:
				if element in available_celllines_dr_response[j]:
					available_celllines.append(element)
		else:
			
			with open(output_string, 'w', encoding='utf-8') as output:
				output.write(available_celllines_dr_response[j][0])
		
		print(len(available_celllines))
		with open(output_string, 'w', encoding='utf-8') as output:
			for cl in available_celllines:
				
				output.write(cl + "\n")




def available_cell_lines_old_drug_data_mut(drug_name, output_file, gdsc1_or_2):
	
	available_celllines = []
	available_celllines_mutation = available_cl_mut()
	available_celllines_dr_response = available_cl_drug(drug_name, gdsc1_or_2)

	
	#print(available_celllines_cnv)
	#print(available_celllines_dr_response)
	#print(available_celllines_ge)
	
	duplicate_measurements = len(available_celllines_dr_response)
	
	for j in range(0, duplicate_measurements):
		if "/" in drug_name:
			drug_name = drug_name.replace("/", "--")
		output_string = output_file + drug_name + "_" + str(j) + ".txt"
		available_celllines = []#empty cell lines for each duplicate measurement
		if not len(available_celllines_dr_response[j]) == 1:
			
			for element in available_celllines_mutation:
				if element in available_celllines_dr_response[j]:
					available_celllines.append(element)
		else:
			
			with open(output_string, 'w', encoding='utf-8') as output:
				output.write(available_celllines_dr_response[j][0])
		
		print(len(available_celllines))
		with open(output_string, 'w', encoding='utf-8') as output:
			for cl in available_celllines:
				
				output.write(cl + "\n")

def available_cell_lines_old_drug_data_ge(drug_name, output_file, gdsc1_or_2):
	
	available_celllines = []
	available_celllines_dr_response = available_cl_drug(drug_name, gdsc1_or_2)
	available_celllines_ge = available_cl_ge()
	
	#print(available_celllines_cnv)
	#print(available_celllines_dr_response)
	#print(available_celllines_ge)
	
	duplicate_measurements = len(available_celllines_dr_response)
	print(str(len(available_celllines_dr_response)))
	
	for j in range(0, duplicate_measurements):
		if "/" in drug_name:
			drug_name = drug_name.replace("/", "--")
		output_string = output_file  + drug_name + "_" + str(j) + ".txt"
		available_celllines = []#empty cell lines for each duplicate measurement
		if not len(available_celllines_dr_response[j]) == 1:
			
			for element in available_celllines_ge:
				if element in available_celllines_dr_response[j]:
					available_celllines.append(element)
		else:
			
			with open(output_string, 'w', encoding='utf-8') as output:
				output.write(available_celllines_dr_response[j][0])
		
		print(len(available_celllines))
		with open(output_string, 'w', encoding='utf-8') as output:
			for cl in available_celllines:
				
				output.write(cl + "\n")


			
################Main function ##################################

def main(drug_name, gene_expression_incl, mut_incl, cnv_incl, output_file, gdsc1_or_2):
	
	
	if gene_expression_incl == "yes" and mut_incl == "yes" and cnv_incl == "yes":
		available_cell_lines_old_drug_data_cnv_mut_ge(drug_name, output_file, gdsc1_or_2)
		
	elif gene_expression_incl == "yes" and mut_incl == "yes":
		available_cell_lines_old_drug_data_mut_ge(drug_name, output_file, gdsc1_or_2)
		
	elif gene_expression_incl == "yes" and cnv_incl == "yes":
		available_cell_lines_old_drug_data_cnv_ge(drug_name, output_file, gdsc1_or_2)
		
	elif mut_incl == "yes" and cnv_incl == "yes":
		available_cell_lines_old_drug_data_cnv_mut(drug_name, output_file, gdsc1_or_2)
		
	elif gene_expression_incl == "yes":
		available_cell_lines_old_drug_data_ge(drug_name, output_file, gdsc1_or_2)
	elif mut_incl == "yes":
		available_cell_lines_old_drug_data_mut(drug_name, output_file, gdsc1_or_2)
	elif cnv_incl == "yes":
		available_cell_lines_old_drug_data_cnv(drug_name, output_file, gdsc1_or_2)
	return



#################################################################
if __name__ == "__main__":
	if len(sys.argv) < 7:
			sys.exit("This program needs the following arguments:\
					\n- drug name\
					\n- should GDSC1 or GDSC2 be used ('GDSC1' or 'GDSC2')\
					\n- yes/no (shall gene expression be included)\
					\n- yes/no (shall mutations be included)\
					\n- yes/no (shall CNVs be included)\
					\n- output file name (including path but without ending)") 

	main(sys.argv[1], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[2]) 	
	
