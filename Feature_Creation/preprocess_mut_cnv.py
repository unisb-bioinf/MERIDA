#Author:Kerstin Lenhof
#Date:11.12.2018
#Script for (point) mutation and copy number data matrix from GDSC

import sys
sys.path.insert(0, '<path_to_Mapping_Facility>')
sys.path.insert(0, '<path_to_Oncogenicity_Annotation>')
sys.path.insert(0, '<path_to_Sensitivity_Annotation>')
from cnv_mutation_annotator import *
from alteration_matrix_container import *
from query_oncogenicity_db import *
from query_sensitivity_db import *


copy_number_data = "<path_to_copy_number_data_on_gene_level>"
mutation_data ="<path_to_mutation_file_from_GDSC>"


def preprocess_copy_number_data_for_matrix():
	
	#samples with attached alterations in genes 
	sample_to_gene_to_alterations_dict ={}
	with open(copy_number_data, 'r', encoding='utf-8') as cnv_data:
		first_line = cnv_data.readline()
		
		for line in cnv_data:
			sline = line.strip('\n').split('\t')
			
			gene_name = sline[1]
			gain_or_loss = sline[2]
			sample_id = sline[0]
			
			if not sample_id in sample_to_gene_to_alterations_dict:
				sample_to_gene_to_alterations_dict[sample_id] = {gene_name: gain_or_loss}
			else:
				if not gene_name in sample_to_gene_to_alterations_dict[sample_id]:
					sample_to_gene_to_alterations_dict[sample_id][gene_name] = gain_or_loss
				else:
					print("Please check file, double occurrence of gene name in a sample. " + "\n" + "Sample: " + sample_id + "\n" + "Gene: " + gene_name + "\n\n")
			
	return sample_to_gene_to_alterations_dict

def preprocess_mutation_data_for_matrix():
	
	#samples with attached alterations in genes
	sample_to_gene_to_alteration_dict ={}
	sample_to_gene_to_truncating_dict = {}
	
	with open(mutation_data, 'r', encoding='utf-8') as mut_data:
		
		first_line = mut_data.readline()
		
		for line in mut_data:
			sline = line.strip('\n').split('\t')
			
			cosmic_id = sline[1].strip()
			gene_name = sline[3].strip()
			amino_acid_change_withp = sline[6].strip()
			is_truncating_info = sline[7].strip()
			
			is_truncating = ""
			if is_truncating_info in ["ess_splice", "frameshift", "nonsense"]:
				is_truncating = True
			else:
				is_truncating = False
			amino_acid_change = amino_acid_change_withp.replace('p.', '')
			
			if not cosmic_id in sample_to_gene_to_alteration_dict:
				
				sample_to_gene_to_alteration_dict[cosmic_id] = {gene_name:[amino_acid_change]}
				
			else:
				if not gene_name in sample_to_gene_to_alteration_dict[cosmic_id]:
					sample_to_gene_to_alteration_dict[cosmic_id][gene_name] = [amino_acid_change]
					
				else:
					sample_to_gene_to_alteration_dict[cosmic_id][gene_name].append(amino_acid_change)
			
			
			if not cosmic_id in sample_to_gene_to_truncating_dict:
				
				sample_to_gene_to_truncating_dict[cosmic_id] = {gene_name:[is_truncating]}
				
			else:
				if not gene_name in sample_to_gene_to_truncating_dict[cosmic_id]:
					sample_to_gene_to_truncating_dict[cosmic_id][gene_name] = [is_truncating]
					
				else:
					sample_to_gene_to_truncating_dict[cosmic_id][gene_name].append(is_truncating)
	
	return [sample_to_gene_to_alteration_dict, sample_to_gene_to_truncating_dict]



def get_available_cl_list(available_cl_fn):
	available_cl = []
	with open(available_cl_fn, 'r', encoding='utf-8') as a_cl_file:
		for line in a_cl_file:
			cl = line.strip()
			if not cl == "":
				available_cl.append(cl)
				
	print(available_cl)
	return available_cl


def get_gene_list(gene_set):
	
	gene_list = []
	
	with open(gene_set, 'r', encoding='utf-8') as gene_set_file:
		
		for line in gene_set_file:
			gene_name = line.strip("\n")
			
			if not gene_name in gene_list:
				gene_list.append(gene_name)
	
	return gene_list

def prepare_matrix_mut_cnv(drug_name, gene_set_mut, gene_set_cnv, available_cl, output_file, gene_onco_mode):
	
	print("Entered the preparation method")
	#List with available cell lines
	available_cl = get_available_cl_list(available_cl)
	
	genes_to_be_considered_mutations = []
	genes_to_be_considered_cnvs = []
	
	if not gene_set_mut == "all":
		genes_to_be_considered_mutations = get_gene_list(gene_set_mut)
		
	if not gene_set_cnv == "all":
		genes_to_be_considered_cnvs = get_gene_list(gene_set_cnv)
	
	
	#GDSC CNV data 
	cnv_sample_to_gene_to_loss_or_gain = preprocess_copy_number_data_for_matrix()
	#GDSC Mutation data
	mut_sample_to_gene_to_alterations_and_truncation = preprocess_mutation_data_for_matrix()
	mut_sample_to_gene_to_alterations = mut_sample_to_gene_to_alterations_and_truncation[0]
	mut_sample_to_gene_to_truncating = mut_sample_to_gene_to_alterations_and_truncation[1]
	
	
	#Initialize the annotator
	annotator = CNV_Mutation_Annotator(drug_name, gene_onco_mode)
	
	print("Loaded all data sets")
	
	#Initialize the container
	container = Alteration_Matrix_Container(available_cl)
	
	#Annotate all available mutations and gather information in matrix container class
	for sample in available_cl:
		
		for gene in mut_sample_to_gene_to_alterations[sample]:
			
			if gene in genes_to_be_considered_mutations or len(genes_to_be_considered_mutations) == 0:
				#For truncating information
				alteration_counter = 0
				for alteration in mut_sample_to_gene_to_alterations[sample][gene]:
					
					annotations = annotator.annotate_mutation(gene, alteration, mut_sample_to_gene_to_truncating[sample][gene][alteration_counter])
					
					if not annotations:
						print("Sample " + sample + " Gene " + gene + " Alteration " + alteration)
						print(annotations)
					
					if annotations[2] in ["Likely Loss-of-function", "Loss-of-function"]:
						
						if annotations[3] == "":
							container.register_mutation(sample, annotations[0] + "_" + "Loss-of-function", "Loss-of-function", "Unknown")
						else:
							container.register_mutation(sample, annotations[0] + "_" + "Loss-of-function", "Loss-of-function", annotations[3])

					elif annotations[2] in ["Likely Gain-of-function", "Gain-of-function"]:
						
						if annotations[3] == "":
							container.register_mutation(sample, annotations[0] + "_" + "Gain-of-function", "Gain-of-function", "Unknown")
						else:
							container.register_mutation(sample, annotations[0] + "_" + "Gain-of-function", "Gain-of-function", annotations[3])
								
		
					elif annotations[2] in ["Likely Neutral", "Neutral"]:
						
						if annotations[3] == "":
							container.register_mutation(sample, annotations[0] + "_" + "Neutral" , "Neutral", "Unknown")
						else:
							container.register_mutation(sample, annotations[0] + "_" + "Neutral" , "Neutral", annotations[3])
					elif annotations[2]== "oncogenic":
						
						print("oncogenic match found")
						if annotations[3] == "":
							container.register_mutation(sample, annotations[0] + "_" + annotations[2].replace(" ", "-"), annotations[2], "Unknown")
						else:
							container.register_mutation(sample, annotations[0] + "_" + annotations[2].replace(" ", "-"), annotations[2], annotations[3])
							

					elif annotations[2] == "Unknown":
						container.register_mutation(sample, annotations[0] + "_" + annotations[2].replace(" ", "-"), annotations[2], "Unknown")
					else:
						if annotations[2] =="":
							container.register_mutation(sample, annotations[0] + "_" + annotations[1], "Unknown", annotations[3])
						else:
							container.register_mutation(sample, annotations[0] + "_" + annotations[1], annotations[2], annotations[3])
					
					alteration_counter = alteration_counter + 1
				
				
	#Annotate copy number variations
	for sample in available_cl:
		
		for gene in cnv_sample_to_gene_to_loss_or_gain[sample]:
			
			if gene in genes_to_be_considered_cnvs or len(genes_to_be_considered_cnvs) == 0:
				annotations = annotator.annotate_cnv(gene, cnv_sample_to_gene_to_loss_or_gain[sample][gene])
				
				#print(annotations)
				if not annotations[0] == "The given alteration was not suitable for testing, please refer to the method description":
					
					if annotations[3] == "":
						container.register_cnv(sample, annotations[0] + "_" + annotations[1], annotations[2], "Unknown")
					else:
						container.register_cnv(sample, annotations[0] + "_" + annotations[1], annotations[2], annotations[3])

				else:
					print(annotations[0])

	#print results into file :)
	
	
	container.print_matrix(output_file)
	

	
	
						
################Main function ##################################

def main(drug_name, gene_set_mut, gene_set_cnv, available_cl, output_file, gene_onco_mode):
	
	prepare_matrix_mut_cnv(drug_name, gene_set_mut, gene_set_cnv, available_cl, output_file, gene_onco_mode)
	return



#################################################################
if __name__ == "__main__":
	if len(sys.argv) < 5:
			sys.exit("This program needs the following arguments:\
					\n- drug name\
					\n- path to a file with a gene set that shall be considered for mutations (can also be 'all' ---> not recommended)\
					\n- path to a file with a gene set that shall be considered for copy number variations (can also be 'all' --> not recommended)\
					\n- name of a file with list of available cell lines\
					\n- output file name (including path)\
					\n- mode for oncogenicity search in genes (gene_onco (strict mode), gene_onco_less_strict (less strict mode))") 

	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6]) 	
	 
