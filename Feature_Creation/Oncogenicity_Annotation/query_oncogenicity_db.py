#Author: Kerstin Lenhof
#Date: 9.11.2018
#Query to get all entries that have something to do with a certain gene in the 'oncogenicity biomarker' database

output_onco1 =  "<path_to_cgi_file>" 
output_onco2 = "<path_to_CIViC_file>"
output_onco3 = "<path_to_oncoKB_file>"


import sys
sys.path.insert(0, '<path_to_Mapping_Facility>')

from useful_tools import * 
from mapping_reader import *

def get_list_of_genes():
	all_genes = []
	onco_db_list = [output_onco1,output_onco2, output_onco3]
	
	for db_file in onco_db_list: 

		with open(db_file, 'r', encoding='utf-8') as current_db:
			
			first_line = current_db.readline()
		
		
			for line in current_db:
				
				sline = line.replace('\n', '').split('\t')
				
				gene_name = sline[0]
				
				if '--' in gene_name:
					two_names = gene_name.split('--')
					
					if not two_names[0].strip() in all_genes:
						all_genes.append(two_names[0].strip())
					if not two_names[1].strip() in all_genes:
						all_genes.append(two_names[1].strip())
						
				elif ';' in gene_name:
					two_names = gene_name.split(';')
					if not two_names[0].strip() in all_genes:
						all_genes.append(two_names[0].strip())
					if not two_names[1].strip() in all_genes:
						all_genes.append(two_names[1].strip())
				else:	
					if not sline[0].strip() in all_genes:
						all_genes.append(sline[0].strip())
	
	return all_genes


def get_onco_db(query_gene_name):
	
	return query_onco_db(query_gene_name, "no_output_file")


def get_onco_db_file(output_file_name):
	
	all_genes = []
	onco_db_list = [output_onco1, output_onco2, output_onco3]
	
	captions = ""
	entries = []
	for db_file in onco_db_list: 

		with open(db_file, 'r', encoding='utf-8') as current_db:
			
			captions = current_db.readline()	
			
			for line in current_db:
				
				entries.append(line.strip('\n'))
				
	with open(output_file_name, 'w', encoding='utf-8') as output:
		
		output.write(captions)
		for entry in entries:
			
			output.write(entry + "\n")
	

def query_onco_db(query_gene_name, output_file):
	
	#print(query_gene_name)
	official_gene_symbols = load_gene_symbols()
	alias_to_symbol = normalized_alias_map()


	db_list =[output_onco1, output_onco2, output_onco3]
	
	relevant_entries = []
	for db in db_list:
		
		with open(db, 'r', encoding='utf-8') as current_db:
			
			current_db.readline()
			
			for line in current_db:
				
				sline = line.split("\t")
				
				gene_name = sline[0]
				
				
				
				several_gene_names = gene_name.split("--")
				
				for gn in several_gene_names:
					
					if normalize_caseless(query_gene_name.strip()) == normalize_caseless(gn.strip()):
						
						relevant_entries.append(line)
						
					elif normalize_caseless(query_gene_name.strip()) in alias_to_symbol:
						
						norm_common_name_of_query = alias_to_symbol[normalize_caseless(query_gene_name.strip())]
						
						
						if not len(norm_common_name_of_query) == 1:
							print('Could not find a unique identifier for the query gene name. Name: ' + query_gene_name)
						else:
							if norm_common_name_of_query[0].strip() == normalize_caseless(gn.strip()):
								relevant_entries.append(line)
					
	
	
	if output_file == "no_output_file":
		
		return relevant_entries
	
	else:
		with open(output_file, 'w', encoding='utf-8') as relevant_output:
			
			relevant_output.write("gene_name" + "\t" + "HGNC_ID"  + "\t" +"alteration" + "\t" + "alteration_type" + "\t"  + "mutation_effect" + "\t" + "association_implication" + "\t" + "primary_tumour_type" + "\t" + "tsg_onco_annotation" + "\t" + "source" + "\t" +  "source_database" + "\t" + "evidence_level" + "\t" + "evidence_statement" + "\t" + "cds_like_mutation" + "\t" + "genome_coordinates" + "\t" + "reference_build" + "\n")
			
			for element in relevant_entries:
				relevant_output.write(element)
		
def query_filtered_onco_db(query_gene_name, output_file, mode):
	
	if mode in ["gene_onco", "gene_neutral"]:
		return query_filtered_onco_db_strict(query_gene_name, output_file, mode)
		
	elif mode == "gene_onco_less_strict":
		return query_filtered_onco_db_less_strict(query_gene_name, output_file, mode)
	
def query_filtered_onco_db_strict(query_gene_name, output_file, mode):

	all_entries = get_onco_db(query_gene_name)
	
	filtered_entries_onco = {}
	filtered_entries_inconclusive  = {}
	filtered_entries_neutral = {}
	
	
	for entry in all_entries:
		
		splitted_entry = entry.strip('\n').split('\t')
		
		alteration = splitted_entry[2].strip()
		mutation_effect = splitted_entry[4].strip()
		association_implication = splitted_entry[5].strip()
		source_database = splitted_entry[9].strip()
		
		if association_implication in  ["validated oncogenic", "Predisposing:Pathogenic", "Predisposing:Likely Pathogenic", "Prognostic:Poor Outcome", "Diagnostic:Positive", "Oncogenic", "Likely Oncogenic"]:
			#only useful entries for comparison with gene expression data in later steps + deletion/amplification etc is removed, which does not matter as there are two columns for deletion and amplification (loss and gain) for CNV anyway
			if mutation_effect in ["Likely Loss-of-function", "Loss-of-function", "Gain-of-function", "Likely Gain-of-function"] or "fs" in alteration:
				gene_name = splitted_entry[0].strip()
				#alteration = splitted_entry[2].strip()
				hash_ident = gene_name + alteration 
				
				
				if not hash_ident in filtered_entries_onco:
					filtered_entries_onco[hash_ident] = entry
				
				elif source_database in ["oncoKB", "CIViC"]:
					filtered_entries_onco[hash_ident] = entry
				
		elif association_implication == "Inconclusive":
			gene_name = splitted_entry[0].strip()
			#alteration = splitted_entry[2].strip()
			hash_ident = gene_name + alteration 
			
			if not hash_ident in filtered_entries_onco:
				filtered_entries_inconclusive[hash_ident] = entry
				
		elif association_implication in ["Neutral", "Likely Neutral"]:
			gene_name = splitted_entry[0].strip()
			#alteration = splitted_entry[2].strip()
			hash_ident = gene_name + alteration 
			
			#print("neutral ???")
			#print(alteration)
			if not hash_ident in filtered_entries_onco:
				if not hash_ident in filtered_entries_inconclusive:
					#print("neutral found")
					filtered_entries_neutral[hash_ident] = entry
				elif source_database in ["oncoKB", "CIViC"]:
					filtered_entries_neutral[hash_ident] = entry
	
	if mode == "gene_onco":
		
		if output_file == "no_output_file":
			return filtered_entries_onco
		
		else:
			with open(output_file, 'w', encoding='utf-8') as relevant_output:
			
				relevant_output.write("gene_name" + "\t" + "HGNC_ID"  + "\t" +"alteration" + "\t" + "alteration_type" + "\t"  + "mutation_effect" + "\t" + "association_implication" + "\t" + "primary_tumour_type" + "\t" + "tsg_onco_annotation" + "\t" + "source" + "\t" +  "source_database" + "\t" + "evidence_level" + "\t" + "evidence_statement" + "\t" + "cds_like_mutation" + "\t" + "genome_coordinates" + "\t" + "reference_build" + "\n")
				
				for element in filtered_entries_onco:
					relevant_output.write(filtered_entries_onco[element])
				
	elif mode == "gene_neutral":
		
		if output_file == "no_output_file":
			return filtered_entries_neutral
		
		else:
			with open(output_file, 'w', encoding='utf-8') as relevant_output:
		
				relevant_output.write("gene_name" + "\t" + "HGNC_ID"  + "\t" +"alteration" + "\t" + "alteration_type" + "\t"  + "mutation_effect" + "\t" + "association_implication" + "\t" + "primary_tumour_type" + "\t" + "tsg_onco_annotation" + "\t" + "source" + "\t" +  "source_database" + "\t" + "evidence_level" + "\t" + "evidence_statement" + "\t" + "cds_like_mutation" + "\t" + "genome_coordinates" + "\t" + "reference_build" + "\n")
				
				for element in filtered_entries_neutral:
					relevant_output.write(filtered_entries_neutral[element])
	
	
def query_filtered_onco_db_less_strict(query_gene_name, output_file, mode):

	all_entries = get_onco_db(query_gene_name)
	
	filtered_entries_onco = {}
	filtered_entries_inconclusive  = {}
	filtered_entries_neutral = {}
	
	
	for entry in all_entries:
		
		splitted_entry = entry.strip('\n').split('\t')
		
		alteration = splitted_entry[2].strip()
		mutation_effect = splitted_entry[4].strip()
		association_implication = splitted_entry[5].strip()
		source_database = splitted_entry[9].strip()
		
		#print(association_implication)
		if association_implication in  ["validated oncogenic", "Predisposing:Pathogenic", "Predisposing:Likely Pathogenic", "Prognostic:Poor Outcome", "Diagnostic:Positive", "Oncogenic", "Likely Oncogenic"]:
			#only useful entries for comparison with gene expression data in later steps + deletion/amplification etc is removed, which does not matter as there are two columns for deletion and amplification (loss and gain) for CNV anyway
			
			gene_name = splitted_entry[0].strip()
			#alteration = splitted_entry[2].strip()
			hash_ident = gene_name + alteration 
			
			
			if not hash_ident in filtered_entries_onco:
				filtered_entries_onco[hash_ident] = entry
			
			elif source_database == "oncoKB":
				filtered_entries_onco[hash_ident] = entry
				
		elif association_implication == "Inconclusive":
			gene_name = splitted_entry[0].strip()
			#alteration = splitted_entry[2].strip()
			hash_ident = gene_name + alteration 
			
			if not hash_ident in filtered_entries_onco:
				filtered_entries_inconclusive[hash_ident] = entry
				
		elif association_implication in ["Neutral", "Likely Neutral"]:
			gene_name = splitted_entry[0].strip()
			#alteration = splitted_entry[2].strip()
			hash_ident = gene_name + alteration 
			
			#print("neutral ???")
			#print(alteration)
			if not hash_ident in filtered_entries_onco:
				if not hash_ident in filtered_entries_inconclusive:
					#print("neutral found")
					filtered_entries_neutral[hash_ident] = entry
				elif source_database in ["oncoKB", "CIViC"]:
					filtered_entries_neutral[hash_ident] = entry
	
	if mode == "gene_onco_less_strict":
		
		if output_file == "no_output_file":
			return filtered_entries_onco
		
		else:
			with open(output_file, 'w', encoding='utf-8') as relevant_output:
			
				relevant_output.write("gene_name" + "\t" + "HGNC_ID"  + "\t" +"alteration" + "\t" + "alteration_type" + "\t"  + "mutation_effect" + "\t" + "association_implication" + "\t" + "primary_tumour_type" + "\t" + "tsg_onco_annotation" + "\t" + "source" + "\t" +  "source_database" + "\t" + "evidence_level" + "\t" + "evidence_statement" + "\t" + "cds_like_mutation" + "\t" + "genome_coordinates" + "\t" + "reference_build" + "\n")
				
				for element in filtered_entries_onco:
					relevant_output.write(filtered_entries_onco[element])
				
	if mode == "gene_neutral":
		
		if output_file == "no_output_file":
			return filtered_entries_neutral
		
		else:
			with open(output_file, 'w', encoding='utf-8') as relevant_output:
		
				relevant_output.write("gene_name" + "\t" + "HGNC_ID"  + "\t" +"alteration" + "\t" + "alteration_type" + "\t"  + "mutation_effect" + "\t" + "association_implication" + "\t" + "primary_tumour_type" + "\t" + "tsg_onco_annotation" + "\t" + "source" + "\t" +  "source_database" + "\t" + "evidence_level" + "\t" + "evidence_statement" + "\t" + "cds_like_mutation" + "\t" + "genome_coordinates" + "\t" + "reference_build" + "\n")
				
				for element in filtered_entries_neutral:
					relevant_output.write(filtered_entries_neutral[element])

################Main function ##################################

def main(gene_name, output_file, mode):
	
	print(mode)
	if mode == "all":
		get_onco_db_file(output_file)
	elif mode == "gene":
		query_onco_db(gene_name, output_file)
	elif mode == "gene_onco" or mode == "gene_neutral":
		query_filtered_onco_db_strict(gene_name, output_file, mode)
	elif mode == "gene_onco_less_strict":
		query_filtered_onco_db_less_strict(gene_name, output_file, mode)
	else:
		print("Unsupported mode given:" + mode)
	return



#################################################################
if __name__ == "__main__":
	if len(sys.argv) < 3:
			sys.exit("This program needs the following arguments:\
					\n- mode\
					\n- gene name (can be 'none' if mode is equal to 'all')\
					\n- output file name (including path)\
					\n Supported modes are: all, gene, gene_onco, gene_neutral") 

	main(sys.argv[2], sys.argv[3], sys.argv[1]) 	
	
