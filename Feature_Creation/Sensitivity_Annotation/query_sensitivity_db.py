#Author: Kerstin Lenhof
#Date: 8.11.2018
#Query to get all entries that have something to do with a certain drug or with a certain gene in the "sensitivity biomarker" database
output1 = "<path_to_CGI>"
output2 = "<path_to_CIViC>"
output3 = "<path_to_oncoKB>"
output4 = "<path_to_COSMIC>"

import sys
sys.path.insert(0, '<path_to_Mapping_Facility>')

from useful_tools import * 
from mapping_reader import *
from query_oncogenicity_db import *

def get_drug_db(query_drug_name):
	
	return query_drug_db(query_drug_name, "no_output_file")
	
def get_gene_db(query_gene_name):
	
	return query_gene_db(query_gene_name, "no_output_file")

def query_drug_db(query_drug_name, output_file):
	
	normalized_syndrug_to_common_dict = load_drugsynnormalized_2_commonnormalized_map()
	normalized_to_common_dict = load_drugnormalized_2_commonname_map()
	
	db_list =[output1,output2, output3, output4]
	
	relevant_entries = []
	for db in db_list:
		
		with open(db, 'r', encoding='utf-8') as current_db:
			
			current_db.readline()
			
			for line in current_db:
				
				sline = line.split("\t")
				
				drug_name = sline[0]
				
				
				several_drug_names = drug_name.split("+")
				
				for dr in several_drug_names:
					
					if normalize_caseless(query_drug_name.strip()) == normalize_caseless(dr.strip()):
						
						relevant_entries.append(line)
						
					elif normalize_caseless(query_drug_name.strip()) in normalized_syndrug_to_common_dict:
						
						norm_common_name_of_query = normalized_syndrug_to_common_dict[normalize_caseless(query_drug_name.strip())]
						
						
						if norm_common_name_of_query.strip() == normalize_caseless(dr.strip()):
							relevant_entries.append(line)
							
						else:
							norm_name_of_current_drug = normalize_caseless(dr.strip())
							
							if norm_name_of_current_drug in normalized_syndrug_to_common_dict:
								norm_common_name_of_current_drug = normalized_syndrug_to_common_dict[norm_name_of_current_drug]
								
								if norm_common_name_of_current_drug.strip() == norm_common_name_of_query.strip():
									
									relevant_entries.append(line)
								
					
	
	if output_file == "no_output_file":
		
		return relevant_entries

	else:
	
		with open(output_file, 'w', encoding='utf-8') as relevant_output:
			
			relevant_output.write("drug_name" + "\t" + "DrugBank_ID" + "\t" + "drug_family" + "\t" + "drug_targets" + "\t" + "mutation_effect" + "\t" + "associated_gene" + "\t" + "HGNC_ID_associated_gene" + "\t" + "pathways_associated_gene" + "\t" + "alteration_type" + "\t"+ "alteration" + "\t" + "association_implication" + "\t" + "evidence" + "\t" + "source" + "\t" + "primary_tumour_type" + "\t" + "tsg_onco_annotation" + "\t" + "source_database" + "\t" + "evidence_statement" + "\t" + "cds_like_mutation" + "\t" + "genome_coordinates" + "\t" + "reference_build" + "\n")
			
			for element in relevant_entries:
				relevant_output.write(element)

def get_elements_from_onco_db(gene_name, alteration):

	onco_db_entries = get_onco_db(gene_name)
	
	relevant_entries = []
	for entry in onco_db_entries:
		splitted_entry = entry.strip('\n').split('\t')
		mutation_effect = splitted_entry[4].strip()
		association_implication = splitted_entry[5].strip()
		
		if association_implication == "Inconclusive" or association_implication == "Likely Neutral" or association_implication == "Neutral":
			continue
		else:
			if alteration == "LOSS-OF-FUNCTION":
				if mutation_effect == "Loss-of-function" or mutation_effect == "Likely Loss-of-function":
					relevant_entries.append(entry)
			elif alteration == "GAIN-OF-FUNCTION":
				if mutation_effect == "Gain-of-function" or mutation_effect == "Likely Gain-of-function":
					relevant_entries.append(entry)
			elif alteration == "ONCOGENIC MUTATION" or alteration == "Oncogenic mutations":
				
				if association_implication == "Predisposing:Pathogenic" or association_implication == "Predisposing:Likely Pathogenic" or  association_implication == "Prognostic:Poor Outcome" or association_implication == "Diagnostic:Positive" or association_implication == "validated oncogenic" or association_implication == "Oncogenic" or association_implication == "Likely Oncogenic":
					relevant_entries.append(entry)
		
	return relevant_entries

def merge_sens_onco(o_entry, onco_db_entries):
	
	merged_entries = []
	for entry in onco_db_entries:
		
		splitted_entry = entry.strip('\n').split('\t')
		alteration = splitted_entry[2].strip()
		alteration_type = splitted_entry[3].strip()
		mutation_effect = splitted_entry[4].strip()
		evidence = splitted_entry[10].strip()
		evidence_statement = splitted_entry[11].strip()
		source = splitted_entry[8].strip()
		source_database = splitted_entry[9].strip()
		primary_tumour_type = splitted_entry[6].strip()
		association_implication = splitted_entry[5].strip()
		
		o_entry_splitted = o_entry.strip('\n').split('\t')
		
		o_entry_splitted[4] = mutation_effect
		o_entry_splitted[8] = alteration_type 
		o_entry_splitted[9] = alteration
		o_entry_splitted[11] = o_entry_splitted[11].strip() + "|" + evidence
		o_entry_splitted[12] = o_entry_splitted[12].strip() + "|" + source
		o_entry_splitted[13] = o_entry_splitted[13].strip() + "|" + primary_tumour_type
		o_entry_splitted[15] = o_entry_splitted[15].strip() + "|" + source_database
		o_entry_splitted[16] = o_entry_splitted[16].strip() + "|" + evidence_statement
		new_entry = "\t".join(o_entry_splitted)
		merged_entries.append(new_entry)
	
	return merged_entries
	
def compare_association(association_duplicate, association_new):
	
	if association_duplicate == association_new:
		return 0
	elif association_duplicate == "resistant" and association_new == "sensitive":
		return 1
	elif association_duplicate == "sensitive" and association_new == "resistant":
		return 1
	elif association_duplicate == "resistant" and association_new == "not resistant":
		return 2
	elif association_duplicate == "resistant" and association_new == "not sensitive":
		return 2
	elif association_duplicate == "sensitive" and association_new == "not sensitive":
		return 2
	elif association_duplicate == "sensitive" and association_new == "not resistant":
		return 2
	elif association_duplicate == "not resistant" and association_new == "resistant":
		return 3
	elif association_duplicate == "not sensitive" and association_new == "resistant":
		return 3
	elif association_duplicate == "not resistant" and association_new == "sensitive":
		return 3
	elif association_duplicate == "not sensitive" and association_new == "sensitive":
		return 3
	elif association_duplicate == "not resistant" and association_new == "not sensitive":
		return 4
	elif association_duplicate == "not sensitive" and association_new == "not resistant":
		return 4
	else:
		return 50 
	
	
def query_drug_filtered(query_drug_name, output_file):
	
	all_entries = get_drug_db(query_drug_name)
	#print(all_entries)
	
	filtered_entries = {}
	
	for entry in all_entries:
		
		splitted_entry = entry.strip('\n').split('\t')
		len(splitted_entry)
		drug_name = splitted_entry[0].strip()
		gene_name = splitted_entry[5].strip()
		alteration = splitted_entry[9].strip()
		association_implication = splitted_entry[10].strip()
		
		if '+' in drug_name:
			continue
		
		if '--' in gene_name or '+' in gene_name:
			continue
		
		#Could be changed to ignoring if desired, also MUTATION is not interpreted because the meaning is unclear
		if alteration == "ONCOGENIC MUTATION" or alteration == "Oncogenic mutations" or alteration == "LOSS-OF-FUNCTION" or alteration == "GAIN-OF-FUNCTION":
			
			onco_db_entries_certain_gene = get_elements_from_onco_db(gene_name, alteration)
			
			new_merged_entries = merge_sens_onco(entry, onco_db_entries_certain_gene)
			
			for element in new_merged_entries:
				
				alteration_merged_entry = element.strip('\n').split('\t')[9].strip()
				hash_ident = drug_name + gene_name + alteration_merged_entry
				
				if hash_ident in filtered_entries:
					
					duplicate_entry = filtered_entries[hash_ident]
					
					splitted_duplicate_entry = duplicate_entry.strip('\n').split('\t')
					
					association_implication_duplicate_entry = splitted_duplicate_entry[10].strip()
					
					which_case = compare_association(association_implication_duplicate_entry, association_implication)
					if which_case == 1 or which_case == 4:#a contradiction 
						 
						new_entry_splitted = duplicate_entry.strip('\n').split('\t')
						confidence_duplicate = new_entry_splitted[20].strip()
						
						my_new_confidence = ""
						if confidence_duplicate == "regular entry":
							my_new_confidence = "contradictory: regular entry vs inferential entry"	
						elif confidence_duplicate == "inferential entry":
							my_new_confidence = "contradictory: inferential entry vs inferential entry"
						elif confidence_duplicate == "contradictory: regular entry vs inferential entry":
							my_new_confidence = "contradictory: regular entry vs inferential entry"
						elif confidence_duplicate == "contradictory: inferential entry vs inferential entry":
							my_new_confidence = "contradictory: inferential entry vs inferential entry"
						elif confidence_duplicate == "supporting: regular entry and inferential entry":
							my_new_confidence = "contradictory: regular entry vs inferential entry"
						elif confidence_duplicate == "supporting: inferential entry and inferential entry":
							my_new_confidence = "contradictory: inferential entry vs inferential entry"
						else:
							my_new_confidence = confidence_duplicate
						
						new_entry_splitted[20] = my_new_confidence
						new_entry = "\t".join(new_entry_splitted)
						
						filtered_entries[hash_ident] = new_entry
					
					elif which_case == 2 or which_case == 3:#do not say it is a contradiction (although it slightly is)
						
						filtered_entries[hash_ident] = duplicate_entry
						
					else:#supporting entries available
						
						new_entry_splitted = duplicate_entry.strip('\n').split('\t')
						
						#print(len(new_entry_splitted))
						confidence_duplicate = new_entry_splitted[20].strip()
						my_new_confidence = ""
						
						if confidence_duplicate == "regular entry":
							my_new_confidence = "supporting: regular entry and inferential entry" 
						elif confidence_duplicate == "inferential entry":
							my_new_confidence = "supporting: inferential entry and inferential entry"
						else:
							my_new_confidence = confidence_duplicate
						new_entry_splitted[20] = my_new_confidence
						
						new_entry = "\t".join(new_entry_splitted)
						
						filtered_entries[hash_ident] = new_entry
				
				else:
					splitted_element = element.strip('\n').split('\t')
					splitted_element.append('inferential entry')
					filtered_entries[hash_ident] = '\t'.join(splitted_element)
		
		elif alteration == "MUTATION": #unclear what kind of mutation
			continue
		else:
			hash_ident = drug_name + gene_name + alteration
			
			if hash_ident in filtered_entries:
				
				duplicate_entry = filtered_entries[hash_ident]
				splitted_duplicate_entry = duplicate_entry.strip('\n').split('\t')
				association_implication_duplicate_entry = splitted_duplicate_entry[10].strip()
				confidence_duplicate = splitted_duplicate_entry[20].strip()
				my_new_confidence = ""
				
				which_case = compare_association(association_implication_duplicate_entry, association_implication)
				if which_case == 1 or which_case==4:#a contradiction
					
					if confidence_duplicate == "regular entry":
						my_new_confidence = "contradictory: regular entry vs regular entry"
						splitted_duplicate_entry[20] = my_new_confidence
						filtered_entries[hash_ident] = '\t'.join(splitted_duplicate_entry)
					elif confidence_duplicate == "inferential entry":
						my_new_confidence = "contradictory: regular entry vs inferential entry"
						splitted_entry.append(my_new_confidence)
						filtered_entries[hash_ident] = '\t'.join(splitted_entry)
					elif confidence_duplicate == "contradictory: inferential entry vs inferential entry":
						my_new_confidence = "contradictory: regular entry vs inferential entry"
						splitted_entry.append(my_new_confidence)
						filtered_entries[hash_ident] = '\t'.join(splitted_entry)
					elif confidence_duplicate == "contradictory: regular entry vs inferential entry":
						my_new_confidence = "contradictory: regular entry vs regular entry"
						splitted_duplicate_entry[20] = my_new_confidence
						filtered_entries[hash_ident] ='\t'.join( splitted_duplicate_entry)
					elif confidence_duplicate == "supporting: regular entry and inferential entry":
						my_new_confidence = "contradictory: regular entry vs regular entry"
						splitted_duplicate_entry[20] = my_new_confidence
						filtered_entries[hash_ident] = '\t'.join(splitted_duplicate_entry)
					elif confidence_duplicate == "supporting: inferential entry and inferential entry":
						my_new_confidence = "contradictory: regular entry vs inferential entry"
						splitted_entry.append(my_new_confidence)
						filtered_entries[hash_ident] = '\t'.join(splitted_entry)
						
					else:
						filtered_entries[hash_ident] ='\t'.join(splitted_duplicate_entry)
				elif which_case == 2:#do not say it is a contradiction, the duplicate already in list is the one without not 
					
					
					if confidence_duplicate == "inferential entry":
						my_new_confidence == "regular entry"
						splitted_entry.append(my_new_confidence)
						filtered_entries[hash_ident] = '\t'.join(splitted_entry)
					elif confidence_duplicate == "supporting: inferential entry and inferential entry":
						my_new_confidence == "regular entry"
						splitted_entry.append(my_new_confidence)
						filtered_entries[hash_ident] = '\t'.join(splitted_entry)
					elif confidence_duplicate == "contradictory: inferential entry vs inferential entry":
						my_new_confidence == "regular entry"
						splitted_entry.append(my_new_confidence)
						filtered_entries[hash_ident] = '\t'.join(splitted_entry)
						
					
				elif which_case == 3:#do not say it is a contradiction, the new sample is the one without not
					
					if confidence_duplicate == "inferential entry":
						my_new_confidence == "regular entry"
						splitted_entry.append(my_new_confidence)
						filtered_entries[hash_ident] = '\t'.join(splitted_entry)
					elif confidence_duplicate == "supporting: inferential entry and inferential entry":
						my_new_confidence == "regular entry"
						splitted_entry.append(my_new_confidence)
						filtered_entries[hash_ident] = '\t'.join(splitted_entry)
					elif confidence_duplicate == "contradictory: inferential entry vs inferential entry":
						my_new_confidence == "regular entry"
						splitted_entry.append(my_new_confidence)
						filtered_entries[hash_ident] = '\t'.join(splitted_entry)
					elif confidence_duplicate == "contradictory: regular entry vs inferential entry":
						my_new_confidence == "regular entry"
						splitted_entry.append(my_new_confidence)
						filtered_entries[hash_ident] 
					elif confidence_duplicate == "regular entry":
						my_new_confidence == "regular entry"
						splitted_entry.append(my_new_confidence)
						filtered_entries[hash_ident] = '\t'.join(splitted_entry)
	
						
				elif which_case == 50:#unknown association found
					
					print("Unknown association found")
					
				else:#no contradiction
					
					if confidence_duplicate == "inferential entry":
						my_new_confidence = "supporting: regular entry and inferential entry"
						splitted_entry.append(my_new_confidence)
						filtered_entries[hash_ident] = '\t'.join(splitted_entry)
					if confidence_duplicate == "supporting: inferential entry and inferential entry":
						my_new_confidence = "supporting: regular entry and inferential entry"
						splitted_entry.append(my_new_confidence)
						filtered_entries[hash_ident] = '\t'.join(splitted_entry)
			else:
				splitted_entry.append("regular entry")
				filtered_entries[hash_ident] = '\t'.join(splitted_entry)
					
				
	
	#filter out contradictory elements: contradictory:inf vs inf and contradictors:reg vs reg
	keys = list(filtered_entries)
	for key in keys: 
		
		splitted_entry = filtered_entries[key].strip('\n').split('\t')
		
		if splitted_entry[20] == "contradictory: inferential entry vs inferential entry":
			filtered_entries.pop(key)
		elif splitted_entry[20] == "contradictory: regular entry vs regular entry":
			filtered_entries.pop(key)
		
		
	if output_file == "no_output_file":
		return filtered_entries
	
	else:
		with open(output_file, 'w', encoding='utf-8') as relevant_output:
			relevant_output.write("drug_name" + "\t" + "DrugBank_ID" + "\t" + "drug_family" + "\t" + "drug_targets" + "\t" + "mutation_effect" + "\t" + "associated_gene" + "\t" + "HGNC_ID_associated_gene" + "\t" + "pathways_associated_gene" + "\t" + "alteration_type" + "\t"+ "alteration" + "\t" + "association_implication" + "\t" + "evidence" + "\t" + "source" + "\t" + "primary_tumour_type" + "\t" + "tsg_onco_annotation" + "\t" + "source_database" + "\t" + "evidence_statement" + "\t" + "cds_like_mutation" + "\t" + "genome_coordinates" + "\t" + "reference_build" + "\n")
			
			for element in filtered_entries:
				relevant_output.write(filtered_entries[element] + "\n")
				
				
def query_gene_db(query_gene_name, output_file):
	
	official_gene_symbols = load_gene_symbols()
	alias_to_symbol = normalized_alias_map()


	db_list =[output1, output2, output3, output4]
	
	relevant_entries = []
	for db in db_list:
		
		with open(db, 'r', encoding='utf-8') as current_db:
			
			current_db.readline()
			
			for line in current_db:
				
				sline = line.split("\t")
				
				gene_name = sline[5]
				
				if '--' in gene_name and '+' in gene_name:
					print('Cannot handle a case with -- and + in gene name')
					continue
				
				if '--' in gene_name:
				
					several_gene_names = gene_name.split("--")
					
					for gn in several_gene_names:
						
						if normalize_caseless(query_gene_name.strip()) == normalize_caseless(gn.strip()):
							
							if sline[0].strip() == "":
								new_entry = sline[5] + "\t" + sline[6] + "\t" + sline[9] + "\t" + sline[8] + "\t" + sline[4] + "\t" + sline[10] + "\t" + sline[13] + "\t" + sline[14] + "\t" + sline[12] + "\t" + sline[15] + "\t" + sline[11] + "\t" + sline[16] +"\t" +  sline[17] +"\t" +  sline[18] +"\t" +  sline[19] + "\t" + sline[2] 
							else:
								new_entry = sline[5] + "\t" + sline[6] + "\t" + sline[9] + "\t" + sline[8] + "\t" + sline[4] + "\t" + sline[10] + "\t" + sline[13] + "\t" + sline[14] + "\t" + sline[12] + "\t" + sline[15] + "\t" + sline[11] + "\t" + sline[16] +"\t" +  sline[17] +"\t" +  sline[18] +"\t" +  sline[19] + "\t" + sline[0]
							
							relevant_entries.append(new_entry)
							
						elif normalize_caseless(query_gene_name.strip()) in alias_to_symbol:
							
							norm_common_name_of_query = alias_to_symbol[normalize_caseless(query_gene_name.strip())]
							
							
							if not len(norm_common_name_of_query) == 1:
								print('Could not find a unique identifier for the query gene name. Name: ' + query_gene_name)
							else:
								if norm_common_name_of_query[0].strip() == normalize_caseless(gn.strip()):
									
									
									if sline[0].strip()=="":
										new_entry = sline[5] + "\t" + sline[6] + "\t" + sline[9] + "\t" + sline[8] + "\t" + sline[4] + "\t" + sline[10] + "\t" + sline[13] + "\t" + sline[14] + "\t" + sline[12] + "\t" + sline[15] + "\t" + sline[11] + "\t" + sline[16] +"\t" +  sline[17] +"\t" +  sline[18] +"\t" +  sline[19] + "\t" + sline[2]
									else:
										new_entry = sline[5] + "\t" + sline[6] + "\t" + sline[9] + "\t" + sline[8] + "\t" + sline[4] + "\t" + sline[10] + "\t" + sline[13] + "\t" + sline[14] + "\t" + sline[12] + "\t" + sline[15] + "\t" + sline[11] + "\t" + sline[16] +"\t" +  sline[17] +"\t" +  sline[18] +"\t" +  sline[19] + "\t" + sline[0]
									relevant_entries.append(new_entry)
				
				elif '+' in gene_name:
					several_gene_names = gene_name.split("+")
					
					for gn in several_gene_names:
						
						if normalize_caseless(query_gene_name.strip()) == normalize_caseless(gn.strip()):
							
							if sline[0].strip() == "":
								new_entry = sline[5] + "\t" + sline[6] + "\t" + sline[9] + "\t" + sline[8] + "\t" + sline[4] + "\t" + sline[10] + "\t" + sline[13] + "\t" + sline[14] + "\t" + sline[12] + "\t" + sline[15] + "\t" + sline[11] + "\t" + sline[16] +"\t" +  sline[17] +"\t" +  sline[18] +"\t" +  sline[19] + "\t" + sline[2]
							else:
								new_entry = sline[5] + "\t" + sline[6] + "\t" + sline[9] + "\t" + sline[8] + "\t" + sline[4] + "\t" + sline[10] + "\t" + sline[13] + "\t" + sline[14] + "\t" + sline[12] + "\t" + sline[15] + "\t" + sline[11] + "\t" + sline[16] +"\t" +  sline[17] +"\t" +  sline[18] +"\t" +  sline[19] + "\t" + sline[0]
							relevant_entries.append(new_entry)
							
						elif normalize_caseless(query_gene_name.strip()) in alias_to_symbol:
							
							norm_common_name_of_query = alias_to_symbol[normalize_caseless(query_gene_name.strip())]
							
							
							if not len(norm_common_name_of_query) == 1:
								print('Could not find a unique identifier for the query gene name. Name: ' + query_gene_name)
							else:
								if norm_common_name_of_query[0].strip() == normalize_caseless(gn.strip()):
									
									if sline[0].strip() == "":
										new_entry = sline[5] + "\t" + sline[6] + "\t" + sline[9] + "\t" + sline[8] + "\t" + sline[4] + "\t" + sline[10] + "\t" + sline[13] + "\t" + sline[14] + "\t" + sline[12] + "\t" + sline[15] + "\t" + sline[11] + "\t" + sline[16] +"\t" +  sline[17] +"\t" +  sline[18] +"\t" +  sline[19] + "\t" + sline[2]
									else:
										new_entry = sline[5] + "\t" + sline[6] + "\t" + sline[9] + "\t" + sline[8] + "\t" + sline[4] + "\t" + sline[10] + "\t" + sline[13] + "\t" + sline[14] + "\t" + sline[12] + "\t" + sline[15] + "\t" + sline[11] + "\t" + sline[16] +"\t" +  sline[17] +"\t" +  sline[18] +"\t" +  sline[19] + "\t" + sline[0]
									relevant_entries.append(new_entry)
				else:
					several_gene_names = gene_name.split("--")
					
					for gn in several_gene_names:
						
						if normalize_caseless(query_gene_name.strip()) == normalize_caseless(gn.strip()):
							
							if sline[0].strip() == "":
								new_entry = sline[5] + "\t" + sline[6] + "\t" + sline[9] + "\t" + sline[8] + "\t" + sline[4] + "\t" + sline[10] + "\t" + sline[13] + "\t" + sline[14] + "\t" + sline[12] + "\t" + sline[15] + "\t" + sline[11] + "\t" + sline[16] +"\t" +  sline[17] +"\t" +  sline[18] +"\t" +  sline[19] + "\t" + sline[2]
							else:
								new_entry = sline[5] + "\t" + sline[6] + "\t" + sline[9] + "\t" + sline[8] + "\t" + sline[4] + "\t" + sline[10] + "\t" + sline[13] + "\t" + sline[14] + "\t" + sline[12] + "\t" + sline[15] + "\t" + sline[11] + "\t" + sline[16] +"\t" +  sline[17] +"\t" +  sline[18] +"\t" +  sline[19] + "\t" + sline[0]
							relevant_entries.append(new_entry)
							
						elif normalize_caseless(query_gene_name.strip()) in alias_to_symbol:
							
							norm_common_name_of_query = alias_to_symbol[normalize_caseless(query_gene_name.strip())]
							
							
							if not len(norm_common_name_of_query) == 1:
								print('Could not find a unique identifier for the query gene name. Name: ' + query_gene_name)
							else:
								if norm_common_name_of_query[0].strip() == normalize_caseless(gn.strip()):
									
									if sline[0].strip() == "":
										new_entry = sline[5] + "\t" + sline[6] + "\t" + sline[9] + "\t" + sline[8] + "\t" + sline[4] + "\t" + sline[10] + "\t" + sline[13] + "\t" + sline[14] + "\t" + sline[12] + "\t" + sline[15] + "\t" + sline[11] + "\t" + sline[16] +"\t" +  sline[17] +"\t" +  sline[18] +"\t" +  sline[19] + "\t" + sline[2]
									else:
										new_entry = sline[5] + "\t" + sline[6] + "\t" + sline[9] + "\t" + sline[8] + "\t" + sline[4] + "\t" + sline[10] + "\t" + sline[13] + "\t" + sline[14] + "\t" + sline[12] + "\t" + sline[15] + "\t" + sline[11] + "\t" + sline[16] +"\t" +  sline[17] +"\t" +  sline[18] +"\t" +  sline[19] + "\t" + sline[0]
									relevant_entries.append(new_entry)
	
	if output_file == "no_output_file":
		
		return relevant_entries
	
	else:
		with open(output_file, 'w', encoding='utf-8') as relevant_output:
			
			relevant_output.write("gene_name" + "\t" + "HGNC_ID"  + "\t" +"alteration" + "\t" + "alteration_type" + "\t"  + "mutation_effect" + "\t" + "association_implication" + "\t" + "primary_tumour_type" + "\t" + "tsg_onco_annotation" + "\t" + "source" + "\t" +  "source_database" + "\t" + "evidence_level" + "\t" + "evidence_statement" + "\t" + "cds_like_mutation" + "\t" + "genome_coordinates" + "\t" + "reference_build" + "\t" + "drug"+ "\n")
			
			for element in relevant_entries:
				relevant_output.write(element)


################Main function ##################################

def main(mode, name, output_file):
	
	if mode == "drug":
		query_drug_db(name, output_file)
		return
	elif mode == "gene":
		query_gene_db(name, output_file)
		return
	elif mode == "drug_filtered":
		query_drug_filtered(name, output_file)
	else:
		print("Unknown mode for drug sensitivity biomarker database, supported modes are: 'drug', 'drug_filtered', and 'gene'")
		return


#################################################################
if __name__ == "__main__":
	if len(sys.argv) < 4:
			sys.exit("This program needs the following arguments:\
					\n- 'drug', 'drug_filtered' or 'gene' depending on whether one wants to request a drug or gene name\
					\n- drug name\
					\n- output file name (including path)") 

	main(sys.argv[1], sys.argv[2], sys.argv[3]) 	
	
