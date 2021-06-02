from mapping_reader import *
from useful_tools import *

cgi_sens_bio = "<path_to_standardized_cgi_file>"
output1 = "<output_file_for_cgi_file>"

civic_sens_bio = "<path_to_standardized_civic_file>"
output2 = "<output_file_for_civic_file>"

oncokb_sens_bio = "<path_to_standardized_oncokb_file>"
output3 = "<output_file_for_oncokb_file>"

cosmic_sens_bio =  "<path_to_standardized_cosmic_file>"
output4 = "<output_file_for_cosmic_file>"

db_list_bio = [cgi_sens_bio, civic_sens_bio, oncokb_sens_bio, cosmic_sens_bio]
db_list_bio_output = [output1, output2, output3, output4]


gene_symbols = load_gene_symbols()
alias_to_symbol = load_alias_map()
symbol_to_id = load_HGNCsymbol_2_id_map()
drugcommon_to_drugbank_id = load_drugnormalizedname_2_id_map()
normalized_to_name = load_drugnormalized_2_commonname_map()
normalized_syn_to_common = load_drugsynnormalized_2_commonnormalized_map()
gene_to_og_tsg_paper = load_og_tsg_annotation_from_paper()



def map_gene(gene_name):
	gene_name = gene_name.strip()
	if not gene_name in gene_symbols:
			#print('Wrong name')
			
			if not gene_name in alias_to_symbol:
				print("Gene name is not annotated: " + gene_name)
				return gene_name
			else:
				if not len(alias_to_symbol[gene_name])== 1:
					print(alias_to_symbol[gene_name])
					print("Several genes could be meant, checking required: " + gene_name)
					#new_name = ''.join(gene_name)
					return gene_name
				else:
					new_name = str(alias_to_symbol[gene_name][0]).strip()
					return new_name
	else:
		return gene_name
	
counter = 0
for db_file in db_list_bio:
	
	db_entries = []
	with open(db_file, 'r', encoding='utf-8') as db:
		
		print(db_file)
		db.readline() #header
		
		for line in db:
			
			sline = line.split("\t")
			#print(sline)
			gene_name = sline[5].strip()
			new_name = gene_name
			
			if '--' in gene_name and ';' in gene_name:
				print('Cannot handle a case with -- and ; in gene name')
				continue
			#fusions
			if '--' in gene_name: 
				s_gene_name = gene_name.split('--')
				first_gene = s_gene_name[0].strip()
				second_gene = s_gene_name[1].strip()
				
				mapped_gene1 = map_gene(first_gene)
				mapped_gene2 = map_gene(second_gene)
				
				hgnc_1 = ""
				if mapped_gene1.strip() in symbol_to_id:
					hgnc_1 = symbol_to_id[mapped_gene1.strip()]
					
				hgnc_2 = ""
				if mapped_gene2.strip() in symbol_to_id:
					hgnc_2 = symbol_to_id[mapped_gene2.strip()]
					
				sline[6] = ('--'.join([hgnc_1, hgnc_2])).strip()
				
				new_name = str(mapped_gene1).strip() + "--" + str(mapped_gene2).strip()
				#print("'" + new_name + "'")
			
			elif ';' in gene_name: #this is not good, should ask this first and then check for --
				s_gene_name = gene_name.split(';')
				first_gene = s_gene_name[0].strip()
				second_gene = s_gene_name[1].strip()
				
				mapped_gene1 = map_gene(first_gene)
				mapped_gene2 = map_gene(second_gene)
				
				hgnc_1 = ""
				if mapped_gene1.strip() in symbol_to_id:
					hgnc_1 = symbol_to_id[mapped_gene1.strip()]
					
				hgnc_2 = ""
				if mapped_gene2.strip() in symbol_to_id:
					hgnc_2 = symbol_to_id[mapped_gene2.strip()]
					
				sline[6] = ('+'.join([hgnc_1, hgnc_2])).strip()
				
				new_name = str(mapped_gene1).strip() + "+" + str(mapped_gene2).strip()
				#print("'" + new_name + "'")
			else:
				new_name = map_gene(gene_name.strip())
				
				if new_name.strip() in symbol_to_id:
				
					sline[6] = symbol_to_id[new_name.strip()]
				
				else:
					print("Could not find hgnc id for " + new_name.strip())
			new_name = new_name.strip()
			sline[5] = new_name
	
			#print("'" + sline[5] + "'")
			#print("'" + sline[6] + "'")
			
			found_gene = False
			for key in gene_to_og_tsg_paper:
				
				if normalize_caseless(key.strip()) == normalize_caseless(sline[5].strip()):
					
					sline[14] = gene_to_og_tsg_paper[key]
					found_gene = True
			
			if not found_gene:
				print("OG/TSG assignment not available for " + sline[5])
				
			drug_names = sline[0]
			
			#Several drugs administered at once

			s_drug_names = drug_names.split('+')
			common_names_list = []
			db_ids = []
			
			for dr in s_drug_names:
				
				if dr.strip() == "":
					continue
				dr_in_brackets = ""
				if '(' in dr and ')' in dr:
					dr_split_bracket = dr.split('(')
					dr_in_brackets = dr_split_bracket[1].replace(')','').strip()
					dr = dr_split_bracket[0].strip()
					
				normalized_dr_name = normalize_caseless(dr.strip())
				if normalized_dr_name in drugcommon_to_drugbank_id:
					common_names_list.append(normalized_to_name[normalized_dr_name])
					db_ids.append(drugcommon_to_drugbank_id[normalized_dr_name])
				else:
					
					if normalized_dr_name in normalized_syn_to_common:
						
						norm_common = normalized_syn_to_common[normalized_dr_name]
						
						common_names_list.append(normalized_to_name[norm_common])
						db_ids.append(drugcommon_to_drugbank_id[norm_common])
						
					
					else:
						normalized_dr_in_brackets = normalize_caseless(dr_in_brackets)
						if normalized_dr_in_brackets in drugcommon_to_drugbank_id:
							common_names_list.append(normalized_to_name[normalized_dr_in_brackets])
							db_ids.append(drugcommon_to_drugbank_id[normalized_dr_in_brackets])
						else:
				
							if normalized_dr_in_brackets in normalized_syn_to_common:
					
								norm_common = normalized_syn_to_common[normalized_dr_in_brackets]
					
								common_names_list.append(normalized_to_name[norm_common])
								db_ids.append(drugcommon_to_drugbank_id[norm_common])
							
							else:
								print("Could not find the drug in DrugBank: " + dr)
								common_names_list.append(dr)
								db_ids.append(' ')
					
						
			
			sline[0] = "+".join(common_names_list)
			sline[1] = "+".join(db_ids)
				
			
			new_entry = sline[0] + "\t" + sline[1] + "\t" + sline[2] + "\t" + sline[3] + "\t" + sline[4] + "\t" + sline[5] + "\t" + sline[6] + "\t" + sline[7] + "\t" + sline[8] + "\t" + sline[9] + "\t" + sline[10] + "\t" + sline[11] + "\t" + sline[12] + "\t" + sline[13] + "\t" + sline[14] + "\t" + sline[15] + "\t" + sline[16] + "\t" + sline[17] + "\t" + sline[18] + "\t" + sline[19]
			db_entries.append(new_entry)
	
	
	with open(db_list_bio_output[counter], 'w', encoding='utf-8') as output:
		
		output.write("drug_name" + "\t" + "DrugBank_ID" + "\t" + "drug_family" + "\t" + "drug_targets" + "\t" + "mutation_effect" + "\t" + "associated_gene" + "\t" + "HGNC_ID_associated_gene" + "\t" + "pathways_associated_gene" + "\t" + "alteration_type" + "\t"+ "alteration" + "\t" + "association_implication" + "\t" + "evidence" + "\t" + "source" + "\t" + "primary_tumour_type" + "\t" + "tsg_onco_annotation" + "\t" + "source_database" + "\t" + "evidence_statement" + "\t" + "cds_like_mutation" + "\t" + "genome_coordinates" + "\t" + "reference_build" + "\n")
		for element in db_entries:
			output.write(element)
	
	counter = counter +1 
	
#for key in normalized_syn_to_common:
	#print("Alias:" + key + ", Name:" +  normalized_syn_to_common[key])
	#print(normalized_syn_to_common[element])
