	
from mapping_reader import *

cgi_onco = "<path_to_standardized_cgi_file>"
output_onco1 =  "<output_file_for_cgi_file>"

civic_onco = "<path_to_standardized_civic_file>"
output_onco2 = "<output_file_for_civic_file>"

oncokb_onco = "<path_to_standardized_oncokb_file>"
output_onco3 = "<output_file_for_oncokb_file>"

db_list_onco = [cgi_onco, civic_onco, oncokb_onco]
db_list_onco_output = [output_onco1, output_onco2, output_onco3]

#common HGNC symbol + ID mapping
gene_symbols = load_gene_symbols()
alias_to_symbol = load_alias_map()
symbol_to_id = load_HGNCsymbol_2_id_map()
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
					return new_name.strip()
	else:
		return gene_name


counter2 = 0
for db_file in db_list_onco:
	db_entries = []
	with open(db_file, 'r', encoding='utf-8') as db:
		
		print(db_file)
		db.readline() #header
		
		for line in db:
			
			sline = line.split("\t")
			#print(sline)
			gene_name = sline[0].strip()
			new_name = gene_name
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
					
				sline[1] = '--'.join([hgnc_1, hgnc_2])
				
				new_name = str(mapped_gene1).strip() + "--" + str(mapped_gene2).strip()
			
			else:
				new_name = map_gene(gene_name)
				
				if new_name.strip() in symbol_to_id:
				
					sline[1] = symbol_to_id[new_name.strip()]
				
				else:
					print("Could not find hgnc id for " + new_name.strip())
			
			new_name = new_name.strip()
			sline[0] = new_name
			
			found_gene = False
			for key in gene_to_og_tsg_paper:
				
				if normalize_caseless(key.strip()) == normalize_caseless(sline[0].strip()):
					
					sline[7] = gene_to_og_tsg_paper[key]
					found_gene = True
			
			if not found_gene:
				print("OG/TSG assignment not available for " + sline[0])
			
			
			new_entry = sline[0] + "\t" + sline[1] + "\t" + sline[2] + "\t" + sline[3] + "\t" + sline[4] + "\t" + sline[5] + "\t" + sline[6] + "\t" + sline[7] + "\t" + sline[8] + "\t" + sline[9] + "\t" + sline[10] + "\t" + sline[11] + "\t" + sline[12] + "\t" + sline[13] + "\t" + sline[14]
			db_entries.append(new_entry)
	
	
	with open(db_list_onco_output[counter2], 'w', encoding='utf-8') as output:
		
		output.write("gene_name" + "\t" + "HGNC_ID"  + "\t" +"alteration" + "\t" + "alteration_type" + "\t" +  "mutation_effect" + "\t" + "association_implication" + "\t" + "primary_tumour_type" + "\t" + "tsg_onco_annotation" + "\t" + "source" + "\t" +  "source_database" + "\t" + "evidence_level" + "\t" + "evidence_statement" + "\t" + "cds_like_mutation" + "\t" + "genome_coordinates" + "\t" + "reference_build" + "\n")
		for element in db_entries:
			output.write(element)
	
	counter2 = counter2 +1 
