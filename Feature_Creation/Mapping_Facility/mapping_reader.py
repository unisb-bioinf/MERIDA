#Author: Kerstin Lenhof
#Date: 31.10.2018
#Map gene names to official gene symbols, drug to drugs to common names....
from useful_tools import * 

official_gene_symbol_file = "<path_to_official_gene_symbol_file>"
alias_file = "<path_to_gene_alias_to_official_gene_symbol_file>"
drugbank_file = "<path_to_drugbank_vocabulary_file>"
hgnc_file = "<path_to_hgnc_symbol_to_hgnc_id_file>"
og_tcg_file_paper = "<path_to_oncogenic_signaling_pathways_in_tcga_file_with_TSG_and_OG_annotations>"


def load_gene_symbols():
	official_gene_symbols = []
	with open(official_gene_symbol_file, 'r', encoding='utf-8') as official_gene_symbol:
		
		for line in official_gene_symbol:
			if line.startswith("#"):
				continue
			else:
				if line.strip() != "":
					official_gene_symbols.append(line.strip())
	return official_gene_symbols

def load_alias_map():
	alias_to_gene_dict = {}
	with open(alias_file, 'r', encoding='utf-8') as alias:
		
		for line in alias:
			if line.startswith('#'):
				continue
			else:
				alias_symbol = line.split('\t')
				
				if alias_symbol[0].strip() in alias_to_gene_dict:
					alias_to_gene_dict[alias_symbol[0].strip()].append(alias_symbol[1].strip())
					#print("Alias found twice, alias:" + alias_symbol[0])
				else:
					if alias_symbol[0].strip() != "":
						alias_to_gene_dict[alias_symbol[0].strip()] = [alias_symbol[1].strip()]
	
	return alias_to_gene_dict


def normalized_alias_map():
	alias_to_gene_dict = {}
	with open(alias_file, 'r', encoding='utf-8') as alias:
		
		for line in alias:
			if line.startswith('#'):
				continue
			else:
				alias_symbol = line.split('\t')
				
				if normalize_caseless(alias_symbol[0].strip()) in alias_to_gene_dict:
					alias_to_gene_dict[normalize_caseless(alias_symbol[0].strip())].append(alias_symbol[1].strip())
					#print("Alias found twice, alias:" + alias_symbol[0])
				else:
					if alias_symbol[0].strip() != "":
						alias_to_gene_dict[normalize_caseless(alias_symbol[0].strip())] = [alias_symbol[1].strip()]
	
	return alias_to_gene_dict


def load_drugcommonName_2_id_map():
	
	name_2_id_map ={}
	with open(drugbank_file, 'r', encoding='utf-8') as drug_bank:
		
		drug_bank.readline()
		
		for line in drug_bank:
			
			if "\"" in line:
				sline = line.split('\"')
				
				left_sline = sline[0].split(',')
				
				#mid_sline = sline[1].split("|")
				
				drug_bank_id = left_sline[0].strip()
				common_name = left_sline[2].strip()
				
				if common_name == "":
					continue
				if common_name in name_2_id_map:
					print('Double occurrence of ' + common_name)
				else:
					name_2_id_map[common_name] = drug_bank_id
				
			else:
				sline = line.split(',')
				
				drug_bank_id, common_name, synonyms = sline[0].strip(), sline[2].strip(), sline[5].strip()
				
				if common_name == "":
					continue
				if common_name in name_2_id_map:
					print('Double occurrence of ' + common_name)
				else:
					name_2_id_map[common_name] = drug_bank_id
				
						
	return name_2_id_map

def load_drugnormalizedname_2_id_map():
	name_2_id_map ={}
	with open(drugbank_file, 'r', encoding='utf-8') as drug_bank:
		
		drug_bank.readline()
		
		for line in drug_bank:
			
			if "\"" in line:
				sline = line.split('\"')
				
				left_sline = sline[0].split(',')
				
				#mid_sline = sline[1].split("|")
				
				if left_sline[2].strip() == "":
					continue
				drug_bank_id = left_sline[0].strip()
				normalized_common_name = normalize_caseless(left_sline[2].strip())
				
				if normalized_common_name in name_2_id_map:
					print('Double occurrence of ' + normalized_common_name)
				else:
					name_2_id_map[normalized_common_name] = drug_bank_id
			else:
				sline = line.split(',')
				
				drug_bank_id, common_name, synonyms = sline[0].strip(), sline[2].strip(), sline[5].strip()
				
				if common_name == "":
					continue
				
				normalized_common_name = normalize_caseless(common_name)
				if normalized_common_name in name_2_id_map:
					print('Double occurrence of ' + common_name)
				else:
					name_2_id_map[normalized_common_name] = drug_bank_id
				
						
	return name_2_id_map
	
def load_drugnormalized_2_commonname_map():
	name_normalized_map = {}
	
	with open(drugbank_file, 'r', encoding='utf-8') as drug_bank:
		
		drug_bank.readline()
		
		for line in drug_bank:
			
			if "\"" in line:
				sline = line.split('\"')
				
				left_sline = sline[0].split(',')
				
				#mid_sline = sline[1].split("|")
				
				#drug_bank_id = left_sline[0]
				common_name = left_sline[2].strip()
				
				if common_name == "":
					continue
				
				normalized_common_name = normalize_caseless(common_name)
				
				if normalized_common_name in name_normalized_map:
					print('Double occurrence of ' + normalized_common_name)
				else:
					name_normalized_map[normalized_common_name] = common_name
			
			else:
				sline = line.split(',')
				
				common_name =  sline[2].strip()
				
				if common_name == "":
					continue
				
				normalized_name = normalize_caseless(common_name)
				if normalized_name in name_normalized_map:
					print('Double occurrence of ' + normalized_name)
				else:
					name_normalized_map[normalized_name] = common_name
				
						
	return name_normalized_map

def load_drugsyn_2_common_map():
	
	synonym_2_common_map = {}
	with open(drugbank_file, 'r', encoding='utf-8') as drug_bank:
		
		drug_bank.readline()
		
		for line in drug_bank:
			
			if "\"" in line:
				sline = line.split('\"')
				
				left_sline = sline[0].split(',')
				
				common_name = left_sline[2].strip()
				
				if common_name == "":
					continue
				synonyms = sline[1].strip()
				
				if synonyms == "":
					continue
				
				else:
					
					split_synonyms = synonyms.split("|")
				
					for syn in split_synonyms:
						
						if syn.strip() in synonym_2_common_map:
							print('Found the same synonym twice: ' + syn)
							
						else:
							synonym_2_common_map[syn.strip()] = common_name
					
			else:
				sline = line.split(',')
				
				drug_bank_id, common_name, synonyms = sline[0].strip(), sline[2].strip(), sline[5].strip()
				
				if common_name == "":
					continue
					
				if synonyms == "":
					continue
				else:
					split_synonyms = synonyms.split('|')
					
					for syn in split_synonyms:
						
						if syn.strip() in synonym_2_common_map:
							print('Found the same synonym twice: ' + syn)
							
						else:
							synonym_2_common_map[syn.strip()] = common_name
						
	return synonym_2_common_map



def load_drugsynnormalized_2_commonnormalized_map():
	
	synonym_2_common_map = {}
	with open(drugbank_file, 'r', encoding='utf-8') as drug_bank:
		
		drug_bank.readline()
		
		for line in drug_bank:
			
			if "\"" in line:
				sline = line.split('\"')
				
				left_sline = sline[0].split(',')
				
				common_name = left_sline[2].strip()
				
				if common_name == "":
					continue
				synonyms = sline[1].strip()
				
				if synonyms == "":
					continue
				
				else:
					
					split_synonyms = synonyms.split("|")
				
					for syn in split_synonyms:
						
						if normalize_caseless(syn.strip()) in synonym_2_common_map:
							print('Found the same synonym twice: ' + syn)
							
						else:
							synonym_2_common_map[normalize_caseless(syn.strip())] = normalize_caseless(common_name)
			
			else:
				sline = line.split(',')
				
				drug_bank_id, common_name, synonyms = sline[0].strip(), sline[2].strip(), sline[5].strip()
				
				if common_name == "":
					continue
					
				if synonyms == "":
					continue
				else:
					split_synonyms = synonyms.split('|')
					
					for syn in split_synonyms:
						
						if normalize_caseless(syn.strip()) in synonym_2_common_map:
							print('Found the same synonym twice: ' + syn)
							
						else:
							synonym_2_common_map[normalize_caseless(syn.strip())] = normalize_caseless(common_name)
						
	return synonym_2_common_map		
	
def load_HGNCsymbol_2_id_map():
	
	symbol_2_id_map = {}
	
	with open(hgnc_file, 'r', encoding='utf-8') as hgnc_id:
		
		hgnc_id.readline()
		
		for line in hgnc_id:
			sline = line.split('\t')
			
			hgnc_id, symbol = sline[0].strip(), sline[1].strip()
			
			if symbol == "":
				continue
			
			if hgnc_id == "":
				continue
			
			if symbol in symbol_2_id_map:
				print("Duplicate gene found " + symbol)
			else:
				symbol_2_id_map[symbol] = hgnc_id
	
	return symbol_2_id_map

def load_og_tsg_annotation_from_paper():
	
	gene_to_og_tsg = {}
	with open(og_tcg_file_paper, 'r', encoding='utf-8') as og_tsg:
		
		for line in og_tsg:
			
			sline = line.split('\t')
			
			gene, assignment = sline[0], sline[1]
			if assignment.strip() == "" or assignment.strip() == "Unknown" or assignment.strip() =="NA" or gene.strip() == "":
				continue
			else:
				
				if gene.strip() in gene_to_og_tsg:
					print("Gene found twice")
					
					if gene_to_og_tsg[gene.strip()] == assignment.strip():
						print("Same annotation for gene " + gene + "\n")
						
					else:
						print("Different annotation for gene " + gene + "\n")
						
				else:
					gene_to_og_tsg[gene.strip()] = assignment.strip()
				
	return gene_to_og_tsg			
				
		
		
		
		
		
		
		
		
