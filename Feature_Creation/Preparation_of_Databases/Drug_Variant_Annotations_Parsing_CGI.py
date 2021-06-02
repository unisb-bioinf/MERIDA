#Author: Kerstin Lenhof
#Date: 19/10/2018
#Description: Databases with information on drug-variant pairs are parsed to obtain one unified database. 

def extractAlteration(item):
	
	first_whitespace_split = item.strip().split(" ")
	#print(first_whitespace_split)
	#print(first_whitespace_split[1])
	if "oncogenic" == first_whitespace_split[1]:
		return "ONCOGENIC MUTATION"
	elif "overexpression" == first_whitespace_split[1]:
		return "OVEREXPRESSION"
	elif "fusion" == first_whitespace_split[1]:
		return "FUSION"
	elif "deletion" == first_whitespace_split[1]:
		return "DELETION"
	elif "amplification" == first_whitespace_split[1]:
		return "AMPLIFICATION"
	elif "biallelic" == first_whitespace_split[1]:
		return "BIALLELIC INACTIVATION"
	elif "promoters" == first_whitespace_split[1]:
		return "PROMOTERS CORE"
	elif "wildtype" == first_whitespace_split[1]:
		return "WILDTYPE"
	elif "underexpression" == first_whitespace_split[1]:
		return "UNDEREXPRESSION"
	elif "expression" == first_whitespace_split[1]:
		return "NORMAL EXPRESSION"
	elif "frameshift" == first_whitespace_split[1]:
		return "FRAMESHIFT"
	elif "splice" == first_whitespace_split[1]:
		if first_whitespace_split[2] == "donor":
			return "SPLICE DONOR VARIANT"
		elif first_whitespace_split[2] == "acceptor":
			return "SPLICE ACCEPTOR VARIANT"
		else:
			print("Unknown biomarker type found:" + item)
			return "UNKNOWN BIOMARKER TYPE"
	elif "proximal" == first_whitespace_split[1]:
		return "CHECKING REQUIRED 1"
	elif "exon" == first_whitespace_split[1]:
		return "CHECKING REQUIRED 2"
	elif "mutation" == first_whitespace_split[1]:
		return "CHECKING REQUIRED 3"
	elif first_whitespace_split[1].startswith("("):
		return "MUTATION"
	elif "inframe" == first_whitespace_split[1]:
		return "CHECKING REQUIRED 4"
	else:
		print("Unknown biomarker type found:" + item)
		return "UNKNOWN BIOMARKER TYPE"


print("Parsing Cancer Genome Interpreter database")
cancer_genome_interpreter_db = "<path_to_cgi_biomarkers>"

cgi_db_entries = []
with open(cancer_genome_interpreter_db, 'r', encoding='utf-8') as cgi_db:
	
	cgi_db.readline()#header
	
	for line in cgi_db:
		sline = line.strip('\n').split('\t')
		if len(sline) < 20:
			print("Line with incorrect length found, length: " + str(len(sline)) + "\n" + "Line: " + line)
			continue
			 
		biomarker, gene_name, alteration_type, alteration, drug_family, drug, association, evidence_level, source, primary_tumor_type = sline[0].strip(), sline[1].strip(), sline[2].strip(), sline[3].strip(), sline[6].strip(), sline[7].strip(), sline[8].strip(), sline[9].strip(), sline[11].strip(), sline[19].strip()
		
		
		if "+" in biomarker and "inframe deletion" in biomarker: 
			print("+ and inframe deletion")
			print(biomarker)
		if association == "Responsive":
			association = "sensitive"
		elif association == "Resistant":
			association = "resistant"
		elif association == "No Responsive":
			association = "not sensitive"
		else:
			print("Uninterpretable association found, association: " + association)
		
		if ";" in drug:
			drug = drug.replace(';', '+')		
		if "[]" == drug:
			drug = " "
		drug = drug.replace('[', '').replace(']', '')
			
		drug_family = drug_family.replace('[', '').replace(']', '')
		drug_family = drug_family.replace(';', '+')
			
		#are there several alterations that have to co-occur 
		split_plus = biomarker.split('+')
		#print(split_plus)
		split_semi = alteration.split(';')
		
		if len(split_plus) > 1:
			
			alt_types = []
			alts = []
			had_to_break = False
			for count_alteration, element in enumerate(split_plus):


				if ',' in element:
					print("Currently parsing such complex expressions is not possible")
					print(biomarker)
					had_to_break = True
					break
				
				extracted_alteration = extractAlteration(element)
				
				if extracted_alteration == "CHECKING REQUIRED 3":
					conseq_alt = split_semi[count_alteration].split(':')[1]
					alts.append(conseq_alt)
					alt_types.append("MUTATION")
				elif extracted_alteration == "CHECKING REQUIRED 2" or extracted_alteration == "CHECKING REQUIRED 1":
					conseq_alt = split_semi[count_alteration].split('consequence::')
					conseq_alt_splitted = conseq_alt[1].split(':')
					
					if conseq_alt_splitted[0] == "inframe_insertion":
						alt_types.append("INFRAME INSERTION")
					elif conseq_alt_splitted[0] =="inframe_deletion":
						alt_types.append("INFRAME DELETION")
					else:
						print("Unknown biomarker type found: " + element)
						
					alts.append(conseq_alt_splitted[1])
				elif extracted_alteration == "FRAMESHIFT":
					
					whitespace = element.split(' ',1)[1]
					
					rm_everything = whitespace.replace('frameshift variant', '')
					rm_everything = rm_everything.replace('(', '')
					rm_everything = rm_everything.replace(')', '')
					
					alts.append(rm_everything)
					alt_types.append(extracted_alteration)
				elif extracted_alteration == "CHECKING REQUIRED 4":
					
					whitespace = element.split(' ',1)[1]
					if whitespace.startswith('inframe deletion'):
						alt_types.append("INFRAME DELETION")
						only_alter = whitespace.replace('inframe deletion', '').replace('(', '').replace(')', '').strip()
						
						alts.append(only_alter)
					elif whitespace.startswith('inframe insertion'):
						alt_types.append('INFRAME INSERTION')
						only_alter = whitespace.replace('inframe insertion', '').replace('(', '').replace(')','').strip()
						
						alts.append(only_alter)
					else:
						alts.append(alter_elem)
						alt_types.append("INFRAME ?")
				
				elif extracted_alteration == "FUSION":
					alts.append(element)
					alt_types.append(extracted_alteration)
					
					s_gene_name = gene_name.split(';')
					
					s_gene_name[count_alteration] = element.split(' ',1)[0].replace('-', '--')
					gene_name = ';'.join(s_gene_name)
					
				else:
					if "(" in element and ")" in element:
				
						b_split = element.split("(")
						b_split_right = b_split[1]
						b_split_right = b_split_right.replace(")", "")
						
						alts.append(b_split_right)#do not know how to handle this case otherwise ...
						alt_types.append("MUTATION")
					else:	
						alts.append(extracted_alteration)
						alt_types.append(extracted_alteration)
			
			if had_to_break == True:
				continue
			alt_type_entry = "+".join(alt_types)
			alt_entry = "+".join(alts)

			new_entry = drug + "\t" +  " " +  "\t" + drug_family + "\t" + " " + "\t" + " " + "\t" + gene_name + "\t" + " " + "\t" + " " + "\t" + alt_type_entry +  "\t" + alt_entry + "\t" +  association + "\t" + evidence_level + "\t" + source + "\t" + primary_tumor_type + "\t" + " " + "\t" + "cancer_genome_interpreter" + "\t" + " " + "\t" + "" + "\t" + "" + "\t" + "" + "\n" 
			
			cgi_db_entries.append(new_entry)
			
		elif len(split_plus) == 1:
			alt_types = []
			alts = []
			
			#print(biomarker)
			if len(biomarker.split(' ')) <=1 : #TODO:remove this... had to be done because they did not process their data uniformly
			  print("I do not like this database")
			  continue
			
			extracted_alteration = extractAlteration(biomarker)
			
			if extracted_alteration == "MUTATION":
				
				first_whitespace_split = biomarker.split(' ',1)[1]
				first_whitespace_clean = first_whitespace_split.replace('(', '').replace(')', '')
				
				
				split_comma = first_whitespace_clean.split(',')
			
				
				for alter_elem in split_comma:
					alter_elem = alter_elem.strip()
					alts.append(alter_elem)
					alt_types.append("MUTATION")
			
						
			
				for counter,element in enumerate(alts):
					new_entry = drug + "\t" +  " " +  "\t" + drug_family + "\t" + " " + "\t" + " " + "\t" + gene_name + "\t" + " " + "\t" + " " + "\t" + alt_types[counter] +  "\t" + element + "\t" +  association + "\t" + evidence_level + "\t" + source + "\t" + primary_tumor_type + "\t" + " " + "\t" + "cancer_genome_interpreter" + "\t" + " " + "\t" + "" + "\t" + "" + "\t" + "" + "\n" 
				
					cgi_db_entries.append(new_entry)
					
			elif extracted_alteration == "CHECKING REQUIRED 4" or extracted_alteration == "FRAMESHIFT":
				
				
				first_whitespace_split = biomarker.split(' ',1)[1]
				first_whitespace_clean = first_whitespace_split.replace('(', '').replace(')', '')
				
				
				split_comma = first_whitespace_clean.split(',')
				
				for alter_elem in split_comma:
					alter_elem = alter_elem.strip()
					if alter_elem.startswith('inframe deletion'):
						alt_types.append("INFRAME DELETION")
						only_alter = alter_elem.replace('inframe deletion', '').strip()
						
						alts.append(only_alter)
					elif alter_elem.startswith('inframe insertion'):
						alt_types.append('INFRAME INSERTION')
						only_alter = alter_elem.replace('inframe insertion', '')
						
						alts.append(only_alter)
					elif alter_elem.startswith('frameshift variant'):
						alt_types.append("FRAMESHIFT")
						only_alter = alter_elem.replace('frameshift variant', '').strip()
						
						alts.append(only_alter)
					else:
						alts.append(alter_elem)
						alt_types.append("MUTATION")
			
						
			
				for counter,element in enumerate(alts):
					new_entry = drug + "\t" +  " " +  "\t" + drug_family + "\t" + " " + "\t" + " " + "\t" + gene_name + "\t" + " " + "\t" + " " + "\t" + alt_types[counter] +  "\t" + element + "\t" +  association + "\t" + evidence_level + "\t" + source + "\t" + primary_tumor_type + "\t" + " " + "\t" + "cancer_genome_interpreter" + "\t" + " " + "\t" + "" + "\t" + "" + "\t" + "" + "\n" 
				
					cgi_db_entries.append(new_entry)
					
			elif extracted_alteration == "CHECKING REQUIRED 3":
				conseq_alt = split_semi[0].split(':')[1]
				
				comma_conseq = conseq_alt.split(',')
				if len(comma_conseq) > 1:
					
					for el in comma_conseq:
						alts.append(conseq_alt)
						alt_types.append("MUTATION")
				else:
					alts.append(conseq_alt)
					alt_types.append("MUTATION")
				for counter,element in enumerate(alts):
					new_entry = drug + "\t" +  " " +  "\t" + drug_family + "\t" + " " + "\t" + " " + "\t" + gene_name + "\t" + " " + "\t" + " " + "\t" + alt_types[counter] +  "\t" + element + "\t" +  association + "\t" + evidence_level + "\t" + source + "\t" + primary_tumor_type + "\t" + " " + "\t" + "cancer_genome_interpreter" + "\t" + " " + "\t" + "" + "\t" + "" + "\t" + "" + "\n" 
				
					cgi_db_entries.append(new_entry)
				
			elif extracted_alteration == "CHECKING REQUIRED 2" or extracted_alteration == "CHECKING REQUIRED 1":
				conseq_alt = split_semi[0].split('consequence::')
				conseq_alt_splitted = conseq_alt[1].split(':')
				
				if conseq_alt_splitted[0] == "inframe_insertion":
					alt_types.append("INFRAME INSERTION")
				elif conseq_alt_splitted[0] =="inframe_deletion":
					alt_types.append("INFRAME DELETION")
				else:
					print("Unknown biomarker type found: " + biomarker)
					
				#print(split_semi)
				#print(conseq_alt_splitted)
				alts.append(conseq_alt_splitted[1])
				
				for counter,element in enumerate(alts):
					new_entry = drug + "\t" +  " " +  "\t" + drug_family + "\t" + " " + "\t" + " " + "\t" + gene_name + "\t" + " " + "\t" + " " + "\t" + alt_types[counter] +  "\t" + element + "\t" +  association + "\t" + evidence_level + "\t" + source + "\t" + primary_tumor_type + "\t" + " " + "\t" + "cancer_genome_interpreter" + "\t" + " " + "\t" + "" + "\t" + "" + "\t" + "" + "\n" 
				
					cgi_db_entries.append(new_entry)
					
			elif extracted_alteration == "FUSION":
				
				alteration_ = biomarker
				alteration_type_ = extracted_alteration
				gene_name = biomarker.split(' ',1)[0].replace('-', '--')
				
				new_entry = drug + "\t" +  " " +  "\t" + drug_family + "\t" + " " + "\t" + " " + "\t" + gene_name + "\t" + " " + "\t" + " " + "\t" + alteration_type_+  "\t" + alteration + "\t" +  association + "\t" + evidence_level + "\t" + source + "\t" + primary_tumor_type + "\t" + " " + "\t" + "cancer_genome_interpreter" + "\t" + " " + "\t" + "" + "\t" + "" + "\t" + "" + "\n" 
				
				cgi_db_entries.append(new_entry)
				
			else:
				new_entry = drug + "\t" +  " " +  "\t" + drug_family + "\t" + " " + "\t" + " " + "\t" + gene_name + "\t" + " " + "\t" + " " + "\t" + extracted_alteration +  "\t" + extracted_alteration + "\t" +  association + "\t" + evidence_level + "\t" + source + "\t" + primary_tumor_type + "\t" + " " + "\t" + "cancer_genome_interpreter" + "\t" + " " + "\t" + "" + "\t" + "" + "\t" + "" + "\n" 
				
				cgi_db_entries.append(new_entry)
			
print("Writing CGI data base intermediate summary")

output_cgi_db = "<output_path>"
with open(output_cgi_db, "w", encoding='utf-8') as output_cgi:
	
	output_cgi.write("drug_name" + "\t" + "DrugBank_ID" + "\t" + "drug_family" + "\t" + "drug_targets" + "\t" + "mutation_effect" + "\t" + "associated_gene" + "\t" + "HGNC_ID_associated_gene" + "\t" + "pathways_associated_gene" + "\t" + "alteration_type" + "\t"+ "alteration" + "\t" + "association_implication" + "\t" + "evidence" + "\t" + "source" + "\t" + "primary_tumour_type" + "\t" + "tsg_onco_annotation" + "\t" + "source_database" + "\t" + "evidence_statement" + "\t" + "cds_like_mutation" + "\t" + "genome_coordinates" + "\t" + "reference_build" + "\n" )
	for element in cgi_db_entries:
		
		output_cgi.write(element)
	
	

	

		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
			
	
	
