
print("Parsing CIViC database")
civic_db = "<path_to_nightly-ClinicalEvidenceSummaries_from_CIViC>"

civic_db_entries = []
with open(civic_db, 'r', encoding='utf-8') as c_db:
	
	c_db.readline()#header	
	
	for line in c_db:
		sline = line.strip().split('\t')
		
		if len(sline)< 38:
			print('Line with incorrect length found, length: ' + len(sline) + "\n" + "Line: " + line)
			continue 
			
		gene_name, alteration, disease, drug, evidence_type, evidence_direction, evidence_level, clinical_significance, evidence_statement, pubmed_id, evidence_status = sline[0].strip(), sline[2].strip(), sline[3].strip(), sline[6].strip(), sline[7].strip(), sline[8].strip(), sline[9].strip(), sline[10].strip(), sline[11].strip(), sline[12].strip(), sline[15].strip()
		
		chromosome, start, stop, ref, var, reference_build, chromosome2, start2, stop2 = sline[20].strip(), sline[21].strip(), sline[22].strip(), sline[23].strip(), sline[24].strip(), sline[31].strip(), sline[26].strip(), sline[27].strip(), sline[28].strip()
		
		cds_like_mutation = ""
		genome_coordinates = ""
		if chromosome =="":
			if not chromosome2 =="":
				genome_coordinates = chromosome2 + ":" + start2 + ".." + stop2
			#else:
				#print("")
				#print("Chromosome 1: " + chromosome)
				#print("Start:" + start)
				#print("Stop:" + stop)
				#print("Chromosome 2: " + chromosome2)
				#print("Start:" + start)
				#print("Stop:" + stop)
				#print(line)
				
		else:
			if chromosome2 == "":
				genome_coordinates = chromosome + ":" + start + ".." + stop
			else:
				genome_coordinates = chromosome + ":" + start + ".." + stop + "_+_" + chromosome2 + ":" + start2 + ".." + stop2
				
		if ref =="" and var == "":
			
			cds_like_mutation = ""
		else:
			cds_like_mutation = ref + ">" + var 
			
		if "(" in alteration and ")" in alteration: 
			s_alteration = alteration.split("(")
			right_a = s_alteration[1]
			right_a_rm = right_a.replace(")", "")
			
			if "c." in right_a_rm:
				right_a_rm_rm = right_a_rm.replace("c.", "")
				cds_like_mutation = right_a_rm_rm
				alteration = s_alteration[0]
			
		alteration_type = " "
		if '-' in alteration:
			split_alteration = alteration.split('-')
			
			split_alteration_left = split_alteration[0].strip()
			split_alteration_right = split_alteration[1].strip()
			
			if gene_name in split_alteration_left or gene_name in split_alteration_right:
				gene_name = alteration.split(' ',1)[0].replace('-', '--')
				alteration_type = "FUSION"
				
			elif "BCR-ABL" in alteration:#had to do this because I could not distinguish between mutations in the form X-Y and and fusions
				gene_name = alteration.split(' ',1)[0]
				alteration_type = "FUSION"
				
		if evidence_type == "Prognostic":
			
			if evidence_direction == "Supports":
				association = ""
				if clinical_significance == "Poor Outcome":
					association = "Prognostic:Poor Outcome"
				elif clinical_significance == "Better Outcome":
					association = "Prognostic:Better Outcome"
				elif clinical_significance == "Negative":
					association = "Prognostic:Poor Outcome"
				elif clinical_significance == "Positive":
					association = "Prognostic:Better Outcome"
				elif clinical_significance == "N/A":
					association = "Prognostic:N/A"
				else:
					print("Unknown clinical significance for prognostic biomarker found:" + clinical_significance)
					print(line)
					association = "Prognostic:" + clinical_significance
				
				evidence = ""
				if evidence_level == "A":
					evidence = "validated association"
				elif evidence_level == "B":
					evidence = "clinical evidence from clinical trial" 
				elif evidence_level == "C":
					evidence = "case study"
				elif evidence_level == "D":
					evidence = "preclinial evidence by in vivo or in vitro models" 
				elif evidence_level == "E":
					evidence = "inferential association"
				else:
					print("Unknown evidence level detected, level: " + evidence_level)
					
				new_entry = gene_name + "\t" + " " + "\t" + alteration + "\t" + alteration_type + "\t" + " " + "\t" +  association + "\t" +  disease + "\t" + " " + "\t" + "PMID:" + pubmed_id + "\t" + "CIViC" + "\t" + evidence + "\t" +  evidence_statement +  "\t" + cds_like_mutation + "\t" + genome_coordinates + "\t" + reference_build + "\n"
				
				civic_db_entries.append(new_entry)
					
			elif evidence_direction == "Does Not Support":
				association = ""
				if clinical_significance == "Poor Outcome":
					association = "Prognostic:Poor Outcome"
				elif clinical_significance == "Better Outcome":
					association = "Prognostic:Better Outcome"
				elif clinical_significance == "Negative":
					association = "Prognostic:Poor Outcome"
				elif clinical_significance == "Positive":
					association = "Prognostic:Better Outcome"
				elif clinical_significance == "N/A":
					association = "Prognostic:N/A"
				else:
					print("Unknown clinical significance for prognostic biomarker found:" + clinical_significance)
					print(line)
					association = "Prognostic:" + clinical_significance
					
				evidence = ""
				if evidence_level == "A":
					evidence = "validated association"
				elif evidence_level == "B":
					evidence = "clinical evidence from clinical trial" 
				elif evidence_level == "C":
					evidence = "case study"
				elif evidence_level == "D":
					evidence = "preclinial evidence by in vivo or in vitro models" 
				elif evidence_level == "E":
					evidence = "inferential association"
				else:
					print("Unknown evidence level detected, level: " + evidence_level)
					
				new_entry = gene_name + "\t" + " " + "\t" + alteration + "\t" + alteration_type + "\t" + " " + "\t" +  association + "\t" +  disease + "\t" + " " + "\t" + "PMID:" + pubmed_id + "\t" + "CIViC" + "\t" + evidence + "\t" +  evidence_statement +  "\t" + cds_like_mutation + "\t" + genome_coordinates + "\t" + reference_build +  "\n"
				
				civic_db_entries.append(new_entry)
			else:
				continue
				
		elif evidence_type == "Diagnostic":
			
			if evidence_direction == "Supports":
				association = ""
				if clinical_significance == "Positive":
					association = "Diagnostic:Positive"
				elif clinical_significance == "Negative":
					association = "Diagnostic:Negative"
				else:
					print("Unknown clinical significance for diagnostic biomarker found:" + clinical_significance)
					association = "Diagnostic:" + clinical_significance
				
				evidence = ""
				if evidence_level == "A":
					evidence = "validated association"
				elif evidence_level == "B":
					evidence = "clinical evidence from clinical trial" 
				elif evidence_level == "C":
					evidence = "case study"
				elif evidence_level == "D":
					evidence = "preclinial evidence by in vivo or in vitro models" 
				elif evidence_level == "E":
					evidence = "inferential association"
				else:
					print("Unknown evidence level detected, level: " + evidence_level)
					
				new_entry = gene_name + "\t" + " " + "\t" + alteration + "\t" + alteration_type + "\t" +  " " + "\t" + association + "\t" +  disease + "\t" + " " + "\t" + "PMID:" + pubmed_id + "\t" + "CIViC" + "\t" + evidence + "\t" +  evidence_statement +  "\t" + cds_like_mutation + "\t" + genome_coordinates + "\t" + reference_build +  "\n"
				
				civic_db_entries.append(new_entry)
			
			elif evidence_direction == "Does Not Support":
				association = ""
				if clinical_significance == "Positive":
					association = "Diagnostic:Not Positive"
				elif clinical_significance == "Negative":
					association = "Diagnostic:Not Negative"
				else:
					print("Unknown clinical significance for diagnostic biomarker found:" + clinical_significance)
					association = "Diagnostic:" + clinical_significance
				
				evidence = ""
				if evidence_level == "A":
					evidence = "validated association"
				elif evidence_level == "B":
					evidence = "clinical evidence from clinical trial" 
				elif evidence_level == "C":
					evidence = "case study"
				elif evidence_level == "D":
					evidence = "preclinial evidence by in vivo or in vitro models" 
				elif evidence_level == "E":
					evidence = "inferential association"
				else:
					print("Unknown evidence level detected, level: " + evidence_level)
					
				new_entry = gene_name + "\t" + " " + "\t" + alteration + "\t" + alteration_type + "\t" +  " " + "\t" + association + "\t" +  disease + "\t" + " " + "\t" + "PMID:" + pubmed_id + "\t" + "CIViC" + "\t" + evidence + "\t" +  evidence_statement +  "\t" + cds_like_mutation + "\t" + genome_coordinates + "\t" + reference_build +  "\n"
				
				civic_db_entries.append(new_entry)
			else:
				continue
			
		elif evidence_type == "Predisposing":
			if evidence_direction == "Supports":
				association = ""
				if clinical_significance == "Pathogenic":
					association = "Predisposing:Pathogenic"
				elif clinical_significance == "Likely Pathogenic":
					association = "Predisposing:Likely Pathogenic"
				elif clinical_significance == "Positive":
					association = "Predisposing:Positive"
				elif clinical_significance == "Uncertain Significance":
					association = "Predisposing:Uncertain Significance"
				elif clinical_significance == "" or clinical_significance == "N/A":
					association = "Predisposing:No Significance Given"
				else:
					print("Unknown clinical significance for predisposing biomarker found:" + clinical_significance)
					print(line)
					association = "Predisposing:" + clinical_significance
				
				evidence = ""
				if evidence_level == "A":
					evidence = "validated association"
				elif evidence_level == "B":
					evidence = "clinical evidence from clinical trial" 
				elif evidence_level == "C":
					evidence = "case study"
				elif evidence_level == "D":
					evidence = "preclinial evidence by in vivo or in vitro models" 
				elif evidence_level == "E":
					evidence = "inferential association"
				else:
					print("Unknown evidence level detected, level: " + evidence_level)
					
				new_entry = gene_name + "\t" + " " + "\t" + alteration + "\t" + alteration_type + "\t" +  " " + "\t" + association + "\t" +  disease + "\t" + " " + "\t" + "PMID:" + pubmed_id + "\t" + "CIViC" + "\t" + evidence + "\t" +  evidence_statement +  "\t" + cds_like_mutation + "\t" + genome_coordinates + "\t" + reference_build +  "\n"
				
				civic_db_entries.append(new_entry)
			
			elif evidence_direction == "Does Not Support":
				association = ""
				if clinical_significance == "Pathogenic":
					association = "Predisposing:Pathogenic"
				elif clinical_significance == "Likely Pathogenic":
					association = "Predisposing:Likely Pathogenic"
				elif clinical_significance == "Positive":
					association = "Predisposing:Positive"
				elif clinical_significance == "Uncertain Significance":
					association = "Predisposing:Uncertain Significance"
				elif clinical_significance == "" or clinical_significance == "N/A":
					association = "Predisposing:No Significance Given"
				else:
					print("Unknown clinical significance for predisposing biomarker found:" + clinical_significance)
					association = "Predisposing:" + clinical_significance
				
				evidence = ""
				if evidence_level == "A":
					evidence = "validated association"
				elif evidence_level == "B":
					evidence = "clinical evidence from clinical trial" 
				elif evidence_level == "C":
					evidence = "case study"
				elif evidence_level == "D":
					evidence = "preclinial evidence by in vivo or in vitro models" 
				elif evidence_level == "E":
					evidence = "inferential association"
				else:
					print("Unknown evidence level detected, level: " + evidence_level)
					
				new_entry = gene_name + "\t" + " " + "\t" + alteration + "\t" + alteration_type + "\t" + " " + "\t" + association + "\t" +  disease + "\t" + " " + "\t" + "PMID:" + pubmed_id + "\t" + "CIViC" + "\t" + evidence + "\t" +  evidence_statement +  "\t" + cds_like_mutation + "\t" + genome_coordinates + "\t" + reference_build + "\n"
				
				civic_db_entries.append(new_entry)
			
			else:
				continue 
			
		elif evidence_type == "Predictive":
			continue
		else:
			print("Unknown evidence type detected: " + evidence_type)
			continue
		

print("Writing CIViC database intermediate summary")
output_civic_db = "<output_path>"
with open(output_civic_db, "w", encoding='utf-8') as ocdb:
	
	ocdb.write("gene_name" + "\t" + "HGNC_ID"  + "\t" +"alteration" + "\t" + "alteration_type" + "\t"  + "mutation_effect" + "\t" + "association_implication" + "\t" + "primary_tumour_type" + "\t" + "tsg_onco_annotation" + "\t" + "source" + "\t" +  "source_database" + "\t" + "evidence_level" + "\t" + "evidence_statement" + "\t" + "cds_like_mutation" + "\t" + "genome_coordinates" + "\t" + "reference_build" + "\n" )
	for element in civic_db_entries:
		
		ocdb.write(element)
		
		
