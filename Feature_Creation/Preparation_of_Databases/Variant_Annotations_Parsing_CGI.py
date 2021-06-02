
print("Parsing Cancer Genome Interpreter database")
cancer_genome_interpreter_db = "<path_to_catalog_of_validated_oncogenic_mutations_from_CGI>"

cgi_db_entries = []
with open(cancer_genome_interpreter_db, 'r', encoding='utf-8') as cgi_db:
	
	cgi_db.readline()#header
	
	for line in cgi_db:
		sline = line.strip().split('\t')
		if len(sline) < 9:
			print("Line with incorrect length found, length: " + str(len(sline)) + "\n" + "Line: " + line)
			continue
			 
		gene_name, cds_like, protein_change, cancer_type, source_db_cgi, source = sline[0].strip(), sline[1].strip(), sline[2].strip(), sline[6].strip(), sline[7].strip(), sline[8].strip()

		if source_db_cgi == "Biomarker":
			continue
		
		
		cds_like_split = []
		if not cds_like =="":
			cds_like_split = cds_like.split("__")
		
		count = 0
		for element in cds_like_split:
			c_elem = element.replace("chr","")
			c_elem2 = c_elem.replace("g.", "")
			
			if count == 0:
				cds_like = c_elem2
			elif not count == 0:
				cds_like = cds_like + "," + c_elem2
			count = count +1.0
		
		protein_change = protein_change.replace("p.", "")
		
		new_entry = gene_name + "\t" + " " + "\t"  + protein_change  + "\t" +  " "  + "\t" + " " + "\t" +  "validated oncogenic" + "\t" + cancer_type + "\t" + " " + "\t" + source + "\t" + "cancer_genome_interpreter" + "(" + source_db_cgi + ")" + "\t" + "validated" + "\t" + " " + "\t" + cds_like + "\t" + " " + "\t" + "missing reference build" + "\n" 
		cgi_db_entries.append(new_entry)
			
print("Writing CGI data base intermediate summary")

output_cgi_db = "<output_path>"
with open(output_cgi_db, "w", encoding='utf-8') as output_cgi:
	
	output_cgi.write("gene_name" + "\t" + "HGNC_ID"  + "\t" +"alteration" + "\t" + "alteration_type" + "\t"  + "mutation_effect" + "\t" + "association_implication" + "\t" + "primary_tumour_type" + "\t" + "tsg_onco_annotation" + "\t" + "source" + "\t" +  "source_database" + "\t" + "evidence_level" + "\t" + "evidence_statement" + "\t" + "cds_like_mutation" + "\t" + "genome_coordinates" + "\t" + "reference_build" + "\n" )
	for element in cgi_db_entries:
		
		output_cgi.write(element)
	
