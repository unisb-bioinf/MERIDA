
print("COSMIC database")
cosmic_db = "<path_to_CosmicResistanceMutations>"

cosmic_db_entries = []
with open(cosmic_db, 'r', encoding='utf-8') as co_db:
	
	co_db.readline()#header
	
	for line in co_db:
		
		sline = line.split("\t")
		
		if len(sline )< 22:
			print('Line with incorrect length found, length: ' + str(len(sline)) + "\n" + "Line: " + line)
			continue 
		
		gene_name, drug_name, protein_change, histology, pubmed_id, evidence_level = sline[2].strip(), sline[5].strip(), sline[7].strip(), sline[12].strip(), sline[15].strip(), sline[21].strip()
		
		s_gene = gene_name.split('_')
		if len(s_gene) > 1:
			gene_name = s_gene[0]
		cds_mut, genome_coordinates, reference_build = sline[8].strip(), sline[20], "GRCh38"
		
		cds_like_mutation = ""
		if not cds_mut.startswith('c.?'):
			
			cds_like_mutation = cds_mut.replace('c.', '')
		
		if protein_change.startswith('p.'):
			protein_change = protein_change.replace('p.', '')
		
		evidence = ""
		if evidence_level == "1":
			evidence = "strong evidence"
		elif evidence_level == "2":
			evidence = "emerging evidence" 
		else:
			evidence = "not specified by COSMIC"
		new_entry = drug_name + "\t" + " "  + "\t" + " " + "\t" + " " + "\t" + " " + "\t" + gene_name + "\t" + " " + "\t" + " " + "\t" + " " + "\t" + protein_change + "\t" + "resistant" + "\t" + evidence + "\t" + "PMID:" + pubmed_id  + "\t" + histology + "\t" + " " + "\t" + "COSMIC" + "\t" + " " + "\t" + cds_like_mutation + "\t" + genome_coordinates + "\t" + reference_build + "\n"
		
		cosmic_db_entries.append(new_entry)
		

print("Writing COSMIC database intermediate summary")
output_cosmic_db = "<output_file>"
with open(output_cosmic_db, "w", encoding='utf-8') as ocodb:
	
	ocodb.write("drug_name" + "\t" + "DrugBank_ID" + "\t" + "drug_family" + "\t" + "drug_targets" + "\t" + "mutation_effect" + "\t" + "associated_gene" + "\t" + "HGNC_ID_associated_gene" + "\t" + "pathways_associated_gene" + "\t" + "alteration_type" + "\t" + "alteration" + "\t"+ "association_implication" + "\t" + "evidence" + "\t" + "source" + "\t" + "primary_tumour_type" + "\t" + "tsg_onco_annotation" + "\t" + "source_database" + "\t" + "evidence_statement" + "\t" + "cds_like_mutation" + "\t" + "genome_coordinates" + "\t" + "reference_build" + "\n" )
	for element in cosmic_db_entries:
		
		ocodb.write(element)		
	
		
