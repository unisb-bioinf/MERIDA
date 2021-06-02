
print("Parsing oncoKB database")
oncoKB_db = "<path_to_allActionableVariants_oncoKB>"

oncoKB_db_entries = []

with open(oncoKB_db, 'r', encoding='utf-8') as okb_db:
	
	okb_db.readline()#header
	
	
	for line in okb_db:
		
		sline = line.split('\t')
		if len(sline) < 11:
			print("Line with incorrect length found, length: " + str(len(sline)) + "\n" + "Line: " + line)
			continue
			
		gene_name, alteration, protein_change, cancer_type, evidence_plus_association, drug, source = sline[3].strip(), sline[4].strip(), sline[5].strip(), sline[6], sline[7].strip(), sline[8].strip(), sline[9].strip()
		
		several_drugs = []
		if "," in drug:
			several_drugs = drug.split(",")
			
		
		alteration_type = " "
		if 'Fusion' in protein_change:
			if gene_name in protein_change:
				
				gene_name = protein_change.split(' ',1)[0].replace('-', '--')
				alteration_type = "FUSION"
				
		if evidence_plus_association == "R1":
			association = "resistant"
			evidence = "standard care biomarker for FDA approved drug in this indication"
		else:
			association = "sensitive"
			
			if evidence_plus_association == "1":
				evidence = "FDA-recognized biomarker for FDA approved drug in this indication"
				
			else:
				if evidence_plus_association == "2A":
					evidence = "standard-care biomarker for FDA approved drug in this indication"
				else:
					if evidence_plus_association == "2B":
						evidence = "standard-care biomarker for FDA approved drug in other indication"
						
					else:
						if evidence_plus_association == "3A":
							evidence = "clinical evidence for this indication and neither drug nor biomarker are standard care"
						
						else:
							if evidence_plus_association == "3B":
								evidence = "clinical evidence for other indication and neither drug nor biomarker are standard care"
							
							else:
								if evidence_plus_association == "4":
									evidence = "biological evidence and neither drug nor biomarker are standard care"
								else:
									print("Wrong level of evidence found:" + evidence_level)
								
		pmids = source.strip().split(",")
		
		if pmids[0] == "":
			source_pmids = ""
		else:
			source_pmids = "PMID:" + pmids[0]
		
		count = 0
		for pmid in pmids:
			if count == 0 and len(pmids) == 1:
				break
			
			count = count +1.0
			
			source_pmids = source_pmids + "," +  "PMID:" + pmid
			
		for dr in several_drugs:
			oncoKB_db_entries.append(dr.strip() + "\t" + " " + "\t" + " " + "\t" + " " + "\t" + " " + "\t" + gene_name + "\t" + " " + "\t" + " " + "\t" + alteration_type + "\t" + protein_change + "\t" + association + "\t" + evidence + "\t" + source_pmids + "\t" + cancer_type + "\t" + " " + "\t" + "oncoKB" + "\t" + " " + "\t" + " " + "\t" + " " + "\t" + " " + "\n")		
		
	
print("Writing oncoKB data base intermediate summary")
output_oncoKB_db = "<output_file>"
with open(output_oncoKB_db, "w", encoding='utf-8') as output_okb:
	
	output_okb.write("drug_name" + "\t" + "DrugBank_ID" + "\t" + "drug_family" + "\t" + "drug_targets" + "\t" + "mutation_effect" + "\t" + "associated_gene" + "\t" + "HGNC_ID_associated_gene" + "\t" + "pathways_associated_gene" + "\t" + "alteration_type" + "\t" + "alteration" + "\t" + "association_implication" + "\t" + "evidence" + "\t" + "source" + "\t" + "primary_tumour_type" + "\t" + "tsg_onco_annotation" + "\t" + "source_database" + "\t" + "evidence_statement" + "\t" + "cds_like_mutation" + "\t" + "genome_coordinates" + "\t" + "reference_build" + "\n" )
	for element in oncoKB_db_entries:
		
		output_okb.write(element)
		
		
