
print("Parsing oncoKB database")
oncoKB_db = "<path_to_allAnnotatedVariants_from_oncoKB>"

oncoKB_db_entries = []

with open(oncoKB_db, 'r', encoding = "latin-1") as okb_db:
	
	okb_db.readline()#header
	for line in okb_db:
		sline = line.split('\t')
		if len(sline) < 10:
			print("Line with incorrect length found, length: " + str(len(sline)) + "\n" + "Line: " + line)
			continue
			
		gene_name, protein_change, onco, mut_effect, source = sline[3].strip(), sline[5].strip(), sline[6].strip(), sline[7].strip(), sline[8].strip()
					
					
		
		alteration_type = " "
		if 'Fusion' in protein_change:
			if gene_name in protein_change:
				
				gene_name = protein_change.split(' ',1)[0].replace('-', '--')
				alteration_type = "FUSION"
				
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
				
		oncoKB_db_entries.append(gene_name + "\t" + " " +  "\t" + protein_change + "\t" + alteration_type + "\t" + mut_effect + "\t" + onco + "\t" + " " + "\t" + " " + "\t"  + source_pmids + "\t" + "oncoKB" + "\t" + " " + "\t" + " " + "\t" + " " + "\t" + " " + "\t" + " " + "\n")		
		

print("Writing oncoKB data base intermediate summary")
output_oncoKB_db = "<output_path>"
with open(output_oncoKB_db, "w", encoding='utf-8') as output_okb:
	
	output_okb.write("gene_name" + "\t" + "HGNC_ID"  + "\t" +"alteration" + "\t" + "alteration_type" + "\t" + "mutation_effect" + "\t" +  "association_implication" + "\t" + "primary_tumour_type" + "\t" + "tsg_onco_annotation" + "\t" + "source" + "\t" +  "source_database" + "\t" + "evidence_level" + "\t" + "evidence_statement" + "\t" + "cds_like_mutation" + "\t" + "genome_coordinates" + "\t" + "reference_build" + "\n")
	for element in oncoKB_db_entries:
		
		output_okb.write(element)
		
