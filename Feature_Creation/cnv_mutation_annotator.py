import sys
sys.path.insert(0, '/home/student/klenhof/Research_Doktorarbeit/DrugSensitivityPredictionDatabase_WS2018/Mapping_Files')
sys.path.insert(0, '/home/student/klenhof/Research_Doktorarbeit/DrugSensitivityPredictionDatabase_WS2018/Database_Query')

from query_oncogenicity_db import *
from query_sensitivity_db import *

#This class can be used to annotate mutations and copy number variations 
class CNV_Mutation_Annotator:
	
	def __init__(self, drug_name, mode_for_oncosearch):
		
		
		#List of alterations that are implicated in sensitivity towards the queried drug
		self.drug_sens_alterations = list(query_drug_filtered(drug_name, "no_output_file").values())
		#Determine list with genes that have sensitivity/resistance information attached
		self.exact_point_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res = {}
		self.approximate_point_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res = {}
		self.copy_number_annotation_gene_to_loss_or_gain_to_sens_res = {}
		self.frameshift_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res = {}
		self.truncating_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res = {}
		for entry in self.drug_sens_alterations:
			
			s_entry = entry.strip('\n').split('\t')
			gene_name = s_entry[5].strip()
			alteration = s_entry[9].strip()
			mutation_effect = s_entry[4].strip()
			implication = s_entry[10].strip()
			
			if "FS*" in alteration:
				alteration.replace("FS*", "fs*")
			
			if alteration in ["loss", "LOSS", "DELETION", "Deletion", "Copy Number Loss"]:
				alteration = "CNVloss" # nomenclature should be equal for the normal copy number loss and the annotated copy number loss
				
				if not gene_name in self.copy_number_annotation_gene_to_loss_or_gain_to_sens_res:
					self.copy_number_annotation_gene_to_loss_or_gain_to_sens_res[gene_name] = {alteration:implication}
					continue
				else:
					if not alteration in self.copy_number_annotation_gene_to_loss_or_gain_to_sens_res[gene_name]:
						self.copy_number_annotation_gene_to_loss_or_gain_to_sens_res[gene_name][alteration] = implication
						continue
					else:
						continue
						
				
			elif alteration in ["AMPLIFICATION", "Amplification"]:
				alteration = "CNVgain"
				
				if not gene_name in self.copy_number_annotation_gene_to_loss_or_gain_to_sens_res:
					self.copy_number_annotation_gene_to_loss_or_gain_to_sens_res[gene_name] = {alteration:implication}
					continue
				else:
					if not alteration in self.copy_number_annotation_gene_to_loss_or_gain_to_sens_res[gene_name]:
						self.copy_number_annotation_gene_to_loss_or_gain_to_sens_res[gene_name][alteration] = implication
						continue
					else:
						continue
					
			if alteration in ["FRAMESHIFT TRUNCATION", "FRAMESHIFT", "frameshift"]:
				if not gene_name in self.frameshift_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res:
					self.frameshift_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res[gene_name] = {"general_frameshift":implication}
					continue
				else:
					if not alteration in self.frameshift_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res[gene_name]:
						self.frameshift_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res[gene_name]["general_frameshift"] = implication
						continue
					else:
						continue
			
			if alteration in ["Truncating Mutations"]:
				if not gene_name in self.truncating_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res:
					self.truncating_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res[gene_name] = {"truncating_mutation":implication}
					continue
				else:
					if not alteration in self.truncating_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res[gene_name]:
						self.truncating_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res[gene_name]["truncating_mutation"] = implication
						continue
					else:
						continue
			
			if 'fs*' in alteration or not '*' in alteration:	
			
				if not gene_name in self.exact_point_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res:
					self.exact_point_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res[gene_name] = {alteration: [mutation_effect, implication]}
				else:
					if not alteration in self.exact_point_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res[gene_name]:
						self.exact_point_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res[gene_name][alteration] = [mutation_effect, implication]
				
					else:
						print("5 Found an alteration in a gene twice: " + gene_name + " " + alteration )
						continue
			
			else:
				#now '*' is in alteration
				
				if not gene_name in self.approximate_point_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res:
					self.approximate_point_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res[gene_name] = {alteration: [mutation_effect, implication]}
				else:
					if not alteration in self.approximate_point_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res[gene_name]:
						self.approximate_point_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res[gene_name][alteration] = [mutation_effect, implication]
						
					else:
						print("6 Found an alteration in a gene twice: " + gene_name + " " + alteration )
						continue


		#List with all genes with reported oncogenic/neutral/inconclusive mutations
		self.all_onco_info_genes= get_list_of_genes()
		#oncogenic alterations
		self.onco_alterations_all_genes = {}
		#neutral alterations
		self.neutral_alterations_all_genes = {}
		
		#print(mode_for_oncosearch)
		#print(len(self.all_onco_info_genes))
		for gene in self.all_onco_info_genes:
			
			print(gene)
			onco_info_one_gene = query_filtered_onco_db(gene, "no_output_file", mode_for_oncosearch)
		
			neutral_info_one_gene = query_filtered_onco_db(gene, "no_output_file", "gene_neutral")
			
			if onco_info_one_gene:
				self.onco_alterations_all_genes[gene] = list(onco_info_one_gene.values())
			if neutral_info_one_gene:
				self.neutral_alterations_all_genes[gene] = list(neutral_info_one_gene.values())
		
		#exact oncogenic point mutations
		self.exact_point_mutations_onco_gene_to_alteration_to_mutation_effect = {}
		
		#approximate oncogenic point mutations
		self.approximate_point_mutations_onco_gene_to_alteration_to_mutation_effect = {}
		
		#frameshift & truncating mutations 
		self.frameshift_mutations_onco_gene_to_alteration_to_mutation_effect = {}
		self.truncating_mutations_onco_gene_to_alteration_to_mutation_effect = {}
		
		for gene in self.all_onco_info_genes:
			
			if gene in self.onco_alterations_all_genes:
				
				for entry in self.onco_alterations_all_genes[gene]:
					sentry = entry.strip('\n').split('\t')
					alteration = sentry[2].strip() #should be the alteration
					mutation_effect = sentry[4].strip()#should be gain/loss of function annotation
					
					if "FS*" in alteration:
						alteration.replace("FS*", "fs*")
					
					
					if alteration in ["FRAMESHIFT TRUNCATION", "FRAMESHIFT", "frameshift"]:
						if not gene in self.frameshift_mutations_onco_gene_to_alteration_to_mutation_effect:
							self.frameshift_mutations_onco_gene_to_alteration_to_mutation_effect[gene] = {"general_frameshift":mutation_effect}
							continue
						else:
							continue
						#else:#TODO:check whether this makes sense ---> solution: should not be neccessary
						#	if not alteration in self.frameshift_mutations_onco_gene_to_alteration_to_mutation_effect[gene]:
						#		self.frameshift_mutations_onco_gene_to_alteration_to_mutation_effect[gene]["general_frameshift"] = mutation_effect
						#		continue
				
					if alteration in ["Truncating Mutations"]:
						if not gene in self.truncating_mutations_onco_gene_to_alteration_to_mutation_effect:
							self.truncating_mutations_onco_gene_to_alteration_to_mutation_effect[gene] = {"truncating_mutation":mutation_effect}
							continue
						else:
							continue
						#else:#TODO:check whether this makes sense ---> solution: should not be neccessary
						#	if not alteration in self.truncating_mutations_onco_gene_to_alteration_to_mutation_effect[gene]:
						#		self.truncating_mutations_onco_gene_to_alteration_to_mutation_effect[gene]["truncating_mutation"] = mutation_effect
						#		continue
				
					
					
					
					if not alteration in ["loss", "LOSS", "DELETION", "Deletion", "Copy Number Loss", "AMPLIFICATION", "Amplification"]:#It could be important to know whether a copy number gain or loss is oncogenic, for now I ignore this
						
						if 'fs*' in alteration or not '*' in alteration:
							
							if not gene in self.exact_point_mutations_onco_gene_to_alteration_to_mutation_effect:
								self.exact_point_mutations_onco_gene_to_alteration_to_mutation_effect[gene] = {alteration:mutation_effect}
							else:
								if not alteration in self.exact_point_mutations_onco_gene_to_alteration_to_mutation_effect[gene]:
									self.exact_point_mutations_onco_gene_to_alteration_to_mutation_effect[gene][alteration] = mutation_effect
								else:
									print(" 1 Found an alteration in a gene twice: " + gene + " " + alteration )
									continue
						else:#now '*' is in string
							
							if not gene in self.approximate_point_mutations_onco_gene_to_alteration_to_mutation_effect:
								self.approximate_point_mutations_onco_gene_to_alteration_to_mutation_effect[gene] = {alteration:mutation_effect}
							else:
								if not alteration in self.approximate_point_mutations_onco_gene_to_alteration_to_mutation_effect[gene]:
									self.approximate_point_mutations_onco_gene_to_alteration_to_mutation_effect[gene][alteration] = mutation_effect
								else:
									print("2 Found an alteration in a gene twice: " + gene + " " + alteration )
									continue
									
		#neutral alterations: exact point mutations
		self.exact_point_mutations_neutral_gene_to_alteration_to_mutation_effect = {}
		
		#neutral alterations: approximate point mutations
		self.approximate_point_mutations_neutral_gene_to_alteration_to_mutation_effect = {}
		
		
		for gene in self.all_onco_info_genes:
			
			if gene in self.neutral_alterations_all_genes:
				
				for entry in self.neutral_alterations_all_genes[gene]:
					sentry = entry.strip('\n').split('\t')
					alteration = sentry[2].strip() #should be the alteration
					mutation_effect = sentry[4].strip()#should be gain/loss of function annotation
					
					if "FS*" in alteration:
						alteration.replace("FS*", "fs*")
					
					if not alteration in ["loss", "LOSS", "DELETION", "Deletion", "Copy Number Loss", "AMPLIFICATION", "Amplification"]:#It could be important to know whether a copy number gain or loss is oncogenic, for now I ignore this
						
						if 'fs*' in alteration or not '*' in alteration:
							
							if not gene in self.exact_point_mutations_neutral_gene_to_alteration_to_mutation_effect:
								self.exact_point_mutations_neutral_gene_to_alteration_to_mutation_effect[gene] = {alteration:mutation_effect}
							else:
								if not alteration in self.exact_point_mutations_neutral_gene_to_alteration_to_mutation_effect[gene]:
									self.exact_point_mutations_neutral_gene_to_alteration_to_mutation_effect[gene][alteration] = mutation_effect
								else:
									print("3 Found an alteration in a gene twice: " + gene_name + " " + alteration )
									continue
							
						else:#now '*' is in string
							
							if not gene in self.approximate_point_mutations_neutral_gene_to_alteration_to_mutation_effect:
								self.approximate_point_mutations_neutral_gene_to_alteration_to_mutation_effect[gene] = {alteration:mutation_effect}
							else:
								if not alteration in self.approximate_point_mutations_neutral_gene_to_alteration_to_mutation_effect[gene]:
									self.approximate_point_mutations_neutral_gene_to_alteration_to_mutation_effect[gene][alteration] = mutation_effect
								else:
									print("4 Found an alteration in a gene twice: " + gene_name + " " + alteration )
									continue
									
	
		#with open("Alterations_logged.txt", "w") as alteration_log:
		#	alteration_log.write(self.exact_point_mutations_onco_gene_to_alteration_to_mutation_effect)
		#	alteration_log.write(self.approximate_point_mutations_onco_gene_to_alteration_to_mutation_effect)
		#print(self.exact_point_mutations_onco_gene_to_alteration_to_mutation_effect)
		#print(self.approximate_point_mutations_onco_gene_to_alteration_to_mutation_effect)
		
	#Goes through all possible annotations in a certain order 
	#Returns an 'annotation object': [gene, alteration, mutation_effect, implication]
	#If a certain alteration cannot be found in the database, the annotation object will be of the form [gene, alteration, "", ""]
	def annotate_mutation(self, gene_name, alteration, is_truncating):
		
		norm_gene_name = gene_name.strip()
		norm_alteration = alteration.strip()#need to know whether alteration is given in the form 'p.' or only the amino acid change
		
		#At first: check if there is sensitivity/resistance information attached
		#Start with the alterations that require a perfect match (without *)
		if norm_gene_name in self.exact_point_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res:
			
			if norm_alteration in self.exact_point_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res[norm_gene_name]:
				
				print("Exact Mutation Match found (Sensitivity and Resistance Information)")
				return([gene_name, alteration, self.exact_point_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res[norm_gene_name][norm_alteration][0],self.exact_point_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res[norm_gene_name][norm_alteration][1]])
				
		#Go on with those that do not have a perfect match (with *)
		if norm_gene_name in self.approximate_point_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res:
			
			for approximate_alteration in self.approximate_point_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res[norm_gene_name]:
				
				prefix_approximate = approximate_alteration.split('*')[0] # * will be in it because the dictionary was constructed that way, will automatically fail if the string startswith *
				
				if prefix_approximate == "":#just in case
					continue
				if prefix_approximate in norm_alteration:
					
					split_norm_alteration = norm_alteration.split(prefix_approximate)[1] #relevant part of alteration when splitted at common prefix
					
					
					if split_norm_alteration in ["", "*", "."] or not split_norm_alteration[0].isdigit():
					
						print("Approximate Match found (Sensitivity and Resistance Information)")
						return([gene_name, approximate_alteration, self.approximate_point_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res[norm_gene_name][approximate_alteration][0],  self.approximate_point_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res[norm_gene_name][approximate_alteration][1]])
			
		#Now check oncogenic mutations 
		if norm_gene_name in self.exact_point_mutations_onco_gene_to_alteration_to_mutation_effect:
			
			if norm_alteration in self.exact_point_mutations_onco_gene_to_alteration_to_mutation_effect[norm_gene_name]:
				
				mutation_effect = self.exact_point_mutations_onco_gene_to_alteration_to_mutation_effect[norm_gene_name][norm_alteration]
				
				if mutation_effect in ["Likely Loss-of-function", "Loss-of-function", "Likely Gain-of-function", "Gain-of-function"]: 
					print("Exact Mutation Match found (Oncogenic Mutation 1)")
					return([gene_name, alteration, self.exact_point_mutations_onco_gene_to_alteration_to_mutation_effect[norm_gene_name][norm_alteration], ""])
				else:
					print("Exact Mutation Match found (Oncogenic Mutation 2)")
					return([gene_name, alteration, "oncogenic", "" ])#TODO:check if returning oncogenic is correct
				
		
		if norm_gene_name in self.approximate_point_mutations_onco_gene_to_alteration_to_mutation_effect:
			
			for approximate_alteration in self.approximate_point_mutations_onco_gene_to_alteration_to_mutation_effect[norm_gene_name]:
				prefix_approximate = approximate_alteration.split('*')[0] # * will be in it because the dictionary was constructed that way, will automatically fail if the string startswith *
				
				if prefix_approximate == "":
					continue
				
				if prefix_approximate in norm_alteration:
					
					
					mutation_effect = self.approximate_point_mutations_onco_gene_to_alteration_to_mutation_effect[norm_gene_name][approximate_alteration]
					
					if mutation_effect in ["Likely Loss-of-function", "Loss-of-function", "Likely Gain-of-function", "Gain-of-function"]:
						
						split_norm_alteration = norm_alteration.split(prefix_approximate)[1] #relevant part of alteration when splitted at common prefix
					
					
						if split_norm_alteration in ["", "*", "."] or not split_norm_alteration[0].isdigit():
							print("Approximate Match found (Oncogenic Mutation 1)")
							return([gene_name, approximate_alteration, self.approximate_point_mutations_onco_gene_to_alteration_to_mutation_effect[norm_gene_name][approximate_alteration], ""])
						
					else:
						split_norm_alteration = norm_alteration.split(prefix_approximate)[1] #relevant part of alteration when splitted at common prefix
					
					
						if split_norm_alteration in ["", "*", "."] or not split_norm_alteration[0].isdigit():
						
							print("Approximate Match found (Oncogenic Mutation 2)")
							return([gene_name, alteration, "oncogenic", "" ])#TODO
							
			
		#Now check neutral mutations
		
		if norm_gene_name in self.exact_point_mutations_neutral_gene_to_alteration_to_mutation_effect:
			
			if norm_alteration in self.exact_point_mutations_neutral_gene_to_alteration_to_mutation_effect[norm_gene_name]:
				print("Exact Match found (Neutral Mutations)")
				return([gene_name, alteration, self.exact_point_mutations_neutral_gene_to_alteration_to_mutation_effect[norm_gene_name][norm_alteration], ""])
			
		if norm_gene_name in self.approximate_point_mutations_neutral_gene_to_alteration_to_mutation_effect:
			
			for approximate_alteration in self.approximate_point_mutations_neutral_gene_to_alteration_to_mutation_effect[norm_gene_name]:
				prefix_approximate = approximate_alteration.split('*')[0] # * will be in it because the dictionary was constructed that way, will automatically fail if the string startswith *
				
				if prefix_approximate == "":
					continue
				
				if prefix_approximate in norm_alteration:
					
					split_norm_alteration = norm_alteration.split(prefix_approximate)[1] #relevant part of alteration when splitted at common prefix
					
					
					if split_norm_alteration in ["", "*", "."] or not split_norm_alteration[0].isdigit():
						print("Approximate Match found (Neutral Mutations)")
						return([gene_name, approximate_alteration, self.approximate_point_mutations_neutral_gene_to_alteration_to_mutation_effect[norm_gene_name][approximate_alteration], ""])	
		
		#Now check mutations with no further genomic position specified 
		if norm_gene_name in self.frameshift_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res:
			
			if "fs" in norm_alteration:
				print("Frameshift Match found (Sensitivity and Resistance Information)")
				return([gene_name, "general_frameshift", "", self.frameshift_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res[norm_gene_name]["general_frameshift"]])
			
			
		if norm_gene_name in self.truncating_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res:
			
			if is_truncating:
				print("Truncating Match found (Sensitivity and Resistance Information)")
				return([gene_name, "truncating_mutation", "", self.truncating_mutations_all_gene_name_to_aberration_to_aberration_effect_and_sens_res[norm_gene_name]["truncating_mutation"]])
			
		if norm_gene_name in self.frameshift_mutations_onco_gene_to_alteration_to_mutation_effect:
			
			if "fs"  in norm_alteration:
				
				print(self.frameshift_mutations_onco_gene_to_alteration_to_mutation_effect)
				mutation_effect = self.frameshift_mutations_onco_gene_to_alteration_to_mutation_effect[norm_gene_name]["general_frameshift"]
				
				if mutation_effect in ["Likely Loss-of-function", "Loss-of-function", "Likely Gain-of-function", "Gain-of-function"]:
					print("Frameshift Match found (Oncogenic Alterations 1)")
					return([gene_name, "general_frameshift", self.frameshift_mutations_onco_gene_to_alteration_to_mutation_effect[norm_gene_name]["general_frameshift"], ""])
				else:
					print("Frameshift Match found (Oncogenic Alterations 2)")
					return([gene_name, "general_frameshift", "oncogenic", "" ])#TODO
			
			
		if norm_gene_name in self.truncating_mutations_onco_gene_to_alteration_to_mutation_effect:
			
			if is_truncating:
				if norm_gene_name == "TSC1":
					print("TSC1 truncating mutation successfully found and matched")
				
				mutation_effect = self.truncating_mutations_onco_gene_to_alteration_to_mutation_effect[norm_gene_name]["truncating_mutation"]
				
				if mutation_effect in ["Likely Loss-of-function", "Loss-of-function", "Likely Gain-of-function", "Gain-of-function"]:
					print("Truncating Match found (Oncogenic Alterations 1)")
					return([gene_name, "truncating_mutation", self.truncating_mutations_onco_gene_to_alteration_to_mutation_effect[norm_gene_name]["truncating_mutation"], ""])
				else:
					print("Truncating Match found (Oncogenic Alterations 2)")
					return([gene_name, "truncating_mutation", "oncogenic", "" ])#TODO
				
		
		#Everything is checked, we found no match in the database:
		
		return([gene_name, alteration, "Unknown", ""])
		
		             
		
		
	#Goes through all possible annotations 
	#Returns an 'annotation object': [gene, alteration, mutation_effect (=""), implication]
	#If a certain alteration cannot be found in the database, the annotation object will be of the form [gene, alteration, "", ""]
	#Currently only alterations named "loss" and "gain" will be recognized
	def annotate_cnv(self, gene_name, alteration):
		
		
		#alterations can only be "loss" or "gain" anything else will not be recognized
		norm_gene_name = gene_name.strip()
		norm_alteration = alteration.strip()
		
		if norm_gene_name in self.copy_number_annotation_gene_to_loss_or_gain_to_sens_res:
			
			if norm_alteration == "loss":
				
				if "CNVloss" in self.copy_number_annotation_gene_to_loss_or_gain_to_sens_res[norm_gene_name]:
					return([gene_name, "CNVloss", "CNVloss", self.copy_number_annotation_gene_to_loss_or_gain_to_sens_res[norm_gene_name]["CNVloss"]])
				else:
					return([gene_name, "CNVloss", "CNVloss", ""])
			if norm_alteration == "gain":
				if "CNVgain" in self.copy_number_annotation_gene_to_loss_or_gain_to_sens_res[norm_gene_name]:
					return([gene_name, "CNVgain", "CNVgain", self.copy_number_annotation_gene_to_loss_or_gain_to_sens_res[norm_gene_name]["CNVgain"]])
				else:
					return([gene_name, "CNVgain", "CNVgain", ""])
			else:
				return(["The given alteration was not suitable for testing, please refer to the method description", "" , "" ,""])
		
		else:
			if norm_alteration == "loss":
				return([gene_name, "CNVloss", "CNVloss", ""])
			elif norm_alteration == "gain":
				return([gene_name, "CNVgain", "CNVgain", ""])
			
			else:
				return(["The given alteration was not suitable for testing, please refer to the method description", "" , "" ,""])
		
		
		
		
		
		
		
