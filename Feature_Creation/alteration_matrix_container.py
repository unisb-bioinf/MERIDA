import sys

class Alteration_Matrix_Container:

	def __init__(self, available_cl_list):
		
		self.available_cl = available_cl_list
		self.all_mutations = []
		self.all_cnvs = []
		self.sample_to_mutation_to_effect_and_implication = {}
		self.sample_to_cnv_to_effect_and_implication = {}
		
	#Registers mutation in a sample
	#Alteration_name consists of gene name and the alteration, mutation_effect and implication can be empty strings
	def register_mutation(self, sample, alteration_name, mutation_effect, implication):
		
		if not sample in self.sample_to_mutation_to_effect_and_implication:
			
			self.sample_to_mutation_to_effect_and_implication[sample] = {alteration_name:[mutation_effect, implication]}
			
		else:
			
			if not alteration_name in self.sample_to_mutation_to_effect_and_implication[sample]:
				
				self.sample_to_mutation_to_effect_and_implication[sample][alteration_name] = [mutation_effect, implication]
		
		
		already_in = False
		for alteration_element in self.all_mutations:
			if alteration_name == alteration_element[0]:
				already_in = True
				
		if not already_in:
			self.all_mutations.append([alteration_name, mutation_effect, implication])
	
	
	#Register copy number variation in a sample
	#Alteration_name consists of gene name and the alteration, mutation_effect and implication can be empty strings
	def register_cnv(self, sample, alteration_name, mutation_effect, implication):
		
		
		if not sample in self.sample_to_cnv_to_effect_and_implication:
			
			self.sample_to_cnv_to_effect_and_implication[sample] = {alteration_name:[mutation_effect, implication]}
			
		else:
			
			if not alteration_name in self.sample_to_cnv_to_effect_and_implication[sample]:
				
				self.sample_to_cnv_to_effect_and_implication[sample][alteration_name] = [mutation_effect, implication]
		
		already_in = False
		for alteration_element in self.all_cnvs:
			if alteration_name == alteration_element[0]:
				already_in = True
				
		if not already_in:
			self.all_cnvs.append([alteration_name, mutation_effect, implication])
	
	def __print_header(self, output):
		first_row = ["sample_id"]
		second_row = ["mutation_effect"]
		third_row = ["implication"]
		for alteration in self.all_mutations:
			first_row.append(alteration[0])
			second_row.append(alteration[1])
			third_row.append(alteration[2])
		
		for alteration in self.all_cnvs:
			first_row.append(alteration[0])
			second_row.append(alteration[1])
			third_row.append(alteration[2])
			
		output.write("\t".join(first_row) + "\n")
		print(len(first_row))
		output.write("\t".join(second_row) + "\n")
		print(len(second_row))
		output.write("\t".join(third_row) + "\n")
		print(len(third_row))
		
	def print_matrix(self, output_file_name):
		
		#print("There were " + str(len(self.all_mutations)) + " different mutations registered.")
		#print("There were " + str(len(self.all_cnvs)) + " different copy number variations registered.")
		all_rows = []
		
		for sample in self.available_cl:
			one_row = [sample]
			
			for alteration in self.all_mutations:
				
				if not sample in self.sample_to_mutation_to_effect_and_implication:
					one_row.append('0')
					
				else:
					if not alteration[0] in self.sample_to_mutation_to_effect_and_implication[sample]:
						one_row.append('0')
						
					else:
						one_row.append('1')
						
						
			for alteration in self.all_cnvs:
				
				if not sample in self.sample_to_cnv_to_effect_and_implication:
					one_row.append('0')
					
				else:
					if not alteration[0] in self.sample_to_cnv_to_effect_and_implication[sample]:
						one_row.append('0')
						
					else:
						one_row.append('1')
						
			all_rows.append(one_row)
			
			
			
		with open(output_file_name, 'w',  encoding='utf-8')	as output:
			
			self.__print_header(output)
		
			for row in all_rows:
				
				#print(len(row))
				output.write("\t".join(row) + "\n")
				
				
				
				