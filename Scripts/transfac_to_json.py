# 82transfac.py by Kenta Hsu 
import sys
import re
import gzip
import json

# things to store
# prob best store this all in a dictionary
# key: (ie. ID, etc) 

# AC: 1 string (should make this the dictionary key)
# ID: 1 string 
# DE: tuple(list w/ 2 strings), 1 string)
# PO: dictionary w/ keys: number of nts in motif
	# value: dictionary w/ keys: A, C, G, T
		# value: prob each base
# tax_group: 1 string
# tf_family: 1 string
# pubmed: 1 string 
# uniprot..: 1 string 
# data_type: 1 string 

# use 24lines.transfac (first 24 lines of the transfac file)
# read each line at a time 
# startswith()?
# after three 'X' starts create new dictionary

path = sys.argv[1]
with gzip.open(path, 'rt') as fp:
	transfac_list = [] # list of dict for each id
	current_id = '' # track ac
	numbers = False

	for line in fp:
		line = line.rstrip()
		tag = line[:2]
		
		if tag == 'XX':
			numbers = False
		if tag == 'ID':
			current_id = line[3:]
			transfac_list.append({'id': current_id})
		if tag == 'PO':
			transfac_list[-1] = {'id': current_id, 'pwm': []}
			numbers = True
			continue # PWM values starting on next line
		
		# when numbers=True
		if numbers == True:
			# store dictionary for current id temporarily
			temp_dict = transfac_list[-1] 

			# add PWM values to a new dictionary
			line_as_list = line.split()
			new_dict = {
				'A': line_as_list[1], 
				'C': line_as_list[2],
				'G': line_as_list[3],
				'T': line_as_list[4]
			}

			# append the pwm list with each nt in consensus
			temp_dict['pwm'].append(new_dict)

			# change actual dictionary 
			transfac_list[-1] = temp_dict



print(json.dumps(transfac_list, indent=4))


