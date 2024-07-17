# Creates log odd score matrix as a dictionary 
# input: fasta file and JSON formattted PWM data for BOTH splice acceptor and donor site

'''
How to Prepare the JSON file 
	1. should first produce PWM for a genome using "weighted_pwn.py" in TRANSFAC
	2. Store output of step 1 in file, then feed output file into 82transfac.py through 
	the command line argument, running the program will convert it to JSON and store
	output in a file (this is the file you should feed into this program)
'''

import mcb185
import sys
import json
import math
import gzip

if len(sys.argv) == 4:
	fasta = sys.argv[1]
	gff = sys.argv[2]
	pwm_json = sys.argv[3]
if len(sys.argv) == 6:
	fasta = sys.argv[1]
	test_fasta = sys.argv[2] # fake fasta to test window log score function
	test_gff = sys.argv[3]
	gff = sys.argv[4]
	pwm_json = sys.argv[5]


# initializes dictionary to store frequencies for each base 
def make_dict(file):
	base_count_dict = {
		'A': 0, 'C': 0, 'G': 0, 'T': 0
	}

	# return the dictionary 
	return base_count_dict

# go through the fasta file line by line
# stores occurance of each base in the whole genome
# returns dictionary
def genome_count(fasta):
	# initialize dictionary
	dictionary = make_dict(fasta)

	# iterate through each seq portion of the fasta file
	for defline, seq in mcb185.read_fasta(fasta):
		for nt in seq:
			nt = nt.upper()
			if nt == 'A':
				dictionary["A"] += 1
			if nt == 'C':
				dictionary["C"] += 1
			if nt == 'G':
				dictionary["G"] += 1
			if nt == 'T':
				dictionary["T"] += 1
	
	# dictionary should be updated now with all the nucleotide counts

	# return the dictionary 
	return dictionary

# converts json files into python data type
def json_reader(json_file):
	with open(json_file) as file:
		# load JSON data from file 
		data = json.load(file)

	# return 
	# data type is a bit confusing 
	# see what it prints for more info
	return data


# function stores log-odds ratio using genome base frequency calculated
# in genome_count() function 
# it also takes in the json path to get the pwm frequency as well.

def log_odds_matrix(fasta, pwm_json):

	# calculate the background frequency as counts
	nt_count_dict = genome_count(fasta)

	# sum of the counts of all nucleotides
	count_sum = nt_count_dict['A'] + nt_count_dict['C'] + nt_count_dict['G'] + nt_count_dict['T']

	# frequency in fraction as dictionary
	# use the background frequencies here for the log odds calculation
	nt_frequency_dict = {
		'A': nt_count_dict['A']/count_sum, 
		'C': nt_count_dict['C']/count_sum, 
		'G': nt_count_dict['G']/count_sum, 
		'T': nt_count_dict['T']/count_sum
	}


	# convert the JSON file to python data type 
	pwm_json_acceptor_list, pwm_json_donor_list = json_reader(pwm_json)


	# extract just the pwm dictionary from the JSON
	# donor and acceptor site data mashed up? should check 
	pwm_json_acceptor = pwm_json_acceptor_list['pwm']
	pwm_json_donor = pwm_json_donor_list['pwm']

	# there will be 4 log-odds scores associated with a single position in the sequence 
	# 4 for each base ATGC

	# need both overall genome nt frequency AND PWM frequencies of splice 
	# donor and acceptors

	# we therefore, need two separate dictionaries.
	# one for splice donor and one for splice acceptor

	# initialize log odds dictionary for donor
	log_odds_donor_dict = {}
	# initialize for the acceptor as well
	log_odds_acceptor_dict = {}

	# for loop to create the log odds dictionaries
	# can probably abstract (don't need two separate for loops for donor and acceptor)

	# get the length of the pwm list to find how many bases the PWM tracks
	donor_len = len(pwm_json_donor)
	acceptor_len = len(pwm_json_acceptor)
	
	# donor
	for i in range(donor_len):
		# create a dictionary for the base position
		log_odds_donor_dict[i+1] = {}
		
		# use in for loop
		bases = 'ACGT' 
		base_sum = 0

		# base occurance sum for base position i 
		for base in bases:
			base_sum += int(pwm_json_donor[i][base])

		# second for loop because we need to add 4 key-values(for every base
		for base in bases:
			# convert # of occurance of base -> frequency of base
			freq = int(pwm_json_donor[i][base]) / base_sum

			# make new key value pair
			log_odds_donor_dict[i+1][base] = math.log10(freq / nt_frequency_dict[base])
		
	# do this for the acceptor too
	for i in range(acceptor_len):
		# create a dictionary for the base position
		log_odds_acceptor_dict[i+1] = {}
		
		# use in for loop
		bases = 'ACGT' 
		base_sum = 0

		# base occurance sum for base position i 
		for base in bases:
			base_sum += int(pwm_json_acceptor[i][base])

		# second for loop because we need to add 4 key-values(for every base
		for base in bases:
			# convert # of occurance of base -> frequency of base
			freq = int(pwm_json_acceptor[i][base]) / base_sum

			# make new key value pair
			log_odds_acceptor_dict[i+1][base] = math.log10(freq / nt_frequency_dict[base])
	
	# return log odds matrices
	return log_odds_donor_dict, log_odds_acceptor_dict


# produce log odds matrix for donor and acceptor
log_odds_donor, log_odds_acceptor = log_odds_matrix(fasta, pwm_json)


# window through the whole genome for windows
# can do each chromosome at a time 
# need fasta genome file, log_odds scoring matrices for donor+acceptor
def window_log_score(fasta, log_odds_donor, log_odds_acceptor):
	# dictionary to store all of the log odds scores of all windows 
	log_odds_all_window_scores_donor = {}
	log_odds_all_window_scores_acceptor = {}

	for defline, seq in mcb185.read_fasta(fasta):
		# add a new dictionary for each chromosome/defline 
		# the key should be the defline, and value is empty dictionary
		log_odds_all_window_scores_donor[defline] = {}
		log_odds_all_window_scores_acceptor[defline] = {}

		# for each fasta seq...
		# window
		w_donor = len(log_odds_donor)
		w_acceptor = len(log_odds_acceptor)

		# score donor log odds
		# window for all possible donor sites
		for i in range(len(seq) -w_donor +1):
			# extract current window seq
			window_seq = seq[i:i+w_donor].upper()

			# using log_odds_donor dictionary to score the window
			# initalize score var
			log_odds_score = 0
			
			# for each nt in the window_seq, score it, then add it to the sum score
			for j, base in enumerate(window_seq):
				# i is the base position (0 through 5)
				# calling the dictionary with key = [base] will
				# call the appropriate score
				log_odds_score += log_odds_donor[j+1][base]

			# we need a way to store all of the log odds score
			# can create a dictionary that has index of the first base (of the window)
			# as the key
			# the value will be the log odds score for the window
			log_odds_all_window_scores_donor[defline][i+1] = log_odds_score

		# score acceptor log odds
		# window for all possible acceptor sites
		for i in range(len(seq) -w_acceptor +1):
			# extract current window seq
			window_seq = seq[i:i+w_acceptor].upper()

			# using log_odds_donor dictionary to score the window
			# initalize score var
			log_odds_score = 0
			
			# for each nt in the window_seq, score it, then add it to the sum score
			for j, base in enumerate(window_seq):
				# i is the base position (0 through 5)
				# calling the dictionary with key = [base] will
				# call the appropriate score
				log_odds_score += log_odds_acceptor[j+1][base]

			# we need a way to store all of the log odds score
			# can create a dictionary that has index of the first base (of the window)
			# as the key
			# the value will be the log odds score for the window
			log_odds_all_window_scores_acceptor[defline][i+1] = log_odds_score

	# return 
	return log_odds_all_window_scores_donor, log_odds_all_window_scores_acceptor


# store all donor and acceptor scores into variables (as a dictionary)

# WHEN RUNNING FOR 1PCT GENOME
donor_scores, acceptor_scores = window_log_score(fasta, log_odds_donor, log_odds_acceptor)

# TESTING STAGE
# donor_scores, acceptor_scores = window_log_score(test_fasta, log_odds_donor, log_odds_acceptor)
#print(donor_scores, acceptor_scores, sep='\n\n\n')

# Now that we have log odds scores for all of the possible windows...
# We need to store all of the coordinates of the windows that are the acceptor and donors

# iterate through gff, find introns, store donor and acceptor
# acceptor_base_length is the # of bases you are taking the PWM for the acceptor site
# can limit which chrom it extracts from with 'lim_chrom' (ex. ['I', 'III'])
# can also limit total # extracted by specifying range 'lim_num' (ex. ['1 500', '1 1000'])
	# use a space char to separate start and end coords (1=based index)
# length of chrom and abs# list MUST MATCH
def extract_positive_coordinates(fasta, gff, acceptor_base_length, d_lim_chrom, d_lim_num, a_lim_chrom, a_lim_num):
	# initialize donor and acceptor coordinate dictionaries
	donor_coords = {}
	acceptor_coords = {}

	# the log odd scores for all windows are stored as 1-based index
		# NOTE: log odds scores are all on the positive strand AND they are stored in
		# different dictionaries for every chromosome

	# create separate lists for every chromosome/defline
	for defline, seq in mcb185.read_fasta(fasta):
		defline_words = defline.split()
		donor_coords[defline_words[0]] = []
		acceptor_coords[defline_words[0]] = []

	# iterate through gff to store
	with gzip.open(gff, 'rt') as file:
		for line in file:
			words = line.split()
			if words[1] != 'RNASeq_splice': continue
			# extract intron information 
			chrom = words[0]
			start = int(words[3])
			end = int(words[4])
			strand = words[6]

			# find window index for donor and acceptor
			donor_index = start
			acceptor_index = end - acceptor_base_length + 1

			# only collecting + strand for now
			if strand != '+': continue

			# for loop to store in the correct dictionary
			for key in donor_coords.keys():
				if key == chrom:
					donor_coords[key].append(donor_index)
					acceptor_coords[key].append(acceptor_index)

	# TESTING
	'''
	print('before limiting')
	for key, num in zip(lim_chrom, lim_num):
		print(key)
		print(len(donor_coords[key]))
		print(len(acceptor_coords[key]))
		print()
	'''

	# create the chrom limited d and a pos dict if requested
	# donor
	if d_lim_chrom != 'NA':
		limited_pos_donor_coords = {}

		for defline, seq in mcb185.read_fasta(fasta):
			defline_words = defline.split()
			limited_pos_donor_coords[defline_words[0]] = []
		# donor
		for key in donor_coords.keys():
			if key not in d_lim_chrom: continue
			# only take the coords if it is the chrom requested
			limited_pos_donor_coords[key] = donor_coords[key]
	# acceptor
	if a_lim_chrom != 'NA':
		limited_pos_acceptor_coords = {}

		for defline, seq in mcb185.read_fasta(fasta):
			defline_words = defline.split()
			limited_pos_acceptor_coords[defline_words[0]] = []
		# acceptor 
		for key in acceptor_coords.keys():
			if key not in a_lim_chrom: continue
			# only take the coords if it is the chrom requested
			limited_pos_acceptor_coords[key] = acceptor_coords[key]

	# if abs number of coord limitation is also requested...
	# donor
	if d_lim_num != 'NA':
		# donor
		for key, num in zip(d_lim_chrom, d_lim_num):
			# first need to get rid of repeats
			limited_pos_donor_coords[key] = list(set(limited_pos_donor_coords[key]))


			# slice appropriately
			words = num.split()
			# when need only one range
			if len(words) == 2:
				start = int(words[0]) - 1 
				end = int(words[1])
				limited_pos_donor_coords[key] = limited_pos_donor_coords[key][start:end]
			
			# when need two ranges 
			if len(words) == 4:
				start1 = int(words[0]) - 1
				end1 = int(words[1])
				start2 = int(words[2]) - 1
				end2 = int(words[3])
				# coord list 1
				coordlist1 = limited_pos_donor_coords[key][start1:end1]
				# coord list 2
				coordlist2 = limited_pos_donor_coords[key][start2:end2]

				# join both lists
				limited_pos_donor_coords[key] = coordlist1 + coordlist2

	# if abs number of coord limitation is also requested...
	# acceptor
	if a_lim_num != 'NA':
		# acceptor
		for key, num in zip(a_lim_chrom, a_lim_num):
			# first need to get rid of repeats
			limited_pos_acceptor_coords[key] = list(set(limited_pos_acceptor_coords[key]))


			# slice appropriately
			words = num.split()
			# when need only one range
			if len(words) == 2:
				start = int(words[0]) - 1 
				end = int(words[1])
				limited_pos_acceptor_coords[key] = limited_pos_acceptor_coords[key][start:end]
			
			# when need two ranges 
			if len(words) == 4:
				start1 = int(words[0]) - 1
				end1 = int(words[1])
				start2 = int(words[2]) - 1
				end2 = int(words[3])
				# coord list 1
				coordlist1 = limited_pos_acceptor_coords[key][start1:end1]
				# coord list 2
				coordlist2 = limited_pos_acceptor_coords[key][start2:end2]

				# join both lists
				limited_pos_acceptor_coords[key] = coordlist1 + coordlist2


	# TESTING
	'''
	for key, num in zip(lim_chrom, lim_num):
		print(key)
		print(len(limited_pos_donor_coords[key]))
		print(len(limited_pos_acceptor_coords[key]))
		print()
	'''

	# return limited versions
	if d_lim_chrom != 'NA' and a_lim_chrom != 'NA':
		return limited_pos_donor_coords, limited_pos_acceptor_coords

	if d_lim_chrom == 'NA' and a_lim_chrom != 'NA':
		return donor_coords, limited_pos_acceptor_coords

	if d_lim_chrom != 'NA' and a_lim_chrom == 'NA':
		return limited_pos_donor_coords, acceptor_coords

	return donor_coords, acceptor_coords

# FOR REAL GENOME
'''
d_lim_chrom = ['I', 'II', 'III', 'IV', 'V', 'X']
d_lim_num = ['1 256 385 512', '1 254 382 508', '1 120 181 240', '1 194 292 388', '1 270 406 540', '1 174 262 348']
a_lim_chrom = ['I', 'II', 'III', 'IV', 'V', 'X']
a_lim_num = ['1 238 358 476', '1 254 382 508', '1 122 184 244', '1 174 262 348', '1 244 367 488', '1 160 241 320']
'''


d_lim_chrom = 'NA'
d_lim_num = 'NA'
a_lim_chrom = 'NA'
a_lim_num = 'NA'


donor_coordinates, acceptor_coordinates = extract_positive_coordinates(fasta, gff, 25, d_lim_chrom, d_lim_num, a_lim_chrom, a_lim_num)


# TESTING
'''
print('positive coords')
for key in ['I', 'II', 'III', 'IV', 'V']:
	print(f'donor {key} {len(donor_coordinates[key])}')
	print(f' acceptor {key} {len(acceptor_coordinates[key])}')
print()
print()
'''




# FOR TESTING
# donor_coordinates, acceptor_coordinates = extract_positive_coordinates(test_fasta, test_gff, 25)
#print(donor_coordinates, acceptor_coordinates, sep='\n\n\n')

# Now we want to collect all of the negatives
# use the positive acceptor and donor coordinates to look at all windows that ARE NOT
# acceptor and donor sites -> if those windows have a GT or AG, then add it to the negative
# coords list

# NOTE: We have to make sure that we check for the 9 bases and 25 bases for both sites

# inputs: donor and acceptor coord dictionary of lists, fasta, and
# donor and acceptor base length for the PWM
# can limit which chrom it extracts from with 'lim_chrom' (ex. ['I', 'III'])
# can also limit total # extracted by specifying range 'lim_num' (ex. ['1 500', '1 1000'])
	# use a space char to separate start and end coords (1=based index)
# length of chrom and abs# list MUST MATCH
def extract_negative_coords(donor_coordinates, acceptor_coordinates, fasta, donor_base_length, acceptor_base_length, d_lim_chrom, d_lim_num, a_lim_chrom, a_lim_num):
	# initialize donor and acceptor dictionaries of lists 
	negative_donor_coords = {}
	negative_acceptor_coords = {} 

	# create separate lists for every chromosome/defline (ie. I, MtDNA, etc)
	for defline, seq in mcb185.read_fasta(fasta):
		defline_words = defline.split()
		negative_donor_coords[defline_words[0]] = []
		negative_acceptor_coords[defline_words[0]] = []

	# iterate one chromosome at a time 
	for defline, seq in mcb185.read_fasta(fasta):
		defline_words = defline.split()
		chrom = defline_words[0]

		# extract the appropriate donor and acceptor coordinates list for the current
		# chromosome
		# Note: these coordinates are 1-based index 
		current_donor_list = donor_coordinates[chrom]
		current_acceptor_list = acceptor_coordinates[chrom]

		# window through each base in the seq
		# donor window first
		for i in range(len(seq) - donor_base_length +1):
			# fix to 1-based index
			base_number = i + 1

			# if base number matches with any of the donor coordinates, skip
			if base_number in current_donor_list: continue

			# if not in donor list, check if first 2 bases are GT
			# if yes, store as negative coord

			if seq[base_number - 1:base_number + 1].upper() == 'GT':
				negative_donor_coords[chrom].append(base_number)

		# after finished with windowing through the donor for the particular chrom,
		# window through the negative
		for i in range(len(seq) - acceptor_base_length +1):
			if seq[i+acceptor_base_length-2:i+acceptor_base_length] == 'AG':
				if i+1 not in current_acceptor_list:
					negative_acceptor_coords[chrom].append(i+1)

		
	# create the chrom limited d and a negative dict if requested
	# donor
	if d_lim_chrom != 'NA':
		limited_neg_donor_coords = {}

		for defline, seq in mcb185.read_fasta(fasta):
			defline_words = defline.split()
			limited_neg_donor_coords[defline_words[0]] = []
		for key in negative_donor_coords.keys():
			if key not in d_lim_chrom: continue
			# only take the coords if it is the chrom requested
			limited_neg_donor_coords[key] = negative_donor_coords[key]

	# acceptor
	if a_lim_chrom != 'NA':
		limited_neg_acceptor_coords = {}

		for defline, seq in mcb185.read_fasta(fasta):
			defline_words = defline.split()
			limited_neg_acceptor_coords[defline_words[0]] = []
		# acceptor 
		for key in negative_acceptor_coords.keys():
			if key not in a_lim_chrom: continue
			# only take the coords if it is the chrom requested
			limited_neg_acceptor_coords[key] = negative_acceptor_coords[key]

	# if abs number of coord limitation is also requested...
	# donor
	if d_lim_num != 'NA':
		# donor
		for key, num in zip(d_lim_chrom, d_lim_num):
			# first need to get rid of repeats
			limited_neg_donor_coords[key] = list(set(limited_neg_donor_coords[key]))


			# slice appropriately
			words = num.split()
			# when need only one range
			if len(words) == 2:
				start = int(words[0]) - 1 
				end = int(words[1])
				limited_neg_donor_coords[key] = limited_neg_donor_coords[key][start:end]
			
			# when need two ranges 
			if len(words) == 4:
				start1 = int(words[0]) - 1
				end1 = int(words[1])
				start2 = int(words[2]) - 1
				end2 = int(words[3])
				# coord list 1
				coordlist1 = limited_neg_donor_coords[key][start1:end1]
				# coord list 2
				coordlist2 = limited_neg_donor_coords[key][start2:end2]

				# join both lists
				limited_neg_donor_coords[key] = coordlist1 + coordlist2

	# if abs number of coord limitation is also requested...
	# acceptor
	if a_lim_num != 'NA':
		# acceptor
		for key, num in zip(a_lim_chrom, a_lim_num):
			# first need to get rid of repeats
			limited_neg_acceptor_coords[key] = list(set(limited_neg_acceptor_coords[key]))


			# slice appropriately
			words = num.split()
			# when need only one range
			if len(words) == 2:
				start = int(words[0]) - 1 
				end = int(words[1])
				limited_neg_acceptor_coords[key] = limited_neg_acceptor_coords[key][start:end]
			
			# when need two ranges 
			if len(words) == 4:
				start1 = int(words[0]) - 1
				end1 = int(words[1])
				start2 = int(words[2]) - 1
				end2 = int(words[3])
				# coord list 1
				coordlist1 = limited_neg_acceptor_coords[key][start1:end1]
				# coord list 2
				coordlist2 = limited_neg_acceptor_coords[key][start2:end2]

				# join both lists
				limited_neg_acceptor_coords[key] = coordlist1 + coordlist2


	# return limited versions
	if d_lim_chrom != 'NA' and a_lim_chrom != 'NA':
		return limited_neg_donor_coords, limited_neg_acceptor_coords

	if d_lim_chrom == 'NA' and a_lim_chrom != 'NA':
		return negative_donor_coords, limited_neg_acceptor_coords

	if d_lim_chrom != 'NA' and a_lim_chrom == 'NA':
		return limited_neg_donor_coords, negative_acceptor_coords

	# return un-limited version otherwise
	return negative_donor_coords, negative_acceptor_coords

# FOR REAL GENOME (with chrom+abs # limits)

'''
d_lim_chrom = ['I', 'II', 'III', 'IV', 'V', 'X']
d_lim_num = ['129 512', '128 508', '61 240', '98 388', '136 540', '88 348']
a_lim_chrom = ['I', 'II', 'III', 'IV', 'V', 'X']
a_lim_num = ['120 476', '128 508', '62 244', '88 348', '123 488', '81 320']
'''

neg_donor_coordinates, neg_acceptor_coordinates =  extract_negative_coords(donor_coordinates, acceptor_coordinates, fasta, 9, 25, d_lim_chrom, d_lim_num, a_lim_chrom, a_lim_num)


# TESTING
'''
print('neg coords')
for key in ['I', 'II', 'III', 'IV', 'V']:
	print(f'donor {key} {len(neg_donor_coordinates[key])}')
	print(f' acceptor {key} {len(neg_acceptor_coordinates[key])}')
'''

# FOR TESTING 
# LIMITING NEGATIVES TO CHROM 'I'
# neg_donor_coordinates, neg_acceptor_coordinates =  extract_negative_coords(donor_coordinates, acceptor_coordinates, test_fasta, 9, 25, 'I')
#print(neg_donor_coordinates, neg_acceptor_coordinates, sep='\n\n\n')


# Next, take the dictionary of lists of positive and negative coords for donor and acceptor
# regions, and combine the data with the the log odds dictionary

# NOTE: The log odds dictionary has the entire defline as the dictionary key, so watch out
# the postive and negative coordinate dictinaries have only the chromosome roman numeral as
# the key

# function combines the positive and negative data with the log odds data
# stores it as a dictionary
def combine_log_odds_and_pos_and_neg(log_odds_donor_dict, log_odds_accept_dict, neg_donor_coords, neg_accept_coords, pos_donor_coords, pos_accept_coords):
	# initialize new dictionary
	combined_donor_dict = {}
	combined_accept_dict = {}

	# for every key in the log odds dict, create a key in the combined_dict with
	# a dictionary for the value

	# donor
	for key in log_odds_donor_dict.keys():
		combined_donor_dict[key] = {}

	# acceptor
	for key in log_odds_accept_dict.keys():
		combined_accept_dict[key] = {}

	# transfer over the log odds data first
	# donor
	for key in combined_donor_dict.keys():
		# call the relevant dictionary the current chromosome
		current_log_odds_dict = log_odds_donor_dict[key]

		# iterate through the current chrom log odds dict
		for window in current_log_odds_dict.keys():
			log_odds = log_odds_donor_dict[key][window]
			combined_donor_dict[key][window] = {
				'log_odds': log_odds,
				'class': 'NA'
			}

	# next transfer over acceptor
	for key in combined_accept_dict.keys():
		# call the relevant dictionary the current chromosome
		current_log_odds_dict = log_odds_accept_dict[key]

		# iterate through the current chrom log odds dict
		for window in current_log_odds_dict.keys():
			log_odds = log_odds_accept_dict[key][window]
			combined_accept_dict[key][window] = {
				'log_odds': log_odds,
				'class': 'NA'
			}

	# next transfer over the positive and negative data 
	# the new dictionary is basically set, just need to change the value of the 'class' key
	# value pair

	# start with the donor (positive and negative)
	# iterate through combined donor dictionary
	for key in combined_donor_dict.keys():
		key_words = key.split()
		chrom = key_words[0]

		# store current chromosome dictionary
		current_dict = combined_donor_dict[key]

		# iterate through all of the dict keys (the window numbers)
		for window in current_dict.keys():
			# we need to check if the current window number is contained in the
			# pos_donor_coords
			if window in pos_donor_coords[chrom]:
				# if window is contained, label the window as '+'
				combined_donor_dict[key][window]['class'] = '+'

			# label the negatives too
			if window in neg_donor_coords[chrom]:
				# label as '-'
				combined_donor_dict[key][window]['class'] = '-'

	# then do the acceptor (pos and neg)
	for key in combined_accept_dict.keys():
		key_words = key.split()
		chrom = key_words[0]

		# store current chromosome dictionary
		current_dict = combined_accept_dict[key]

		# iterate through all of the dict keys (the window numbers)
		for window in current_dict.keys():
			# we need to check if the current window number is contained in the
			# pos_accept_coords
			if window in pos_accept_coords[chrom]:
				# if window is contained, label the window as '+'
				combined_accept_dict[key][window]['class'] = '+'

			# label the negatives too
			if window in neg_accept_coords[chrom]:
				# label as '-'
				combined_accept_dict[key][window]['class'] = '-'

	# return 
	return combined_donor_dict, combined_accept_dict

# FOR REAL GENOME
donor, acceptor = combine_log_odds_and_pos_and_neg(donor_scores, acceptor_scores, neg_donor_coordinates, neg_acceptor_coordinates, donor_coordinates, acceptor_coordinates)

# convert to JSON
# first put dictionaries into list
json_array = [donor, acceptor]

# print JSON
print(json.dumps(json_array, indent=4))


'''
# CHECKING ACCURACY OF DONOR POSITIVES
# just loop through all of the donor positives and check if it actually is positive

# checking all donor positives in first chromosome 

donor_list = donor_coordinates['I']
pos_count = 0
neg_count = 0
for coord in donor_list:
	if donor['I 1..150725'][coord]['class'] == '+':
		pos_count += 1
	else: neg_count += 1

print(pos_count) # 907
print(neg_count) # 0



# checking donor positives in second chromosome
donor_list = donor_coordinates['II']
pos_count = 0
neg_count = 0
for coord in donor_list:
	if donor['II 1..152795'][coord]['class'] == '+':
		pos_count += 1
	else: neg_count += 1

print(pos_count) # 853
print(neg_count) # 0



# checking donor positives in chromosome V
# checking donor positives in second chromosome
donor_list = donor_coordinates['V']
pos_count = 0
neg_count = 0
for coord in donor_list:
	if donor['V 1..209242'][coord]['class'] == '+':
		pos_count += 1
	else: neg_count += 1

print(pos_count) # 872
print(neg_count) # 0


# CHECKING ACCURACY OF DONOR NEGATIVES

# checking chromosome I donor negatives

neg_donor_list = neg_donor_coordinates['I']
pos_count = 0
neg_count = 0
for coord in neg_donor_list:
	if donor['I 1..150725'][coord]['class'] == '+':
		pos_count += 1
	if donor['I 1..150725'][coord]['class'] == '-': 
		neg_count += 1

print(pos_count) # 0 
print(neg_count) # 6451


# checking chromosome II donor negatives

neg_donor_list = neg_donor_coordinates['II']
pos_count = 0
neg_count = 0
for coord in neg_donor_list:
	if donor['II 1..152795'][coord]['class'] == '+':
		pos_count += 1
	if donor['II 1..152795'][coord]['class'] == '-': 
		neg_count += 1

print(pos_count) # 0 
print(neg_count) # 7202


# checking chromosome V donor negatives
neg_donor_list = neg_donor_coordinates['V']
pos_count = 0
neg_count = 0
for coord in neg_donor_list:
	if donor['V 1..209242'][coord]['class'] == '+':
		pos_count += 1
	if donor['V 1..209242'][coord]['class'] == '-': 
		neg_count += 1

print(pos_count) # 0
print(neg_count) # 8774

# CHECKING ACCURACY FOR ACCEPTOR POSITIVE

# check chromosome I

accept_list = acceptor_coordinates['I']
pos_count = 0
neg_count = 0
for coord in accept_list:
	if acceptor['I 1..150725'][coord]['class'] == '+':
		pos_count += 1
	if acceptor['I 1..150725'][coord]['class'] == '-': 
		neg_count += 1

print(pos_count) # 907
print(neg_count) # 0


# check chromosome II

accept_list = acceptor_coordinates['II']
pos_count = 0
neg_count = 0
for coord in accept_list:
	if acceptor['II 1..152795'][coord]['class'] == '+':
		pos_count += 1
	if acceptor['II 1..152795'][coord]['class'] == '-': 
		neg_count += 1

print(pos_count) # 853
print(neg_count) # 0


# check chromosome V

accept_list = acceptor_coordinates['V']
pos_count = 0
neg_count = 0
for coord in accept_list:
	if acceptor['V 1..209242'][coord]['class'] == '+':
		pos_count += 1
	if acceptor['V 1..209242'][coord]['class'] == '-': 
		neg_count += 1

print(pos_count) # 872
print(neg_count) # 0


# CHECKING ACCURACY FOR NEGATIVE ACCEPTORS

# check chromosome I

neg_accept_list = neg_acceptor_coordinates['I']
pos_count = 0
neg_count = 0
for coord in neg_accept_list:
	if acceptor['I 1..150725'][coord]['class'] == '+':
		pos_count += 1
	if acceptor['I 1..150725'][coord]['class'] == '-': 
		neg_count += 1

print(pos_count) # 0
print(neg_count) # 6939


# check chromosome II

neg_accept_list = neg_acceptor_coordinates['II']
pos_count = 0
neg_count = 0
for coord in neg_accept_list:
	if acceptor['II 1..152795'][coord]['class'] == '+':
		pos_count += 1
	if acceptor['II 1..152795'][coord]['class'] == '-': 
		neg_count += 1

print(pos_count) # 0
print(neg_count) # 7272


# check chromosome V

neg_accept_list = neg_acceptor_coordinates['V']
pos_count = 0
neg_count = 0
for coord in neg_accept_list:
	if acceptor['V 1..209242'][coord]['class'] == '+':
		pos_count += 1
	if acceptor['V 1..209242'][coord]['class'] == '-': 
		neg_count += 1

print(pos_count) # 0
print(neg_count) # 9564



# Is there overlap in coordinates in the pos and negative acceptor coord lists?

match_count = 0
for coords in acceptor_coordinates['I']:
	for coordinate in neg_acceptor_coordinates['I']:
		if coords == coordinate:
			match_count += 1

print(match_count) # 821


# take a look at which sequences matched

match_count = 0
for coords in acceptor_coordinates['I']:
	for coordinate in neg_acceptor_coordinates['I']:
		if coords == coordinate:
			match_count += 1
			print(coords)


# print some of these sequences
# 3640 was a match on chromosome I

for defline, seq in mcb185.read_fasta(fasta):
	defline_words = defline.split()
	chrom = defline_words[0]
	if chrom == 'I':
		print(seq[3639:3639+26]) # ACTCCCCCTTTAACAACCACCCGAGG

# sequence does end in a 'AG' so qualifies for negative
# however, since it is also a pos, it should not MATCH
# maybe problem with 0-based vs. 1 based indexing?

'''

