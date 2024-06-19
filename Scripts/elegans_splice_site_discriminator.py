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

fasta = sys.argv[1]
pwm_json = sys.argv[2]

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
	# starting out with the donor (can probably abstract this into a function
	# later)
	
	for i in range(6):
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
	for i in range(7):
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


print(log_odds_donor)
print()
print(log_odds_acceptor)
print()


# window through the whole genome for 6 and 7 base windows
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
		w_donor = 6
		w_acceptor = 7

		# score donor log odds
		# window for all possible donor sites
		for i in range(len(seq) -w_donor +1):
			# extract current window seq
			window_seq = seq[i:i+w_donor]

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
			window_seq = seq[i:i+w_acceptor]

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
donor_scores, acceptor_scores = window_log_score(fasta, log_odds_donor, log_odds_acceptor)

'''
print(donor_scores)
print()
print(acceptor_scores)
'''


# want to check if output is correct
# manually calculate first number 
	# need to go to fasta file, look at first 6 bases, 
	# print out the log odds score matrix for the donor
		# based on it score the first 6 bases and compare with result


