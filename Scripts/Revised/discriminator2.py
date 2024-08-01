import sys
import mcb185
import json
import math
import random
import numpy as np
import matplotlib.pyplot as plt

genome_fasta = sys.argv[1]
d_fasta = sys.argv[2]
a_fasta = sys.argv[3]
n_d_fasta = sys.argv[4]
n_a_fasta = sys.argv[5]
pwm_json = sys.argv[6]

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

# produce log odds scoring matrix for donor and acceptor
log_odds_donor, log_odds_acceptor = log_odds_matrix(genome_fasta, pwm_json)


class Seq:
	def __init__(self, sequence, log_odd, classification, d_or_n):
		self.sequence = sequence
		self.log_odd = log_odd
		self.classification = classification # '+' or '-'
		self.d_or_n = d_or_n # 'd' or 'n'

	# returns TP, TN, FP, FN
	def discriminate(self, threshold):
		if self.log_odd >= threshold:
			if self.classification == '+': return 'TP'
			else: return 'FP'
		if self.log_odd < threshold:
			if self.classification == '-': return 'TN'
			else: return 'FN'


# classification = '+' or '-'
# d_or_n = 'd' or 'n'
# fasta = d_fasta, a_fasta, etc...
def fasta_to_list(fasta, classification, d_or_n):
	the_list = []
	for defline, seq in mcb185.read_fasta(fasta):
		log_odd = 0
		for i, base in enumerate(seq):
			if d_or_n == 'd':
				log_odd += log_odds_donor[i+1][base]
			if d_or_n == 'a':
				log_odd += log_odds_acceptor[i+1][base]
		the_list.append(Seq(seq, log_odd, classification, d_or_n))
	return the_list

all_d = fasta_to_list(d_fasta, '+', 'd') + fasta_to_list(n_d_fasta, '-', 'd')
all_a = fasta_to_list(a_fasta, '+', 'a') + fasta_to_list(n_a_fasta, '-', 'a')

random.shuffle(all_d)
random.shuffle(all_a)

# split data 
training_ratio = 0.7 
validation_ratio = 0.15
test_ratio = 0.15

d_train_size = int(len(all_d) * training_ratio)
d_validation_size = int(len(all_d) * validation_ratio)
d_test_size = len(all_d) - d_train_size - d_validation_size

a_train_size = int(len(all_a) * training_ratio)
a_validation_size = int(len(all_a) * validation_ratio)
a_test_size = len(all_a) - a_train_size - a_validation_size

# the split data
d_train = all_d[:d_train_size]
d_validation = all_d[d_train_size:d_train_size+d_validation_size]
d_test = all_d[d_train_size+d_validation_size:]

a_train = all_a[:a_train_size]
a_validation = all_a[a_train_size:a_train_size+a_validation_size]
a_test = all_a[a_train_size+a_validation_size:]



# the_list = training, validation, or test list of Seq objects
# return TP, TN, FP, FN 
def discriminator(the_list, threshold):
	tp = 0
	tn = 0
	fp = 0
	fn = 0
	epsilon = 1e-10

	for seq in the_list:
		if seq.discriminate(threshold) == 'TP': tp += 1
		if seq.discriminate(threshold) == 'TN': tn += 1
		if seq.discriminate(threshold) == 'FP': fp += 1
		if seq.discriminate(threshold) == 'FN': fn += 1

	accuracy = (tp+tn) / (tp+tn+fp+fn+epsilon)

	return accuracy

def graph(list1, list2, upper, lower, step, title):
	if list2 == None:
		thres_list = []
		acc_list = []
		for i in np.arange(lower, upper, step):
			acc_list.append(discriminator(list1, i))
			thres_list.append(i)

		plt.plot(thres_list, acc_list, marker='o')
		plt.title(title)
		plt.xlabel('Threshold')
		plt.ylabel('Accuracy')
		plt.ylim([0,1])
		plt.show()
	else:
		t_list = []
		acc_list1 = []
		acc_list2 = []

		for i in np.arange(lower, upper, step):
			acc_list1.append(discriminator(list1, i))
			acc_list2.append(discriminator(list2, i))
			t_list.append(i)

		plt.plot(t_list, acc_list1, marker='o', label='Training')
		plt.plot(t_list, acc_list2, marker='o', label='Validation')

		plt.title(title)
		plt.xlabel("Threshold")
		plt.ylabel("Accuracy")
		plt.ylim([0,1])
		plt.legend()
		plt.show()


def best_thres(the_list, upper, lower, step):
	best_accuracy = 0
	best_thres = 0
	for i in np.arange(lower, upper, step):
		if discriminator(the_list, i) > best_accuracy:
			best_accuracy = discriminator(the_list, i)
			best_thres = i
	return best_thres, best_accuracy


# print(best_thres(d_train, 10, -10, 0.01))

# graph(d_train, 10, -10, 0.1, "Donor Training Set Accuracy")

graph(d_train, d_validation, 10, -10, 0.1, "Donor Train and Validation Set Accuracy Over Thresholds")





