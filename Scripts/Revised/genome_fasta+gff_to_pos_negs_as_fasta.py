import sys
import gzip
import mcb185
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from Bio import SeqIO

fasta = sys.argv[1]
gff = sys.argv[2]
pick = sys.argv[3] # d, a, n_d, or n_a

# program splits available pos and neg coordinates into train, validation, and test

def extract_positive_coordinates(fasta, gff, acceptor_base_length, d_lim_chrom, d_lim_num, a_lim_chrom, a_lim_num):
	# initialize donor and acceptor coordinate dictionaries
	donor_coords = {}
	acceptor_coords = {}

	# the log odd scores for all windows are stored as 1-based index
		# NOTE: log odds scores are all on the positive strand AND they are stored in
		# different dictionaries for every chromosome

	# create separate lists for every chromosome/defline
	# may want to consider doing gzip.open instead since we don't need use 'seq'
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
			# get rid of repeats
			limited_pos_donor_coords[key] = list(set(limited_pos_donor_coords[key]))

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
			limited_pos_acceptor_coords[key] = list(set(limited_pos_acceptor_coords[key]))


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

def extract_negative_coords(donor_coordinates, acceptor_coordinates, fasta, donor_base_length, acceptor_base_length, d_lim_num, a_lim_num):
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
		if chrom == 'MtDNA': break

		# extract the appropriate donor and acceptor coordinates list for the current
		# chromosome
		# Note: these coordinates are 1-based index 
		current_donor_list = donor_coordinates[chrom]
		current_acceptor_list = acceptor_coordinates[chrom]

		# window through each base in the seq
		# donor window first
		count = 0
		for i in range(len(seq) - donor_base_length +1):
			if count == int(d_lim_num / 6): break
			# fix to 1-based index
			base_number = i + 1
			# if not in donor list, check if first 2 bases are GT
			# if yes, store as negative coord

			if seq[base_number - 1:base_number + 1].upper() == 'GT':
				if base_number in current_donor_list: continue
				count += 1
				negative_donor_coords[chrom].append(base_number)

		# after finished with windowing through the donor for the particular chrom,
		# window through the acceptor
		count = 0
		for i in range(len(seq) - acceptor_base_length +1):
			if count == int(a_lim_num / 6): break
			if seq[i+acceptor_base_length-2:i+acceptor_base_length].upper() == 'AG':
				if i+1 not in current_acceptor_list:
					count += 1
					negative_acceptor_coords[chrom].append(i+1)

	return negative_donor_coords, negative_acceptor_coords

d_lim_chrom = ['I', 'II', 'III', 'IV', 'V', 'X']
d_lim_num = 'NA' 
a_lim_chrom = ['I', 'II', 'III', 'IV', 'V', 'X']
a_lim_num = 'NA' 

donor_coordinates, acceptor_coordinates = extract_positive_coordinates(fasta, gff, 25, d_lim_chrom, d_lim_num, a_lim_chrom, a_lim_num)

# for 1pct genome

'''
d_lim_num = 2547
a_lim_num = 2394
'''

# for 100pct genome
# will take around 70 minutes for either donor or acceptor
d_lim_num = 263477
a_lim_num = 251983

neg_donor_coordinates, neg_acceptor_coordinates =  extract_negative_coords(donor_coordinates, acceptor_coordinates, fasta, 9, 25, d_lim_num, a_lim_num)

# open genome as dictionary of SeqIO sequence objects

with gzip.open(fasta, 'rt') as handle:
	    genome_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))# Function to get sequence


# convert dictionary of donor or acceptor coord list into list of tuples
def to_tuples(coord_dictionary, window_len):
	# create list
	tuple_list = []
	for key in coord_dictionary.keys():
		for coord in coord_dictionary[key]:
			tuple_list.append((key, coord, coord+window_len-1))
	return tuple_list


# Function to get sequence
# ASSUMES ALL COORDS ARE ON POS STRAND
# 'chrom' is just the roman numeral 
def get_sequence(chrom, start, end):
	return str(genome_dict[chrom].seq[start-1:end])

all_coords_list = [donor_coordinates, acceptor_coordinates, neg_donor_coordinates, neg_acceptor_coordinates]
d, a, n_d, n_a = [], [], [], []
lists = [d, a, n_d, n_a]

for i, coord_dict in enumerate(all_coords_list):
	if i % 2 == 0: lists[i] = [get_sequence(chrom, start, end) for chrom, start, end in to_tuples(coord_dict, 9)]
	if i % 2 != 0: lists[i] = [get_sequence(chrom, start, end) for chrom, start, end in to_tuples(coord_dict, 25)]

# convert into FASTA format
# print FASTA for d, n_d, a, or n_a
def print_fasta(pick):
	if pick == 'd':
		for i, seq in enumerate(lists[0]):
			print(f'>d_seq_{i}')
			print(seq)
	if pick == 'a':
		for i, seq in enumerate(lists[1]):
			print(f'>a_seq_{i}')
			print(seq)
	if pick == 'n_d':
		for i, seq in enumerate(lists[2]):
			print(f'>n_d_seq_{i}')
			print(seq)
	if pick == 'n_a':
		for i, seq in enumerate(lists[3]):
			print(f'>n_a_seq{i}')
			print(seq)

print_fasta(pick)




