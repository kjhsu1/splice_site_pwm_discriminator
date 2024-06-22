# Creates the weighted pwm from C.elegans gff and genommic fasta data

# input: gff and fasta
# output: Weighted PWM in JSON format

# gff has feature name "intron" so can extract the first 7 bases
# (donor) and the last 6 bases of the intron region (acceptor) 

# 84splicesite.py already does this, but it only counts each entry once. 

# the gff data is gathered from when they did a RNAseq of the whole cell
# one entry, denotes the specific intron that they observed
# from the experiment 
# In each entry, it contains the # number of occurances of the specific
# intron as well. 
	# 84splicesite.py doesn't use this data to create a weighted-by-
	# occurance PWM 

# In other words, we just have to make sure that program doesn't just
# +1 for each "intron" entry" it sees, and instead, +occurance


###############
# Code Below


'''
- read gff file
- find all the introns

- first 6nt == splice donor 
- last 7nt == splice acceptor 
- ABOVE WAS TRUE, BUT NOW CAN CHANGE # BASES TO EXTRACT

- create 2 dictionaries
	- one for splice acceptor 
	- one for splice donor 
- go through gff file linearly, for each 'intron'
+= 1 (or +=# of occurance to appropriate bases for each position 
- at the end should have filled dictionary
'''

import sys 
import mcb185
import gzip
import json
import dogma

gff = sys.argv[1]
fa = sys.argv[2]


def po_dict_maker(site_length):
	po = {}
	for i in range(1, site_length+1):
		po[i] = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
	return po

def printer(org, a_or_d, pwm):
	print(f'AC {org}') # org = organism name 
	print('XX')
	print(f'ID {org}')
	print('XX')
	print(f'DE {a_or_d}')
	print(f'PO{"A":>15}{"C":>15}{"G":>15}{"T":>15}')
	# this goes through each key in acceptor dict
	for i, nt in enumerate(pwm.keys()): 
		i += 1
		a = pwm[nt]['A']
		c = pwm[nt]['C']
		g = pwm[nt]['G']
		t = pwm[nt]['T']
		print(f'{i}{a:>15}{c:>15}{g:>15}{t:>15}')
	print('XX')
	print('//')

def splice_site_pwm(gff, fa, donor_bases, acceptor_bases, org):
	introns = []
	with gzip.open(gff, 'rt') as fp:
		for line in fp:
			words = line.split()
			if words[1] != 'RNASeq_splice': continue
			chrom = words[0]
			start = int(words[3]) - 1
			end = int(words[4]) - 1
			n = words[5]
			strand = words[6]
			introns.append((chrom, start, end, n, strand))

	donor = po_dict_maker(donor_bases)
	acceptor = po_dict_maker(acceptor_bases)

	for defline, seq in mcb185.read_fasta(fa):
		defline_words = defline.split()
		for intron in introns:
			# FOR WEIGHTED PMW
			# store the # occurances of each intron
			n = int(float(intron[3]))
			# if correct chromosome
			if defline_words[0] == intron[0]:
				# if intron is '+' strand
				if intron[4] == '+':
					intron_seq = seq[intron[1]:intron[2]+1]
				# intron has to be at least 60 nt long
				#if len(intron_seq) < 60: continue
				# if minus strand, do reverse comp 
				elif intron[4] == '-':
					intron_seq = mcb185.anti_seq(seq[intron[1]:intron[2]+1])
				# donor and acceptor seq
				d_seq = intron_seq[:donor_bases]
				a_seq = intron_seq[-acceptor_bases:]
				# ONLY PART MODIFIED FOR WEIGHTED VERSION
				# instead of adding only += 1 for each nucleotide in 
				# donor and acceptor seq, add += 1 * number of occurance
				# of the particular intron that the donor and acceptor
				# was extracted from
					
				# donor
				for i, nt in enumerate(d_seq):
					index = i + 1
					donor[index][nt] += 1 * n 
				# acceptor 
				for i, nt in enumerate(a_seq):
					index = i + 1
					acceptor[index][nt] += 1 * n 
	# print 
	printer(org, 'splice acceptor', acceptor)
	printer(org, 'splice donor', donor)


# running for 9 bases for donor and 25 bases for acceptor
splice_site_pwm(gff, fa, 9, 25, 'C.elegans')




