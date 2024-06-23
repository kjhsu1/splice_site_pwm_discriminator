# Program takes in labeled log odds dictionary, and threshold log odds score, and spits out assessment

import json
import sys
import matplotlib.pyplot as plt
import numpy as np

json_file = sys.argv[1] # should be JSON formatted labeled log odds dictionary
upper_threshold = float(sys.argv[2])
lower_threshold = float(sys.argv[3])

# convert json into python data type

def json_reader(json_file):
	with open(json_file, 'r') as file:
		# load JSON data from file 
		data = json.load(file)

	# return 
	# data type is a bit confusing 
	# see what it prints for more info
	return data

labeled_donor_dict, labeled_acceptor_dict = json_reader(json_file)

# function takes in labeled donor and acceptor dict, and donor and acceptor thresholds
# returns list with TP, FP, TN, FN for both donor and acceptor
def discriminator(labeled_donor_dict, labeled_acceptor_dict, d_threshold, a_threshold):
	# create 4 dictionaries for sorting (pos donor, neg donor, pos acceptor, neg acceptor)

	pos_donor_dict = {}
	neg_donor_dict = {}
	pos_accept_dict = {}
	neg_accept_dict = {}

	# donor first
	# iterate through all chromosomes
	for key in labeled_donor_dict.keys():
		# add defline/chrom as key
		pos_donor_dict[key] = {
			'true+': [],
			'false+': []
		}

		neg_donor_dict[key] = {
			'true-': [],
			'false-': []
		}

		# iterate through all windows in the chrom
		for window in labeled_donor_dict[key].keys():
			# skip window if unlabeled
			if labeled_donor_dict[key][window]['class'] == 'NA': continue

			# check log odds score of window

			# if window has score greater or equal than threshold
			if d_threshold <= labeled_donor_dict[key][window]['log_odds']:
				# if it's also an actual pos
				if labeled_donor_dict[key][window]['class'] == '+':
					pos_donor_dict[key]['true+'].append(window)
				# if it's an actual neg
				if labeled_donor_dict[key][window]['class'] == '-':
					pos_donor_dict[key]['false+'].append(window)

			# if window has score less than threshold
			if d_threshold > labeled_donor_dict[key][window]['log_odds']:

				# if it is an actual neg
				if labeled_donor_dict[key][window]['class'] == '-':
					neg_donor_dict[key]['true-'].append(window)

				# if it is a actual pos
				if labeled_donor_dict[key][window]['class'] == '+':
					neg_donor_dict[key]['false-'].append(window)

	# do the acceptor as well
	for key in labeled_acceptor_dict.keys():
		# add defline/chrom as key
		pos_accept_dict[key] = {
			'true+': [],
			'false+': []
		}

		neg_accept_dict[key] = {
			'true-': [],
			'false-': []
		}

		# iterate through all windows in the chrom
		for window in labeled_acceptor_dict[key].keys():
			# skip window if unlabeled
			if labeled_acceptor_dict[key][window]['class'] == 'NA': continue

			# check log odds score of window

			# if window has score greater or equal than threshold
			if a_threshold <= labeled_acceptor_dict[key][window]['log_odds']:
				# if it's also an actual pos
				if labeled_acceptor_dict[key][window]['class'] == '+':
					pos_accept_dict[key]['true+'].append(window)
				# if it's an actual neg
				if labeled_acceptor_dict[key][window]['class'] == '-':
					pos_accept_dict[key]['false+'].append(window)

			# if window has score less than threshold
			if a_threshold > labeled_acceptor_dict[key][window]['log_odds']:

				# if it is an actual neg
				if labeled_acceptor_dict[key][window]['class'] == '-':
					neg_accept_dict[key]['true-'].append(window)

				# if it is a actual pos
				if labeled_acceptor_dict[key][window]['class'] == '+':
					neg_accept_dict[key]['false-'].append(window)

	# retrieve TP, FP, and so on for both donor and acceptor
	# donor first
	d_TP = 0
	d_FP = 0
	d_TN = 0
	d_FN = 0

	# iterate through all chroms
	# can just use the same key for both
	for key in pos_donor_dict.keys():
		d_TP += len(pos_donor_dict[key]['true+'])
		d_FP += len(pos_donor_dict[key]['false+'])
		d_TN += len(neg_donor_dict[key]['true-'])
		d_FN += len(neg_donor_dict[key]['false-'])

	# then acceptor
	a_TP = 0
	a_FP = 0
	a_TN = 0
	a_FN = 0

	for key in pos_accept_dict.keys():
		a_TP += len(pos_accept_dict[key]['true+'])
		a_FP += len(pos_accept_dict[key]['false+'])
		a_TN += len(neg_accept_dict[key]['true-'])
		a_FN += len(neg_accept_dict[key]['false-'])
	
	
	'''
	print(pos_donor_dict)
	print()
	print(neg_donor_dict)
	print() 
	print(pos_accept_dict)
	print()
	print(neg_donor_dict)
	'''

	# return 
	return [d_TP, d_FP, d_TN, d_FN], [a_TP, a_FP, a_TN, a_FN]


# testing multiple thresholds and making precision recall curve
def precision_recall_curve(upper_thresh, lower_thresh, step, donor_or_accept):
	# make 3 separate
	dTP = []
	dFP = []
	dFN = [] 

	aTP = []
	aFP = []
	aFN = []

	for i in np.arange(lower_thresh, upper_thresh, step):
		donor_data, acceptor_data = discriminator(labeled_donor_dict, labeled_acceptor_dict, i, i)

		# data is organized in the order: TP FP TN FN
		# add info to the list
		dTP.append(donor_data[0])
		aTP.append(donor_data[0])

		dFP.append(donor_data[1])
		aFP.append(donor_data[1])

		dFN.append(donor_data[3])
		aFN.append(donor_data[3])

	# precision and recall
	epsilon = 1e-10

	d_precision = [tp / (tp + fp + epsilon) for tp, fp in zip(dTP, dFP)]
	d_recall = [tp / (tp + fn + epsilon) for tp, fn in zip(dTP, dFN)]

	a_precision = [tp / (tp + fp + epsilon) for tp, fp in zip(aTP, aFP)]
	a_recall = [tp / (tp + fn + epsilon) for tp, fn in zip(aTP, aFN)]

	print(d_precision)
	print(a_precision)

	# if you want to plot donor
	if donor_or_accept == 'donor':
		# Plotting the precision-recall graph
		plt.plot(d_recall, d_precision, marker='o')

		# Adding title and labels
		plt.title('Donor Site Precision-Recall Graph')
		plt.xlabel('Recall')
		plt.ylabel('Precision')

		# setting the limit for x and y axis
		plt.xlim([0, 1])
		plt.ylim([0, 1])

		# Display the graph
		plt.show()

	if donor_or_accept == 'acceptor':
		# Plotting the precision-recall graph
		plt.plot(a_recall, a_precision, marker='o')

		# Adding title and labels
		plt.title('Acceptor Site Precision-Recall Graph')
		plt.xlabel('Recall')
		plt.ylabel('Precision')

		# setting the limit for x and y axis
		plt.xlim([0, 1])
		plt.ylim([0, 1])

		# Display the graph
		plt.show()


precision_recall_curve(upper_threshold, lower_threshold, 0.1, 'donor')