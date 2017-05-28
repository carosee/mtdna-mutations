import numpy as np
import matplotlib.pyplot as plt
import csv
from matplotlib import cm
from scipy.stats import linregress
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Nexus import Nexus
from collections import Counter



# clean up data

def get_counts(data):
	counts = []
	for loc in data.T:
		# print(loc)
		counts.append(Counter(loc))
	return counts

def clean_counts(counts): #remove 'N' and '-' from counts. maybe better to assign them the majority?
	for d in counts:
		# d.pop('N', None)
		d.pop('-', None)
		d.pop('N', None)
	return counts

def clean_data(data):
	# replace '-' with N
	data[data == '-'] = 'N'
	data[data == '?'] = 'N'
	return data




# change between dna - codon - aa arrays

def get_codon_array(sequence_array): # input: 8000 x 80 array; return 2899 x 80 array
	codon_array = []
	for species in sequence_array:
		s = ''.join(species)
		codon_lst = [(s[i:i + 3]) for i in range(0, len(species), 3)]
		codon_array.append(codon_lst)
	return np.array(codon_array)


def get_aa_array(data):
	# data = a numpy array of species x nucleotide. must use clean_data first
	# return: a numpy array of species x translated proteins
	sequences = []
	for row in data:
		s = ''.join(row)
		seq = Seq(s, generic_dna)
		sequences.append(seq.translate(table="Vertebrate Mitochondrial"))
	# print(sequences)
	split = []
	for t in sequences:
		split.append(list(t))
	return np.array(split)


def base_to_codon(base_index):
	return base_index // 3

def codon_to_base(codon_index):
	return (codon_index*3)+2




# filter rates 


def filter_4fold(nex_data):
	translated = get_aa_array(nex_data) # first dimension is species, second dim is protein locus

	counts = get_counts(translated) # a list of counters
	majority = [max(count, key=count.get) for count in counts]

	degenerate = ['A', 'G', 'L', 'P', 'R', 'S', 'T', 'V']
	codon_indices = [i for i,x in enumerate(majority) if x in degenerate]
	return codon_indices


def filter_AG(nex_data):
	counts = clean_counts(get_counts(nex_data))

	filtered_base_indices = [i for i in range(len(counts)) if i % 3 == 2 and counts[i].keys() == ['A', 'G']] # only 3rd codon positions
	# print(filtered_base_indices)
	codon_indices = [base_to_codon(i) for i in filtered_base_indices]
	return list(set(codon_indices)) # remove duplicates, will return unordered

def filter_AGtest(nex_data):
	counts = clean_counts(get_counts(nex_data))

	filtered_base_indices = [i for i in range(len(counts)) if i % 3 == 2 and counts[i].keys() == ['A', 'G']] # only 3rd codon positions
	# print(filtered_base_indices)
	ag_ratios = [counts[i]['A'] for i in filtered_base_indices]
	codon_indices = [base_to_codon(i) for i in filtered_base_indices]
	return list(set(codon_indices)) # remove duplicates, will return unordered

def filter_CT(nex_data):
	counts = clean_counts(get_counts(nex_data))

	filtered_base_indices = [i for i in range(len(counts)) if i % 3 == 2 and counts[i].keys() == ['C', 'T']] # only 3rd codon positions
	# print(filtered_base_indices)
	codon_indices = [base_to_codon(i) for i in filtered_base_indices]
	return list(set(codon_indices)) # remove duplicates, will return unordered


def filter_selection(csv_data, pos_tolerance=0.9, neg_tolerance=0.9):
	# pos_tolerance: exclude points with larger probs of positive selection (a smaller value excludes more points)
	# neg_tolerance: exclude points with larger probs of negative selection
	codon_indices = [row for row in range(len(csv_data)) if csv_data[row][5] < pos_tolerance and csv_data[row][6] < neg_tolerance]
	# 5 is prob positive selection, 6 is prob negative selection
	return codon_indices

def filter_outliers(csv_data, tolerance=3):
	codon_indices = [row for row in range(len(csv_data)) if csv_data[row][2] < tolerance]
	return codon_indices




# plot rates

def regression(x,y):
	slope, intercept, r_value, p_value, std_err = linregress(x,y)
	print('slope', slope)
	x = np.array(range(max(x)))
	formula = lambda x: slope * x + intercept
	y = formula(x)
	plt.plot(x, y, 'r-')
	return slope


def plot_all(indices, title="Substitution rates per codon", file="plot_all.png"):
	csv_data = np.genfromtxt('mergedcsv2.csv', delimiter=',')
	# x = [i for i in range(len(csv_data)) if i in indices]
	x = indices
	y = [csv_data[i][2] for i in x]
	z = [csv_data[i][0] for i in x]
	plt.scatter(x, y, c=z, marker = 'o', cmap = cm.flag)
	axes = plt.gca()
	axes.set_xlim([0,2900])
	axes.set_ylim([0,11])
	# plt.plot(x,y,'bo')
	slope = regression(x,y)
	plt.title(title)
	plt.text(100, 10, 'slope: ' + str(slope))
	plt.savefig(file)
	plt.show()


def plot_avg(n, indices, title="Average substitution rates per n codons", file="plot_avg.png"):
	# plot avg value every n indices
	csv_data = np.genfromtxt('mergedcsv2.csv', delimiter=',')
	vals = range(0, max(indices), n)
	x = []
	y = []
	z = []
	for i in vals:
		interval = [elem for elem in indices if elem > i and elem < i+n]
		if interval:
			rates = [csv_data[i][2] for i in interval]
			avg = float(sum(rates))/float(len(interval))
			x.append(i)
			y.append(avg)
			z.append(csv_data[i][0])
	plt.scatter(x, y, c=z, marker = 'o', cmap = cm.flag, s=50)
	axes = plt.gca()
	axes.set_xlim([0,3600])
	axes.set_ylim([0,2])
	slope = regression(x,y)
	plt.title(title)
	plt.text(100, 1.9, 'slope: ' + str(slope))
	plt.savefig(file)
	plt.show()


def plot_logprod(n, indices):
	# plot log product every n indices
	csv_data = np.genfromtxt('mergedcsv2.csv', delimiter=',')
	vals = range(0, max(indices), n)
	x = []
	y = []
	for i in vals:
		interval = [elem for elem in indices if elem > i and elem < i+n]
		if interval:
			rates = [csv_data[i][2] for i in interval]
			logprod = np.prod(np.log(rates))
			x.append(i)
			y.append(logprod)
	print(len(x))
	slope = regression(x,y)
	plt.plot(x,y,'bo')
	plt.text(100, 0.1, 'slope: ' + str(slope))
	plt.show()



# run actual code

handle = open('nex/mtcombined.nex', "rU")
data_lst = list(SeqIO.parse(handle, "nexus"))
nex_data = clean_data(np.array(data_lst)) # first dimension (row) is species, second dim is base locus
csv_data = np.genfromtxt('mergedcsv2.csv', delimiter=',')

# print(nex_data)
# codons = get_codon_array(nex_data)
# print(codons)
# print(get_counts(codons)[200])

# plot all points
# plot_all(range(len(csv_data)), title="Substitution rates per codon", file="plot_all.png")

# filter points for outliers
# plot_avg(30, filter_outliers(csv_data), title="Substitution rates per 30 codons, excluded rates > 3", file="plot_avg_filter_outliers.png")

# filter points for selection
selection_pts = filter_selection(csv_data)
print(len(selection_pts))

# plot_all(filter_selection(csv_data), title="Substitution rates per codon, filtered for selection", file="plot_all_filter_selection.png")
# plot_avg(30, filter_selection(csv_data), title="Average substitution rates per 30 codons, filtered for selection", file="plot_avg_filter_selection.png")

# filter 4 fold degenerate codons
degenerate_pts = filter_4fold(nex_data)
print(len(degenerate_pts))
# plot_all(filter_4fold(nex_data), title="Substitution rates of 4 fold degenerate codons", file="plot_all_4fold.png")


intersect = list(set(degenerate_pts) & set(selection_pts))
# print(intersect)
# print(len(intersect))
plot_all(intersect, title="subsitution rates of intersection of 4fold and selection filtered", file="plot_all_4fold_selection")


# filter A-G substitutions
# plot_all(filter_AG(nex_data), title="Substitution rates of codons with only 3rd position AG substitutions", file="plot_all_AG.png")

# filter C-T substitutions
# plot_all(filter_CT(nex_data), title="Substitution rates of codons with only 3rd position CT substitutions", file="plot_all_CT.png")

# a = filter_4fold(nex_data)
# b = filter_AG(nex_data)
# c = filter_selection(csv_data)

# plot_all(list(set(a) & set(b) & set(c)), title="Substitution rates of 4fold && AG only, filtered for selection", file="plot_all_4fold_AG_selection.png")






