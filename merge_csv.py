import csv

csv_list = ['CO1_FUBAR_report.csv', 'CO2_FUBAR_report.csv', 'ATP8_FUBAR_report.csv', 'ATP6_FUBAR_report.csv', 'CO3_FUBAR_report.csv',
			 'ND3_FUBAR_report.csv', 'ND4L_FUBAR_report.csv', 'ND4_FUBAR_report.csv', 'ND5_FUBAR_report.csv', 'Cytb_FUBAR_report.csv',
			  'ND1_FUBAR_report.csv', 'ND2_FUBAR_report.csv']

reverse_csv = 'ND6_FUBAR_report.csv'


start_dict = { 'ND1' : 3358,
				'ND2' : 4882,
				'CO1' : 82,
				'CO2' : 1231,
				'ATP8' : 1730,
				'ATP6' : 1828,
				'CO3' : 2281,
				'ND3' : 2844,
				'ND4L' : 3115,
				'ND4' : 3308,
				'ND5' : 4327,
				'ND6_r' : 5502+166, # beginning + length; starting at the ND6 start codon and counting backwards
				'Cytb' : 5880,
			}

def reverse_ND6(input_file, output_file):
	fout = open(output_file,"a")
	f = open(input_file)
	f.next()
	gene = 'ND6_r'
	distance = start_dict[gene]
	for line in f:
		fout.write( gene + "," + str(distance) + "," + line ) # add name of gene & distance
		distance -= 2
	f.close()
	fout.close()

def ND1_model2(input_file, output_file):
	fout = open(output_file,"a")
	f = open(input_file)
	f.next()
	gene = 'ND1'
	distance = 9086
	for line in f:
		fout.write( gene + "," + str(distance) + "," + line ) # add name of gene & distance
		distance += 2
	f.close()
	fout.close()


def merge_csv(input_list, output_file): # seems to be working fine now
	fout = open(output_file,"a")
	for file in input_list:
		f = open(file)
		f.next() # skip the header
		gene = file[:-17]
		distance = start_dict[gene]
		for line in f:
			fout.write( gene + "," + str(distance) + "," + line ) # add name of gene & distance
			distance += 2
		f.close() # not really needed
	fout.close()


# reverse_ND6(reverse_csv, 'ND6_r_distances.csv')
# merge_csv(csv_list, 'mt_distances.csv')
ND1_model2('ND1_FUBAR_report.csv', 'ND1_M2_distances.csv')

