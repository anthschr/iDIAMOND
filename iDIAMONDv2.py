# python3 iDIAMONDv2.py batch-0.tsv batch-4-1.tsv 29833021 7931228 


# USAGE: ./python3 iDiamond.py fasta1 fasta2
# OR
# USAGE: ./python3 iDiamond.py out1 out2
print("To run iDIAMOND, you need the following")
print("2 output files, the number of letters in each output file")
print("USAGE: python3 iDIAMOND.py out1 out2 DB1_letter_count DB2_letter_count")

import sys
import csv
import copy

def tsv_to_2d_list(filename):
	out=[]
	with open(filename) as file:
		tsv_file = csv.reader(file,delimiter="\t")
		for line in tsv_file:
			out.append(line)
	return out

def print_2d_list(list):
	for item in list:
		N = len(item)
		for i in range(N):
			if i < N-1:
				if type(item[i]) == type(3.4) and item[i] < 1:
					print("{:.5e}\t".format(item[i]),end='')
				else:
					print("{}\t".format(item[i]),end='')
			else:
				print(item[i])
	return

def S2_rescale(results1,results2, DB1_length, DB2_length):
	# ASSUMPTION: that DIAMOND returns as many results for as many sequences in the database
	# that it doesn't discard certain results below a certain E value
	# Take the E value from running a section of the database and multiply it by (n_total / n_part)
	# must rescale evalue from both results
	# rescale from first set of results first
	DB1_length=int(DB1_length)
	DB2_length=int(DB2_length)
	DB_total_length = DB1_length + DB2_length
	scaled_results1 = rescale(results1,DB1_length,DB_total_length)
	scaled_results2 = rescale(results2,DB2_length,DB_total_length)
	return [scaled_results1, scaled_results2]

def rescale(results, results_length, total_length):
	results = copy.deepcopy(results)
	for row in results:
		e_value = float(row[DIAMOND_COLUMNS["E-value"]])
		e_value = e_value * total_length / results_length
		row[DIAMOND_COLUMNS["E-value"]] = e_value
	return results

def S3_merge_results(R1,R2):
	R1=copy.deepcopy(R1)
	R2=copy.deepcopy(R2)
	n=0;m=0
	e_value_1=R1[m][DIAMOND_COLUMNS["E-value"]]; score1=R1[m][DIAMOND_COLUMNS["Bit_score"]]
	e_value_2=R2[n][DIAMOND_COLUMNS["E-value"]]; score2=R2[n][DIAMOND_COLUMNS["Bit_score"]]
	merged_results = []
	num_hits = len(R1) + len(R2)
	print("Num Hits:%s" % num_hits)
	for i in range(num_hits):
		if ((e_value_1 < e_value_2) or (e_value_1 == e_value_2 and score1 > score2)) and (m < (len(R1)-1)):
			merged_results.append(R1[m])
			m+=1
			e_value_1=R1[m][DIAMOND_COLUMNS["E-value"]]; score1=R1[m][DIAMOND_COLUMNS["Bit_score"]]
		elif n < (len(R2)-1):
			merged_results.append(R2[n])
			n+=1
			e_value_2=R2[n][DIAMOND_COLUMNS["E-value"]]; score2=R2[n][DIAMOND_COLUMNS["Bit_score"]]
	print("i={}".format(i))
	print("n={}".format(n))
	print("m={}".format(m))
	print("Length of R1: {}".format(len(R1)))
	print("Length of R2: {}".format(len(R2)))
	return merged_results

DIAMOND_COLUMNS = {
# Format for old output batches: --outfmt/-f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen
# Format for new output batches: --outfmt/-f 6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore 
	"Query_Seq_ID":0,
	"Query_Seq_Len":1,
	"Subject_Seq_ID":2,
	"Subject_Seq_Len":3,
	"Sequence_identity":4,
	"Length":5,
	"Mismatches":6,
	"Gap_openings":7,
	"Query_start":8,
	"Query_end":9,
	"Target_start":10,
	"Target_end":11,
	"E-value":12,
	"Bit_score":13,
	"Sequence_Length":14}

def make_numbers(results):
	results = copy.deepcopy(results)
	for i in range(len(results)):
		if type(results[i]) == type(list()):
			#print("Type:List")
			results[i] = make_numbers(results[i])
		else:
			try: 
				results[i] = int(results[i])
			except:
				try:
					results[i] = float(results[i])
				except:
					pass

	return results


# Main Code
num_args = len(sys.argv)
argv=sys.argv
# print("Num args {}".format(num_args))
if num_args != 5:
	print("4 arguments required\nUSAGE: python3 iDIAMOND.py out1 out2 DB1_letter_count DB2_letter_count")
	sys.exit(2)

out1_filename=argv[1]
out2_filename=argv[2]
DB1_length=argv[3]
DB2_length=argv[4]
out1=[]
out2=[]
out1=make_numbers(tsv_to_2d_list(out1_filename))
out2=make_numbers(tsv_to_2d_list(out2_filename))

scaled_results = S2_rescale(out1,out2, DB1_length, DB2_length)
merged_results = S3_merge_results(scaled_results[0],scaled_results[1])
with open("merged_output.tsv",'w') as fp:
	for row in merged_results:
		for i in range(len(row)):
			if i < len(row)-1:
				fp.write("{}\t".format(row[i]))
			else:
				fp.write("{}\n".format(row[i]))



