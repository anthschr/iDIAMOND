import sys
import csv
import copy
import pandas as pd
import matplotlib.pyplot as plt
from tabulate import tabulate

def to_fwf(df,fname):
	content=tabulate(df.values.tolist(),list(df.columns),tablefmt="plain")
	open(fname,"w").write(content)

def arrangeDataFrames(df):
	df[0] = df[0] + "->" + df[2]
	df.drop([1,2,3,4,5,6,7,8,9,10,11],inplace=True,axis=1)
	df.columns=['QuerytoTarget','EVALUE','BitScore']

# Main Code
print("2 arguments required\nUSAGE: python3 compareDIAMOND.py Original.tsv Merged.tsv")
num_args = len(sys.argv)
argv=sys.argv
# print("Num args {}".format(num_args))
if num_args != 3:
	print("2 arguments required\nUSAGE: python3 compareDIAMOND.py Original.tsv Merged.tsv")
	sys.exit(2)

original=argv[1]
merged=argv[2]
print("Reading Data")
original_data=pd.read_csv(original,sep="\t", header=None)
merged_data=pd.read_csv(merged,sep="\t", header=None)
print("Arranging Dataframes")
arrangeDataFrames(original_data)
arrangeDataFrames(merged_data)
print("Joining Dataframes")
LeftJoin = pd.merge(original_data, merged_data,on='QuerytoTarget',how='left')
LeftJoin.rename(columns={'QuerytoTarget':'QuerytoTarget','EVALUE_x':'Orig_Evalue','BitScore_x':'Orig_BitScore','EVALUE_y':'iDIAMOND_Evalue','BitScore_y':'iDIAMOND_BitScore'},inplace=True)
LeftJoin['Perc_Diff'] = 100*(LeftJoin['iDIAMOND_Evalue'] - LeftJoin['Orig_Evalue'])/LeftJoin['Orig_Evalue']
print("Saving to TSV")
# Use tabulate hear
to_fwf(LeftJoin,"FWF.txt")
# LeftJoin.to_csv(path_or_buf='LeftJoin.tsv', sep = '\t',float_format="%e.6")
# print("Creating Histogram")
# perc_diff_nums = LeftJoin['Perc_Diff'].tolist()
# plt.hist(perc_diff_nums,bins=100)
# plt.title("Distribution of Percent Error")
# plt.xlabel("Percent Difference from DIAMOND")
# plt.ylabel("Frequency")
# plt.show()
print("Done")

