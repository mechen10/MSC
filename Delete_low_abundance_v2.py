#!/bin/bash/python


import argparse
import os, os.path
import sys
import subprocess

#=========================================================

parser = argparse.ArgumentParser(
	description="OTU table; in biom format")
parser.add_argument(
	'-i',
	'--otu_table',
	help = "Biom format",
	required = True)
parser.add_argument(
	'-t',
	'--threshold_abs',
	help = 'Threshold for minimum number of reads for each OTU to be included in a sample. Only one of -t and -p should be included.',
	required = False,
	default = False)
parser.add_argument(
	'-p',
	'--threshold_perc',
	help = "Threshold for minimum relative abundance for each OTU to be included in a sample. Only one of -t and -p should be included.",
	required = False,
	default = False)
parser.add_argument(
	'-T',
	'--threshold_abs_overall',
	help = "Threshold for minimum number of reads for each OTU in entire OTU table for OTU to be included in table. Only one of -T and -P should be included.",
	required = False,
	default = False)
parser.add_argument(
	'-n',
	'--minimum_per_sample',
	help = "Threshold for minimum number of reads per sample for sample to be included in final table.",
	required = False,
	default = False)
parser.add_argument(
	'-m',
	'--maximum_per_sample',
	help = "Threshold for maximum number of reads per sample for sample to be included in final table.",
	required = False,
	default = False)
parser.add_argument(
	'-o',
	'--output',
	help = 'Output name for new OTU table',
	required = True)


#=========================================================
# SET ARGS

args = parser.parse_args()

otu_table = args.otu_table
threshold_abs = args.threshold_abs
threshold_perc = args.threshold_perc
threshold_abs_overall = args.threshold_abs_overall
min_sample = args.minimum_per_sample
max_sample = args.maximum_per_sample
output = args.output

# FOR TROUBLESHOOTING

# otu_table='/Users/melissachen/Documents/Masters/Project_Environmental/FERDOUS_filtered_otu_tables/18S/OTU_Table_wtaxa.biom'
# threshold_abs=3
# threshold_perc=False
# threshold_abs_overall=1
# min_sample=1
# max_sample=False
# output='otu_table_wtaxa_filtered_t3T1'
#=========================================================

# Functions

def loadOTUTable(OTUFP): # Loads and converts to relative abundance
	global threshold_abs_overall
	# Make OTU Table
	os.system('biom convert -i ' + OTUFP + ' --to-tsv --header-key taxonomy -o OTU_Table_text_FORFILT.txt')
	# Load OTU Table
	biomOpen = open("OTU_Table_text_FORFILT.txt", 'r') # This is in Unix format
	biomTemp = []
	for i in biomOpen:
		tempLine = i.strip()
		tempLineSplit = tempLine.split("\t")
		biomTemp.append(tempLineSplit)
	biomOpen.close()
	# Get rid of last column
	OTUTable = {} # Now make dictionary
	taxaIDs = {} # Make taxonomy reference 
	for lineN in range(len(biomTemp)):
		if lineN == 0: # This is first line; skip this
			pass
		elif lineN == 1: # This is the header line
			headers = biomTemp[lineN][1:len(biomTemp[lineN])-1]
		else:
			OTUTable[str(biomTemp[lineN][0])] = {}
			taxaIDs[str(biomTemp[lineN][0])] = biomTemp[lineN][len(biomTemp[lineN])-1]
			for abund in range(len(biomTemp[lineN][1:])-1):
				OTUTable[str(biomTemp[lineN][0])][headers[abund]] = float(biomTemp[lineN][1:][abund])
	# If you entered a threshold for overall, delete those here
	if threshold_abs_overall != bool:
		todelete = []
		for otu in OTUTable:
			totalreads = sum([OTUTable[otu][i] for i in OTUTable[otu]])
			if totalreads < float(threshold_abs_overall):
				todelete.append(otu)
		# get rid of those reads
		for key in todelete:
			del OTUTable[key]
	# Make a reads version of OTU Table
	OTUTable_abs = OTUTable.copy()
	# Get total counts for each site
	totalCounts = {}
	for h in range(len(headers)):
		totalCounts[headers[h]] = 0
	for OTU in OTUTable:
		for h in range(len(headers)):
			totalCounts[headers[h]] += float(OTUTable[OTU][str(headers[h])])
	# Convert to relative abundance
	for OTU in OTUTable.keys():
		tempOTUlist = OTUTable[OTU].copy()
		for sites in OTUTable[OTU].keys():
			tempOTUlist[sites] = float(OTUTable[OTU][sites])/float(totalCounts[sites])
		OTUTable[OTU] = tempOTUlist
	os.system("rm OTU_Table_text_FORFILT.txt")
	return OTUTable, OTUTable_abs, taxaIDs # Output is a 2-layer dictionary; first is OTU IDs and second is samples. Also, one-layer dict with taxa IDs
	
	
def switchMatDirections(table):
	switched = {}
	first = True
	for otu in table.keys():
		for sample in table[otu].keys():
			if first:
				switched[sample] = {otu: table[otu][sample]}
			else:
				switched[sample][otu] = table[otu][sample]
		first = False
	return switched
	

	
#=========================================================

# Load OTU table
OTUTable,OTUTable_abs, taxaIDs = loadOTUTable(otu_table)

# Switch table around
otu_switch = switchMatDirections(OTUTable)
otu_abs_switch = switchMatDirections(OTUTable_abs)
	
# DO SAMPLE-BY-SAMPLE 
# For relative abundance threshold:
if type(threshold_perc) != bool and type(threshold_abs) != bool:
	print "Please choose ONE: percent threshold or absolute threshold. Not both."
elif type(threshold_perc) != bool:
	todelete = {}
	for sample in otu_switch:
		todelete[sample] = []
		for otu in otu_switch[sample]:
			if otu_switch[sample][otu] < float(threshold_perc):
				todelete[sample].append(otu)
# For absolute abundance threshold:
elif type(threshold_abs) != bool:
	todelete = {}
	for sample in otu_abs_switch:
		todelete[sample] = []
		for otu in otu_abs_switch[sample]:
			if otu_abs_switch[sample][otu] < float(threshold_abs):
				todelete[sample].append(otu)
else:
	print "Your threshold for perc or abs is not in the correct format. Please enter a whole number for perc and a decimal for abs."
			

# Now, delete the ones that marked
OTUTable_filtered = OTUTable_abs.copy()
for otu in OTUTable_filtered:
	for sample in OTUTable_filtered[otu]:
		if otu in todelete[sample]:
			OTUTable_filtered[otu][sample] = 0

# Delete samples below a certain threshold; 
# First, switch over matrix directions
OTUTable_filtered_switched = switchMatDirections(OTUTable_filtered)

# Then,figure out min and max sample abundances
if type(min_sample) == bool:
	min_sample = 0
if type(max_sample) == bool:
	max_sample = sum([OTUTable_filtered_switched[sample][otu] for otu in OTUTable_filtered_switched[sample] for sample in OTUTable_filtered_switched])

# now filter
todelete = []
for sample in OTUTable_filtered_switched:
	sampSum = sum([OTUTable_filtered_switched[sample][i] for i in OTUTable_filtered_switched[sample]])
	if sampSum < float(min_sample) or sampSum > float(max_sample):
		todelete.append(sample)
for sample in todelete:
	del OTUTable_filtered_switched[sample]
	
OTUTable_final = switchMatDirections(OTUTable_filtered_switched)
	

# Print OTU table
first = True
toWrite = "#OTU ID"
for row in OTUTable_final:
	if first:
		for sample in OTUTable_final[row].keys():
			toWrite += "\t" + sample 
		toWrite += "\t"+ "taxonomy" + "\n"
		first = False
	toWrite += row
	for abund in OTUTable_final[row]:
		toWrite += "\t" + str(OTUTable_filtered[row][abund])
	toWrite += "\t" + taxaIDs[row] + "\n"

open(output+".txt", 'w').write(toWrite)
os.system('biom convert -i ' + output + ".txt" + ' --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy -o '+ output + ".biom")
