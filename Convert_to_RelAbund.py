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
	'-o',
	'--output',
	help = 'Output name for new OTU table',
	required = True)


#=========================================================
# SET ARGS

args = parser.parse_args()

otu_table = args.otu_table
output = args.output

#=========================================================

# Functions

def loadOTUTable(OTUFP): # Loads and converts to relative abundance
	# Make OTU Table
	os.system('biom convert -i ' + OTUFP + ' --to-tsv --header-key taxonomy -o OTU_Table_text_FORRELABUND.txt')
	# Load OTU Table
	biomOpen = open("OTU_Table_text_FORRELABUND.txt", 'r') # This is in Unix format
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
	os.system("rm OTU_Table_text_FORRELABUND.txt")
	return OTUTable, OTUTable_abs, taxaIDs # Output is a 2-layer dictionary; first is OTU IDs and second is samples. Also, one-layer dict with taxa IDs


#=========================================================

# Load OTU table
OTUTable,OTUTable_abs, taxaIDs = loadOTUTable(otu_table)

# Print OTU table
first = True
toWrite = "#OTU ID"
for row in OTUTable:
	if first:
		for sample in OTUTable[row].keys():
			toWrite += "\t" + sample 
		toWrite += "\t"+ "taxonomy" + "\n"
		first = False
	toWrite += row
	for abund in OTUTable[row]:
		toWrite += "\t" + str(OTUTable[row][abund])
	toWrite += "\t" + taxaIDs[row] + "\n"

open(output+".txt", 'w').write(toWrite)
os.system('biom convert -i ' + output + ".txt" + ' --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy -o '+ output + ".biom")
