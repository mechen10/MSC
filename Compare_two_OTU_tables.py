#!/bin/bash

#!/bin/bash/python
#### OTUs IN COMMON ###

# Built under QIIME 1.9.0 and python 2.7.3
# Made 13Feb2017
# Not test yet; but should work?

# Script for determining shared microbiomes in different treatment groups all at once. 
import argparse
import os
import sys
import copy


#################################

parser = argparse.ArgumentParser(
	description="Takes multiple OTU tables and merges them. Uses assign_taxonomy output to determine matching OTUs. **Uses OTU IDs from first OTU table**")
parser.add_argument(
	'-i',
	'--OTU_table',
	help = "Comma separated list of OTU table file paths from QIIME containing all samples",
	required = True,)
parser.add_argument(
	'-m',
	'--metadata',
	help = 'Comma separated mapping files for each OTU table; must be in same order as OTU_table',
	required = True)
parser.add_argument(
	'--align_log',
	help = 'The log file from align_seqs; or, can be any tab delimited file with the following headers: "candidate sequence ID", "template ID", "BLAST percent identity to template". (The last one is used to tell whether the file has ended; just need something in that column ',
	required = True)
parser.add_argument(
	'-o',
	'--outputFolder',
	help = 'Output folder [Default: batch_core_analysis]',
	required = False,
	default = 'batch_core_analysis')
	
	
args = parser.parse_args()

OTUFP = args.OTU_table
OTU1,OTU2 = OTUFP.split(",")
metadataFP = args.metadata
MF1,MF2 = metadataFP.split(",")
align_log = args.align_log
outputFolder = args.outputFolder


#################################
# FUNCTIONS

# LOAD DATA

def loadOTUTableRaw(OTUFP): # Loads and DOES NOT convert to relative abundance
	# Make OTU Table
	os.system('biom convert -i ' + OTUFP + ' --to-tsv --header-key taxonomy -o OTU_Table_text.txt')
	# Load OTU Table
	biomOpen = open("OTU_Table_text.txt", 'U') # This is in Unix format
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
				OTUTable[str(biomTemp[lineN][0])][headers[abund]] = biomTemp[lineN][1:][abund]
	os.system("rm OTU_Table_text.txt")
	return OTUTable,taxaIDs # Output is a 2-layer dictionary; first is OTU IDs and second is samples. Also, one-layer dict with taxa IDs
	
def loadMetadata(metadataFP):
	metadataOpen = open(metadataFP, 'U') # U is for 'Universal read'-- automatically turns into Unix LF
	metadataTemp = []
	for i in metadataOpen:
		lineTemp = i.strip()
		lineTempSplit = lineTemp.split("\t")
		metadataTemp.append(lineTempSplit)
	metadataOpen.close()
	metadata = {}
	metadataSites = []
	for lineN in range(len(metadataTemp)):
		if lineN == 0:
			headerList = metadataTemp[lineN]
			for headerName in metadataTemp[lineN]:
				metadata[headerName] = {}
		else:
			for i in range(1,len(metadataTemp[lineN])):
				metadataSites.append(metadataTemp[lineN][0])
				sortHeader = headerList[i]
				metadata[sortHeader][metadataTemp[lineN][0]] = metadataTemp[lineN][i]
	return metadata # output is 2-layer dictionary: first is Metadata and second is samples

def loadAligmentLog(align_log):
	# Open file
	AL = open(align_log,'U')
	# Save line by line
	alignedData = []
	first = True
	for line in AL:
		templine = line.strip().split("\t")
		# Save number of columns
		if first:
			ncolumns = len(templine)
			alignedData.append(templine)
		# if number of columns == ncolumn, then save
		elif len(templine) == ncolumns:
			alignedData.append(templine)
		# number of clumns != ncolumn, then add None until it does
		else:
			for i in range(ncolumns-len(templine)):
				templine.append(None)
			alignedData.append(templine)
		first = False
	AL.close()
	# Take header and find index of things we want
	otu1header = alignedData[0].index('candidate sequence ID')
	otu2header = alignedData[0].index('template ID')
	percIDheader = alignedData[0].index('BLAST percent identity to template')
	# Now, for all matches (with percIDheader), make it into a two column reference table
	first = True
	otumap = {'otu1_first':{}, 'otu2_first':{}}
	for line in alignedData:
		if not first and line[percIDheader]:
			otumap['otu1_first'][line[otu1header]] = line[otu2header]
		if not first and line[percIDheader]:
			otumap['otu2_first'][line[otu2header]] = line[otu1header]
		first = False
	return(otumap)		

	
def printTableFromDictionary(dictionary, output):
	toPrint = ''
	first = True
	for row in dictionary:
		if first == True:
			for column in dictionary[row]:
				toPrint += "\t" + str(column)
				first = False
			toPrint += "\n"
		toPrint += str(row)
		for column in dictionary[row]:
			toPrint += "\t" + str(dictionary[row][column])
		toPrint += "\n"
	open(str(output)+".txt", 'w').write(toPrint)
	print "DONE"
	
def printBiomFromDictionary(otutable, taxaIDs, output): # prints a biom file and text file from dictionary
	toPrint = '#OTU ID'
	first = True
	for row in otutable:
		if first == True:
			for column in otutable[row]:
				toPrint += "\t" + str(column)
			first = False
			toPrint += "\ttaxonomy\n"
		toPrint += str(row)
		for column in otutable[row]:
			toPrint += "\t" + str(otutable[row][column])
		toPrint += "\t" + taxaIDs[row] + "\n"
	open(str(output)+".txt", 'w').write(toPrint)
	os.system('biom convert -i ' + str(output) + '.txt --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy -o ' + str(output) + ".biom")
	print "DONE"
				

#################################
align_log = '/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/Feb_Neighbours_work_with3analysis/pynast_aligned/MC_neighbours_OTUreps_log.txt'
OTU1 = '/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/Feb_Neighbours_work_with3analysis/OTU_Table_nochlpmito_m800.biom'
OTU2 = ''
MF1 = '/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/Feb_Neighbours_work_with3analysis/OTU_Table_nochlpmito_m800.biom /Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/Feb_Neighbours_work_with3analysis/Merged_mapping_file.txt'
MF2 = ''
outputFolder = "TEST"

os.system(str("mkdir " + outputFolder))

# Load OTU tables
OTUTable1,taxaIDs1 = loadOTUTableRaw(OTU1)
OTUTable2,taxaIDs2 = loadOTUTableRaw(OTU2)
print "DONE LOADING OTU TABLE"

otumap = loadAligmentLog(align_log)
megaOTU = {}
# Keep track of how many seqs are lost
lostOTU1 = 0
lostOTU2 = 0
for OTU in otumap['otu1_first'].keys():
	metaOTU[OTU] = {}
	otu1ID = OTU
	otu2ID = otumap['otu1_first'][OTU]
	if otu1ID in OTUTable1.keys():
		megaOTU[OTU] = OTUTable1[otu1ID]
	else:
		lostOTU1 += 1
	if otu2ID in OTUTable2.keys():
		megaOTU[OTU].update(OTUTable2[otu2ID])
	else:
		lostOTU2 += 1
		
# Print OTU table
printBiomFromDictionary(megaOTURaw, taxaIDs1, str(outputFolder + "/merged_otutable"))
		
#################################
