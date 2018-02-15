#!/bin/bash

#!/bin/bash/python
#### OTUs IN COMMON ###

# Built under QIIME 1.9.0 and python 2.7.3
# Working as of 13Feb2017; only one trial

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
	'--align_log',
	help = 'The log file from align_seqs; or, can be any tab delimited file with the following headers: "candidate sequence ID", "template ID", "BLAST percent identity to template". (The last one is used to tell whether the file has ended; just need something in that column ',
	required = True)
parser.add_argument(
	'-k',
	'--keep_unique_otus',
	help = 'If included, output will INCLUDE OTUs that are NOT in common between datasets. Be aware that this means you cannot differentiate between zero OTUs, and OTUs that were simply not in dataset.',
	action = 'store_true')
parser.add_argument(
	'-o',
	'--outputFolder',
	help = 'Output folder [Default: batch_core_analysis]',
	required = False,
	default = 'batch_core_analysis')
	
args = parser.parse_args()

OTUFP = args.OTU_table
OTU1,OTU2 = OTUFP.split(",")
align_log = args.align_log
keep = args.keep_unique_otus
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
			if line[otu2header] not in otumap['otu2_first'].keys():
				otumap['otu2_first'][line[otu2header]] = [line[otu1header]]
			else:
				otumap['otu2_first'][line[otu2header]].append(line[otu1header])
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
	
def printBiomFromDictionary(otutable, taxaIDs2, taxaIDs1, output): # prints a biom file and text file from dictionary
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
		newrow = row.replace('.1','')
		newrow = newrow.replace('.2','')
		if newrow in taxaIDs2.keys():
			toPrint += "\t" + taxaIDs2[newrow] + "\n"
		elif newrow in taxaIDs1.keys():
			toPrint += "\t" + taxaIDs1[newrow] + "\n"
		else:
			toPrint += "\t\n"
	open(str(output)+".txt", 'w').write(toPrint)
	os.system('biom convert -i ' + str(output) + '.txt --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy -o ' + str(output) + ".biom")
	print "DONE"
				

#################################

os.system(str("mkdir " + outputFolder))

# Load OTU tables
OTUTable1,taxaIDs1 = loadOTUTableRaw(OTU1)
OTUTable2,taxaIDs2 = loadOTUTableRaw(OTU2)
print "DONE LOADING OTU TABLE"

otumap = loadAligmentLog(align_log)
megaOTU = {}
# Find headers for both tables so we can put 0 OTUs for missing ones
table1Samples = OTUTable1[OTUTable1.keys()[2]].keys() # 2 is arbirary; just use any number < n
table2Samples = OTUTable2[OTUTable2.keys()[2]].keys()
zeros1 = dict.fromkeys(table1Samples,0)
zeros2 = dict.fromkeys(table2Samples,0)

if keep: # If you want to keep all OTUs, even ones that are not in common
	for OTU in otumap['otu2_first'].keys():
		usedseq1 = 0
		otu2ID = OTU
		otu1ID = otumap['otu2_first'][OTU]
		if otu2ID in OTUTable2.keys(): 
			megaOTU[OTU] = OTUTable2[otu2ID].copy() # copies table 2's numbers
		else:
			megaOTU[OTU] = zeros2
		sumseqvar1 = zeros1
		for seqvar in otu1ID: # now add all sequence variant counts	
			if seqvar in OTUTable1.keys():
				usedseq1 += 1
				sumseqvar1 = { k: float(sumseqvar1.get(k)) + float(OTUTable1[seqvar].get(k)) for k in set(sumseqvar1) | set(OTUTable1[seqvar]) }
		# Now, add those summed counts to the megaOTU.
		megaOTU[OTU].update(sumseqvar1)
	for OTU in OTUTable1.keys():
		if OTU not in otumap['otu1_first'].keys():
			megaOTU[OTU+'.1'] = OTUTable1[OTU].copy()
			megaOTU[OTU+'.1'].update(zeros2)
	for OTU in OTUTable2.keys():
		if OTU not in otumap['otu2_first'].keys():
			megaOTU[OTU+'.2'] = OTUTable2[OTU].copy()
			megaOTU[OTU+'.2'].update(zeros1)
else: # If you do NOT want to keep all OTUs; just the ones that are in common
	for OTU in otumap['otu2_first'].keys():
		usedseq1 = 0
		otu2ID = OTU
		otu1ID = otumap['otu2_first'][OTU]
		if otu2ID in OTUTable2.keys(): 
			megaOTU[OTU] = OTUTable2[otu2ID].copy() # copies table 2's numbers
		else:
			megaOTU[OTU] = zeros2
		sumseqvar1 = zeros1
		for seqvar in otu1ID: # now add all sequence variant counts	
			if seqvar in OTUTable1.keys():
				usedseq1 += 1
				sumseqvar1 = { k: float(sumseqvar1.get(k)) + float(OTUTable1[seqvar].get(k)) for k in set(sumseqvar1) | set(OTUTable1[seqvar]) }
		# Now, add those summed counts to the megaOTU.
		megaOTU[OTU].update(sumseqvar1)
	lostOTU1 = len(OTUTable1.keys()) - usedseq1
	lostOTU2 = len(OTUTable2.keys()) - len(megaOTU.keys())

if keep:
	open(outputFolder+"/LOG.txt", 'w').write('All OTUs kept; cannot differentiate between different zeros')
else:
	open(outputFolder+"/LOG.txt", 'w').write('OTUS lost from Table 1:'+str(lostOTU1)+'\nOTUs lost from Table 2:'+str(lostOTU2))

# Print OTU table
printBiomFromDictionary(megaOTU, taxaIDs2, taxaIDs1, str(outputFolder + "/merged_otutable"))
		
#################################
