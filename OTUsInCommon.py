#!/bin/bash/python
#### OTUs IN COMMON ###

# Built under QIIME 1.9.0 and python 2.7.3
# Works as of 13Feb2017

# Script for determining shared microbiomes in different treatment groups all at once. 



import argparse
import os
import sys
import copy


#################################

parser = argparse.ArgumentParser(
	description="Takes OTU table, metadata, and column:treatment to find shared OTUs in each treatment. ")
parser.add_argument(
	'-i',
	'--OTU_table',
	help = "OTU table from QIIME containing all samples",
	required = True,)
parser.add_argument(
	'-m',
	'--metadata',
	help = 'Mapping file for OTU table',
	required = True)
parser.add_argument(
	'-c',
	'--column_name',
	help = 'Column name for treatment groups in metadata',
	type = str,
	required = True)
parser.add_argument(
	'--groups',
	help = 'Comma separated list of treatment groups to include: if not included, will do all treatment groups',
	required = False,
	default = 'False')
parser.add_argument(
	'-o',
	'--outputFolder',
	help = 'Output folder [Default: batch_core_analysis]',
	required = False,
	default = 'batch_core_analysis')
	
	
args = parser.parse_args()

OTUFP = args.OTU_table
metadataFP = args.metadata
columnName = args.column_name
groups = args.groups
outputFolder = args.outputFolder


#################################
# FUNCTIONS

# LOAD DATA

def loadOTUTable(OTUFP): # Loads and converts to relative abundance
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
	
def makeOTUandMetaSame(OTUtable, Metadata): # filters out samples that are not in respective files
	otuSamplesToDel = []
	metaSamplesToDel = []
	otutable = OTUtable.copy()
	metadata = Metadata.copy()
	for sample in otutable[otutable.keys()[1]]:
		if sample not in metadata[metadata.keys()[1]].keys():
			otuSamplesToDel.append(sample)
	for sample in metadata[metadata.keys()[1]].keys():
		if sample not in otutable[otutable.keys()[1]].keys():
			metaSamplesToDel.append(sample)
	for sample in otuSamplesToDel:
		for OTU in otutable.keys():
			del otutable[OTU][sample]
	for sample in metaSamplesToDel:
		for environ in metadata.keys():
			if environ != '#SampleID':
				del metadata[environ][sample]
	return(otutable, metadata)
	
def printTableFromDictionary(dictionary, output):
	toPrint = ''
	first = True
	for row in sorted(dictionary):
		if first == True:
			for column in sorted(dictionary[row]):
				toPrint += "\t" + str(column)
				first = False
			toPrint += "\n"
		toPrint += str(row)
		for column in sorted(dictionary[row]):
			toPrint += "\t" + str(dictionary[row][column])
		toPrint += "\n"
	open(str(output)+".txt", 'w').write(toPrint)
	print "DONE"

def getColNames(metadata, columnName):
	allTreatmentList = {}
# 	tempNames = []
	for sample in metadata[columnName].keys():
# 		tempNames.append(metadata[columnName][sample])
		if not metadata[columnName][sample] in allTreatmentList:
			allTreatmentList[metadata[columnName][sample]] = []
		allTreatmentList[metadata[columnName][sample]].append(sample)
	return(allTreatmentList)
	
def getColNames2(metadata, columnName, listNames):
	allTreatmentList = {}
	for treatment in listNames:
		allTreatmentList[treatment] = []
		for sample in metadata[columnName].keys():
			if metadata[columnName][sample] == treatment:
				allTreatmentList[treatment].append(sample)
	return allTreatmentList
	
def getShared(OTUTable, metadata, columnName, Treatment1, Treatment2):
	# list of shared OTU for each pair of samples
	# Get Treatment
	AllColumn = metadata[columnName]
	# Get all sample names in desired treatment
	sampleList1 = []
	for i in AllColumn:
		if AllColumn[i] == Treatment1:
			sampleList1.append(i)
	sampleList2 = []
	for i in AllColumn:
		if AllColumn[i] == Treatment2:
			sampleList2.append(i)
	# Get these samples from the OTU table
	OTUShare = []
	for OTU in OTUTable.keys():
		shareTest1 = 0
		shareTest2 = 0
		for sample in sampleList1:
			if OTUTable[OTU][sample] > 0:
				shareTest1 += 1
		for sample in sampleList2:
			if OTUTable[OTU][sample] > 0:
				shareTest2 += 1
		if shareTest1 >= 1 and shareTest2 >=1:
			OTUShare.append(OTU)
	return(OTUShare)
	
def getTotalNumOTUs(OTUTable, metadata, columnName, Treatment):
	# Get total number of OTUs in sample
	# Get Treatment
	AllColumn = metadata[columnName]
	# Get all sample names in desired treatment
	sampleList = []
	for i in AllColumn:
		if AllColumn[i] == Treatment:
			sampleList.append(i)
	# Get these samples from the OTU table
	OTUtotal = []
	for OTU in OTUTable.keys():
		presence = 0
		counter = 0
		while presence == 0 and counter < len(sampleList):
			sample = sampleList[counter]
			if OTUTable[OTU][sample] > 0:
				presence += 1
			counter += 1
		if presence > 0:
			OTUtotal.append(OTU)
	return(OTUtotal)

def getAveAbundInGroup(OTUTable, metadata, columnName, Treatment):
	AllColumn = metadata[columnName]
	# Get all sample names in desired treatment
	sampleList = []
	for i in AllColumn:
		if AllColumn[i] == Treatment:
			sampleList.append(i)
	# Get these samples from the OTU table
	SampleOTUCount = {}
	for OTU in OTUTable.keys():
		Readsubtotal = 0
		for sample in sampleList:
			Readsubtotal += OTUTable[OTU][sample]
		SampleOTUCount[OTU] = Readsubtotal/(len(sampleList))
	ReadTotal = 0
	for OTU in SampleOTUCount.keys():
		ReadTotal += SampleOTUCount[OTU]
	return(SampleOTUCount,ReadTotal)

#################################


os.system(str("mkdir " + outputFolder))

OTUTableFull,taxaIDs = loadOTUTable(OTUFP)
print "DONE LOADING OTU TABLE"
metadataRaw = loadMetadata(metadataFP)
print "DONE LOADING METADATA"

OTUTable,metadata = makeOTUandMetaSame(OTUTableFull, metadataRaw)

# Get list of groups
if groups == 'False':
	allTreatmentList = getColNames(metadata, columnName)
else:
	listNames = groups.split(",")
	allTreatmentList = getColNames2(metadata, columnName, listNames)

# Calculate shared OTUs for each

ALLSHARED = {}
for Treatment1 in allTreatmentList:
	ALLSHARED[Treatment1] = {}
	for Treatment2 in allTreatmentList:
		ALLSHARED[Treatment1][Treatment2] = getShared(OTUTable, metadata, columnName, Treatment1, Treatment2)

# Calculate total number of OTUs each
TOTALOTUS = {}
for Treatment in allTreatmentList:
	TOTALOTUS[Treatment] = getTotalNumOTUs(OTUTable, metadata, columnName, Treatment)
	
#################################

# Sample matrix with raw OTU numbers
sharedOTUsmatrix = {}
for Treatment1 in allTreatmentList:
	sharedOTUsmatrix[str(Treatment1)+'('+str(len(TOTALOTUS[Treatment1]))+')'] = {}
	for Treatment2 in allTreatmentList:
		sharedOTUsmatrix[str(Treatment1)+'('+str(len(TOTALOTUS[Treatment1]))+')'][str(Treatment2)+'('+str(len(TOTALOTUS[Treatment2]))+')'] = len(ALLSHARED[Treatment1][Treatment2])
	
# Sample matrix with OTUs as percentage of horizontal row
sharedOTUspercent = {}
for Treatment1 in allTreatmentList:
	sharedOTUspercent[str(Treatment1)+'('+str(len(TOTALOTUS[Treatment1]))+')'] = {}
	for Treatment2 in allTreatmentList:
		sharedOTUspercent[str(Treatment1)+'('+str(len(TOTALOTUS[Treatment1]))+')'][Treatment2] = len(ALLSHARED[Treatment1][Treatment2])/float(len(TOTALOTUS[Treatment1]))

# By horizontal row: what amount is represented by the shared OTUs
sharedOTUsproportion = {}
for Treatment1 in allTreatmentList:
	aveAbundTemp,total = getAveAbundInGroup(OTUTable, metadata, columnName, Treatment1)
	sharedOTUsproportion[str(Treatment1)+'('+str(len(TOTALOTUS[Treatment1]))+')'] = {}
	for Treatment2 in allTreatmentList:
		totalPerc = 0
		for OTU in ALLSHARED[Treatment1][Treatment2]:
			totalcount = aveAbundTemp[OTU]
			totalPerc += totalcount/total
		sharedOTUsproportion[str(Treatment1)+'('+str(len(TOTALOTUS[Treatment1]))+')'][Treatment2] = totalPerc

		
# Print all of these	
printTableFromDictionary(sharedOTUsmatrix, str(outputFolder + "/SHAREDMATRIX_rawnumbers"))
printTableFromDictionary(sharedOTUspercent, str(outputFolder + "/SHAREDMATRIX_percentByRows"))
printTableFromDictionary(sharedOTUsproportion, str(outputFolder + "/SHAREDMATRIX_proportionByRows"))

#################################

# Print taxaIDLegend

toWrite = "OTUID\ttaxonomy\n"
for i in taxaIDs.keys():
	toWrite += i + "\t" + str(taxaIDs[i]) + "\n"
open(str(outputFolder + "/taxaIDLegend.txt"), 'w').write(toWrite)

#################################

# Optional: print OTUtable in python formatted dictionary?

