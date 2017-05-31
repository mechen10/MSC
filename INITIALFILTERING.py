#!/bin/bash

import argparse
import subprocess
import os



#=========================================================

parser = argparse.ArgumentParser(
	description="Initial analyzing and processing steps after creating OTU table")
parser.add_argument(
	'-i',
	'--input_biom',
	help = "Biom file for input",
	required = True)
# parser.add_argument(
# 	'-o',
# 	'--output_tag',
# 	help = 'Tag to add onto the end of all outputs',
# 	required = True)
parser.add_argument(
	'-t',
	'--tree_fp',
	help = 'Tree for phylogenetic methods. [Default: None]',
	required = False,
	default = None)
parser.add_argument(
	'-A',
	'--Alpha_parameter_fp',
	help = 'Parameter file for alpha, if desired. [Default: None]',
	default = None)
parser.add_argument(
	'-B',
	'--Beta_parameter_fp',
	help = 'Parameter file for beta, if desired. [Default: None]',
	default = None)
parser.add_argument(
	'-s',
	'--samples_tokeep',
	help = " use formatting, 'factor:type' to designate which controls to delete",
	default = None)
parser.add_argument(
	'-m',
	'--metadata',
	help = 'metadata file to delete controls',
	required = True)
parser.add_argument(
	'-n',
	'--min_persample',
	help = 'minimum number of outs in sample for it to be retained',
	required = False,
	default = 1000)
parser.add_argument(
	'-c',
	'--collapse_fields',
	help = 'fields to collapse on; name at least one header. [Default:ColRep]',
	default = 'ColRep',
	required = False)
	
args = parser.parse_args()

input_biom = args.input_biom
# output_tag = args.output_tag
tree_fp = args.tree_fp
Alpha_parameter_fp = args.Alpha_parameter_fp
Beta_parameter_fp = args.Beta_parameter_fp
samples_tokeep = args.samples_tokeep
metadata = args.metadata
min_persample = args.min_persample
collapse_fields = args.collapse_fields


#=========================================================

# First thing we do is delete the samples that we don't want; aka, chloroplasts, mitochondria, etc

print '\nFiltering Chloroplast, Mitocondria out of data...\n'

subprocess.call(['biom','convert','-i',input_biom,'--header-key','taxonomy','--to-tsv','-oTODELETE.txt'])
os.system('grep -v "Chloroplast" TODELETE.txt | grep -v "Mitochondria" >> TOKEEP.txt')
subprocess.call(['filter_otus_from_otu_table.py','-i',input_biom,'-oOTU_Table_nochlpmito.biom','-eTOKEEP.txt','--negate_ids_to_exclude'])
os.system('rm TODELETE.txt')

# Then we get rid of controls

print 'Filtering controls and others from data...\n'
subprocess.call(['filter_samples_from_otu_table.py','-iOTU_Table_nochlpmito.biom','-oOTU_Table_nochlpmito_noCon_m1000.biom','-m',metadata,'-s',samples_tokeep,'--output_mapping_fp','MF_All_m1000.txt','-n',min_persample])
os.system('rm OTU_Table_nochlpmito.biom')


# Now, we start doing OTU-table making.

print 'Making new OTUtables with filters...\n'
os.system('mkdir OTU_Tables_and_MF')

os.system('mv OTU_Table_nochlpmito_noCon_m1000.biom ./OTU_Tables_and_MF')
os.system('mv MF_All_m1000.txt ./OTU_Tables_and_MF')

subprocess.call(['collapse_samples.py','-b./OTU_Tables_and_MF/OTU_Table_nochlpmito_noCon_m1000.biom','--output_biom_fp','./OTU_Tables_and_MF/OTU_Table_nochlpmito_noCon_m1000_COLL.biom','-m./OTU_Tables_and_MF/MF_All_m1000.txt','--output_mapping_fp','./OTU_Tables_and_MF/MF_All_m1000_COLL.txt','--collapse_fields',collapse_fields])

subprocess.call(['touch',Alpha_parameter_fp])

print 'Running alpha_rarefaction.py...\n'
subprocess.call(['alpha_rarefaction.py','-i./OTU_Tables_and_MF/OTU_Table_nochlpmito_noCon_m1000.biom','-oalpha_rare','-p',Alpha_parameter_fp,'-m./OTU_Tables_and_MF/MF_All_m1000.txt','-t',tree_fp])

print 'DONE'
	
	
