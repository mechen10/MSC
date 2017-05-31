#!/macqiime/bin/python

# V2017-05-31

# Filters OTU for potential contaminants; produces an output file of contaminantes
# Does this by comparing NEG OTUs with sample OTUs; if an OTU is at least 0.5x the amount of otus in the largest sample OTU read count, it is marked

