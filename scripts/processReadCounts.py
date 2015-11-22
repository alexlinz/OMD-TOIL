###############################################################################
# processReadCounts.py
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Process read count data from mapping of OMD-TOIL MT reads to our reference
# genomes.
################################################################################

#%%#############################################################################
### Import packages
################################################################################

import os
import pandas as pd
import re

#%%#############################################################################
### Static folder structure
################################################################################
# Define fixed input and output files
pathToReads = '../data/readCounts'

#%%#############################################################################
### Step 0 - Populate MT read frame. Create empty dataframes.
################################################################################

# Read in MT reads
mtReads = pd.read_csv('../data/miscData/totalReads.csv', index_col=0, sep=',')

# Create empty dataframe for genome read counts
alignedMatrix = pd.DataFrame(columns=[mtReads.index])
unalignedMatrix = pd.DataFrame(columns=[mtReads.index])

for genomeFile in os.listdir(pathToReads):
    if genomeFile.endswith('.out'): 
        MT = re.search('(.+)-(.+).out', genomeFile).group(1)
        genome = re.search('(.+)-(.+).out', genomeFile).group(2)
# Some .out files may be empty, so this section needs to be embedded with try/except loop
        try:
            genomeReads = pd.read_csv(pathToReads+'/'+genomeFile, index_col=0, sep='\t')
# Grab the 'unaliged' row and add to the 'unaligned' DF
# Grab the 'total reads' from the 'mtReads' DF and subtract the 'unaligned.' Add to the 'aligned' DF.
# Normalize all DFs
            unalignedMatrix.loc[genome, MT] = (genomeReads.loc['__not_aligned'][0]) / mtReads.loc[MT][0].astype(float)
            alignedMatrix.loc[genome, MT] = (mtReads.loc[MT][0] - genomeReads.loc['__not_aligned'][0]) / mtReads.loc[MT][0].astype(float)
        except:
            unalignedMatrix.loc[genome, MT] = mtReads.loc[MT][0] / mtReads.loc[MT][0].astype(float)
            alignedMatrix.loc[genome, MT] = 0

# Convert to a percent
alignedMatrix = alignedMatrix*100
unalignedMatrix = unalignedMatrix*100
# Write to CSV file
alignedMatrix.to_csv('../data/miscData/mappedReads.csv', sep=',')
