# coding: utf-8

# TODO: Why is there no ".loc" at turbineSeparation?

# Import modules
import glaes as gl
import pandas as pd

# Load input files
parameterFile = pd.read_csv(snakemake.input[0], index_col='name')['value']
eligibilityFile = snakemake.input[1]
regionFile = snakemake.input[2]

# Set parameters according to parameter file
regionID = int(parameterFile.loc['regionID'])
turbineSeparation = float(parameterFile['turbineSeparation'])

# Initiate exclusion calculator object
ec = gl.ExclusionCalculator(regionFile, where=regionID)

# Load land eligibility result
ec.excludeRasterType(eligibilityFile, value=0)

# Execute placement algorithm and save output file
items = ec.distributeItems(turbineSeparation, output=snakemake.output[0])
