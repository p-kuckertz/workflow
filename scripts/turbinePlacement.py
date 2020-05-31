# coding: utf-8

# TODO: Why is the substitution of the reagion and eligibility file paths not working?
#       Why is there no ".loc" at turbineSeparation?

# Import modules
import glaes as gl
import pandas as pd

# Load input files
parameterFile = pd.read_csv(snakemake.input[0], index_col='name')['value']
regionFile = "output_data/gadm36_DEU_1.shp"
eligibilityFile = "output_data/eligibility_result.tif"

# Set parameters according to parameter file
regionID = int(parameterFile.loc['regionID'])
turbineSeparation = float(parameterFile['turbineSeparation'])

# Initiate exclusion calculator object
ec = gl.ExclusionCalculator(regionFile, where=regionID)

# Load land eligibility result
ec.excludeRasterType(eligibilityFile, value=0)

# Execute placement algorithm and save output file
items = ec.distributeItems(turbineSeparation, output=snakemake.output[0])
