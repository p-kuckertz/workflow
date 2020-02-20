# coding: utf-8

# Import modules
import glaes as gl
import pandas as pd
PARAMS = pd.read_csv('input_data/landEligibilityParameters.csv', index_col='name')['value']
gl.Priors.loadDirectory("input_data")

# INPUTS
regionFile = "data/gadm36_DEU_1.shp"
regionID = int(PARAMS.loc['regionID'])
settlement_proximity = int(PARAMS.loc['exc_settlement_proximity'])

# OUTPUTS
outputFile = "output_data/eligibility_result.tif"

# Initiate Exclusion Calculator object
ec = gl.ExclusionCalculator(regionFile, where=regionID)

# Do a simple exclusion
ec.excludePrior("settlement_proximity", value=(0,settlement_proximity)) # Exclude urban areas

# TODO: Add into inputs
ec.excludePrior("roads_main_proximity", value=(0,200)) # Exclude urban areas
ec.excludePrior("roads_secondary_proximity", value=(0,200)) # Exclude urban areas
ec.excludePrior("protected_landscape_proximity", value=(0,1000)) # Exclude urban areas

# Save result
ec.save(outputFile)
