# coding: utf-8

# Import modules
import glaes as gl
import pandas as pd
PARAMS = pd.read_csv( "input_data/turbinePlacementParameters.csv", index_col='name')['value']

# INPUT FILES
regionFile = "output_data/gadm36_DEU_1.shp"
regionID = int(PARAMS.loc['regionID'])
eligibilityFile = "output_data/eligibility_result.tif"
turbineSeparation = float(PARAMS['turbineSeparation'])

# OUTPUTS
outputFile = "output_data/placements.shp"

# Initiate Exclusion Calculator object
ec = gl.ExclusionCalculator(regionFile, where=regionID)

# load eligibility result from RULE 1
ec.excludeRasterType(eligibilityFile, value=0)

# Do Placement algorithm
items = ec.distributeItems(turbineSeparation, output=outputFile)
