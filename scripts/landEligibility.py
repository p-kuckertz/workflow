# coding: utf-8

# TODO: The params file should only use the parameters of the glaes interface.
#       So we can reference the glaes interface documentation to explain parameter meaning and usage.
#       "proximity_start" is temporarily added to avoid hardcoding 0s.
#		Remove path from "loadDirectory"

# Import modules
import glaes as gl
import pandas as pd

# Load input files
parameterFile = pd.read_csv(snakemake.input[0], index_col='name')['value']
regionFile = snakemake.input[1]
# Load directory of prior dataset
gl.Priors.loadDirectory("input_data/priorDataSet")

# Set parameters according to parameter file
regionID = int(parameterFile.loc['regionID'])
proximity_start = int(parameterFile.loc['proximity_start'])
settlement_proximity = int(parameterFile.loc['settlement_proximity'])
roads_main_proximity = int(parameterFile.loc['roads_main_proximity'])
roads_secondary_proximity = int(parameterFile.loc['roads_secondary_proximity'])
protected_landscape_proximity = int(parameterFile.loc['protected_landscape_proximity'])

# Initiate exclusion calculator object
ec = gl.ExclusionCalculator(regionFile, where=regionID)

# Exclude areas
# Exclude urban areas
ec.excludePrior("settlement_proximity", value=(proximity_start,settlement_proximity))
# Exclude urban areas
ec.excludePrior("roads_main_proximity", value=(proximity_start,roads_main_proximity))
# Exclude urban areas
ec.excludePrior("roads_secondary_proximity", value=(proximity_start,roads_secondary_proximity))
# Exclude urban areas
ec.excludePrior("protected_landscape_proximity", value=(proximity_start,protected_landscape_proximity))

# Save result
ec.save(snakemake.output[0])
