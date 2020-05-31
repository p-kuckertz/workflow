# coding: utf-8

# TODO: Substitute input file paths.
#       Get wind data from somewhere.
#       Do not synthesize (and hardcode) anything.

# Import modules
import geokit as gk
from reskit import windpower
import pandas as pd
import numpy as np

# Load input files
parameterFile = pd.read_csv(snakemake.input[0], index_col='name')['value']
gwaFile = "/data/s-ryberg/data/geography/global_wind_atlas/v3/gwa3_250_wind-speed_100m.tif"
placementFile = "output_data/placements.shp"

# Set parameters according to parameter file
turbineDesign = windpower.TurbineLibrary.loc[parameterFile['turbineDesign']]

# Get Placements
placements = gk.vector.extractFeatures(placementFile)

# Get mean windspeed for each placement
placements['ws100m'] = gk.raster.interpolateValues( gwaFile, placements.geom )

# Make synthetic wind speed data
locs = gk.LocationSet(placements.geom)
np.random.seed(0)
windspeedValues = pd.DataFrame(np.random.normal(placements['ws100m'], placements['ws100m']/4, (8760, placements.shape[0], )), columns=locs)

# Wind turbine simulation
generation = windpower.simulateTurbine(windspeedValues, powerCurve=turbineDesign.PowerCurve)

# Save result
generation.to_csv(snakemake.output[0])
