# coding: utf-8

# Import modules
import geokit as gk
from reskit import windpower
import pandas as pd
import numpy as np
PARAMS = pd.read_csv( "input_data/turbineSimulationParameters.csv", index_col='name')['value']

# INPUT FILES
gwaFile = "/data/s-ryberg/data/geography/global_wind_atlas/v3/gwa3_250_wind-speed_100m.tif"
placementFile = "output_data/placements.shp"
turbineDesign = windpower.TurbineLibrary.loc[PARAMS['turbineDesign']]

# OUTPUTS
outputFile = "output_data/generation.csv"

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
generation.to_csv(outputFile)
