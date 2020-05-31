# coding: utf-8

# TODO: FINE does not accept file ending for result file name. So "snakemake.output[0]" cannot be used. How to circumvent?

# Import modules
import os
import FINE as fn 

# Load input files
mergedEnergySystemScenarioDataFile = snakemake.input[0]

# Run energy sistem optimization "from Excel"
esM = fn.energySystemModelRunFromExcel(mergedEnergySystemScenarioDataFile)
# Remove unwanted output file (standard name and location do not suit this setup)
if os.path.exists("scenarioResults.xlsx"):
    os.remove("scenarioResults.xlsx")

# Save result
fn.writeOptimizationOutputToExcel(esM, outputFileName='output_data/energySystemOptimizationResults')
