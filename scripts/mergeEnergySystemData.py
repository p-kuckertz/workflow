# coding: utf-8

# TODO: Get scenarioInput.xlsx from GitHub repository.
#       Is there a more elegant way?

# Import modules
import pandas as pd

# Load input files
turbineSimulationResultFile = snakemake.input[0]
energySystemOptimizationInputData = snakemake.input[1]

# Load data
# Load energy system optimization input data 
file = pd.ExcelFile(energySystemOptimizationInputData)
# Set indices of spreadsheets (to avoid empty index cells caused by merged cells within Excel)
df0 = pd.read_excel(file, sheet_name='EnergySystemModel', index_col=[0])
df1 = pd.read_excel(file, sheet_name='Misc', index_col=[0,1,2])
df2 = pd.read_excel(file, sheet_name='Conversion', index_col=[0])
df3 = pd.read_excel(file, sheet_name='Sink', index_col=[0])
df4 = pd.read_excel(file, sheet_name='SinkTimeSeries', index_col=[0,1,2])
df5 = pd.read_excel(file, sheet_name='Source', index_col=[0])
df6 = pd.read_excel(file, sheet_name='SourceLocSpecs', index_col=[0,1,2])
df7 = pd.read_excel(file, sheet_name='SourceTimeSeries', index_col=[0,1,2])
df8 = pd.read_excel(file, sheet_name='Storage', index_col=[0])
df9 = pd.read_excel(file, sheet_name='StorageLocSpecs', index_col=[0,1,2])
df10 = pd.read_excel(file, sheet_name='Transmission', index_col=[0])
df11 = pd.read_excel(file, sheet_name='TransmissionLocSpecs', index_col=[0,1,2])
# Load turbine simulation results
df = pd.read_csv(turbineSimulationResultFile)

# Merge data
# Turbine simulation result data is put into energy system optimization input data:
# Within the sheet "SourceTimeSeries" the line specified by "Wind (onshore)", "operationRateMax" and "cluster_0" is overwritten (completely from row index 0 to row index 8759)
df7.iloc[0,:] = df.iloc[:, 1]

# Save result
with pd.ExcelWriter(snakemake.output[0]) as fileout:
    df0.to_excel(fileout, sheet_name='EnergySystemModel')
    df1.to_excel(fileout, sheet_name='Misc')
    df2.to_excel(fileout, sheet_name='Conversion')
    df3.to_excel(fileout, sheet_name='Sink')
    df4.to_excel(fileout, sheet_name='SinkTimeSeries')
    df5.to_excel(fileout, sheet_name='Source')
    df6.to_excel(fileout, sheet_name='SourceLocSpecs')
    df7.to_excel(fileout, sheet_name='SourceTimeSeries')
    df8.to_excel(fileout, sheet_name='Storage')
    df9.to_excel(fileout, sheet_name='StorageLocSpecs')
    df10.to_excel(fileout, sheet_name='Transmission')
    df11.to_excel(fileout, sheet_name='TransmissionLocSpecs')
