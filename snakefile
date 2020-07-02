from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

rule all:
    input: "output_data/energySystemOptimizationResults.xlsx"

rule fetchRegionData:
    input:
        region=HTTP.remote("biogeo.ucdavis.edu/data/gadm3.6/shp/gadm36_DEU_shp.zip", keep_local=False)
    output:
        # Shapefile
        "output_data/gadm36_DEU_1.shp",
        "output_data/gadm36_DEU_1.shx",
        "output_data/gadm36_DEU_1.prj",
        "output_data/gadm36_DEU_1.dbf",
        "output_data/gadm36_DEU_1.cpg",
    shell:
        ("unzip -o {input.region} -d output_data"
        " && rm output_data/gadm36_DEU_0.*"
        " && rm output_data/gadm36_DEU_2.*"
        " && rm output_data/gadm36_DEU_3.*"
        " && rm output_data/gadm36_DEU_4.*")
    
rule landEligibility:
    input:
        # Parameters
        "input_data/landEligibilityParameters.csv",
        # Shapefile
        "output_data/gadm36_DEU_1.shp",
        "output_data/gadm36_DEU_1.shx",
        "output_data/gadm36_DEU_1.prj",
        "output_data/gadm36_DEU_1.dbf",
        "output_data/gadm36_DEU_1.cpg",
        # Geospatial constraints
        "input_data/priorDataSet/airport_proximity.1497341471.tif",
        "input_data/priorDataSet/protected_biosphere_proximity.1497677528.tif",
        "input_data/priorDataSet/protected_wilderness_proximity.1497687503.tif",
    output:
        # Land eligibility result
        "output_data/eligibility_result.tif"
    conda:
        "envs/renewables-env.yml"
    script:
        "scripts/landEligibility.py"
    
rule turbinePlacement:
    input:
        # Parameters
        "input_data/turbinePlacementParameters.csv",
        # Land eligibility result
        "output_data/eligibility_result.tif",
        # Shapefile
        "output_data/gadm36_DEU_1.shp",
        "output_data/gadm36_DEU_1.shx",
        "output_data/gadm36_DEU_1.prj",
        "output_data/gadm36_DEU_1.dbf",
        "output_data/gadm36_DEU_1.cpg",
    output:
        # Turbine placement Result
        "output_data/placements.shp",
        "output_data/placements.shx",
        "output_data/placements.prj",
        "output_data/placements.dbf",
    conda:
        "envs/renewables-env.yml"
    script:
        "scripts/turbinePlacement.py"

rule fetchWeatherData:
    input:
        weatherData=HTTP.remote("https://globalwindatlas3.s3-eu-west-1.amazonaws.com/country_tifs/DEU_wind-speed_100m.tif", keep_local=False)
    output:
        "output_data/DEU_wind-speed_100m.tif",
    shell:
        ("cp {input.weatherData} output_data")

rule turbineSimulation:
    input:
        # Parameters
        "input_data/turbineSimulationParameters.csv",
        # Data Files
        "output_data/DEU_wind-speed_100m.tif",
        # Turbine placement Result
        "output_data/placements.shp",
        "output_data/placements.shx",
        "output_data/placements.prj",
        "output_data/placements.dbf",        
    output:
        # Turbine simulation result
        "output_data/generation.csv",
    conda:
        "envs/renewables-env.yml"
    script:
        "scripts/turbineSimulation.py"

rule fetchScenarioData:
    input:
        # TODO: use FINE 2.x instead of 1.0.2 because of blanks in file path. Ask developers to tag a corresponding commit as v.2
        scenarioData=HTTP.remote("github.com/FZJ-IEK3-VSA/FINE/raw/v.1.0.2/examples/Model%20Run%20from%20Excel/scenarioInput.xlsx", keep_local=False)
    output:
        "output_data/scenarioInput.xlsx",
    shell:
        ("cp {input.scenarioData} output_data")

rule mergeEnergySystemData:
    input:
        # Turbine simulation result
        "output_data/generation.csv",
        # Energy system scenario data
        "output_data/scenarioInput.xlsx"
    output:
        "output_data/mergedEnergySystemScenarioData.xlsx"
    conda:
        "envs/renewables-env.yml"
    script:
        "scripts/mergeEnergySystemData.py"
        
rule energySystemOptimization:
    input:
        "output_data/mergedEnergySystemScenarioData.xlsx"
    output:
        "output_data/energySystemOptimizationResults.xlsx"
    conda:
        "envs/fine-env.yml"
    script:
        "scripts/energySystemOptimization.py"
