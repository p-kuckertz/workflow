from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

rule all:
    input: "output_data/generation.csv"

rule fetchResources:
    input:
        region=HTTP.remote("biogeo.ucdavis.edu/data/gadm3.6/shp/gadm36_DEU_shp.zip", keep_local=True)
    output:
        # Shapefile
        "output_data/gadm36_DEU_1.shp",
        "output_data/gadm36_DEU_1.shx",
        "output_data/gadm36_DEU_1.prj",
        "output_data/gadm36_DEU_1.dbf",
        "output_data/gadm36_DEU_1.cpg",
    shell:
        "unzip -o {input.region} -d output_data"
    
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
        "input_data/airport_proximity.1497341471.tif",
        "input_data/protected_biosphere_proximity.1497677528.tif",
        "input_data/protected_wilderness_proximity.1497687503.tif",
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
        # Shapefile
        "output_data/gadm36_DEU_1.shp",
        "output_data/gadm36_DEU_1.shx",
        "output_data/gadm36_DEU_1.prj",
        "output_data/gadm36_DEU_1.dbf",
        "output_data/gadm36_DEU_1.cpg",
        # Land eligibility result
        "output_data/eligibility_result.tif",
    output:
        # Turbine placement Result
        "output_data/placements.shp",
    conda:
        "envs/renewables-env.yml"
    script:
        "scripts/turbinePlacement.py"
        
rule turbineSimulation:
    input:
	    # Parameters
        "input_data/turbineSimulationParameters.csv",
        # Data Files
        "/data/s-ryberg/data/geography/global_wind_atlas/v3/gwa3_250_wind-speed_100m.tif",
        # Turbine placement Result
        "output_data/placements.shp",
    output:
        # Turbine simulation result
        "output_data/generation.csv",
    conda:
        "envs/renewables-env.yml"
    script:
        "scripts/turbineSimulation.py"