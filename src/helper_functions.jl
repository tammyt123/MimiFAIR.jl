# #----------------------------------------------------------------------------------------------------------------------
# #----------------------------------------------------------------------------------------------------------------------
# This file contains functions and other snippets of code that are used in various calculations for Mimi-FAIR.
# #----------------------------------------------------------------------------------------------------------------------
# #----------------------------------------------------------------------------------------------------------------------


#######################################################################################################################
# CALCULATE SCALE FACTOR FOR CONSISTENT CO₂ RADIATIVE FORCING
########################################################################################################################
# Description: For Etminan et al. radiative forcing equations, this function calculates a scaling factor to ensure
#              modeled carbon dioxide radiative forcing remains consistent with the user-specified forcing from a doubling
#              of carbon dioxide parameter (F2x).
#
# Function Arguments:
#
#       F2x:   Radiative forcing from a doubling of carbon dioxide concentrations (Wm⁻²).
#       CO₂_0: Pre-industrial carbon dioxide concentrations (ppm).
#       N₂O_0: Pre-industrial nitrous oxide concentrations (ppb).
#----------------------------------------------------------------------------------------------------------------------

function co2_rf_scale(F2x::Float64, CO₂_0::Float64, N₂O_0::Float64)
    # Calcualte forcing from doubling of CO₂ (following approach from original FAIR Python code).
    co2_doubling = (-2.4e-7 * CO₂_0^2 + 7.2e-4 * CO₂_0 - 2.1e-4 * N₂O_0 + 5.36) * log(2)
    # Calculate scaling factor, given user supplied F2x parameter.
    CO₂_scale = F2x / co2_doubling
    return CO₂_scale
end



#######################################################################################################################
# CALCULATE RADIATIVE FORCING COEFFICIENTS GIVEN USER PARAMETER SETTINGS
########################################################################################################################
# Description: Given user parameter settings governing climate system behavior, calculate the two radiative forcing "q"
#              coefficients, which correspond to [1] thermal equilibration of deep ocean & [2] thermal admustment of
#              upper ocean (K W⁻¹m²).
#
# Function Arguments:
#
#       tcr:     Transient climate response (K).
#       ecs:     Equilibrium climate sensitivity (K).
#       d:       Two-element array of coefficients governing slow and fast thermal response times for upper and deep oceans (years).
#       f2x:     Effective radiative forcing from a doubling of carbon dioxide concentrations (Wm⁻²).
#       tcr_dbl: Time to a doubling of carbon dioxide concentrations under 1% per year carbon dioxide increase (years).
#----------------------------------------------------------------------------------------------------------------------

function calculate_q(tcr, ecs, d, f2x; tcr_dbl=69.661)

    k1 = 1.0 - (d[1]/tcr_dbl)*(1.0 - exp(-tcr_dbl/d[1]))
    k2 = 1.0 - (d[2]/tcr_dbl)*(1.0 - exp(-tcr_dbl/d[2]))

    q1 = (1.0/f2x) * (1.0 / (k1-k2)) * (tcr - ecs * k2)
    q2 = (1.0/f2x) * (1.0 / (k1-k2)) * (ecs * k1 - tcr)

    return [q1, q2]
end



#######################################################################################################################
# CALCULATE TROPOSPHERIC OZONE RADIATIVE FORCING TEMPERATURE FEEDBACK
########################################################################################################################
# Description: Calclate a small, negative temperature feedback on tropospheric ozone radiative forcing using a curve
#              fit to year 2000, 20230, and 2100 temperature changes under RCP8.5 in Stevenson et al. (2013) (10.5194/acp-13-3063-2013)
#
# Function Arguments:
#
#       temperature: Global mean surface temperature anomaly (K).
#----------------------------------------------------------------------------------------------------------------------

function temperature_feedback(temperature::Float64)
    # Assume no temperature feedback if temp < 0.
    if  temperature <= 0.0
        temperature_feedback = 0.0
    else
        temperature_feedback = 0.03189267 * exp(-1.34966941 * temperature) - 0.03214807
    end

    return temperature_feedback
end



#######################################################################################################################
# LOAD ALL DATA AND PARAMETER VALUES NEEDED FOR Mimi-FAIR
########################################################################################################################
# Description: This function loads and cleans up the Mimi-FAIR data so it can be more easily incorporated into the model.
#
# Function Arguments:
#
#       start_year:      First year to run the model.
#       end_year:        Final year to run the model.
#       rcp_scenario:    A string indicating which RCP scenario to use ("RCP26", "RCP45", "RCP60", & "RCP85").
#       usg_scenario:    A string indicating which USG scenario to use ("USG1", "USG2", "USG3", "USG4", & "USG5").
#----------------------------------------------------------------------------------------------------------------------

function load_fair_data(start_year::Int, end_year::Int, rcp_scenario::String, usg_scenario::String)

    # Calculate indicies to extract RCP data (RCP data spans 1765-2500)
    start_index, end_index = findall((in)([start_year, end_year]), collect(1765:2500))

    # Create vector of names for minor greenhouse gases to loop over.
    other_ghg_names = ["CF4", "C2F6", "C6F14", "HFC23", "HFC32", "HFC43_10", "HFC125", "HFC134a", "HFC143a", "HFC227ea", "HFC245fa", "SF6", "CFC_11", "CFC_12", "CFC_113", "CFC_114", "CFC_115", "CARB_TET", "MCF", "HCFC_22", "HCFC_141B", "HCFC_142B", "HALON1211", "HALON1202", "HALON1301", "HALON2402", "CH3BR", "CH3CL"]

    #---------------------------------------
    # Read in relevant data sets.
    #---------------------------------------

    # RCP scenario emissions.
    rcp_emissions_raw = DataFrame(load(joinpath(@__DIR__, "..", "data", "model_data", rcp_scenario*"_EMISSIONS.csv"), skiplines_begin=36))
    # rcp_emissions_raw = DataFrame(load(joinpath(@__DIR__, "data", "model_data", rcp_scenario*"_EMISSIONS.csv"), skiplines_begin=36)) # personal path
    # Natural emissions for methane and nitrous oxide as estimated by FAIR team.
    natural_emissions_raw = DataFrame(load(joinpath(@__DIR__, "..", "data", "model_data", "natural_emissions.csv"), skiplines_begin=3))
    # natural_emissions_raw = DataFrame(load(joinpath(@__DIR__, "data", "model_data", "natural_emissions.csv"), skiplines_begin=3)) # personal path
    # CMIP6 Solar forcing.
    cmip6_solar_forcing = DataFrame(load(joinpath(@__DIR__, "..", "data", "model_data", "cmip6_solar.csv"), skiplines_begin=6))[start_index:end_index, Symbol("Radiative forcing")]
    # CMIP6 volcanic forcing.
    cmip6_volcano_forcing = DataFrame(load(joinpath(@__DIR__, "..", "data", "model_data", "cmip6_volcanic.csv"), skiplines_begin=8))[start_index:end_index, Symbol("Volcanic ERF")]
    # Fraction of NOx emissions attibutable to aviation (for contrail RF).
    aviation_fraction_raw = DataFrame(load(joinpath(@__DIR__, "..", "data", "model_data", "aviNOx_fraction.csv"), skiplines_begin=4))[:,Symbol(rcp_scenario)]
    # Time-varying shares of anthropogenic methane attribuatable to fossil sources.
    ch4_fossil_frac_raw = DataFrame(load(joinpath(@__DIR__, "..", "data", "model_data", "fossilCH4_fraction.csv"), skiplines_begin=4))[:,Symbol(rcp_scenario)]
    # Information on various gas specieis.
    gas_data = DataFrame(load(joinpath(@__DIR__, "..", "data", "model_data", "fair_ghg_species_data.csv"), skiplines_begin=10))

    ## EMF CO2 emissions data
    
    # PAGE INPUT DATA (from Charles)
    EMF_CO2_data = readxlsx(joinpath(@__DIR__, "..", "data", "model_data", "FAIR-NCEE", "EMF22 ref data - with extrapolations for PAGE09, 9-1-11.xlsx"))
    # EMF_CO2_data = readxlsx(joinpath(@__DIR__, "data", "model_data", "FAIR-NCEE", "EMF22 ref data - with extrapolations for PAGE09, 9-1-11.xlsx")) # personal path
    EMF_CO2 = EMF_CO2_data["GLOBAL extrapolations"]["C11:AZ41"] # 2000 to 2300, decadal time step
    EMF_CO2[ismissing.(EMF_CO2)] .= 0 # replace missing values with 0 

    gtco2_to_gtc = 12/44 # FAIR CO2 emissions are in GtC

    annual_industrial_CO2_emissions = Dict()
    annual_industrial_CO2_emissions["USG1"] = EMF_CO2[:,6] * gtco2_to_gtc
    annual_industrial_CO2_emissions["USG2"] = EMF_CO2[:,16] * gtco2_to_gtc 
    annual_industrial_CO2_emissions["USG3"] = EMF_CO2[:,26] * gtco2_to_gtc
    annual_industrial_CO2_emissions["USG4"] = EMF_CO2[:,36] * gtco2_to_gtc
    annual_industrial_CO2_emissions["USG5"] = EMF_CO2[:,46] * gtco2_to_gtc

    annual_land_CO2_emissions = Dict()
    annual_land_CO2_emissions["USG1"] = EMF_CO2[:,10] * gtco2_to_gtc
    annual_land_CO2_emissions["USG2"] = EMF_CO2[:,20] * gtco2_to_gtc
    annual_land_CO2_emissions["USG3"] = EMF_CO2[:,30] * gtco2_to_gtc
    annual_land_CO2_emissions["USG4"] = EMF_CO2[:,40] * gtco2_to_gtc
    annual_land_CO2_emissions["USG5"] = EMF_CO2[:,50] * gtco2_to_gtc

    # annual_CO2_emissions = Dict()
    # annual_CO2_emissions["USG1"] = EMF_CO2[:,6] + EMF_CO2[:,10] # sum of land and industrial CO2, GtCO2/year
    # annual_CO2_emissions["USG2"] = EMF_CO2[:,16] + EMF_CO2[:,20]
    # annual_CO2_emissions["USG3"] = EMF_CO2[:,26] + EMF_CO2[:,30]
    # annual_CO2_emissions["USG4"] = EMF_CO2[:,36] + EMF_CO2[:,40]
    # annual_CO2_emissions["USG5"] = EMF_CO2[:,46] + EMF_CO2[:,50]
    
    old_years = collect(2000:10:2300) # EMF emissions years
    annual_years = collect(2000:1:2300)

    for scen in ["USG1", "USG2", "USG3", "USG4", "USG5"]
        annual_industrial_CO2_emissions[scen] = MimiIWG._interpolate(annual_industrial_CO2_emissions[scen], old_years, annual_years)
        annual_land_CO2_emissions[scen] = MimiIWG._interpolate(annual_land_CO2_emissions[scen], old_years, annual_years)
    end

    # merge with RCP CO2 emissions for pre-2000: RCP8.5 with USG1-4 and RCP4.5 with USG5

    # RCP8.5 scenario emissions.
    rcp85_emissions_raw = DataFrame(load(joinpath(@__DIR__, "..", "data", "model_data", "RCP85_EMISSIONS.csv"), skiplines_begin=36))
    # rcp85_emissions_raw = DataFrame(load(joinpath(@__DIR__, "data", "model_data", "RCP85_EMISSIONS.csv"), skiplines_begin=36)) # personal path

    # rcp85_emissions_raw[rcp85_emissions_raw[:,1].<2000,:] # pre-2000
    # rcp85_emissions_raw[rcp85_emissions_raw[:,1].<2000,2] # pre-2000 fossil (industrial) CO2
    # rcp85_emissions_raw[rcp85_emissions_raw[:,1].<2000,3] # pre-2000 other (land) CO2

    # RCP4.5 scenario emissions.
    rcp45_emissions_raw = DataFrame(load(joinpath(@__DIR__, "..", "data", "model_data", "RCP45_EMISSIONS.csv"), skiplines_begin=36))
    # rcp45_emissions_raw = DataFrame(load(joinpath(@__DIR__, "data", "model_data", "RCP45_EMISSIONS.csv"), skiplines_begin=36)) # personal path

    # rcp45_emissions_raw[rcp45_emissions_raw[:,1].<2000,:] # pre-2000
    # rcp45_emissions_raw[rcp45_emissions_raw[:,1].<2000,2] # pre-2000 fossil (industrial) CO2
    # rcp45_emissions_raw[rcp45_emissions_raw[:,1].<2000,3] # pre-2000 other (land) CO2

    # merge industrial CO2
    annual_industrial_CO2_emissions["USG1"] = append!(rcp85_emissions_raw.FossilCO2[rcp85_emissions_raw[:,1].<2000], annual_industrial_CO2_emissions["USG1"])
    annual_industrial_CO2_emissions["USG2"] = append!(rcp85_emissions_raw.FossilCO2[rcp85_emissions_raw[:,1].<2000], annual_industrial_CO2_emissions["USG2"])
    annual_industrial_CO2_emissions["USG3"] = append!(rcp85_emissions_raw.FossilCO2[rcp85_emissions_raw[:,1].<2000], annual_industrial_CO2_emissions["USG3"])
    annual_industrial_CO2_emissions["USG4"] = append!(rcp85_emissions_raw.FossilCO2[rcp85_emissions_raw[:,1].<2000], annual_industrial_CO2_emissions["USG4"])
    annual_industrial_CO2_emissions["USG5"] = append!(rcp45_emissions_raw.FossilCO2[rcp45_emissions_raw[:,1].<2000], annual_industrial_CO2_emissions["USG5"]) # merge USG5 with RCP4.5

    # annual_industrial_CO2_emissions

    # merge land CO2 (column 3 is "other CO2")
    annual_land_CO2_emissions["USG1"] = append!(rcp85_emissions_raw[rcp85_emissions_raw[:,1].<2000,3], annual_land_CO2_emissions["USG1"])
    annual_land_CO2_emissions["USG2"] = append!(rcp85_emissions_raw[rcp85_emissions_raw[:,1].<2000,3], annual_land_CO2_emissions["USG2"])
    annual_land_CO2_emissions["USG3"] = append!(rcp85_emissions_raw[rcp85_emissions_raw[:,1].<2000,3], annual_land_CO2_emissions["USG3"])
    annual_land_CO2_emissions["USG4"] = append!(rcp85_emissions_raw[rcp85_emissions_raw[:,1].<2000,3], annual_land_CO2_emissions["USG4"])
    annual_land_CO2_emissions["USG5"] = append!(rcp45_emissions_raw[rcp45_emissions_raw[:,1].<2000,3], annual_land_CO2_emissions["USG5"]) # merge USG5 with RCP4.5

    # annual_land_CO2_emissions

    ## EMF N2O and CH4 emissions data
    EMF_N2O_CH4 = readxlsx(joinpath(@__DIR__, "..", "data", "model_data", "FAIR-NCEE", "CH4N20emissions_annualversion.xlsx"))
    # EMF_N2O_CH4 = readxlsx(joinpath(@__DIR__, "data", "model_data", "FAIR-NCEE", "CH4N20emissions_annualversion.xlsx")) # personal path
    
    # N2O
    EMF_N2O = EMF_N2O_CH4["N20annual"]["B2:F297"] # 2005 to 2300, annual
    annual_N2O_emissions = Dict()
    annual_N2O_emissions["USG1"] = EMF_N2O[:,1]
    annual_N2O_emissions["USG2"] = EMF_N2O[:,2]
    annual_N2O_emissions["USG3"] = EMF_N2O[:,3]
    annual_N2O_emissions["USG4"] = EMF_N2O[:,4]
    annual_N2O_emissions["USG5"] = EMF_N2O[:,5]
    
    # CH4
    EMF_CH4 = EMF_N2O_CH4["CH4annual"]["B2:F297"] # 2005 to 2300, annual
    annual_CH4_emissions = Dict()
    annual_CH4_emissions["USG1"] = EMF_CH4[:,1]
    annual_CH4_emissions["USG2"] = EMF_CH4[:,2]
    annual_CH4_emissions["USG3"] = EMF_CH4[:,3]
    annual_CH4_emissions["USG4"] = EMF_CH4[:,4]
    annual_CH4_emissions["USG5"] = EMF_CH4[:,5]

    # merge with RCP N2O emissions for pre-2005

    rcp85_emissions_raw[rcp85_emissions_raw[:,1].<2005,5] # RCP8.5 pre-2005 N2O
    rcp45_emissions_raw[rcp45_emissions_raw[:,1].<2005,5] # RCP4.5 pre-2005 N2O

    annual_N2O_emissions["USG1"] = append!(rcp85_emissions_raw[rcp85_emissions_raw[:,1].<2005,5], annual_N2O_emissions["USG1"])
    annual_N2O_emissions["USG2"] = append!(rcp85_emissions_raw[rcp85_emissions_raw[:,1].<2005,5], annual_N2O_emissions["USG2"])
    annual_N2O_emissions["USG3"] = append!(rcp85_emissions_raw[rcp85_emissions_raw[:,1].<2005,5], annual_N2O_emissions["USG3"])
    annual_N2O_emissions["USG4"] = append!(rcp85_emissions_raw[rcp85_emissions_raw[:,1].<2005,5], annual_N2O_emissions["USG4"])
    annual_N2O_emissions["USG5"] = append!(rcp45_emissions_raw[rcp45_emissions_raw[:,1].<2005,5], annual_N2O_emissions["USG5"]) # merge USG5 with RCP4.5

    # merge with RCP CH4 emissions for pre-2005

    rcp85_emissions_raw[rcp85_emissions_raw[:,1].<2005,4] # RCP8.5 pre-2005 CH4
    rcp45_emissions_raw[rcp45_emissions_raw[:,1].<2005,4] # RCP4.5 pre-2005 CH4

    annual_CH4_emissions["USG1"] = append!(rcp85_emissions_raw[rcp85_emissions_raw[:,1].<2005,4], annual_CH4_emissions["USG1"])
    annual_CH4_emissions["USG2"] = append!(rcp85_emissions_raw[rcp85_emissions_raw[:,1].<2005,4], annual_CH4_emissions["USG2"])
    annual_CH4_emissions["USG3"] = append!(rcp85_emissions_raw[rcp85_emissions_raw[:,1].<2005,4], annual_CH4_emissions["USG3"])
    annual_CH4_emissions["USG4"] = append!(rcp85_emissions_raw[rcp85_emissions_raw[:,1].<2005,4], annual_CH4_emissions["USG4"])
    annual_CH4_emissions["USG5"] = append!(rcp45_emissions_raw[rcp45_emissions_raw[:,1].<2005,5], annual_CH4_emissions["USG5"]) # merge USG5 with RCP4.5

    # Calculate new index for years 1765:2300
    new_start_index, new_end_index = findall((in)([start_year, end_year]), collect(1765:2300))

    #---------------------------------------
    # Emissions
    #---------------------------------------
    emissions = DataFrame()

    # emissions.FossilCO2   = rcp_emissions_raw[start_index:end_index, :FossilCO2]
    # emissions.OtherCO2    = rcp_emissions_raw[start_index:end_index, :OtherCO2]
    emissions.FossilCO2     = annual_industrial_CO2_emissions[usg_scenario][new_start_index:new_end_index] 
    emissions.OtherCO2      = annual_land_CO2_emissions[usg_scenario][new_start_index:new_end_index] 
    # emissions.CH4         = rcp_emissions_raw[start_index:end_index, :CH4]
    emissions.CH4         = annual_CH4_emissions[usg_scenario][new_start_index:new_end_index]
    emissions.NaturalCH4  = natural_emissions_raw[start_index:end_index, :ch4]
    # emissions.N2O         = rcp_emissions_raw[start_index:end_index, :N2O]
    emissions.N2O         = annual_N2O_emissions[usg_scenario][new_start_index:new_end_index]
    emissions.NaturalN2O  = natural_emissions_raw[start_index:end_index, :n2o]
    emissions.NMVOC       = rcp_emissions_raw[start_index:end_index, :NMVOC]
    emissions.CO          = rcp_emissions_raw[start_index:end_index, :CO]
    emissions.NOx         = rcp_emissions_raw[start_index:end_index, :NOx]
    emissions.SOx         = rcp_emissions_raw[start_index:end_index, :SOx]
    emissions.BC          = rcp_emissions_raw[start_index:end_index, :BC]
    emissions.OC          = rcp_emissions_raw[start_index:end_index, :OC]
    emissions.NH3         = rcp_emissions_raw[start_index:end_index, :NH3]

    # Other greenhouse gases
    for i in other_ghg_names
        emissions[!,Symbol(i)] = rcp_emissions_raw[start_index:end_index, Symbol(i)]
    end

    #---------------------------------------
    # Gas Fractions
    #---------------------------------------
    gas_fractions = DataFrame()

    gas_fractions.nox_aviation = aviation_fraction_raw[start_index:end_index]
    gas_fractions.ch4_fossil = ch4_fossil_frac_raw[start_index:end_index]

    #---------------------------------------
    # Emission to Concentration Conversions
    #---------------------------------------
    emiss_conversions = DataFrame()

    # Set names for concentration conversions.
    emiss_conversions.gases = vcat("CO2", "CH4", "N2O", other_ghg_names)

    # Mass of atmosphere (kg).
    mass_atmos = 5.1352e18

    # Molecular weights of air, CO₂, N₂, CH₄, and other greenhouse gases.
    mol_wt_air    = gas_data[gas_data.gas .== "AIR", :mol_weight][1]
    mol_wt_carbon = gas_data[gas_data.gas .== "C", :mol_weight][1]
    mol_wt_n2     = gas_data[gas_data.gas .== "N2", :mol_weight][1]
    mol_wt_ch4    = gas_data[gas_data.gas .== "CH4", :mol_weight][1]
    mol_wt_others = gas_data[findall((in)(other_ghg_names), gas_data.gas), :mol_weight]

    # Calculate CO₂ conversion from GtC to ppm.
    emiss2conc_carbon = (mass_atmos / 1.0e18) * (mol_wt_carbon / mol_wt_air)

    # Emission to concentration for CH₄.
    emiss2conc_ch4 = (mass_atmos / 1.0e18) * (mol_wt_ch4 / mol_wt_air)

    # Emission to concentration for N₂O.
    # Note: Use N₂ for N₂O conversion, from FAIR: "Funny units for nitrogen emissions - N₂O is expressed in N₂ equivalent."
    emiss2conc_n2o= (mass_atmos / 1.0e18) * (mol_wt_n2 / mol_wt_air)

    # Conversion factors for other Kyoto Protocol or ozone-depleting gases.
    emiss2conc_others = (mass_atmos / 1.0e18) .* (mol_wt_others ./ mol_wt_air)

    # Combine values conversion values into a single data frame.
    emiss_conversions.emiss2conc = vcat(emiss2conc_carbon, emiss2conc_ch4, emiss2conc_n2o, emiss2conc_others)

    return emissions, cmip6_volcano_forcing, cmip6_solar_forcing, gas_data, gas_fractions, emiss_conversions
end