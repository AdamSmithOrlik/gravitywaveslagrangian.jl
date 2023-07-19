###############################
# AUTHOR: Adam Smith-Orlik
# DATE: March 2023
# DESCRIPTION: Module that runs the binary evolution using the input parameters in "run_dictionary"
# CONTACT: asorlik@yorku.ca
###############################

module gravitywaveslagrangian

include("partition.jl")

export run

###############################################################
#               DICTIONARIES DEFINING RUNS                    #  
#   MUST PASS THE DICTIONARY NAME TO THE run() FUNCTION       #  
#    AND THE DICTIONARY MUST ALREADY BE DEFINED HERE...       #       
###############################################################

# INSTRUCTIONS:
# Copy defualt_dictionary and edit the parameters as needed. Change the dict name and add it to 
# the run function check to avoid errors. 
# Run the code using "gravitywaveslagrangian.run(<your_dict_name>)"

# test dictionary to make sure code is working... 
defualt_dictionary = Dict(
    # Name of hdf5 produced 
    :Run_name=>"10_3_testing",  

    # Beginning and end radius of inpiral 
    :a_initial=>10.0,  # initial radius in units of r_isco
    :a_final=>3.0,      # final radius in units of r_isco (where integration stops)
        
    # System Parameters 
    :Mbh=>1e5,  # Black hole [Msol] 
    :Mc=>1e3,   # Compact object [Msol]
    :D=>1e5,    # Lumiosity distance 
    :e0=>0.0,   # Ellipticity 
    :phi0=>0.0, # Inital angle 
    :lisa_bandwidth=> [1e-5, 1e0], # frequency range for lisa 

    # time parameters 
    :time_resolution=>1e-8, # resolution of the solution arrays and the strain, recommended 1e-9 for optimal GW strain resolution 
    :t_start=>0.0,          # [years]
    :t_max=>1e6,            # [years] upper limit of the time allowed for the evolution to run without terminating 
    :dt=>0.01,              # time steps between t_start and t_max, must be < t_max. dt > 0.1 uses 30+GBs of RAM. 
    # NOTE: dt ~ 0.05 work for longer inspirals. If you want small inspirals (e.g. 10-8 risco) you may need to make this value smaller.

    # orbital data 
    :save_orbital_data=>false,  # save the orbital arrays, if false only the strain is saved
    :orbital_resolution=>1e-5, # used to sample the full resolution arrays by steps round(orbit_res/time_res). Must be greater than time resolution
    
    # Forces included in insprial 
    :gwEmission=>true,
    :dynamicalFriction=>true,
    :postNewtonian=>true,
    
    # DM halo model flags 
    :dmModel=>"Relativistic", # "Newtonian" or "Relativistic"
    :dfModel=>"Relativistic", # "Newtonian" or "Relativistic"
    :plot_halo_model=>true,
    # Inital DM halo density and redshift from Eda et al. (arxiv:1408.3534)
    # :rho0=> 3.8e-22 .* 7.071e8, # [g/cm^3 to 1/pc^2]
    # low rho0->high density
    :rho0=> 8.91e-25 .* 7.071e8, # [g/cm^3 to 1/pc^2] -> Equivalent to 0.5 GeV/cm^3

    # redshift taken from arxiv:1408.3534
    :z_f => 20,

    # coulomb log 
    :coulombLog=>3.0
)

# test dictionary to make sure code is working... 
eda_dictionary = Dict(
    # Name of hdf5 produced 
    :Run_name=>"100_3_Eda_noDM_test",  

    # Beginning and end radius of inpiral 
    :a_initial=>100.0,  # initial radius in units of r_isco
    :a_final=>3.0,      # final radius in units of r_isco (where integration stops)
        
    # System Parameters 
    :Mbh=>1e5,  # Black hole [Msol] 
    :Mc=>1e3,   # Compact object [Msol]
    :D=>1e5,    # Lumiosity distance 
    :e0=>0.0,   # Ellipticity 
    :phi0=>0.0, # Inital angle 
    :lisa_bandwidth=> [1e-5, 1e0], # frequency range for lisa 

    # time parameters 
    :time_resolution=>1e-8, # resolution of the solution arrays and the strain, recommended 1e-8 for optimal GW strain resolution 
    :t_start=>0.0,          # [years]
    :t_max=>1e6,            # [years] upper limit of the time allowed for the evolution to run without terminating 
    :dt=>0.05,              # time steps between t_start and t_max, must be < t_max. dt > 0.1 uses 30+GBs of RAM

    # orbital data 
    :save_orbital_data=>true,  # save the orbital arrays, if false only the strain is saved
    :orbital_resolution=>1e-5, # used to sample the full resolution arrays by steps round(orbit_res/time_res). Must be greater than time resolution
    
    # Forces included in insprial 
    :gwEmission=>true,
    :dynamicalFriction=>false, # if false no DM is included
    :postNewtonian=>true,
    
    # DM halo model flags 
    :dmModel=>"Relativistic", # "Newtonian" or "Relativistic"
    :dfModel=>"Relativistic", # "Newtonian" or "Relativistic"
    :plot_halo_model=>true,
    # Inital DM halo density and redshift from Eda et al. (arxiv:1408.3534)
    :rho0=> 3.8e-22 .* 7.071e8, # [g/cm^3 to 1/pc^2]
    # redshift taken from arxiv:1408.3534
    :z_f => 20,

    # coulomb log 
    :coulombLog=>3.0
)

run0_dictionary = Dict(
    # Name of hdf5 produced 
    :Run_name=>"high_res_100_3_risco_Mbh_1e5_Mc_1e3_rho_0_5GeV_dm_full_rel",  
    # :Run_name=>"10_6_test",  

    # Beginning and end radius of inpiral 
    :a_initial=>100.0,  # initial radius in units of r_isco
    :a_final=>3.0,      # final radius in units of r_isco (where integration stops)
        
    # System Parameters 
    :Mbh=>1e5,  # Black hole [Msol] 
    :Mc=>1e3,   # Compact object [Msol]
    :D=>1e5,    # Lumiosity distance 
    :e0=>0.0,   # Ellipticity 
    :phi0=>0.0, # Inital angle 
    :lisa_bandwidth=> [1e-5, 1e0], # frequency range for lisa 

    # time parameters 
    :time_resolution=>1e-11, # resolution of the solution arrays and the strain, recommended 1e-8 for optimal GW strain resolution 
    :t_start=>0.0,          # [years]
    :t_max=>1e6,            # [years] upper limit of the time allowed for the evolution to run without terminating 
    :dt=>0.01,              # time steps between t_start and t_max, must be < t_max. dt > 0.1 uses 30+GBs of RAM

    # orbital data 
    :save_orbital_data=>false,  # save the orbital arrays, if false only the strain is saved
    :orbital_resolution=>1e-5, # used to sample the full resolution arrays by steps round(orbit_res/time_res). Must be greater than time resolution
    
    # Forces included in insprial 
    :gwEmission=>true,
    :dynamicalFriction=>true,
    :postNewtonian=>true,
    
    # DM halo model flags 
    :dmModel=>"Relativistic", # "Newtonian" or "Relativistic"
    :dfModel=>"Relativistic", # "Newtonian" or "Relativistic"
    :plot_halo_model=>true,
    # Inital DM halo density and redshift from Eda et al. (arxiv:1408.3534)
    # :rho0=> 3.8e-22 .* 7.071e8, # [g/cm^3 to 1/pc^2]
    :rho0=> 8.91e-25 .* 7.071e8, # [g/cm^3 to 1/pc^2] -> Equivalent to 0.5 GeV/cm^3
  
    # redshift taken from arxiv:1408.3534
    :z_f => 20,

    # coulomb log 
    :coulombLog=>3.0
)

# 0.1 GeV from 100-3 risco Mc = 1e3
run1_dictionary = Dict(
    # Name of hdf5 produced 
    :Run_name=>"10_3_risco_Mbh_1e5_Mc_1e3_rho_0_5GeV_full_newtonian",  

    # Beginning and end radius of inpiral 
    :a_initial=>10.0,  # initial radius in units of r_isco
    :a_final=>3.0,      # final radius in units of r_isco (where integration stops)
        
    # System Parameters 
    :Mbh=>1e5,  # Black hole [Msol] 
    :Mc=>1e3,   # Compact object [Msol]
    :D=>1e5,    # Lumiosity distance 
    :e0=>0.0,   # Ellipticity 
    :phi0=>0.0, # Inital angle 
    :lisa_bandwidth=> [1e-5, 1e0], # frequency range for lisa 

    # time parameters 
    :time_resolution=>1e-8, # resolution of the solution arrays and the strain, recommended 1e-8 for optimal GW strain resolution 
    :t_start=>0.0,          # [years]
    :t_max=>1e6,            # [years] upper limit of the time allowed for the evolution to run without terminating 
    :dt=>0.01,              # time steps between t_start and t_max, must be < t_max. dt > 0.1 uses 30+GBs of RAM

    # orbital data 
    :save_orbital_data=>true,  # save the orbital arrays, if false only the strain is saved
    :orbital_resolution=>1e-5, # used to sample the full resolution arrays by steps round(orbit_res/time_res). Must be greater than time resolution
    
    # Forces included in insprial 
    :gwEmission=>true,
    :dynamicalFriction=>true,
    :postNewtonian=>true,
    
    # DM halo model flags 
    :dmModel=>"Newtonian", # "Newtonian" or "Relativistic"
    :dfModel=>"Newtonian", # "Newtonian" or "Relativistic"
    :plot_halo_model=>true,
    # Inital DM halo density and redshift from Eda et al. (arxiv:1408.3534)
    # :rho0=> 3.8e-22 .* 7.071e8, # [g/cm^3 to 1/pc^2]
    
    # 0.1 GeV
    # :rho0=> 1.78e-25 .* 7.071e8, # [g/cm^3 to 1/pc^2] -> Equivalent to 0.1 GeV/cm^3
    :rho0=> 8.91e-25 .* 7.071e8, # [g/cm^3 to 1/pc^2] -> Equivalent to 0.5 GeV/cm^3
    # redshift taken from arxiv:1408.3534
    :z_f => 20,

    # coulomb log 
    :coulombLog=>3.0
)

# rho 0.5 GeV from 100-3 risco Mc = 1e2
run2_dictionary = Dict(
    :Run_name=>"10_5_risco_Mbh_1e5_Mc_1e3_rho_0_5GeV_Newtonian",  

    # Beginning and end radius of inpiral 
    :a_initial=>10.0,  # initial radius in units of r_isco
    :a_final=>8.0,      # final radius in units of r_isco (where integration stops)
        
    # System Parameters 
    :Mbh=>1e5,  # Black hole [Msol] 
    :Mc=>1e3,   # Compact object [Msol]
    :D=>1e5,    # Lumiosity distance 
    :e0=>0.0,   # Ellipticity 
    :phi0=>0.0, # Inital angle 
    :lisa_bandwidth=> [1e-5, 1e0], # frequency range for lisa 

    # time parameters 
    :time_resolution=>1e-9, # resolution of the solution arrays and the strain, recommended 1e-8 for optimal GW strain resolution 
    :t_start=>0.0,          # [years]
    :t_max=>1e6,            # [years] upper limit of the time allowed for the evolution to run without terminating 
    :dt=>0.05,              # time steps between t_start and t_max, must be < t_max. dt > 0.1 uses 30+GBs of RAM

    # orbital data 
    :save_orbital_data=>true,  # save the orbital arrays, if false only the strain is saved
    :orbital_resolution=>1e-5, # used to sample the full resolution arrays by steps round(orbit_res/time_res). Must be greater than time resolution
    
    # Forces included in insprial 
    :gwEmission=>true,
    :dynamicalFriction=>true,
    :postNewtonian=>true,
    
    # DM halo model flags 
    :dmModel=>"Newtonian", # "Newtonian" or "Relativistic"
    :dfModel=>"Newtonian", # "Newtonian" or "Relativistic"
    :plot_halo_model=>true,
    
    # Inital DM halo density and redshift from Eda et al. (arxiv:1408.3534)
    :rho0=> 3.8e-22 .* 7.071e8, # [g/cm^3 to 1/pc^2]
    
    # :rho0=> 1.78e-25 .* 7.071e8, # [g/cm^3 to 1/pc^2] -> Equivalent to 0.1 GeV/cm^3
    # :rho0=> 8.91e-25 .* 7.071e8, # [g/cm^3 to 1/pc^2] -> Equivalent to 0.5 GeV/cm^3
    # redshift taken from arxiv:1408.3534
    :z_f => 20,

    # coulomb log 
    :coulombLog=>3.0
)

###############################################################
#                     RUN FUNCTION                            #        
###############################################################
function run(dictname::Symbol)
    if !(dictname in (:defualt_dictionary, :eda_dictionary, :run0_dictionary, :run1_dictionary, :run2_dictionary, :run3_dictionary, :run4_dictionary, :run5_dictionary))
        error("Dictionary $dictname is not defined in the gravitywaveslagrangian module. Please define it and try again.")
    else
        dict = eval(dictname)
        run_partition(; dict...)
    end
end


end # module