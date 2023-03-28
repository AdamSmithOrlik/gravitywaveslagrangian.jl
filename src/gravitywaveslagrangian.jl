###############################
# AUTHOR: Adam Smith-Orlik
# DATE: March 2023
# DESCRIPTION: Module that runs the binary evolution using the input parameters in "run_dictionary"
# CONTACT: asorlik@yorku.ca
###############################

module gravitywaveslagrangian

include("partition.jl")

export run

run_dictionary = Dict(
    # Name of hdf5 produced 
    :Run_name=>"test_strain",  

    # Beginning and end radius of inpiral 
    :a_initial=>8.0,  # initial radius in units of r_isco
    :a_final=>5.0,      # final radius in units of r_isco (where integration stops)
        
    # System Parameters 
    :Mbh=>1e5, # Black hole [Msol] 
    :Mc=>1e3, # Compact object [Msol]
    :D=>1e5, # Lumiosity distance 
    :e0=>0.0, # Ellipticity 
    :phi0=>0.0, # Inital angle 
    :lisa_bandwidth=> [1e-5, 1e0], # frequency range for lisa 

    # time parameters 
    :time_resolution=>1e-6, # resolution of the outoput solution arrays, recommended 1e-8 for GW strain resolution 
    :t_start=>0.0, # [years]
    :t_max=>1e6, # [years] upper limit of the time allowed for the evolution to run without terminating 
    :dt=>0.005, # time steps between t_start and t_max, must be < t_max

    # orbital data 
    :save_orbital_data=>true, # save the orbital arrays, if false only the strain is saved
    :orbital_resolution=>1e-6, # used to sample the full resolution arrays by steps round(orbit_res/time_res). Must be greater than time resolution
    
    # Forces included in insprial 
    :gwEmission=>true,
    :dynamicalFriction=>true,
    :postNewtonian=>true,
    
    # DM halo model flags 
    :dmModel=>"Relativistic", # "Newtonian" or "Relativistic"
    :dfModel=>"Relativistic", # "Newtonian" or "Relativistic"
    
    # Inital DM halo density and redshift from Eda et al. (arxiv:1408.3534)
    :rho0=> 3.8e-22 .* 7.071e8, # [g/cm^3 to 1/pc^2]
    # redshift taken from arxiv:1408.3534
    :z_f => 20,

    # coulomb log 
    :coulombLog=>3.0
)

function run()
    run_partition(; run_dictionary...)
end


end # module
