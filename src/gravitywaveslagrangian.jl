
module gravitywaveslagrangian

include("partition.jl")

export run

run_dictionary = Dict(
    # Name of hdf5 produced 
    :Run_name=>"julia_test",  

    # Beginning and end radius of inpiral 
    :a_initial=>8.0,  # initial radius in units of r_isco
    :a_final=>3.0,      # final radius in units of r_isco (where integration stops)
        
    # System Parameters 
    :Mbh=>1e5, # Black hole [Msol] 
    :Mc=>1e3, # Compact object [Msol]
    :D=>1e5, # Lumiosity distance 
    :e0=>0.0, # Ellipticity 
    :phi0=>0.0, # Inital angle 
    
    # Evolution total orbits and orbital resolution, i.e. time steps per evolution
    :norbits=>1e6,   # upper limit on number of orbits for total evolution 
    :t_orbital=>20,   # number of time steps per orbit, i.e. time resolution per orbit 
    :t_start=>0.0,
    
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
