###############################
# AUTHOR: Adam Smith-Orlik
# DATE: March 2023
# DESCRIPTION: Defines the cosmological parameters used in halo.jl
# CONTACT: asorlik@yorku.ca
###############################

# export constant and functions 
# export hubble, Omega_m_0, Omega_L_0, hubbleParameter, criticalDensity, Omega_m

# define constants
const hubble = 2.3e-10 # [1/pc]
const Omega_m_0 = 0.3111
const Omega_L_0 = 0.6889

# functions
function hubbleParameter(z::Integer)
    return hubble * sqrt(Omega_L_0) / tanh(asinh(sqrt(Omega_L_0/Omega_m_0/(1+z)^3))) 
end

function criticalDensity(z)
    H = hubbleParameter(z)
    return 3.0 * H^2 / 8.0 / Float64(pi)
end

function Omega_m(z)
    rho_c = criticalDensity(z)
    return 1.0 - 3*hubble^2 * Omega_L_0 / (8.0 * Float64(pi) * rho_c)
end

