###############################
# AUTHOR: Adam Smith-Orlik
# DATE: March 2023
# DESCRIPTION: Defines various parameters used to compute the evolution 
# CONTACT: asorlik@yorku.ca
###############################

const solar_mass_to_pc = 4.8e-14

# system parameters

# radius of the innermost stable circular orbit
function r_isco(; kwargs...)    
    """
    Radius of the Innermost Stable Circular Orbit (ISCO) of a Schwarzschild black hole with mass m
    """
    
    Mbh = kwargs[:Mbh]
    Mbh *= solar_mass_to_pc
    return 6.0*float(Mbh)
end

function a_0(; kwargs...)
    """
    a_initial * r_isco
    """
    return kwargs[:a_initial]*r_isco(; kwargs...)
    
end

# initial arrays

# reduced mass of the bianry system 
function mass_A(; kwargs...)
    Mbh = kwargs[:Mbh] * solar_mass_to_pc
    Mc = kwargs[:Mc] * solar_mass_to_pc
    
    m_total = (Mbh .+ Mc) 
    
    mu = (Mbh .* Mc) ./ (Mbh .+ Mc)
    
    
    return [Mbh, Mc, m_total, mu]
end

# start and stop times for each block
function time_params(mass; kwargs...)  
    
    t_start = kwargs[:t_start]
    a0 = a_0(; kwargs...)
    F0 = sqrt(mass[3]/a0^3)/2/pi
    norbits = kwargs[:norbits]
    t_end = norbits / F0
    
    return [t_start, t_end]

end

function initial_conditions(mass; kwargs...)
    
    """
    Calculate the initial conditions for a Keplerian orbit with parameters a, e
    """
    mt = mass[3]
    a0 = a_0(; kwargs...)
    e0 = kwargs[:e0]
    phi0 = kwargs[:phi0]

    r0 = a0 * (1. - e0^2) / (1. + e0 * cos(phi0))
    dphi0 = sqrt(mt * a0 * (1. - e0^2)) / r0^2
    dr0 = a0* (1. - e0^2) / (1 + e0*cos(phi0))^2 * e0 *sin(phi0)*dphi0

    out = [r0, phi0, dr0, dphi0]

end