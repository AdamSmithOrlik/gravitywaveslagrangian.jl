###############################
# AUTHOR: Adam Smith-Orlik
# DATE: March 1st
# DESCRIPTION: Defines the halo model/s for the dark matter spike 
# in the binary orbit black hole system 
# CONTACT: asorlik@yorku.ca
###############################

using Roots
using QuadGK
using Plots
using Measures

include("cosmology.jl")

const solar_mass_to_pc = 4.8e-14

# Parameters for defining different spiked profiles
spike_params = Dict(
    ## Parameters for the spiked DM profile fit in arxiv:2204.12508

    # NFW parameters
    :alpha=>0.331,
    :beta=>-1.66,
    :gamma=>0.32,
    :delta=>-0.000282,
    
    # NFW relativistic parameters
    :eta=>1,
    :A=>6.42,
    :w=>1.82,
    :q=>1.91,
    
    # spike param
    :alpha_ini=>1.0
)

function mass_concentration_relation(M_vir, z_f)
    """
    Takes
        M_vir: Virial mass of the dark matter halo in solar masses
        z_f: redshift, 0 being the present 
    Returns
        rho_s: scale density in units 1/pc^2
        r_s: scale radius in units pc 
    Equation for mass concentarion relation taked from Eda et al. (arxiv:1408.3534)
    """

    M_vir = M_vir * solar_mass_to_pc
    A_200 = 5.71; B_200 = -0.084; C_200 = -0.47; M_piv = 1.0/0.7 *1e14* solar_mass_to_pc 
    
    c_200 = A_200 * (M_vir/M_piv)^B_200 * (1.0 + z_f)^C_200
    rho_crit_0 = 3.0* hubble^2 / 8.0 / Float64(pi)
    Omegam = Omega_m(z_f)
    rho_m = Omega_m_0*rho_crit_0 * (1.0 +z_f)^3.0
    Delta_vir = 18.0 * Float64(pi)^2 *(1.0 + 0.4093 * (1/Omegam - 1.0)^0.9052)
    r_vir = (3.0 *M_vir / (4.0 * Float64(pi) * Delta_vir * rho_m))^(1.0/3.0)
    r_s = r_vir/c_200
    f = log(1.0 + c_200) - c_200/(1.0 + c_200)
    rho_s = 1.0/3.0 / f * Delta_vir * rho_m * c_200^3.0
    
    return rho_s, r_s
    
end

function rs_from_rhos(rhos, z_f)
    
    function f(expt)
        Mvir = 10.0^expt
        LHS = rhos / solar_mass_to_pc # scales results
        RHS, na = mass_concentration_relation(Mvir, z_f)
        RHS = RHS / solar_mass_to_pc
        return LHS - RHS
    end
    
    # find the exponent on M_vir by numerically finding the root of the function
    # mass_concentration_relation
    # res = find_zero(f, (1,14))
    res = find_zero(f, 4)

    Mvir = 10.0^res
    rho_s, r_s =  mass_concentration_relation(Mvir, z_f)
    
    return r_s
end

function get_Rs(; kwargs...) 
    """
    Takes: None
    Args: None
    Kwargs: 
         Mass of black hole 'Mbh' in units [Msol]
    
    Returns: Schwartzchild radius for black hole of mass Mbh in units [pc]
    """
    # conversion factors
    
    Mbh = kwargs[:Mbh]
    
    G = 4.3e-3  # units [pc/Msol * (km/s)^2]
    c = 3e5 # units [km/s]
    Rs = 2 .* G .* Mbh ./ c^2 
    
    return Rs
end

function rho_nfw(r; kwargs...)
    """
    Takes: 
        Radius 'r' in [pc]
    Args: 
        rhos in units [1/pc^2]
        rs in [pc] 
        (arxiv:1408.3534)
    Kwargs: None
    
    Returns: Dark matter NFW density at r [pc] in units [1/pc^2]
    c=G=h=1 units
    """
    rhos, z = kwargs[:rho0], kwargs[:z_f] 
    rs = rs_from_rhos(rhos, z)
    
    x = r./rs
    rho = rhos./ (x .* (1 .+ x).^2)

    return rho
end

function rho_generalized(r, alpha=1, beta=3, gamma=1; kwargs...)
    """
    Takes: 
        Radius 'r' in [pc]
        alpha, beta, gamme: determines behavior at small r, large r, and transion regions.
        (1,3,1) corresponds to NFW and (1,4,1) to Hernquist.
    Args: 
        rhos in units [1/pc^2]
        rs in [pc] 
        (arxiv:2204.12508 eqn 1)
    Kwargs: None
    
    Returns: Dark matter NFW density at r [pc] in units [1/pc^2]
    c=G=h=1 units
    """
    rhos, z = kwargs[:rho0], kwargs[:z_f] 
    rs = rs_from_rhos(rhos, z)
    
    x = r./rs
    exp = (gamma-beta)/alpha
    rho = rhos.*x^(-gamma)*(1+x^alpha)^exp
    
    return rho
end

function rho_spike(r; rsp=0.54, kwargs...)
    """
    Takes: 
        Radius 'r' in [pc]
    Args: 
        Radius of spike 'rsp' calculated numerically using rsp = 0.2*rh from arxiv:1408.3534
    see function get_rsp()
    Kwargs:
        Mass of black hole 'Mbh' in units [Msol]
        alpha initial 'alpha_ini' dimensionless in range (0,2)
    
    Returns: Dark matter NFW density at r [pc] in units [1/pc^2]
    c=G=h=1 units
    """
    alpha_ini = spike_params[:alpha_ini]
    
    rhosp = rho_nfw(rsp; kwargs...)
    alpha = (9.0 .-(2.0 .*alpha_ini))/(4.0 .-alpha_ini)
    Rs = get_Rs(; kwargs...)
    g_alpha = (1 .- 4 .*(Rs ./ r)).^3

    rho = rhosp .* g_alpha .* (rsp ./ r).^alpha # [GeV/cm^3]

    return ifelse.(r .< rsp, rho, rho_nfw(r; kwargs...))
end

function get_rsp(; kwargs...)
    """
    Takes: 
        None
    Args: 
        None
    Kwargs:
        Mass of black hole 'Mbh' in units [Msol]
    
    Returns: r spike 'rsp' calculate via rsp ~ 0.2*rh where rh is defined in arxiv:1408.3534
    as: M(<rh) = int_0^rh 4pi rho_dm(r)r^2 = 2M_bh
    
    NB:
        rhosp is defined as rho_nfw(rsp) and sets the normalization for rho_spike. Important
    to get this right as it determines the magnitude of density at low radii. 
    c=G=h=1 units
    """
    
    Mbh = kwargs[:Mbh]
    
    function f(rh)
        rsp = 0.2*rh
        I(r) = 4 * pi .* rho_spike(r; rsp=rsp, kwargs...) .* r.^2 
        LHS, err = quadgk(I, 1e-6, rh) 
        RHS = 2 * Mbh .* solar_mass_to_pc 
        return LHS .- RHS
    end
    
    rh = find_zero(f, 1.0)
    return 0.2.*rh # rsp 
end

function rho_relativistic(r; kwargs...)
    """
    Takes: 
        Radius 'r' in [pc]
    Args: 
        Mass of black hole 'Mbh' in units [Msol]
    Kwargs: 
        rho0 and a are scale paramters form arxiv:1408.3534
        alpha, beta, gamma and delta are scale fit parameters for Eq 7. from arXiv:2204.12508v1
        A, w, q, eta are fit parameters from Table 1. from arXiv:2204.12508v1
        
    
    Returns: Dark matter effective density in [1/pc^2] as a function of r, valid in the range
    r << a, a >= 100 pc, or r << 100 pc
    c=G=h=1 units
    """
    # unit conversion 
    unity_pc2_to_Gev_cm3 = 7.9338e14
    pc_to_kpc = 1e-3
    
    Mbh = kwargs[:Mbh]
    # scale params
    rho0, z = kwargs[:rho0] , kwargs[:z_f] 
    # "a" is the scale radius "rs". Kept as "a" to match Eq (7) arXiv:2204.12508v1  
    a = rs_from_rhos(rho0, z)
    
    # convert to the units used in arXiv:2204.12508v1
    rho0 *= unity_pc2_to_Gev_cm3
    a *= pc_to_kpc
    
    # profile params
    alpha, beta, gamma, delta = spike_params[:alpha], spike_params[:beta], spike_params[:gamma], spike_params[:delta]
    # fit params 
    A, w, q, eta= spike_params[:A], spike_params[:w], spike_params[:q], spike_params[:eta]
   
    A *= 1e-43 / (solar_mass_to_pc)^2 #[ Msol^-2] to [pc^-2]
    x_tilde = r / (Mbh * solar_mass_to_pc)
    
    # rho_bar is the only dimensionful quantity ([rho_bar]=1/pc^2) converted so that the result is in untis 1/pc^2 
    rho_bar = A*(1.0 - 4.0 * (eta / x_tilde))^w *(4.17 * (1e11/x_tilde))^q
    
    rho = rho_bar*(10^delta)*((rho0/0.3)^alpha)*((Mbh/1e6)^beta)*((a/20.0)^gamma)

    return rho
end

# Plotting the halo model 
function plot_dm_spike(; kwargs...)
    dmModel = kwargs[:dmModel]
    runname = kwargs[:Run_name]
    aini = kwargs[:a_initial]
    afin = kwargs[:a_final]
    risco = r_isco(; kwargs...)
    
    inspiral_range = [aini, afin]
    x = collect(range(inspiral_range[1],inspiral_range[2],3))
    ym = zeros(length(x))
    yt = zeros(length(x))
    
    r_list = collect(10 .^ range(log10(risco), log10(1e8*risco), 10000)) # log spaced range 
    # generate nfw density as reference density
    rhonfw = [rho_nfw(r; kwargs...) for r in r_list]
    # generate the spike model density
    if dmModel == "Newtonian"
        rsp = get_rsp(; kwargs...)
        rho = [rho_spike(r; rsp=rsp, kwargs...) for r in r_list]
        model = "GS"
    elseif dmModel == "Relativistic"
        rho = [rho_relativistic(r; kwargs...) for r in r_list]
        model = "Relativistic"
    else
        throw(DomainError(dmModel, "must be either 'Newtonian' or 'Relativistic'."))
    end

    # remove all negative values 
    rho_length = length(rho)
    rho = filter(x -> x >= 0, rho)
    rho_length_filtered = length(rho)
    start = (rho_length-rho_length_filtered) + 1
    
   # plot the results 
    plot(r_list[start:end] ./risco, rho, xscale=:log10, yscale=:log10, left_margin=10mm, right_margin=10mm, xlabel="r/risco", ylabel="Density [1/pc^2]",labelfontsize=6,framestyle=:box, label="Spike Model")
    plot!(r_list ./risco, rhonfw, xscale=:log10, yscale=:log10, label="NFW profile")
    vline!(inspiral_range, linecolor=:grey, linestyle=:dash, label=nothing)
    ymin, ymax = ylims()
    ym .= ymin
    yt .= ymax
    plot!(x, ym, fillrange = yt, fillalpha = 0.1, c =:black, label="Inspiral Range")
    title!("$model Spike Model for $runname",  titlefont=font(8))

    data_directory = pwd() * "/data/"
     isdir(data_directory) ? nothing : mkdir(data_directory)
 
     filename = data_directory * runname * ".png"
    savefig(filename)
end