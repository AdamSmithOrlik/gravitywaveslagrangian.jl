###############################
# AUTHOR: Adam Smith-Orlik
# DATE: March 2023
# DESCRIPTION: Defines the evolve function that numerically integrates to solve the orbital evolution of a binary system
# CONTACT: asorlik@yorku.ca
###############################

using DifferentialEquations

include("halo.jl")
include("runparameters.jl")

function Evolve(u0, initial, final, resolution; kwargs...)
    """
    Evolve the system of differential equations from t_start to t_end with nSteps with initial conditions u0
    """
    m = mass_A(; kwargs...)
    tspan = (initial, final)    
    
    # Toggle forces
    dynamicalFriction = kwargs[:dynamicalFriction]
    postNewtonian = kwargs[:postNewtonian]
    gwEmission = kwargs[:gwEmission]
    coulombLog = kwargs[:coulombLog]
    
    # Model flags
    dmModel = kwargs[:dmModel]
    dfModel = kwargs[:dfModel]
    
    #---------------------------Derivative Function---------------------------------------
    
    function orbit!(du, u, m, t)
        
        # unpack input arrays 
        r, phi, dr, dphi = u
        m1, m2, mt, mu = m 
        eta = mu/mt

        # define velocity
        v = sqrt(dr^2 + (r^2)*(dphi^2))
        
        # choose dm model 
        if dmModel == "Newtonian"
            rho = rho_spike(r; kwargs...)
        elseif dmModel == "Relativistic"
            rho = rho_relativistic(r; kwargs...)
        else
            throw(DomainError(dmModel, "must be either 'Newtonian' or 'Relativistic'."))
        end
        
        # choose df model 
        if dfModel == "Newtonian"
            zeta = 0
        elseif dfModel == "Relativistic"
            zeta = 1
        else
            throw(DomainError(dfModel, "must be either 'Newtonian' or 'Relativistic'."))
        end
        
        # relativistic correction term 
        gamma = sqrt(1 .- (zeta .*v^2))
        correction = gamma^2 .* (1 .+ (zeta .* v^2))^2
        
        # forces in terms of power defined in arxiv:XXXX.XXXX (our paper)
        P_df     = dynamicalFriction ? 4*pi*(m2^2 / mu) .* rho .* correction  .* (coulombLog/v^3) : 0.
        
        P_pn_r   = postNewtonian ? (mt/r^2)*( (4 + 2*eta)*(mt/r) - (3*eta +1)*(r^2)*(dphi^2) + (3 - 7*eta/2)*dr^2) : 0.
        
        P_pn_phi = postNewtonian ? (mt/r^2)*(4 - 2*eta)*dr*dphi : 0.
        
        P_gw_r   = gwEmission ? ( (8/5)*mu*(mt/r^3)*(2*v^2 + (8/3)*(mt/r) )) : 0.
        
        P_gw_phi = gwEmission ? ( (8/5)*mu*(mt/r^3)*(v^2 + 3*(mt/r) )) : 0. 
        
        # differential equations 
        du[1] = dr  
        du[2] = dphi 
        du[3] = -(-P_gw_r + P_df)*(dr) + P_pn_r - (mt/r^2) + r*(dphi^2)  # ddr
        du[4] = -(P_gw_phi + P_df)*(dphi) + P_pn_phi - (2*dr*dphi/r)          # ddphi
    end
    
    #---------------------Call-back functions---------------------------------------------
    
    function terminate_condition(u,t,integrator)       # condition at which the integration terminates
        u[1]< kwargs[:a_final]*r_isco(; kwargs...)
    end
    
    function terminate_affect!(integrator)
        terminate!(integrator)
    end
    
    terminate_cb = DiscreteCallback(terminate_condition,terminate_affect!)
 
    
    #-------------------------Calling the solver--------------------------------------------
    
    # call-backs do not work with lsoda() but they do work with other algorithms
    
    prob = ODEProblem(orbit!, u0, tspan, m, callback=terminate_cb)  # exclude the callback part when using lsoda()
    alg = VCABM() #lsoda() #AN5() #VCABM5() #Tsit5()  #DP5() #AutoVern7(Rodas5())  #Vern9(lazy=false) #Feagin12() #Vern7() 
    @time sol = solve(prob, alg, abstol=1e-14, reltol=1e-12, saveat=collect(range(initial, final, step=resolution))[1:end-1], dense=false, maxiters=Int(1e9))    # For predefined time-steps 
    #@time sol = solve(prob, alg, abstol=1e-14, reltol=1e-12, dense=true, maxiters=Int(1e9))                   # For adaptive time-steps
    return sol

end
