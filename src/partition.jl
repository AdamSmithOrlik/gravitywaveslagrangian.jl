###############################
# AUTHOR: Adam Smith-Orlik
# DATE: March 2023
# DESCRIPTION: Runs the evolution in n partitions to optimise memmory usuage. Saves all orbital data to a structured hdf5 file
# with the name specified in the dictionary as "run_name" located in the data/ directory of this package. 
# CONTACT: asorlik@yorku.ca
###############################

using HDF5

include("evolution.jl")

const hz_to_invpc = 1.029e8

function run_partition(; kwargs...)
    
    """
    Inputs
        n: Number of partitions that the code will attempt to split the evolution from a_initial to a_final into. 
        The size of n is determined by estimating the total number of time steps in the entire evolution by taking 
        the product of the number of orbits and the time steps per orbit (norbits and t_orbital respectivley) from 
        the input dictionary and dividing it by n. 
    
        If norbits is not large enough for the evolution to terminate at a_final the code will terminate early.
    
        kwagrs: Input dictionary including mass of binary object, luminosity distance, conditions for evolutions,
        flags for adding or removing forces, properties of the dark matter halo including rho_0 and relativisitc 
        corrections, etc. 
    
    Ouput
        hdf5 file with folder to save all the input parameters as attributes of a folder "Input Parameters" as well
        as the t, r, dr, phi, dphi array solution from the evolve function for each block in the format: parameter->
        partition number->data array
    
    """
    # target partition size and maximum partitions 
    partition_size = 25e3 #25k time elements per partition keeps RAM pressure low
    
    # mass 
    mass = mass_A(; kwargs...)
    mt = mass[3]
    Mbh = kwargs[:Mbh]
    
    # time parameters 
    t_start, t_end = time_params(mass; kwargs...) # t start and approximate t end
    norbits = kwargs[:norbits] # upper limit of orbits before code terminates prematurely 
    t_orbital = kwargs[:t_orbital] # number of time steps per orbit
    nSteps = Int(norbits * t_orbital) # upper limit on time steps for total evolution 
    step_size = (t_end - t_start) / (nSteps-1) # the time resolution of the evolution 
    
    # define n, the number of maximum partitions         
    frac_n_elements = partition_size
    n = Int(round((nSteps+1)/frac_n_elements)) # maximum partitions
    println("Maximum partitions: " * string(n))
    
    # directory 
    data_directory = pwd() * "/data/"
    run_name = kwargs[:Run_name]
    filename = data_directory * run_name * "_insprial_data.h5"
    
    # risco 
    risco = r_isco(; kwargs...)
    
    # Limits of integration in multiples of r_isco
    a_initial = kwargs[:a_initial]
    a_final = kwargs[:a_final]
    
    # termination condition 
    r_terminate = a_final*risco
    
    initial = t_start
    final = initial
    steps = 0 # tracks total time steps 
    
    # Useful quantities
    D = kwargs[:D]
    
    # evolution flag 
    evolution_complete = false
    nmax = n 
    
    # Partitioning 
    for i in 1:n
        
        if i == 1
            
            println("BEGIN PARTITIONING")
            println("Evolution begins at " * string(a_initial) * " risco and is set to terminate at " * string(a_final) * " risco" )
            
            # creates an hdf5 with the input parametes for the run
            # "w" will overwrite any hdf5 with the same filename 
            h5open(filename, "w") do f 
                # initilize folders for later use 
                create_group(f, "t")
                create_group(f, "r")
                create_group(f, "dr")
                create_group(f, "phi")
                create_group(f, "dphi")
                create_group(f, "final conditions")
                create_group(f, "run report")
                # create a folder to hold all input parameters as attributes
                g = create_group(f, "Input parameters")
                for (k,v) in kwargs
                    attrs(g)[string(k)]=v
                end # for loop for attributes
                println("hdf5 data file initiliazed...")
            end # hdf5 writer
            
        end # if statement
        
        # Create the time array for block i 
        final = frac_n_elements * step_size
        t_f = (i-1) * final

        # Load initial conditions from the mass array 
        if i == 1
            # if first block, calculate initial conditions and create a folder to store them
            u0 = initial_conditions(mass; kwargs...) 
            # write intial conditions to hdf5
            # "cw" ensures we dont overwrite existing hdf5 file
        else
            # read in the initial conditions 
            f = h5open(filename,"r")
            u0 = read(f["final conditions"][string(i-1)]) # i-1 becasue it uses final conditions of the previous block as the new initial conditions 
            close(f)
        end # if statement to load IC
        
        # Evolution
        solution = Evolve(u0, initial, final, t_f, step_size; kwargs...)
        
        # Saving results 
        t, r, phi, dr, dphi = [ solution.t,  solution'[:,1], solution'[:,2], solution'[:,3], solution'[:,4] ] 
        
        # saving memory by erasing future unused arrays
        solution = nothing
        
        t_elemenets = length(t)
        dt_total = t[2] - t[1]
        steps += length(t)
        
        # erasing the last step which was just saved (to avoid overlap in FFT)
        t = t[1:end-1]
        r = r[1:end-1]
        phi = phi[1:end-1]
        dr = dr[1:end-1]
        dphi = dphi[1:end-1]
        #saving the new final conditions
        u_new = [r[end], phi[end], dr[end], dphi[end]]
        
        a_end = r[end]/risco
        
        println("Partition " * string(i) * " evolution complete at " * string(steps) * " total time steps")
        println("Radius of inspiral has reached " * string(a_end) * " risco")
    
        if r[end] <= r_terminate
            evolution_complete = true
        end
        
        # Save data 
        h5open(filename, "cw") do f             
            # write datasets for each partition i 
            write(create_dataset(f["t"], string(i), Float64, (length(t),)), t) 
            write(create_dataset(f["r"], string(i), Float64, (length(r),)), r)
            write(create_dataset(f["dr"], string(i), Float64, (length(dr),)), dr)
            write(create_dataset(f["phi"], string(i), Float64, (length(phi),)), phi)
            write(create_dataset(f["dphi"], string(i), Float64, (length(dphi),)), dphi)
            write(create_dataset(f["final conditions"], string(i), Float64, (length(u_new),)), u_new)
        end # hdf5
        
                
        # Check if evolution is complete  
        if evolution_complete
#             println("Block number " * string(i) *" completed with "* string(t_elemenets) *" time elements at "* string(dt_total ) *" resolution")
            println("EVOLUTION COMPLETE")
            
            # run report 
            h5open(filename, "cw") do f 
                rep = f["run report"]
                attrs(rep)["Evolution"]="Success"
                attrs(rep)["Partitions"]=i
                attrs(rep)["Total steps"]=steps
                attrs(rep)["Time resolution"]=step_size
                attrs(rep)["Final radius [r_isco]"]=a_end
            end # writing
            
            break # exit loop if evolution is complete 
            
        else
            println("CONTINUING EVOLUTION")
            println("------------------------------------------------------------")
            
            # run report 
            if i == nmax
                println("NOT ENOUGH ORBITS TO COMPLETE EVOLUTION.")
                h5open(filename, "cw") do f 
                    rep = f["run report"]
                    attrs(rep)["Evolution"]="Failed"
                    attrs(rep)["Partitions"]=i
                    attrs(rep)["Total steps"]=steps
                    attrs(rep)["Time resolution"]=step_size
                    attrs(rep)["Final radius [r_isco]"]=a_end
                end # writing 
            end # if condition 
            
        end # evolution complete condition
        
            
    end # for loop 
    
end # function end 