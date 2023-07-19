###############################
# AUTHOR: Adam Smith-Orlik
# DATE: March 2023
# DESCRIPTION: Runs the evolution in n partitions to optimise memmory usuage. Saves all orbital data to a structured hdf5 file
# with the name specified in the dictionary as "run_name" located in the data/ directory of this package. 
# CONTACT: asorlik@yorku.ca
###############################

using HDF5

include("evolution.jl")
include("strain.jl")

function run_partition(; kwargs...)
    
    """
    Inputs
        positional args: none.
    
        kwagrs: Input dictionary including mass of binary object, luminosity distance, conditions for evolutions,
        flags for adding or removing forces, properties of the dark matter halo including rho_0 and relativisitc 
        corrections, etc. 

    Ouputs
        An h5 file is created where the input dictionary is saved in the Input parameters folder as attributes. A Run report folder will contain the 
        status of the run, along with key information like the total time of evolution, the partitions, the resolution and so on. Finally the gravity
        waves folder contains the frequency and the strain of the total evolution, whether complete or continued. 

    Functionality

        PARTITIONING
        The run partition function splits up the integration domains for the evolution of the binary system to reduce the RAM pressure of the
        orbital computation in the Evolve function. It does so by defining the number of partitions, n, to split the computation into as the number
        of steps, dt, between t_start and t_end. 

        dt determines the time equivalent size of each partition. If this value is greater than 0.1 the RAM becomes greater than 30GBs. 

        The code will then partition the runs into blocks of size dt and solve given initial conditions from the block prior, adding to the 
        characteristic strain function and saving orbital information if requested by setting save_orbital_data to true in the dictionary. 

        NB: saving orbital data becomes a file size issue for evolutions that persist longer than 1e4 years. If one needs orbital data for 
        long evolutions, it is suggested to set the orbital_resolution parameter to a larger value i.e. 1e-3.

        CONTINUED CODE
        If the code does not complete the evolution in the timespan entered, i.e. the termination condition is not met in the number of partitions
        allocated, all of the data is saved into the h5 file, and if the evolution is run again the code will start where it left off. 
    
        NB: Currently (April 2023) the continued functionality does not work with the strain. This is due to the fact that it is
        difficult to save the strain information between runs. This will be added in a later version. For now runs need to be given enough time to evolve in one run.
 
    
    """
#------------------------------------------#
#           LOAD INPUT PARAMETERS          #
#------------------------------------------#
    println("Loading input parameters...")
    hz_to_invpc = 1.029e8
    year_to_pc = 0.3064
    D = kwargs[:D]
    lisa_bandwidth = kwargs[:lisa_bandwidth].*hz_to_invpc
    
    # mass 
    mass = mass_A(; kwargs...)
    mt = mass[3]
    Mbh = kwargs[:Mbh]
    
    # time parameters 
    t_start = kwargs[:t_start]
    t_endyr = kwargs[:t_max]  # maximum time set for evolution in years
    t_end = t_endyr * year_to_pc
    dt = kwargs[:dt] #0.005 # dt > 0.1 uses >30GB of memory allocation 
    
    # time step array to pass to evolve 
    partitions = range(t_start, t_end, step=dt)
    n = length(partitions)
    start = 1
    println("Maximum partitions: " * string(n))
    
    # resolution of the solution output arrays 
    resolution = kwargs[:time_resolution]
    orbital_res = kwargs[:orbital_resolution]
    orbital_data_flag = kwargs[:save_orbital_data]
    
    # used to sample the output arrays, saving memory 
    res_ratio = orbital_res >= resolution ? Int(round(orbital_res/resolution)) : throw(DomainError(orbital_res, "Orbital resolution cannot be less than the time resolution."))
    resolution >= dt ? throw(DomainError(resolution, "Must be less than $dt")) : nothing
 
     # directory 
     data_directory = pwd() * "/data/"
     isdir(data_directory) ? nothing : mkdir(data_directory)
 
     run_name = kwargs[:Run_name]
     filename = data_directory * run_name * "_insprial_data.h5"
    
    # risco 
    risco = r_isco(; kwargs...)
    
    # Limits of integration in multiples of r_isco
    a_initial = kwargs[:a_initial]
    a_final = kwargs[:a_final]

#------------------------------------------#
#              PLOT HALO MODEL             #
#------------------------------------------#
plot_halo_model = kwargs[:plot_halo_model] 
if plot_halo_model
    plot_dm_spike(; kwargs...)
end

#------------------------------------------#
#       CHECK IF RUN IS CONTINUED          #
#------------------------------------------#
    continued_run = false
    if isfile(filename)
        a = h5open(filename)
        try
            read_attribute(a["Run report"],"Evolution")
        catch e
            println("HDF5 is not a gravity wave dataset.")
        end
        # check if run finished
        if read_attribute(a["Run report"],"Evolution") == "Success"
            println("$filename already exists and succeded evolving. Please change the filename and re-run.")
            return nothing
        elseif read_attribute(a["Run report"],"Evolution") == "Failed"
            println("$filename aleady exists, and will continue running for $t_endyr more years.")
            continued_run = true
        end # if 
        close(a) # close hdf5 
        println("access #1: ", a)
    end # check existing file
    

    
    if continued_run
        # load all input parameters
        a = h5open(filename) 
# TEST
            # nparts = length(a["t"])
            # if orbital data is not saved the above will not work...
            # can retreive the npartis from the partitions attribute 
            nparts = read_attribute(a["Run report"],"Partitions")
            
            last = string(nparts)
            start = nparts + 1
        
            t_max = read_attribute( a["Input parameters"], "t_max" ) * year_to_pc
            t_start = read( a["t"][last] )[end]
            t_end = t_start + t_max

            partitions = range(t_start, t_end, step=dt)
            n += length(partitions)
        
            a_start = round(read( a["r"][last] )[end]/ risco, digits=2)
        close(a)
        println("access #2: ", a)
        println("CONTINUING PARTIONING")
        println("Starting at $a_start risco and set to terminate at $a_final")
        
    end
    
#------------------------------------------#
#               PARTITIONING               #
#------------------------------------------#
    nmax = n - 1
    println("Will start at n=", start, " and go until n=", nmax)
    t_step = 0
    evolution_complete = false
    steps = 0 # tracks total time steps 
    r_terminate = a_final*risco
    strain = 0
    lengths = 0
    freq = 0
    final_strain = 0
    size_of_strain_array = 0
    
    # Partitioning 
    for i in start:nmax

        t_step += 1
        
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
                create_group(f, "gravity waves")
                create_group(f, "final conditions")
                create_group(f, "Run report")
                # create a folder to hold all input parameters as attributes
                g = create_group(f, "Input parameters")
                for (k,v) in kwargs
                    attrs(g)[string(k)]=v
                end # for loop for attributes
                println("hdf5 data file initiliazed...")
            end # hdf5 writer
            
        end # if statement

        # Load initial conditions from the mass array 
        if i == 1
            # if first block, calculate initial conditions and create a folder to store them
            u0 = initial_conditions(; kwargs...) 
        else
            # read in the initial conditions 
            a = h5open(filename)
                u0 = read(a["final conditions"][string(i-1)]) # i-1 becasue it uses final conditions of the previous block as the new initial conditions 
            close(a)
        end # if statement to load IC
        
        # Evolution
        initial = partitions[t_step]
        final = partitions[t_step+1]
 
#         println(initial)
#         println(final)
        solution = Evolve(u0, initial, final, resolution; kwargs...)
        
        # Saving results 
        t, r, phi, dr, dphi = [ solution.t,  solution'[:,1], solution'[:,2], solution'[:,3], solution'[:,4] ] 
        
        # saving memory by erasing future unused arrays
        solution = nothing
        
        t_elemenets = length(t)
        dt_res = t[2] - t[1]
        steps += length(t)
        
        #saving the new final conditions
        u_new = [r[end], phi[end], dr[end], dphi[end]]
        
        a_end = round(r[end]/risco, digits=2)
        time_elapsed = round(t[end]/year_to_pc, digits=2)

        println("Partition $i evolution complete at $steps total time steps")
        println("Radius of inspiral has reached $a_end risco")
        println("Time elapsed since beginning of evolution: $time_elapsed")
    
        # check if evolution have complete 
        if r[end] <= r_terminate
            evolution_complete = true
        end
        
        # calculate the gravity wave strain 
        if i == 1
            ## ADD compatibility with continued runs...
            waveform = strain2d(u_new, t, r, phi; kwargs...)
            freq, s  = strainFFT2d(t, waveform, lisa_bandwidth)
            strain = zeros(ComplexF64, length(s))
            size_of_strain_array = length(strain)
            # size_s = length(s)
            # println("Size of strain array is now: $size_s")
            strain .+= s
            
            # first element of the lengths array
            lengths = [length(waveform)]
            
        else
            if !evolution_complete
                waveform = strain2d(u_new, t, r, phi; kwargs...)
                f, s = strainFFT2d(t, waveform, lisa_bandwidth)
                # add to strain as a superpositon of waveforms 
                # size_s = length(s)
                # println("Size of strain array is now: $size_s")
                if length(s) != size_of_strain_array
                    s = s[1:size_of_strain_array]
                end
                size_s = length(s)
                # println("Size after check: $size_s")
                strain .+= s
                # add length of waveform to lengths array
                push!(lengths, length(waveform))
                
            else # if evolution is complete the final array will be of a different length than the preceeding arrays, and must be interpolated prior to being added to the total strain
                waveform = strain2d(u_new, t, r, phi; kwargs...)
                f, s = strainFFT2d(t, waveform, lisa_bandwidth)
                push!(lengths, length(waveform))
                
                # last step needs to be interpolated becasue it need not be the same size as previous partitions 
                f = real(f)
                freq = real(freq)
                
                conditions = (freq .> f[1]) .&& (freq .< f[end])
                freq = freq[conditions]
                # interpolate to match size 
                last_strain = imaginary_interpd(f, s, freq)
                # define normalization constants 
                c1, c2 = lengths[1]/sum(lengths[1:end]), lengths[end]/sum(lengths[1:end-1])
                # add interpolated final step to total strain 
                final_strain = 2 * freq * c1 .* abs.( strain[conditions] .+ (c2 .* last_strain) )
            end
            
        end
        
        # sample the output arrays by the ratio of input resolutions to save space
        t = t[1:res_ratio:end]
        r = r[1:res_ratio:end]
        dr = dr[1:res_ratio:end]
        phi = phi[1:res_ratio:end]
        dphi = dphi[1:res_ratio:end]
        
        
        # Save data 
        w = h5open(filename, "cw")   
            if orbital_data_flag
                # write datasets for each partition i 
                write(create_dataset(w["t"], string(i), Float64, (length(t),)), t) 
                write(create_dataset(w["r"], string(i), Float64, (length(r),)), r)
                write(create_dataset(w["dr"], string(i), Float64, (length(dr),)), dr)
                write(create_dataset(w["phi"], string(i), Float64, (length(phi),)), phi)
                write(create_dataset(w["dphi"], string(i), Float64, (length(dphi),)), dphi)
            end    
            write(create_dataset(w["final conditions"], string(i), Float64, (length(u_new),)), u_new)
        close(w) # hdf5
        
                
        # Check if evolution is complete  
        if evolution_complete
#             println("Block number " * string(i) *" completed with "* string(t_elemenets) *" time elements at "* string(dt_total ) *" resolution")
            println("EVOLUTION COMPLETE")
            
            # run report 
            w = h5open(filename, "cw") 
                rep = w["Run report"]
                attrs(rep)["Evolution"]="Success"
                attrs(rep)["Partitions"]=i
                attrs(rep)["Time step"]=dt
                attrs(rep)["Total steps"]=steps
                attrs(rep)["Time resolution"]=resolution
                attrs(rep)["Final radius [r_isco]"]=a_end
            
                write(create_dataset(w["gravity waves"], "frequency", Float64, (length(freq),)), freq)
                write(create_dataset(w["gravity waves"], "strain", Float64, (length(final_strain),)), final_strain)
            close(w) # writing
            
            break # exit loop if evolution is complete 
            
        else
            println("CONTINUING EVOLUTION")
            println("------------------------------------------------------------")
            
            # run report 
            if i == nmax
                t_end = round(t[end]/year_to_pc, digits=2)
                println("EVOLUTION NOT ABLE TO COMPLETE IN $t_end YEARS" )
                w = h5open(filename, "cw")
                    rep = w["Run report"]
                    attrs(rep)["Evolution"]="Failed"
                    attrs(rep)["Partitions"]=i
                    attrs(rep)["Time step"]=dt
                    attrs(rep)["Total steps"]=steps
                    attrs(rep)["Time resolution"]=resolution
                    attrs(rep)["Final radius [r_isco]"]=a_end
                
                    write(create_dataset(w["gravity waves"], "frequency", Float64, (length(freq),)), freq)
                    write(create_dataset(w["gravity waves"], "strain", ComplexF64, (length(strain),)), strain)
                close(w) # writing 
            end # if condition 
            
        end # evolution complete condition
        
            
    end # for loop 
    
end # function end 

