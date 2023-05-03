# gravitywaveslagrangian.jl

# Gravitational wave signal from intermediate mass-ration inspiral of a black hole and a compact object using a novel Lagrangian method
Using the lagrangian frameowrk outlined in (arxiv:...) [1] we model the orbital inspiral of an intermediate mass-ratio black hole binary system consisting of a black hole (BH) and a compact object (CO) where the mass ratio BH/CO is >1 and <1e5 in the presence of a dark matter (DM) spike enveloping the system. In the frameowrk we consider the frist two post-Newtonian corrections and dynamical friction between the inspiraling bodies and the DM and calculate the characteristic strain observable by modern and next generation gravity wave (GW) detectors. Our analysis can be performed in a Newtonian or fully relativistic framework.

# How to use the julia package 
Usage:
The package contains one module--gravitywavelagrangian--contianed in the file gravitywavelagrangian.jl with one method--run(:defualt_dictionary)--that runs the inspiral with the input parameters defined in the defualt_dictionary also defined in gravitywavelagrangian.jl. 

There are two ways to run the code. 

1)In the Julia REPL: In the terminal navigate to the root directory of the package and type $ julia. Inside the julia REPL type $ ] which changes to the package manager. Inside package manager do $ activate . and then hit backspace to go back to julia. Now compile the package with $ using gravitywaveslagrangian. After compilation you can do $ gravitywaveslagrangian.run(:your dictionary name) to run the code. The printed output will be displayed in the terminal and the outputted data will be saved to a file called <Run_name>.h5 in the data folder where 'Run name'
is defined in the input dictionary.

2)Using a bash script: In the terminal navigate to the root directory of the package and locate the file run.sh. To run this file execute $ ./run.sh and the above steps will be run internally, with the printed output being piped to a file "output.txt" and the data saved to a file called <Run_name>.h5 in the data folder where 'Run name' is defined in the input dictionary.

NOTE: To make your own run, copy the default_dictionary and change the name along with the parameters of your choice. Then add the dictionaries 
name to the list in the run() function. If you are running your own dictionary, be sure to change the command in run.sh to "gravitywaveslagrangian.run(:your dictionary name)". The code will then evaluation the dictionary name and cross-check that it is defined in the gravitywaveslagrangian scope and then evaluate it as a normal dictionary. 

Input Parameters:
Inside gravitywavelagrangian.jl there is a dictionary called "default_dictionary" including all the default parameters for the run. The input
parameters are defined in the default_dictionary, with the important ones listed below. 

Run_name - the name assigned to the hdf5 file that is created for the run
a_initial - The initial separation between the binary objects in multiples of r_isco (r_isco = 6 * Mbh[pc] O(1e-14))
a_final - the target minimum radius at which point the code will terminate and save the insprial data

Mbh - Mass of the black hole
Mc - Mass of the compact object where Mc < Mbh
time_resolution - Resolution passed to the differential solves, recommended 1e-8
t_start - Beginnning of evolution, keep at 0
t_max - Upper limit of years to inspiral. If system does not reach target inspiral in the time alloted that indicates a problem.
dt - The time steps for the total evolution time alloted. If this value is > 0.1 RAM usage is >30GB. Keep low to reduce computational demand. 

save_orbital_data - Flag to save orbital data. If the evolution is more than 10 r_isco this takes up a lot of memory (O(GB)) 
orbital_resolution - This will sample the outputted solution array which are by default computed with the time_resolution. If you want the highest resolution then set orbital_resolution == time_resolution. If the orbital data can be saved with less resolution then set orbital_resolution < time_resoltuion. orbital_resolution > time_resoltuion will throw an error. 

gwEmission - Flag for the 2.5 PN correction, or the gravity wave emissions
dynamicalFriction - Flag for the dynamical friction term
postNewtonian - Flag for the first order PN correction 

dmModel - Either Newtonian (NFW+spike) or Relativistic (fully-relativistic spike model) [2]
dfModel - Either Newtonian or Relativistic [3,5]

rho0 - Normalization density for the DM model [3]

Output:
The run() method prints periodic updates on the progress of the insprial including the radius--distant between bianries--as a function of the evolution time and whether or not the evolution has completed, i.e. the binaries have merged. 

Moreover, and hdf5 file will be created in the data folder that contains the following folders: Input parameters, containing all inputted parameters from the input dictionary as attributes of the folder; Run report, which includes information about whether the evolution was successful, the start and stop radii, and how many partitions the insprial took; orbital parameters r, t, phi and derivitives dr and dphi are saved for each partition IF the save_orbital_data is set to true; final conditions which are passed to the differential solver and saved for each partition; and gravity waves, which include the calculated characteristic strain and the frequency calculated directly from the orbital data.
   
# Analysis 
In addition to the julia files, there is an analysis folder which contains a script, analysis.py, and a notebook, analysis.ipynb. The python script contains a class called 'analysis' that reads the outputted hdf5 file from the julia package and creates an object with all of the information from the hdf5 stored in easily accessible attributes or methods. 

To see a list of all the availbale data stored, first load an object as follows:

from analysis import analysis
data = analysis(filename, verbose=True)

and to see what data is available run

data.info() .

For more information and examples see the analysis notebook. 

References:
[1]: Our paper
[2]: https://arxiv.org/pdf/2204.12508.pdf (Relativistic DM model)
[3]: https://arxiv.org/pdf/1408.3534.pdf (Eda model)
[4]: https://arxiv.org/abs/1102.5192 (PN correction)
[5]: https://arxiv.org/pdf/2106.08280.pdf (DF correction)
