import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
from scipy.fft import fft, fftfreq
import scipy as sc
import h5py
import time as tm
import os 

# change for windows
data_path = "../data/"
# print(data_path)

# constants
hz_to_invpc = 1.029e8
s_to_pc = 9.716e-9
m_to_pc = 3.241e-17
solar_mass_to_pc = 4.8e-14
g_cm3_to_invpc2 = 7.072e8
year_to_pc = 0.3064
yr_to_s = 31556952

# analysis class object
class analysis:
    
    def __init__(self, filename, verbose=False):
        self.filename = filename
###########################################################
#                        LOAD DATA                        # 
###########################################################
        try:
            file = data_path + self.filename
            data = h5py.File(file, 'r')
        except FileNotFoundError:
            print("File '%s' not found" %filename)
        else:
            # check that file succeeded in evolving
            if not data['Run report'].attrs['Evolution']==b"Success":
                raise Exception("The hdf5 file you are attempting to load did not succeed in evolving. Please try another file.")
                
###########################################################
#                   LOAD RUN PARAMETERS                   # 
###########################################################
            
            # load input parameters into a dictionary 
            self.inputs = {}
            for key in data['Input parameters'].attrs.keys():
                self.inputs[key] = data['Input parameters'].attrs[key]
            
            if verbose:
                print("Input parameters loaded into dictionary 'inputs.'")
                print("Inputs = ", self.inputs)
            
            # load the run report
            self.run_report = {}
            for key in data['Run report'].attrs.keys():
                self.run_report[key] = data['Run report'].attrs[key]
                
            if verbose:
                print("Run parameters loaded into dictionary 'run_report.'")
                print("Run report = ", self.run_report)
                
            # other quantities 
            self.risco = r_isco(data['Input parameters'].attrs['Mbh']) # result in pc
            self.Mbh = solar_mass_to_pc * data['Input parameters'].attrs['Mbh']
            self.Mc = solar_mass_to_pc* data['Input parameters'].attrs['Mc']
            self.mt = self.Mbh + self.Mc
            self.mu = (self.Mbh * self.Mc) / self.mt 
            self.D = data['Input parameters'].attrs['D']
            self.lisa_bandwith = data['Input parameters'].attrs['lisa_bandwidth']
            
            f_gw = np.geomspace(self.lisa_bandwith[0]*hz_to_invpc, self.lisa_bandwith[1]*hz_to_invpc, 100)
            self.fgw = f_gw
            self.noise_curve = np.sqrt(f_gw*NoiseSpectralDensity(f_gw))
            
###########################################################
#                   LOAD GW PARAMETERS                    # 
###########################################################
            if verbose:
                print("Loading gravity wave data.")

            self.strain = data['gravity waves/strain'][:]
            self.frequency = data['gravity waves/frequency'][:]
            
###########################################################
#                 LOAD ORBITAL PARAMETERS                 # 
###########################################################
            # Check if orbital parameters are saved
            # number of partitions
            self.partitions = data["Run report"].attrs["Partitions"]
            keys = list(map(str, range(1,self.partitions+1,1)))
            
            if len(data['r'])==0:
                self.orbital_data = False
                if verbose:
                    print("No orbital data saved.")
            else:
                self.orbital_data = True
                print("Loading orbital data...")
                start = tm.time()
                
                # load datasets into arrays 
                self.r = [np.array(data["r"][key]) for key in keys]
                self.dr = [np.array(data["dr"][key]) for key in keys]
                self.phi = [np.array(data["phi"][key]) for key in keys]
                self.dphi = [np.array(data["dphi"][key]) for key in keys]
                self.t = [np.array(data["t"][key]) for key in keys]
                
                # flatten the arrays and redefine them 
                self.r = np.array(np.concatenate(self.r).flat)
                self.dr = np.array(np.concatenate(self.dr).flat)
                self.phi = np.array(np.concatenate(self.phi).flat)
                self.dphi = np.array(np.concatenate(self.dphi).flat)
                self.t = np.array(np.concatenate(self.t).flat)
                
                self.semi_major_axis = semi_major_axis(self.r, self.dr, self.dphi, self.mt)
                
                end = tm.time()
                print("Orbital arrays with %d elements loaded in %d seconds." %(len(self.t), end-start))
                
            
###########################################################
#                CLASS METHODS (dataframes)               # 
###########################################################
        if verbose:
            print("Creating dataframe accesible via 'object.characteristic_strain()'.")
            print("Creating dataframe accesible via 'object.dephasing()'.")
            print("Creating dataframe accesible via 'object.numer_of_cycles()'.")

    def characteristic_strain(self):
        char_strain = np.array(list(zip(list(self.frequency), list(self.strain))))
        return pd.DataFrame(char_strain, columns=['f [Hz]', 'h_c(f)'])
        
    def dephasing(self):
        if not self.orbital_data:
            print("Cannot calcualte dephasing without orbital data.")
            return None
        
        hplus = strain2d(self.r, self.phi, self.t, self.mu, self.D)
        t = self.t / year_to_pc
        
        deph = np.array(list(zip(list(t), list(hplus))))
        
        return pd.DataFrame(deph, columns=['time [yr]', 'h_+(t)'])
    
    def number_of_cycles(self):
        if not self.orbital_data:
            print("Cannot calcualte dephasing without orbital data.")
            return None
        
        frequency = np.sqrt(self.mt/self.semi_major_axis**3)/np.pi
        cycles = Ncycles(self.phi)
        
        cyc = np.array(list(zip(list(frequency), list(cycles))))
        
        return pd.DataFrame(cyc, columns=['frequency [Hz]','Cycles'])
    
###########################################################
#                   ADDITIONAL FUNCTIONS                  # 
###########################################################
    def info(self):
        return [attr for attr in dir(self) if not attr.startswith('__')]
    
    def copy(self):
        newobj = analysis(self.filename)
        return newobj
    
###########################################################
#               ACCESORY ANALYSIS FUNCTIONS               # 
###########################################################
def strain(u0, r, phi, t, mu, D):
    # Observables
    theta_0 = 0    # inclination_angle
    phi_0 = 0      # pericenter_angle
    
    # Scaling 
    r0 = u0[1]
    rs = r/r0
    dt = t[2]-t[1]
    
    # rotating body parametes 
    x = rs*np.cos(phi)
    y = rs*np.sin(phi)
    z = np.zeros(len(x))
    
    # Quadrupole Moment Tensor
    Q  = np.array([ [mu*x*x, mu*x*y, mu*x*z],
                    [mu*y*x, mu*y*y, mu*y*z],
                    [mu*z*x, mu*z*y, mu*z*z] ] )
    
    # Derivitive of Q w.r.t time 
    def Mdt(Q):     
        """
        Returns the first derivative of Quadrupole Moment Tensor with respect to time
        """
        
        dQdt = np.array([ [np.gradient(Q[0][0], dt), np.gradient(Q[0][1], dt), np.gradient(Q[0][2], dt)], 
                          [np.gradient(Q[1][0], dt), np.gradient(Q[1][1], dt), np.gradient(Q[1][2], dt)],
                          [np.gradient(Q[2][0], dt), np.gradient(Q[2][1], dt), np.gradient(Q[2][2], dt)] ] )
        
        return dQdt
    
    # twice differentiate 
    d2Qd2t = np.array(Mdt(Mdt(Q))) * (r0**2)


    h_plus =  (1.0/D) * np.array(   d2Qd2t[0][0]*(np.cos(phi_0)**2 - np.sin(phi_0)**2 * np.cos(theta_0)**2) 
                          + d2Qd2t[1][1]*(np.sin(phi_0)**2 - np.cos(phi_0)**2 * np.cos(theta_0)**2) 
                          - d2Qd2t[2][2]*(np.sin(theta_0)**2) 
                          - d2Qd2t[0][1]*(np.sin(2*phi_0)*(1.0 + np.cos(theta_0)**2))
                          + d2Qd2t[0][2]*(np.sin(phi_0)* np.sin(2*theta_0)) 
                          + d2Qd2t[1][2]*(np.cos(phi_0)*np.sin(2*theta_0))     ) 
    
    h_cross = (1.0/D) * np.array(   (d2Qd2t[0][0]-d2Qd2t[1][1])*np.sin(2*phi_0)*np.cos(theta_0)
                             + 2*d2Qd2t[0][1]*np.cos(2*phi_0)*np.cos(theta_0) 
                             - 2*d2Qd2t[0][2]*np.cos(phi_0)*np.sin(theta_0) 
                             + 2*d2Qd2t[1][2]*np.sin(theta_0)*np.sin(phi_0)    )
    
    return np.array([h_plus, h_cross])

def strain2d(r, phi, t, mu, D):
    # Observables
    theta_0 = 0    # inclination_angle
    phi_0 = 0      # pericenter_angle
    
    # Scaling 
    r0 = r[0]
    rs = r/r0
    dt = t[2]-t[1]
    
    # rotating body parametes 
    x = rs*np.cos(phi)
    y = rs*np.sin(phi)
    
    # Quadrupole Moment Tensor
    Q  = np.array([ [mu*x*x, mu*x*y],
                    [mu*y*x, mu*y*y] ] )
    
    # Derivitive of Q w.r.t time 
    def Mdt(Q):     
        """
        Returns the first derivative of Quadrupole Moment Tensor with respect to time
        """
        
        dQdt = np.array([ [np.gradient(Q[0][0], dt), np.gradient(Q[0][1], dt)], 
                          [np.gradient(Q[1][0], dt), np.gradient(Q[1][1], dt)] ] )
        
        return dQdt
    
        # twice differentiate 
    d2Qd2t = np.array(Mdt(Mdt(Q))) * (r0**2)
    
    h_plus =  (1.0/D) * np.array(   d2Qd2t[0][0]*(np.cos(phi_0)**2) 
                          + d2Qd2t[1][1]*(-np.cos(phi_0)**2 * np.cos(theta_0)**2)  )
    
    return np.array(h_plus)


def strainFFT(t, strain, f_bin):
    N = len(t)
    dt = t[2] - t[1]
    
    h_plus_fft = fft(strain[0,:])/(2*np.pi*N)
    h_cross_fft = fft(strain[1,:])/(2*np.pi*N)
    xf = fftfreq(N, dt)  #[1:Int(N ÷ 2)] 
    
    f_range = np.where((xf > f_bin[0]) & (xf < f_bin[1]))
    
    h_plus_fft = h_plus_fft[f_range]
    h_cross_fft = h_cross_fft[f_range]
    xf = xf[f_range]
    
    return np.array([xf, h_plus_fft, h_cross_fft])

def strainFFT2d(t, strain2d, f_bin):
    N = len(t)
    dt = t[1] - t[0]
    
    h_plus_fft = fft(strain2d[:])/(2*np.pi*N)
    xf = fftfreq(N, dt) #[1:Int(N ÷ 2)]
    
    f_range = np.where((xf > f_bin[0]) & (xf < f_bin[1]))
    
    h_plus_fft = h_plus_fft[f_range]
    xf = np.real(xf[f_range])
    
    return [np.real(xf) , h_plus_fft]

def semi_major_axis(r, dr, dphi, mt):
    v   = np.sqrt(dr**2 + (r**2 * dphi**2))
    return mt*(1/np.abs(v**2 - 2*mt/r))

def Ncycles(phi):
    return phi / (2*np.pi)

def Ncycle_frequency(r, dr, dphi, mt):
    sma = semi_major_axis(r, dr, dphi, mt)
    return np.sqrt(mt/sma**3)/np.pi

def frequency(strainfft):
    return np.abs(strainfft[0]) / 1.029e8 # hz -> 1/pc

def imaginary_interpd(f, h, f_new):
    #f and h are the computed strain and frequency (say 1,000,000 elements) 
    #f_new is the new frequency we'd like to have (say 10,000 elements)  scipy.interpolate import interp1d
    
    real_h = interp1d(f, np.real(h), kind="previous")
    im_h = interp1d(f, np.imag(h), kind="previous")
    
    return real_h(f_new) + 1j*im_h(f_new)

def r_isco(Mbh):    
    """
    Radius of the Innermost Stable Circular Orbit (ISCO) of a Schwarzschild black hole with mass m
    """ 
    return 6.0*Mbh*solar_mass_to_pc
            
def NoiseSpectralDensity(f):
    P_oms = (1.5e-11 * m_to_pc)**2  * (1. + (2e-3 * hz_to_invpc/f)**4) / hz_to_invpc
    P_acc = (3e-15 * m_to_pc / s_to_pc**2)**2 * (1. + (0.4e-3 * hz_to_invpc/f)**2) * (1. + (f/8e-3/hz_to_invpc)**4) / hz_to_invpc
    f_s = 19.09e-3 * hz_to_invpc
    L = 2.5e9 * m_to_pc
    return 10./3./L**2  * (P_oms + 2.* (1. + np.cos(f/f_s)**2 ) * P_acc / (2.*np.pi*f)**4) * (1. + 6./10. * (f/f_s)**2)
            
            
            
            
                
            