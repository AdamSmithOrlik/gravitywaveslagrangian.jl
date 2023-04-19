###############################
# AUTHOR: Adam Smith-Orlik
# DATE: March 2023
# DESCRIPTION: Calculates the gravity wave characteristic strain using FFT's
# CONTACT: asorlik@yorku.ca
###############################

using PyCall

# init method to define a python function in the julia context
function __init__()
    py"""
    import numpy as np

    def grad(A, dx):
        return np.gradient(A, dx)
    """
end

# np = pyimport("numpy")
using SciPy


function test()
    println(py"grad"([1,2,3,4], 0.5))
end


function strain(mass, u0, D, r, phi, t)
    """
    This function takes in the orbit solution arrays (r and phi) and returns the gravitational wave pattern 
    """
    
    #---------Observer parameters--------------
    mu = mass[4]
    theta_o = 0    # inclination_angle
    phi_o = 0      # pericenter_angle
    
    #---------Scaling--------------------------
    r0 = u0[1]
    r /= r0
    T = t[2]-t[1]
    
    #---------Rotating body parameters---------
    x = r.*cos.(phi)
    y = r.*sin.(phi)
    z = zeros(length(x))
    
    #--------Quadrupole Moment Tensor----------
    Q  = [ [mu*x.*x, mu*x.*y, mu*x.*z],
           [mu*y.*x, mu*y.*y, mu*y.*z],
           [mu*z.*x, mu*z.*y, mu*z.*z] ]
    
    #--------Derivative Functions---------------
    
    function Mdt(Q)     
        """
        Returns the first derivative of Quadrupole Moment Tensor with respect to time
        """
        
        dQdt = [ [py"grad"(Q[1][1], T), py"grad"(Q[1][2], T), py"grad"(Q[1][3], T)], 
                 [py"grad"(Q[2][1], T), py"grad"(Q[2][2], T), py"grad"(Q[2][3], T)],
                 [py"grad"(Q[3][1], T), py"grad"(Q[3][2], T), py"grad"(Q[3][3], T)] ]
        
        return dQdt
    end
    
    
    function Mdt2(Q) 
        """
        Returns the second derivative of Quadrupole Moment Tensor with respect to time
        """
        
        dQ2dt2 = Mdt(Mdt(Q))
        
        return dQ2dt2*(r0^2)
    end
    
    
    d2Qd2t = Mdt2(Q)
    
    
    h_plus =  (1.0/D) * (   d2Qd2t[1][1].*(cos.(phi_o).^2 - sin.(phi_o).^2 .*cos.(theta_o).^2) 
                          + d2Qd2t[2][2].*(sin.(phi_o).^2 - cos.(phi_o).^2 .*cos.(theta_o).^2) 
                          - d2Qd2t[3][3].*(sin.(theta_o).^2) 
                          - d2Qd2t[1][2].*(sin.(2*phi_o).*(1.0 .+ cos.(theta_o).^2))
                          + d2Qd2t[1][3].*(sin.(phi_o).*sin.(2*theta_o)) 
                          + d2Qd2t[2][3].*(cos.(phi_o).*sin.(2*theta_o))     ) 
    
    h_cross = (1.0/D) * (   (d2Qd2t[1][1]-d2Qd2t[2][2]).*sin.(2*phi_o).*cos.(theta_o)
                             + 2*d2Qd2t[1][2].*cos.(2*phi_o).*cos.(theta_o) 
                             - 2*d2Qd2t[1][3].*cos.(phi_o).*sin.(theta_o) 
                             + 2*d2Qd2t[2][3].*sin.(theta_o).*sin.(phi_o)    )
    
    return [h_plus, h_cross]
    
end


function strain2d(u0, t, r, phi; kwargs...)
    
    """
    This function takes in the orbit solution arrays (r and phi) and returns the gravitational wave pattern 
    """
    
    #-------------Observer parameters----------
    mass = mass_A(; kwargs...)
    mu = mass[4]
    D = kwargs[:D]
    theta_o = 0    # inclination_angle
    phi_o = 0      # pericenter_angle
    
    #----------------Scaling-------------------
    r0 = u0[1]
    rs = r/r0
    T = t[2]-t[1]
    
    #---------Rotating body parameters---------
    x = rs.*cos.(phi)
    y = rs.*sin.(phi)
    
    #--------Quadrupole Moment Tensor----------
    Q  = [ [mu*x.*x, mu*x.*y],
           [mu*y.*x, mu*y.*y] ]
    
    #--------Derivative Functions---------------

    function Mdt(Q)     
        """
        Returns the first derivative of Quadrupole Moment Tensor with respect to time
        """
        
        dQdt = [ [py"grad"(Q[1][1], T), py"grad"(Q[1][2], T)], 
                 [py"grad"(Q[2][1], T), py"grad"(Q[2][2], T)], ]
        
        return dQdt
    end
    
    
    function Mdt2(Q) 
        """
        Returns the second derivative of Quadrupole Moment Tensor with respect to time
        """
        
        dQ2dt2 = Mdt(Mdt(Q))
        
        return dQ2dt2*(r0^2)
    end
    
    
    d2Qd2t = Mdt2(Q)

    
    
    h_plus =  (1.0/D) * (   d2Qd2t[1][1].*(cos.(phi_o).^2) 
                          + d2Qd2t[2][2].*(-cos.(phi_o).^2 .*cos.(theta_o).^2)  )
    
    return h_plus
    
end

function strainFFT(t, strain, f_bin)
    N = length(t)
    T = t[2] - t[1]
    
    h_plus_fft = SciPy.fft.fft(strain[1,:])/(2*pi*N)
    h_cross_fft = SciPy.fft.fft(strain[2,:])/(2*pi*N)
    xf = SciPy.fft.fftfreq(N, T)  #[1:Int(N ÷ 2)] 
    
    h_plus_fft = h_plus_fft[(xf .> f_bin[1]) .&& (xf .< f_bin[2])]
    h_cross_fft = h_cross_fft[(xf .> f_bin[1]) .&& (xf .< f_bin[2])]
    xf = xf[(xf .> f_bin[1]) .&& (xf .< f_bin[2])]

    return [xf, h_plus_fft, h_cross_fft]
end


function strainFFT2d(t, strain, f_bin)
    N = length(t)
    T = t[2] - t[1]
    
    h_plus_fft = SciPy.fft.fft(strain[:])/(2*pi*N)
    xf = SciPy.fft.fftfreq(N, T) #[1:Int(N ÷ 2)]
    
    h_plus_fft = h_plus_fft[(xf .> f_bin[1]) .&& (xf .< f_bin[2])]
    xf = xf[(xf .> f_bin[1]) .&& (xf .< f_bin[2])]
    
    return [xf , h_plus_fft]
end

function imaginary_interpd(f, h, f_new)
    #f and h are the computed strain and frequency (say 1,000,000 elements) 
    #f_new is the new frequency we'd like to have (say 10,000 elements)  scipy.interpolate import interp1d
    
    real_h = SciPy.interpolate.interp1d(f, real(h), kind="previous")
    im_h = SciPy.interpolate.interp1d(f, imag(h), kind="previous")
    
    return real_h(f_new) + im*im_h(f_new)
end

