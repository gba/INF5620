from numpy import *
import math
from matplotlib.pyplot import *

def solver(v0, a_s, a_q, b, Re_const, T, dt):
    """
    This is a model for vertical motion in a fluid.

    The function solves the ODE v' = a_s*v + b, for 
    Re < 1 and v' = a_q|v|v + b, for Re >= 1.
    
    v(0) = I, t in (0, T], time step dt
    
    Re_const = constant in expression of 
               Reynolds number Re

    The ODE's are solved by the Crank-Nicolson 
    scheme.
    
    """
    
    dt = float(dt) 
    N = int(round(T/dt))      # Number of time steps
    T = N*dt
    v = zeros(N+1)            # Init the vector v

    t = linspace(0, T, N+1)   # Time mesh vector
    
    v[0] = v0                 # Initial condition
    for n in range(0, N):     # n = 0, 1, ..., N-1
        Re = Re_const*v[n]    # Evaluate Re
        
        if abs(Re) < 1:
            #print 'Re = ', abs(Re)
            #     If Re < 1, use Stokes' drag force 
            v[n+1] = ((1.0 - 0.5*a_s*dt)*v[n] + b*dt) \
                /(1.0 + 0.5*a_s*dt)
        else:
            #print 'Re = ', abs(Re)
            #     If Re >= 1, use the quadratic
            #     drag force
            v[n+1] = (v[n] + dt*b)/(1.0 + dt*a_q*abs(v[n]))
    return v, t


def explore(v0, a_s, a_q, b, Re_const, T, dt):
    """
    Do a simulation with solver and plot the result.
    """
    
    #     Numerical solution
    v, t = solver(v0, a_s, a_q, b, Re_const, T, dt)      
    
    #     Make plot
    figure()
    plot(t, v, 'r--')
    xlabel('t')
    ylabel('v')
    title('Vertical motion in viscous fluid')
    savefig('viscous.png')
    show()

def main():
    #     Air
    mu = 1.8E-5     # Units Pa*s
    density_m = 0.79   # Density of medium, kg/m^3

    #     Water
    #mu = 8.9E-4      # Units Pa*s
    #density_m = 1000    # Units kg/m^3
    
    #CD = 0.45        # Sphere
    #CD = 0.42        # Semi-sphere
    #CD = 1.05        # Cube
    #CD = 0.82        # Long cylinder
    #CD = 0.75        # Rocket
    CD = 1.0         # Man
    #CD = 1.3         # Flat plat perpend.
    #CD = 0.04        # Streamlined body
    
    diameter_b = 0.5    # Body diameter, m
    area_b = 0.9      # Area of body, m^2
    volume_b = 0.08    # Volume of body, m^3
    density_b = 1003.0   # Body density, kg/m^3

    #     Proportionality constant for
    #     Reynolds number
    Re_const = diameter_b*density_m/mu
    
    g = 9.81          # Gravitational const.

    v0 = 0.0          # Initial velocity, m/s
    
    dt = 0.1         # Time step, s
    T = 30.0           # Total time, s
    
    a_s = 3*math.pi*diameter_b*mu/(density_b*volume_b)
    a_q = 0.5*CD*density_m*area_b/(density_b*volume_b)
    b = g*(density_m/density_b - 1.0)
    #print 'a_s, a_q, b = ', a_s, a_q, b
    #print 'Re_const = ', Re_const

    explore(v0, a_s, a_q, b, Re_const, T, dt)
    
if __name__ == '__main__':
    main()
