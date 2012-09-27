from numpy import *
import math
from matplotlib.pyplot import *
import vertical_motion as vm

def main():
    #     Air
    mu = 1.8E-5     # Units Pa*s
    density = 0.79   # Density of medium, kg/m^3

    #     Water
    #mu = 8.9E-4      # Units Pa*s
    #density = 1000    # Units kg/m^3
    
    #CD = 0.45        # Sphere
    #CD = 0.42        # Semi-sphere
    #CD = 1.05        # Cube
    #CD = 0.82        # Long cylinder
    #CD = 0.75        # Rocket
    CD = 1.0         # Man
    #CD = 1.3         # Flat plat perpend.
    #CD = 0.04        # Streamlined body
    
    diameter = 0.5    # Body diameter, m
    area = 0.9      # Area of body, m^2
    volume = 0.08    # Volume of body, m^3
    density_b = 1003.0   # Body density, kg/m^3

    #     Proportionality constant for
    #     Reynolds number
    Re_const = diameter*density/mu
    
    g = 9.81          # Gravitational const.

    v0 = 0.0          # Initial velocity, m/s
    
    dt = 0.1         # Time step, s
    T = 30.0           # Total time, s
    
    a_s = 3*math.pi*diameter*mu/(density_b*volume)
    a_q = 0.5*CD*density*area/(density_b*volume)
    b = g*(density/density_b - 1.0)
    #print 'a_s, a_q, b = ', a_s, a_q, b
    #print 'Re_const = ', Re_const

    #     Simulate a free fall of a parachute
    #     jumper
    vm.explore(v0, a_s, a_q, b, Re_const, T, dt)
    
if __name__ == '__main__':
    main()
