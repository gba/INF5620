from matplotlib.pyplot import *
import vertical_motion as vm
from numpy import *
import math


def calculate_forces(v0, mu, density_m, CD, diameter_b, \
                         area_b, volume_b, density_b, \
                         dt, T):
    """
    Calculate the gravity force, buoyancy force and
    the drag force of a falling body in a viscous 
    fluid.
    """
    
    #     Gravitational const. m/s^2
    g = 9.81         
    #     Proportionality constant for
    #     Reynolds number
    Re_const = diameter_b*density_m/mu
    
    a_s = 3*math.pi*diameter_b*mu/(density_b*volume_b)
    a_q = 0.5*CD*density_m*area_b/(density_b*volume_b)
    b = g*(density_m/density_b - 1.0)
    
    #     Numerical solution gives velocity as 
    #     a function of time.
    v, t = vm.solver(v0, a_s, a_q, b, Re_const, T, dt) 

    #     Initialize vectors
    Fg = zeros(len(v))
    Fb = zeros(len(v))
    Fd = zeros(len(v))

    #     Loop over time steps
    for n in range(0, len(v)):
        #     Evaluate Reynolds number
        Re = Re_const*v[n]   
        
        #     Gravity force
        Fg[n] = -density_b*volume_b*g
        #     Bouyancy force
        Fb[n] = density_m*g*volume_b
        
        #     Drag force
        if abs(Re) < 1:
            #     If Re < 1, use Stokes' drag force 
            Fd[n] = -3.0*math.pi*diameter_b*mu*v[n]
        else:
            #     If Re >= 1, use the quadratic
            #     drag force
            Fd[n] = -0.5*CD*density_m*area_b*abs(v[n])*v[n]

    
    return Fg, Fb, Fd, t        


def plot_forces(v0, mu, density_m, CD, diameter_b, \
                         area_b, volume_b, density_b, \
                         dt, T):
    """
    Calculate forces and plot the result.
    """

    #     Calculate the forces
    Fg, Fb, Fd, t = calculate_forces(v0, mu, density_m, \
                                      CD, diameter_b, area_b, \
                                      volume_b, density_b, dt, T)

    #     Plot the forces
    figure()
    plot(t, Fg, 'r--')
    plot(t, Fb, 'b-')
    plot(t, Fd, 'g--')
    legend(['Gravity', 'Bouyancy', 'Drag'])
    xlabel('Time (s)')
    ylabel('Force (N)')
    title('Forces working on a falling body')
    savefig('forces.png')
    show()

def main():
    #     Air
    mu = 1.8E-5     # Units Pa*s
    density_m = 0.79  # Density of medium, kg/m^3

    #     Water
    #mu = 8.9E-4      # Units Pa*s
    #density_m = 1000    # Units kg/m^3
    
    #CD = 0.45        # Sphere
    #CD = 0.42        # Semi-sphere
    #CD = 1.05        # Cube
    #CD = 0.82        # Long cylinder
    #CD = 0.75        # Rocket
    CD = 1.0        # Man
    #CD = 1.3         # Flat plat perpend.
    #CD = 0.04        # Streamlined body
    
    diameter_b = 0.5    # Body diameter, m
    area_b = 0.9      # Area of body, m^2
    volume_b = 0.08    # Volume of body, m^3
    density_b = 1003.0   # Body density, kg/m^3
    
    v0 = 0.0        # Initial velocity, m/s
    
    dt = 0.1        # Time step, s
    T = 30.0        # Total time, s
        
    plot_forces(v0, mu, density_m, CD, diameter_b, \
                         area_b, volume_b, density_b, \
                         dt, T)

if __name__ == '__main__':
    main()
