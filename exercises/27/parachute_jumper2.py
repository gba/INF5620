from matplotlib.pyplot import *
from numpy import *
import math


class Problem:
    #                 height0=5000.0, pressure0=55000.0, \
    #                 velocity0=0.001, mol_mass=0.029, \
    #                 gravity_const=9.807, gas_const=8.314, \
    #                 temp_const=-0.0065, temp0=288, viscosity=1.8E-5, \
    #                 drag_coeff1=1.0, drag_coeff2=1.3, \
    #                 diameter_b=0.5, area_b1=0.9, area_b2=17.65, \
    #                 volume_b=0.08, density_b=1003.0, time_max=250.0, \
    #                 time_release=70.0
    def __init__(self, height0=5000.0, pressure0=55000.0, \
                     velocity0=0.001, mol_mass=0.029, \
                     gravity_const=9.807, gas_const=8.314, \
                     temp_const=-0.0065, temp0=288.0, viscosity=1.8E-5, \
                     drag_coeff1=1.0, drag_coeff2=1.3, \
                     diameter_b=0.5, area_b1=0.9, area_b2=17.65, \
                     volume_b=0.08, density_b=1003.0, time_max=250.0, \
                     time_release=70.0):
        """
        height0 = initial altitude (m)
        pressure0 = initial pressure (Pa)
        velocity0 = initial velocity (m/s)
        
        mol_mass = molar mass of Earth's air (kg/mol)
        gravity_const = acceleration of gravity (m/s^2)
        gas_const = the universal gas constant (Nm/(mol K)
        temp_const = constant in temperature model (K/m)     
        
        temp0 = temperature at sea level (K)
        viscosity = dynamic viscosity of fluid (Pa s)
        
        drag_coeff1 = drag coefficient of body without
                      parachute
        drag_coeff2 = drag coefficient of body with
                      parachute
        diameter_b = body diameter (m)
        area_b1 = body area without parachute (m) 
        area_b2 = body area with parachute (m)
        volume_b = body volume (m^3)
        density_b = body density (kg/m^3)
        
        time_max = total time (s)
        time_release = time when parachute is released (s)
        
        Re_const = proportionally constant
        for Reynolds number
        """
        self.height0 = height0
        self.pressure0 = pressure0
        self.velocity0 = velocity0
        self.mol_mass = mol_mass
        self.gravity_const = gravity_const
        self.gas_const = gas_const
        self.temp_const = temp_const
        self.temp0 = temp0
        self.viscosity = viscosity
        self.drag_coeff1 = drag_coeff1
        self.drag_coeff2 = drag_coeff2
        self.diameter_b = diameter_b
        self.area_b1 = area_b1
        self.area_b2 = area_b2
        self.volume_b = volume_b
        self.density_b = density_b
        self.time_max = time_max
        self.time_release = time_release
        
        
class Solver:
    def __init__(self, problem, time_step=0.01, theta=0.5):
        self.problem = problem
        self.time_step = float(time_step)
        self.theta = float(theta)
        
    def solve(self):
        self.z, self.p, self.v, self.t = \
            self.solver_theta(self.time_step, \
                                  self.theta)
        
      
    def f1(self, velocity_n):
        """
        Function f^1 in the ODE z'(t) = f^1(z(t), p(t), v(t)),
        where z = height, p = pressure, v = velocity, and
        t = time.
        """
        return velocity_n
    
    def f2(self, height_n, pressure_n, velocity_n):
        """
        Function f^2 in the ODE p'(t) = f^2(z(t), p(t), v(t)),
        where z = height, p = pressure, v = velocity, and
        t = time.
        """
        problem = self.problem
        nominator = -problem.mol_mass*problem.gravity_const \
            *pressure_n
        denominator = problem.gas_const \
            *(problem.temp0 +  problem.temp_const*height_n) \
            *velocity_n
        
        funct2 = nominator/denominator
        #funct2 = 0.0
        return funct2
    
    
    def f3(self, height_n, pressure_n, velocity_n, drag_coeff, \
               area_b):
        """
        Function f^3 in the ODE v'(t) = f^3(z(t), p(t), v(t)),
        where z = height, p = pressure, v = velocity, and
        t = time.
        """
        problem = self.problem
        temperature = problem.temp0 + problem.temp_const*height_n
        Reynolds = problem.mol_mass*problem.diameter_b \
            *pressure_n*velocity_n \
            /(problem.gas_const*problem.viscosity*temperature)
        
        if abs(Reynolds) < 1:
            #     If Reynolds number < 1, use Stokes' drag force 
            const1 = -3.0*math.pi*problem.diameter_b \
                *problem.viscosity \
                /(problem.density_b*problem.volume_b)
            const2 = problem.gravity_const*problem.mol_mass \
                /(problem.gas_const*problem.density_b* \
                      temperature)
            const3 = -problem.gravity_const
            
            funct3 = const1*velocity_n + const2*pressure_n + const3
        else:
            #     If Reynolds number >= 1, use the quadratic
            #     drag force
            const1 = -0.5*problem.mol_mass \
                /(problem.gas_const*temperature*problem.density_b \
                      *problem.volume_b)
            const2 = problem.gravity_const*problem.mol_mass \
                /(problem.gas_const*temperature*problem.density_b)
            const3 = -problem.gravity_const
            
            funct3 = const1*drag_coeff*area_b*abs(velocity_n) \
                *velocity_n*pressure_n + const2*pressure_n \
                + const3    
        
        return funct3
    
    def get_drag_coeff(self, time):
        """     
        Get drag coefficient before and after
        the parachute has been released
        """
        if time < self.problem.time_release:
            drag_coeff = self.problem.drag_coeff1
        else:
            drag_coeff = self.problem.drag_coeff2
        return drag_coeff
  
    def get_area_b(self, time):
        """
        Get body area before and after
        the parachute has been released
        """
        if time < self.problem.time_release:
            area_b = self.problem.area_b1
        else:
            area_b = self.problem.area_b2
        return area_b
    
    def solver_theta(self, time_step, theta):
        """
        Solve the ODE system
        
        z'(t) = f^1(z(t), p(t), v(t))
        p'(t) = f^2(z(t), p(t), v(t))
        v'(t) = f^3(z(t), p(t), v(t)),
        
        using the theta-rule. Here t is in the interval 
        (0, time_max] and the time step is time_step.
        
        """
    
        #     Number of time steps
        n_steps = int(round(self.problem.time_max/time_step))  
        time_max = n_steps*time_step
        #     Initialize vectors z, p, and v
        z_ground = 0.0001
        v_ground = 0.0001
        z = z_ground*ones(n_steps+1)
        p = zeros(n_steps+1)
        v = zeros(n_steps+1)
        z_old = zeros(n_steps+1)
        p_old = zeros(n_steps+1)
        v_old = zeros(n_steps+1)
        area_b = zeros(2)
        drag_coeff = zeros(2)
        
        #     Time mesh vector
        t = linspace(0, time_max, n_steps+1)   
        #     Initial conditions
        z[0] = self.problem.height0
        p[0] = self.problem.pressure0
        v[0] = self.problem.velocity0
        tolerance = 1E-3
        
        #     n = 0, 1, ..., N-1
        for n in range(0, n_steps):     
            #print 'n = ', n
            #     Initial guess
            z_old[n+1] = z[n]
            p_old[n+1] = p[n]
            v_old[n+1] = v[n]
            #     Get body area and drag coefficient
            area_b[0] = self.get_area_b(t[n]) 
            area_b[1] = self.get_area_b(t[n+1])
            drag_coeff[0] = self.get_drag_coeff(t[n])
            drag_coeff[1] = self.get_drag_coeff(t[n+1])
            
            difference = 100.0
            #     Selfconsistency loop
            while difference > tolerance:
                #     The theta rule
                z[n+1] = z[n] + time_step*theta*self.f1(v_old[n+1]) \
                    + time_step*(1.0 - theta)*self.f1(v[n])
                p[n+1] = p[n] + time_step*theta \
                    *self.f2(z_old[n+1], p_old[n+1], v_old[n+1]) \
                    + time_step*(1.0 - theta) \
                    *self.f2(z[n], p[n], v[n])
                v[n+1] = v[n] + time_step*theta \
                    *self.f3(z_old[n+1], p_old[n+1], v_old[n+1], \
                                 drag_coeff[1], area_b[1]) \
                                 + time_step*(1.0 - theta) \
                                 *self.f3(z[n], p[n], v[n], \
                                              drag_coeff[0], area_b[0])
                
                #     Relative difference from previous iteration
                diff_z = abs((z[n+1] - z_old[n+1])/z_old[n+1])
                diff_p = abs((p[n+1] - p_old[n+1])/p_old[n+1])
                diff_v = abs((v[n+1] - v_old[n+1])/v_old[n+1])
                difference = max(diff_z, diff_p, diff_v)
                #print 'diff = ', difference
                
                if z[n+1] < 0.0:
                    z[n+1] = z_ground
                    v[n+1] = v_ground
                    difference = -1.0
                #print 'diff = ', difference
                z_old[n+1] = z[n+1]
                p_old[n+1] = p[n+1]
                v_old[n+1] = v[n+1]
                
        return z, p, v, t

class Visualizer:
    def __init__(self, problem, solver):
        self.problem = problem
        self.solver = solver
        
    def plot_velocity(self):
        """
        Plot the velocity as a function of time.
        """
        import scitools.std
        plt = scitools.std
        
        plt.plot(self.solver.t, self.solver.v, '--')
        plt.xlabel('time (s)')
        plt.ylabel('velocity (m/s)')
        plt.title('Velocity of parachute jumper')
        plt.savefig('velocity.png')
        return plt
    
    def plot_altitude(self):
        """
        Plot the altitude as a function of time.
        """
        
        import scitools.std
        plt = scitools.std
        
        plt.plot(self.solver.t, self.solver.z, '--')
        plt.xlabel('time (s)')
        plt.ylabel('altitude (m)')
        plt.title('Altitude of parachute jumper')
        plt.savefig('altitude.png')
        return plt

def main():
    problem = Problem()
    solver = Solver(problem)
    visu = Visualizer(problem, solver)
    
    #     Solve a system of ODE's to model
    #     a parachute jump
    solver.solve()
    #     Plot the velocity and altitude
    #     as functions of time
    plot_v = visu.plot_velocity()
    plot_v.show()
    plot_a = visu.plot_altitude()
    plot_a.show()
    
 
if __name__ == '__main__':
    main()
    
