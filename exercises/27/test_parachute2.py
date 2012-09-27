import parachute_jumper2 as pj
import nose.tools as nt
import numpy as np

def test_linear():
    #
    # Assumption: density of body is very large,
    #            density of medium is very small.
    #            This should give linear velocity
    #            according to v(t) = v0 - g*t.
    #
    
    problem = pj.Problem()              # Create problem object
    problem.density_b = 1E10            # High density of body
    problem.mol_mass = 1E-10            # Low density of medium
    problem.time_max = 10.0             # Total time (s)
    
    problem.height0 = 5000.0
    problem.pressure0 = 55000.0
    problem.velocity0 = 0.001
    problem.gravity_const = 9.807
    problem.gas_const = 8.314
    problem.temp_const = -0.0065
    problem.temp0 = 288
    problem.viscosity = 1.8E-5
    problem.drag_coeff1 = 1.0
    problem.drag_coeff2 = 1.3
    problem.diameter_b = 0.5
    problem.area_b1 = 0.9
    problem.area_b2 = 3.14
    problem.volume_b = 0.08
    problem.time_release = 60.0

    solver = pj.Solver(problem)         # Create solver object
    solver.time_step = 0.001
    solver.theta = 0.5

    #     Solve a system of ODE's to model
    #     a parachute jump
    solver.solve()
    
    v_a = np.zeros(len(solver.t))
    #     The analytic velocity
    for n in range(0, len(solver.t)):
        v_a[n] = problem.velocity0 - problem.gravity_const*solver.t[n]
    #     Difference
    diff = np.abs(solver.v-v_a).max()
    nt.assert_almost_equal(diff, 0, delta=1E-10)
    
def test_terminal_velocity_q():
    #
    # In the quadratic drag regime the final
    # velociy should be v = -sqrt(abs(a_q)/b).
    # Observe that this holds only when the
    # density is constant.
    #
    
    problem = pj.Problem()              # Create problem object
    problem.time_max = 100.0            # Total time (s)
    problem.mol_mass = 0.029            # Zero molar mass => 
                                        # zero density of medium
    
    problem.height0 = 5000.0
    problem.pressure0 = 55000.0
    problem.velocity0 = 0.001
    problem.gravity_const = 9.807
    problem.gas_const = 8.314
    problem.temp_const = 0.0 #-0.0065
    problem.temp0 = 288.0
    problem.viscosity = 1.8E-5
    problem.drag_coeff1 = 1.0
    problem.drag_coeff2 = 1.3
    problem.diameter_b = 0.5
    problem.area_b1 = 0.9
    problem.area_b2 = 3.14
    problem.volume_b = 0.08
    problem.density_b = 1003.0
    problem.time_release = 0.0001
    
    solver = pj.Solver(problem)         # Create solver object
    solver.time_step = 0.001
    solver.theta = 0.5
    
    #     Solve a system of ODE's to model
    #     a parachute jump
    solver.solve()
    
    n_max = len(solver.z)-1
    temperature = problem.temp0 + problem.temp_const \
        *solver.z[n_max]
    density_m = solver.p[n_max]*problem.mol_mass \
        /(problem.gas_const*temperature)
    #     b and a_q are coefficients of the quadratic 
    #     drag force model
    b = problem.gravity_const*(density_m/problem.density_b - 1.0)
    a_q = 0.5*problem.drag_coeff2*density_m*problem.area_b2 \
        /(problem.density_b*problem.volume_b)
    
    print ' '
    print ' ################################################### '
    print ' '
    print '     NB! The function f2 must be set to 0.0 to get'
    print '     the desired precision'
    print ' '
    print ' ################################################### '
    print ' '
    
    #     Analytic value
    vmax = -np.sqrt(abs(b)/a_q)
    #     Difference		
    diff = abs(vmax - solver.v[n_max])
    nt.assert_almost_equal(diff, 0, delta=1E-11)
    
