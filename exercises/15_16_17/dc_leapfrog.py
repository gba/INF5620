from matplotlib.pyplot import *
from numpy import *
#import sys

def a(t):
    return exp(-t)

def b(t):
    return sin(t)

def solver_theta(I, a, b, T, dt, theta):
    """
    Solve the ODE u'(t) = -a*u(t) + b, u(0) = I,
    using the theta-rule.
    
    Time interval: t in (0, T].
    """
    
    dt = float(dt)
    N = int(round(T/dt))      # Number of time steps
    T = N*dt
    u = zeros(N+1)            # Init the vector v
    
    t = linspace(0, T, N+1)   # Time mesh vector
    
    u[0] = I                  # Initial condition
    for n in range(0, N): 
        u[n+1] = ((1-dt*(1.0-theta)*a)*u[n] \
                      + dt*(theta*b+(1-theta)*b)) \
                      /(1.0 + dt*theta*a)
    return u, t
    

def solver_leapfrog(I, a, b, T, dt, theta):
    """
    Solve the ODE u'(t) = -a(t)*u(t) + b(t), u(0) = I,
    using the Leapfrog scheme.
    
    Time interval: t in (0, T].
    
    The value u[1] is solved using the theta-rule.
         theta = 0: Forward Euler
         theta = 1: Backward Euler
         theta = 1/2: Crank-Nicolson
    """
    
    dt = float(dt)
    N = int(round(T/dt))      # Number of time steps
    T = N*dt
    u = zeros(N+1)            # Init the vector v
    
    t = linspace(0, T, N+1)   # Time mesh vector
    
    u[0] = I                  # Initial condition
    
    #     Solve u[1] with the theta-rule
    u01, t01 = solver_theta(I, a(t[0]), b(t[0]), dt, dt, theta) 
    u[1] = u01[1]
    
    #     Solve u[n+1], n = 1, ..., with the Leapfrog
    #     method
    for n in range(1, N):     # n = 1, 2, ..., N-1
        u[n+1] = u[n-1] + 2.0*dt*(-a(t[n])*u[n] + b(t[n]))
    return u, t

def solver_leapfrog_case1(T, dt, theta):
    """
    Solve the ODE u' = -u + 1, u(0) = 0, 
    t in (0, T], using the Leapfrog method.
    """
    
    def a1(t):
        return 1.0
    def b1(t):
        return 1.0
    #     Solve numerically
    u, t = solver_leapfrog(I=0.0, a=a1, b=b1, \
                                  T=T, dt=dt, \
                                  theta=theta)       
    return u, t

def analytical_case1(t):
    """
    Analytical solution of the ODE
    u' = -u + 1, u(0) = 0. 
    """
    return 1.0 - np.exp(-t)

def explore_leapfrog(I, a, b, T, dt, theta):
    """
    Do a general simulation with solver_leapfrog 
    and plot the result.
    """
    
    #     Numerical solution
    u, t = solver_leapfrog(I, a, b, T, dt, theta)
    
    #     Make plot
    figure()
    plot(t, u, 'r--')
    xlabel('t')
    ylabel('u(t)')
    legend(['dt = %g' % dt], loc='upper left')
    title('Du(t) = -a(t)*u(t) + b(t) with the Leapfrog method')
    savefig('leapfrog_dt%g.png' % dt)
    show()

def explore_leapfrog_case1(T, dt, theta, do_plot=True):
    """
    Simulate the ODE u' = -u + 1, u(0) = 0, 
    t in (0, T], with solver_leapfrog and plot 
    the numerical result together with the
    analytical solution u_e = 1 - exp(-t).
    """
    
    #     Numerical calculation
    u, t = solver_leapfrog_case1(T, dt, theta)
    #     Analytical solution
    u_e = analytical_case1(t)
    #     Local error 
    error = u_e - u
    #     Integrated error
    integr_error = sqrt(dt*sum(error**2))
    
    #     Make plot
    if do_plot:
        figure()
        t_e = linspace(0, T, 1001) 
        u_e = analytical_case1(t_e)
        plot(t, u, 'r--x')
        plot(t_e, u_e, 'b-')
        xlabel('t')
        ylabel('u(t)')
        legend(['Numerical', 'Analytical'])
        title('Du(t) = -u(t) + 1, u(0) = 0, with the Leapfrog method')
        savefig('leapfrog_case1.png')
        show()
    return integr_error


def convergence_rates(dt_array, integr_errors):
    """
    Calculate convergence rates r using time step
    values and corresponding integrated errors.
    """
    
    r = zeros(len(dt_array))    
    for i in range(1, len(dt_array)):
        #     Calculate the convergence rate
        r[i-1] = log(integr_errors[i-1]/integr_errors[i]) \
            /log(dt_array[i-1]/dt_array[i])
    return r

def define_command_line_options():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--I', type=float, default=1.0, \
                            help='The initial value u(0)', \
                            metavar='I')
    parser.add_argument('--T', type=float, default=1.0, \
                            help='Total time of simulation', \
                            metavar='T')
    parser.add_argument('--theta', type=float, default=0.0, \
                            help='Theta variable in the theta-rule', \
                            metavar='theta')
    parser.add_argument('--option', type=str, \
                            default='simulations', \
                            help='Option: simulations or test_conv', \
                            metavar='option')
    parser.add_argument('--dt', type=float, default=[1.0], \
                            help='List of time step lengths', \
                            metavar='dt', nargs='+', dest='dt_array')
    return parser

def read_command_line():
    parser = define_command_line_options()
    args = parser.parse_args()
    
    if len(sys.argv) < 10:
        print 'Usage: %s I T theta option dt1 dt2 dt3 ...' % \
            sys.argv[0]; sys.exit(1)
    else:
        print 'I={}, T={}, theta={}, option={}, dt_array={}'.format(
            args.I, args.T, args.theta, args.option, args.dt_array) 
    return args.I, args.T, args.theta, args.option, args.dt_array

def main():
    I, T, theta, option, dt_array = read_command_line()

    #I = 0.8
    #T = 10.0
    #theta = 0.0
    #dt_array = [1.0, 0.5, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001]
    
    #     Choose one of the following options
    #         'simulation': A general Leapfrog calculation
    #         'test_conv': Test convergence rate with a 
    #                      special case ODE
    #
    #option = 'simulations'

    if option == 'simulations':
        #     Do Leapfrog method simulations with different 
        #     time steps dt
        for dt in dt_array:
            explore_leapfrog(I=I, a=a, b=b, T=T, dt=dt, \
                                 theta=theta)
    
    if option == 'test_conv':
        #     Test the convergence rate for a special case ODE
        #
        #     Get error array
        integr_error = zeros(len(dt_array))
        for n in range(0, len(dt_array)):
            integr_error[n] = explore_leapfrog_case1(T, dt_array[n], \
                                                         theta, do_plot=True)
        #     Evaluate convergence rates
        r = convergence_rates(dt_array, integr_error)
        print r[0:len(r)-1]
    
if __name__ == '__main__':
    main()
