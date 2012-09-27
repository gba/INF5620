import dc_leapfrog as lf
from numpy import *
import sys

def setup_experiment(experiment):
    """
    Set up input parameters and functions a(t) and b(t)
    for a given experiment.
    """
    if experiment == 1:
        I = 0.8
        T = 10.0
        theta = 0.0
        dt_array = [1.0, 0.5, 0.2, 0.1, 0.05, 0.01, \
                        0.005, 0.001]
        def a(t):
            return 1.0
        def b(t):
            return 1.0

    if experiment == 2:
        I = 0.8
        T = 5.0
        theta = 0.0
        dt_array = [1.0, 0.5, 0.2, 0.1, 0.05, 0.01, \
                        0.005, 0.001]
        def a(t):
            return t
        def b(t):
            return sin(t)

    if experiment == 3:
        I = 0.8
        T = 20.0
        theta = 0.0
        dt_array = [1.0, 0.5, 0.2, 0.1, 0.05, 0.01]
        def a(t):
            return exp(-t)
        def b(t):
            return 1.0/(t**2+1.0) 

    return I, T, theta, dt_array, a, b

def define_command_line_options():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--exp', type=int, default=1, \
                            help='Choose experiment', \
                            metavar='exp')
    return parser

def read_command_line():
    parser = define_command_line_options()
    args = parser.parse_args()
    
    if len(sys.argv) < 2:
         print 'Usage: %s exp' % \
            sys.argv[0]; sys.exit(1)
    else:
        print 'exp={}'.format(args.exp)
    return args.exp
        
def main():
    experiment = read_command_line()
    #experiment = 3
    
    I, T, theta, dt_array, a, b = setup_experiment(experiment)
    
    #     Do Leapfrog method simulations with different 
    #     time steps dt
    for dt in dt_array:
        lf.explore_leapfrog(I=I, a=a, b=b, T=T, dt=dt, \
                             theta=theta)
        
if __name__ == '__main__':
    main()
