import dc_leapfrog as lf
import nose.tools as nt
import numpy as np
import math as mt

def test_analytic():
    #
    #     Test case: Solve u' = -u + 1, u(0) = 0.
    #
    #     Compare the numerical solution given by
    #     solver_leapfrog() with the analytic
    #     solution u_e(t) = 1 - exp(-t).
    #
    
    def a(t):
        return 1.0
    def b(t):
        return 1.0
    #     Solve numerically
    u, t = lf.solver_leapfrog(I=0.0, a=a, b=b, \
                                  T=1.0, dt=1.0E-6, \
                                  theta=0.0)
    #     Analytic expression
    u_e = 1.0 - np.exp(-t)
    #     Maximum difference
    diff = np.abs(u-u_e).max()
    #     Test approximative equality
    nt.assert_almost_equal(diff, 0, delta=1E-12)
    
def test_linear_solution():
    #
    #     Method of manufactured solution: 
    #     Test solver_leapfrog() with a
    #     linear solution u_e(t) = ct + I.
    #
    
    #     Exact solution
    def linear_solution(t):
        return c*t + I
    
    def a(t):
        return mt.sin(t)      # Arbitrary function
    def b(t):
        return c + a(t)*linear_solution(t)
    
    c = 0.4
    I = 3.0
    #     Solve numerically
    u, t = lf.solver_leapfrog(I=I, a=a, b=b, \
                                  T=1.2, dt=1.0E-6, \
                                  theta=0.0)
    #     Analytic expression
    u_e = linear_solution(t)
    #     Maximum difference
    diff = np.abs(u-u_e).max()
    #     Test approximative equality
    nt.assert_almost_equal(diff, 0, delta=1E-10)   
