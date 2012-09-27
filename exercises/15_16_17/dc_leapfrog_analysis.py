from sympy import *

A = Symbol('A')
a = Symbol('a')
dt = Symbol('dt')

print ' '
print '-----------------------------------------------------------------'
print ' '
print """    We analyse the Leapfrog method applied to
    the ODE u' = -au, u(0) = I. 

    Assume that we can write u^n = A^n, where A 
    is a constant. Substituting u^n = A^n into the  
    numerical equation of the Leapfrog scheme, we 
    get the second order polynomial equation 
    A**2 + 2*a*dt - 1 = 0. """
print ' '

print """    Solve A**2 + 2*a*dt - 1 = 0 w.r.t. A: """
print ' '
A1, A2 = solve(A**2+2*a*dt*A-1, A)
print '         A1 = ', A1
print '         A2 = ', A2
print ' '

x = Symbol('x')

def A11():
    return -x + (x**2 + 1)**Rational(1,2)
def A22():
    return -x - (x**2 + 1)**Rational(1,2)

print '         A11 = ', A11()
print '         A22 = ', A22()
print ' '

#     Get the four first terms of the Taylor
#     expansions of A1 and A2
A1T = A1.series(a, 0, 5)
A2T = A2.series(a, 0, 5)

print '    Taylor expansions: '
print ' '
print '         A1T = ', A1T
print '         A2T = ', A2T
print ' '

print """    Since there are two roots A, the numerical solution 
    must be of the form u^n = C1*A1**n + C2*A2**n, where
    C1 and C2 are real constants determined by the initial
    conditions.

    From the Taylor expansions of A1 and A2 we find
    that abs(A2) > 1 for sufficiently small a*dt.
    
    Since abs(A2) > 1, abs(A2^n) will become very 
    large after a large number of time steps n,
    and the solution will therefore be unstable
    as n increases.

=================================================================
"""
