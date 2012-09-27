import vertical_motion as vm
import nose.tools as nt
import numpy as np

def test_heavy_slow():
    #
    # Assumption: density of body is very large,
    #            density of medium is very small.
    #            This should give linear velocity
    #            according to v(t) = v0 - g*t.
    #
    g = 9.81
    v0 = 0.0
    #     Solve numerically
    v, t = vm.solver(v0, a_s=0.0, a_q=0.0, b=-g, \
                      Re_const=0.0, T=2.0, dt=0.001)
    
    v_a = np.zeros(len(t))
    #     The analytic velocity
    for n in range(0, len(t)):
        v_a[n] = v0 - g*t[n]
        
    diff = np.abs(v-v_a).max()
    nt.assert_almost_equal(diff, 0, delta=1E-12)

def test_terminal_velocity_q():
    #
    # In the quadratic drag regime the final
    # velociy should be v = -sqrt(abs(a_q)/b).
    # 
    a_s = 1E-6
    a_q = 0.0044
    b = -9.80
    Re_const = 21944.0
    T = 100.0
    dt = 0.1
    v0 = 0.0
    v, t = vm.solver(v0, a_s, a_q, b, Re_const, T, dt)

    vmax = -np.sqrt(abs(b)/a_q)
    print 'vmax = ', vmax
    diff = abs(vmax - v[len(v)-1])
    nt.assert_almost_equal(diff, 0, delta=1E-12)
