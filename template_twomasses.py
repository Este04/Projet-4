#
# LEPL1504 - -    MISSION  -   Modelling a two-mass system
#
# @date 2025
# @author Robotran team
# 
# Universite catholique de Louvain


# import useful modules and functions
import matplotlib.pyplot as plt
from math import sin, cos, pi
import numpy as np
from scipy.integrate import solve_ivp


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
class MBSData:
    """ Class to define the parameters of the the double mass-spring model
    and of the problem to solve.
     
    It contains the following fields (in SI units): 
    
    general parameters:
    -------------------
    g:     Gravity
    
    masses:
    -------
    m1:    Unsprung mass
    m2:    Sprung mass
    
    parameters of the tyre:
    -----------------------
    k01:   Stiffness
    d1:    Damping
    z01:   Neutral length
    
    parameters of the suspension:
    -----------------------------
    k02:   Stiffness
    d2:    Damping
    z02:   Neutral length   
    
    parameter of the external force:
    --------------------------------
    Fmax:  Semi-amplitude of the force
    Zmax:  Semi-amplitude of the displacement
    f0:    Frequency at start
    t0:    Time for specifying frequency at start
    f1:    Frequency at end
    t1:    Time for specifying frequency at end
    
    equilibrium positions:
    ----------------------
    q1:    Equilibrium position coordinate of m1
    q2:    Equilibrium position coordinate of m2
    """
    
    def __init__(self):
    # Write your code here
        self.g = 9.81
        self.m1 = 35
        self.m2 = 315
        self.k01 = 190e3
        self.d1 = 107
        self.z01 = 0.375
        self.k02 = 37e3
        self.d2 = 4000
        self.z02 = 0.8
        self.Fmax = 1e4
        self.Zmax = 0.1
        self.f0 = 1
        self.t0 = 0
        self.f1 = 10
        self.t1 = 10
        self.q1 = 0.357445263
        self.q2 = 0.716482432
            
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
def sweep(t, t0, f0, t1, f1, Fmax):
    """ Compute the value of a force sweep function at the given time.
    The sweep function has a sinusoidal shape with constant amplitude 
    and a varying frequency. This function enables to consider linear
    variation of the frequency between f0 and f1 defined at time t0 and t1.

	:param t: the time instant when to compute the function.
	:param t0: the time instant when to specify f0.
	:param f0: the frequency at time t0.
	:param t1: the time instant when to specify f1.
	:param f1: the frequency at time t1.
	:param Fmax: the semi-amplitude of the function.
		
	:return Fext: the value of the sweep function.
    """
    return Fmax*sin(2*pi*t*(f0+((f1-f0)/(t1-t0))*(t/2)))

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *    
def compute_derivatives(t, y, data):
    """ Compute the derivatives yd for a given state y of the system.
        The derivatives are computed at the given time t with
        the parameters values in the given data structure.
        
        It is assumed that the state vector y contains the following states:
          y = [q1, q2, qd1, qd2] with:
             - q1: the mass 1 position
             - q2: the mass 2 position 
             - qd1: the mass 1 velocity
             - qd2: the mass 2 velocity 

        :param  t: the time instant when to compute the derivatives.
        :param  y: the numpy array containing the states 
        :return: yd a numpy array containing the states derivatives  yd = [qd1, qd2, qdd1, qdd2]
        :param data: the MBSData object containing the parameters of the model
    """                 
    # Write your code here
    # TODO   
    # sweep function should be called here: sweep(t, data.t0, data.f0, data.t1, data.f1, data.Fmax)
    Fext = sweep(t, data.t0, data.f0, data.t1, data.f1, data.Fmax)
    yd = np.zeros(4)
    yd[0] = y[2]
    yd[1] = y[3]
    
    A = np.array([[data.m1 + data.m2, data.m2], [data.m2, data.m2]])
    C = np.array([[data.d1, 0], [0, data.d2]])
    K = np.array([[data.k01, 0], [0, data.k02]])
    D = np.array([-data.g*(data.m1+data.m2), -data.g*data.m2 - Fext])

    B = D - K@y[0:2] - C@y[2:4]
    yd[2:4] = np.linalg.solve(A, B)

    return yd



# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
def compute_dynamic_response(data):
    """  Compute the time evolution of the dynamic response of the two-mass system
         for the given data. Initial and final time are determined
         by the t0 and t1 parameter of the parameter data structure.
         Results are saved to three text files named dirdyn_q.res, dirdyn_qd.res and dirdyn_qdd.res
 
        Time evolution is computed using an time integrator (typically Runge-Kutta).
 
       :param data: the MBSData object containing the parameters of the model
     """
    # ### Runge Kutta ###   should be called via solve_ivp()
    # to pass the MBSData object to compute_derivative function in solve_ivp, you may use lambda mechanism:
    #
    #    fprime = lambda t,y: compute_derivatives(t, y, data)
    #
    # fprime can be viewed as a function that takes two arguments: t, y
    # this fprime function can be provided to solve_ivp
    # Note that you can change the tolerances with rtol and atol options (see online solve_iv doc)
    #
    # Write some code here
    # TODO

    fprime = lambda t,y: compute_derivatives(t, y, data)
    t_span = [data.t0, data.t1]
    y0 = [data.q1, data.q2, 0, 0]
    sol = solve_ivp(fprime, t_span, y0, t_eval=np.linspace(data.t0, data.t1, 1000), method='RK45')
    q1 = sol.y[0]
    q2 = sol.y[1]  
    plt.plot(sol.t, q1, label='q1')
    plt.plot(sol.t, q2, label='q2')
    plt.show()



  


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# Main function

if __name__ == '__main__':
    mbs_data = MBSData()
    compute_dynamic_response(mbs_data)  
