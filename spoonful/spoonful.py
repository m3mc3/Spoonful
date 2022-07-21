import matplotlib.pyplot as plt
import numpy as np
import time

# Welcome to the Spoonful Library - an attempt/passion project to combine all
# of the most essential tools in scientific methods!

# This package contains the beginnings of a journey to efficiency in

# - Numerical Integration
# - Root Finding
# - Differential Equations

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# General Integrator Class

class Integrator(object):
    
    """A general class of numerical integration techniques
    
    
    Parameters
    __________
    
    a : float
    
        start of interval
    
    b : float
    
        end of interval
    
    n : int
        
        number of evaluation points
        
    f : function
        
        (single-variable) function to be integrated
    
    
    Returns
    _______
    
    float
        
        numerical result for the area under the function over [a, b]
    """
    
    def __init__(self, a, b, n):
        self.a = a
        self.b = b
        self.n = n
        self.dx = (b - a)/n
    
    def Simpson(self, f):
        
        """Numerical integration via Simpson's 1/3 Rule
        
        Parameters
        __________
        
        f : function
        
            (single-variable) function to be integrated
        
        
        Returns
        _______
        
        float
            
            numerical result for the area under the function over [a, b]
        """
        
        start = time.perf_counter()
        # Evaluation Points
        xvals = np.linspace(self.a + self.dx, self.b - self.dx, self.n)
        # Initial Sum
        S = (self.dx/3)*(f(self.a) + f(self.b)) 
        # Integrated Sum
        for j in xvals[::2]:
            S += (4*self.dx/3)*f(j)
        for k in xvals[1::2]:
            S += (2*self.dx/3)*f(k)
        
        end = time.perf_counter()
        
        print("Finished in {} second(s)".format(end - start))
        return S

    def Trapezoid(self, f):
        
        """Numerical integration via Trapezoidal Sums
        
        Parameters
        __________
        
        f : function
        
            (single-variable) function to be integrated
        
        
        Returns
        _______
        
        float
            
            numerical result for the area under the function over [a, b]
        """
        
        start = time.perf_counter()
        # Evaluation Points 
        # Use Generators since we don't need list operations
        xvals = np.linspace(self.a + self.dx, self.b - self.dx, self.n)
        # Initial Sum
        S = (self.dx/2)*(f(self.a) + f(self.b))
        # Integrated Sum
        for xval in xvals:
            S += (self.dx/2)*f(xval)
        
        end = time.perf_counter()
        
        print("Finished in {} second(s)".format(end - start))
        return S


# General ODE Integration Class

class ODE(object):
    
    """A general class of numerical techniques to integrating ODEs
    
    
    Parameters
    __________
    
        unspecified
    
    
    Results
    _______
    
        unspecified
    """

    def __init__(self, t_0, t_f, n, IC):
        self.t_0 = t_0 
        self.t_f = t_f 
        self.IC = IC 
        self.t = np.linspace(t_0, t_f, n)
        self.h = self.t[1] - self.t[0]
        
    """Runge Kutta Method of 4th Order
    
    
    Parameters
    __________
    
    f : function
        
        dy/dt = f(y, t)
    
    Optional Parameters
    ___________________
    
    plot = True
    
        automatically returns matplotlib graph of the solution over [t_0, t_f]
    
    label = ''
    
        title of auto-generated plot
    
    Returns
    _______
    
    list
    
        list of function evaluations corresponding to the solution of the ODE.
        
    """
    
    def RK4(self, f, plot=False):
        
        # Initialize list-solution
        self.y = np.zeros(len(self.t))
        
        for i in range(len(self.t) - 1):
            K1 = f( self.t[i], self.y[i] )
            K2 = f( self.t[i] + self.h/2, self.y[i] + self.h*K1/2 )
            K3 = f( self.t[i] + self.h/2, self.y[i] + self.h*K2/2 )
            K4 = f( self.t[i] + self.h, self.y[i] + self.h*K3 )
            
            self.y[i + 1] = self.y[i] + (self.h/6)*(K1 + 2*K2 + 2*K3 + K4)
        
        if plot==True:
            plt.plot(self.t, self.y)
            plt.label("{}".format(self.label))
            plt.grid()
            plt.show()
        
        return self.y

def f(y):
    return(1 - 3*y)/(1 + y**2)

print(ODE(0, 5, 1000, 1).RK4(f))















