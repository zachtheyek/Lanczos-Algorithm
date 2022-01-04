import sympy as sympy
import numpy as np
import scipy as scipy
from scipy import integrate
from scipy.misc import derivative
import matplotlib.pyplot as plt 

np.set_printoptions(suppress=True)
np.set_printoptions(precision=5)



# This function will take derivatives for you...

bounds = 16
n=251.

x_low = -1*bounds
x_high = bounds


def DERIV(y,x_low,x_high,n):
    x = np.linspace(x_low,x_high,n)
    
    dy = np.zeros(y.shape,np.float)
    dy[0:-1]=np.diff(y)/np.diff(x)
    dy[-1]=(y[-1]-y[-2])/(x[-1]-x[-2])
#    dy[-100]=0
   
   
    return dy


beta = 0
alpha = 0
m=10
# T:Tridiagonal matrix
T = np.zeros( (m,m), float)
# 1. We define a symbol x, which is a real number.
# 2. We define a function phi(x) = exp(−x2/2)
# 3. We define a function lam_phi(x) = lambdify(x, phi(x), modules=['numpy'])
# 4. We define a vector x_val = np.linspace(x_low, x_high, int(n))
# 5. We define a vector phi = lam_phi(x_val)
# 6. We define a vector phi_sq = phi*phi
# 7. We define a vector int_phi_sq = np.trapz(phi_sq, x_val)
# 8. We define a scalar norm = np.sqrt(int_phi_sq)
# 9. We print the norm of the first vector
x = sympy.Symbol('x', real=True)
#Lets make a hamiltonian operator
sym_phi = sympy.exp((-x**2)/2)
lam_phi = sympy.lambdify(x, sym_phi, modules=['numpy'])

x_val = np.linspace(x_low, x_high, int(n))
phi = lam_phi(x_val)

phi_sq = phi*phi
int_phi_sq = np.trapz(phi_sq, x_val)

norm  =  np.sqrt(int_phi_sq)

print('The norm of the first vector is ',norm)

phi = phi/norm
#plt.scatter(x_val, phi, label= "stars", color= "0",  marker= ".", s=30) 
#plt.show()

d_phi = DERIV(phi,x_low,x_high, int(n))
dd_phi = DERIV(d_phi,x_low,x_high, int(n))
dd_phi = -0.5*dd_phi

sym_pot = 0.5*x**2
lam_pot = sympy.lambdify(x, sym_pot, modules=['numpy'])
pot = lam_pot(x_val)

H_phi = dd_phi + pot*phi
#plt.scatter(x_val, H_phi, label= "stars", color= "0",  marker= ".", s=30) 
#plt.show()

phi_H_phi = phi*H_phi

int_phi_H_phi = np.trapz(phi_H_phi, x_val)
alpha = int_phi_H_phi

print('a0 = ', alpha)
T[0,0 ] = alpha

 
w = H_phi - alpha*phi
w_2 = w*w
int_w_2 = np.trapz(w_2, x_val)
beta = np.sqrt(int_w_2) 
print('Beta1 is ', beta)

phi_prev = phi
for j in range(1,m):
    phi = w/beta
 
    d_phi = DERIV(phi,x_low,x_high, int(n))
    dd_phi = DERIV(d_phi,x_low,x_high, int(n))
    dd_phi = -1.*dd_phi

    H_phi = dd_phi + pot*phi

    w = H_phi - alpha*phi - beta*phi_prev

    phi_H_phi = phi*H_phi

    alpha = np.trapz(phi_H_phi, x_val)
    
    #Build the tridiagonal matrix
    T[j,j  ] = alpha 
    T[j-1,j] = beta 
    T[j,j-1] = beta

    print ('\n')
    print ('The tridiagonal matrix: \n', T)
    print ('\n')
    esT, vsT = np.linalg.eig( T )
    print ("The corresponding eigen-values are:", esT)  
    print ('\n')    
    w_2 = w*w
    int_w_2 = np.trapz(w_2, x_val)
    beta = np.sqrt(int_w_2) 
    
    phi_prev = phi  