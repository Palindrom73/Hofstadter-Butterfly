from __future__ import unicode_literals

import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse 
from scipy.sparse.linalg import eigsh
import matplotlib as mpl

from matplotlib.ticker import MultipleLocator

#plot package settings
mpl.use('GTK3Cairo')
plt.rc('text', usetex = True)
plt.rc('font', family ='serif')


#mathematical requirements

#my trigonometric interpolation polynomial for the upper bound of the spectrum
def p(theta,u,v,w):
    return u+v+np.sqrt(u**2+v**2)+(u+v-2*w-np.sqrt(u**2+v**2))*np.cos(2*np.pi*theta)
#rotation angle on non-commutative torus/ magnetic flux in QHE
def rho(theta):
    return np.exp(2*np.pi*1j*theta)

# amtrix U from non-commutative torus
def MU(z1,rho,q):
    D1 = rho * np.ones(q)
    D2 = np.ones(q).cumsum() -1
    D = D1 ** D2
    M = z1 * np.diag(D)
    return M

# amtrix V from non-commutative torus
def MV(z2,q):
    M = np.diag(np.ones(q-1), k=-1)
    M[0,-1] = 1
    return z2 * M

#studied matrix or operator in the non-commutative torus
def MH(z1,z2,rho,q,u,v):
    A1 = MU(z1,rho,q)
    A2 = MV(z2,q)
    return u*(A1 + A1.conjugate()) + v*(A2 + ((A2.conjugate()).T))

#coprime pairs p,q for theta=p/q up to maximal matrix size n>=q
def get_fractions(n):
    coprime = []
    for q in range (1,n):
        for p in range(0,q):
            if np.gcd(q,p) == 1:
                coprime.append([p,q])
    coprime.append([1,1])
    return coprime

def sparseM(M):
	sparseM = sparse.csr_matrix(M)
	return sparseM

#calculate the maximum norm in resp. to the largest eigenvalue in absolute value of our operator, note that the operators we observing have real eigenvalues obtained by C*-algebra theory!
def maximumnorm(q_max,u , v, w):
    z1, z2 = 1, 1
    Thetas = []
    Lambdas = []
    
    fraction = get_fractions(q_max)
    for f in fraction:
        p = f[0]
        q = f[1]
        theta = p/q
        Rho = rho(theta)
        M= MH(z1, z2, Rho,q,u,v)
        
        if q<=2:
            ew = np.linalg.norm(M, ord=2)-2 * w * np.cos(2 * np.pi * theta)
        else:
            ew = eigsh(M, k=1, which='LA', return_eigenvectors=False)[0] - 2 * w * np.cos(2 * np.pi * theta)
        Thetas.append(theta)
        Lambdas.append(ew)
    Thetas= np.array(Thetas)
    Lambdas= np.array(Lambdas)    
    return Thetas, Lambdas
    
#main function for the plots with given parameters       
def main(q_max, u, v, w):
    Thetas, Lambdas = maximumnorm(q_max,u , v, w)
    x = np.linspace(0,1,1000)
    y = p(x, u, v, w)
    
    fig= plt.figure()
    
    ax = fig.add_subplot(111)
    
    ax.plot(0.1e6,'o',markersize=2, color= 'red', label=r'$\lambda_{\mathrm{max}}$')
    ax.plot(Thetas,Lambdas,'.', markersize=0.1, color='red')
    ax.plot(x, y, '-', color='black', label=r'$p(\theta)$')
    
    ax.set_xlim(0,1)
    ax.set_ylim(0.5,1.1)
    ax.set_xlabel(r'$\theta$', fontsize=20)
    ax.set_ylabel(r'$\lambda_{\mathrm{max}}$', fontsize=20)
    ax.set_title('Upper bound of the Spectral Radius', fontsize=22, pad=30, x= 0.6)
    
    ax.tick_params(axis = 'both', which = 'major', width = 1, length= 6, direction = 'in', labelsize =16, zorder = 5)
    ax.tick_params(axis = 'both', which = 'minor', width = 1, length= 2, direction = 'in', labelsize =4, zorder = 5)
    ax.tick_params(which = 'both', top = True, labeltop = True, right = True, labelright = True, zorder = 5)
    
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(MultipleLocator(0.04))
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.04))
    
    ax.legend(fontsize = 18, loc='upper right')
    
    ax.grid(True, which = 'major', linestyle='-', zorder= 2)
    ax.grid(True, which = 'minor', linestyle=':', zorder= 2, alpha= 0.75)
    
#    plt.savefig('maximumnorm.pdf')
    
    plt.show()
    
if __name__=='__main__':

# maximal matrix size (something >100 is required for a sufficient detailed plot)
    q_max = 120

# parameter for the random walk weights in front of the generators of the discrete Heisenberg group/ mulitplicatives in front of the matrix representation of the non-commmutative torus/some parameters given by the atomic lattice in QHE    
    u = 11/48
    v = 12/48
    w = 1/48
    

    main(q_max,u, v, w)
    
    
    
    
