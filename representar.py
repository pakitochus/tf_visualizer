'''
=========================================
ANÁLISIS DE LA FUNCIÓN DE TRANSFERENCIA
=========================================

Este script genera una animación en la que se observa el valor absoluto de la
función de transferencia |H(s)| en función de los valores real e imaginario 
de s. También anima para ver el corte en 0 y observar el filtro, así como
la proyección vertical donde se observa el diagrama de polos y ceros. 

Created on Mon Apr  6 18:09:49 2020

@author: pakitochus
'''

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
from numpy.polynomial.polynomial import polyval

# Para limpiar las ventanas como debe ser, que para eso estamos confinados. 
plt.close('all')

# Generamos los datos.
xmin, xmax = -5, 5
ymin, ymax = -5, 5

X = np.linspace(xmin, xmax, 401) # real
Y = np.linspace(ymin, ymax, 401) # imag
X, Y = np.meshgrid(X, Y)
s = X + 1j*Y # variable s

b = [1, 0, 1] # c_0 + c_1s + c_2s^2 etc...
a = [1, 1.4, 4]

# Esta función define la función de transferencia
def h(s, b, a):
    z = polyval(s, b[::-1])/polyval(s, a[::-1])
    return z

Z = abs(h(s, b, a))

# Parámetros para la visualización: 
vmin, vmax = 0, 3
# Esto lo hacemos para saturar los valores infinitos (polos)
Z[Z>vmax] = vmax

params = {'inicial': {'elev': 30, 'azim': -60, 'cut': None},
          'corte': {'elev': 30, 'azim': -60, 'cut': 0},
          'frontal': {'elev': 0, 'azim': 0, 'cut': 0},
          'superior': {'elev': 90, 'azim':-90, 'cut': None}}

stage = 'corte'
fig = plt.figure(figsize=(6,5))

for ix, stage in enumerate(params.keys()):
    
    if stage=='superior':
        ax = fig.add_subplot(f'22{ix+1}')
        ax.pcolormesh(X,Y,Z, cmap=cm.coolwarm, vmin=vmin, vmax=vmax)
        from matplotlib.patches import Circle
        from scipy.signal import tf2zpk
        p = Circle((0,0), 1, fill=False)
        ax.add_patch(p)
        z, p, k = tf2zpk(b, a)
        ax.scatter(p.real, p.imag, marker='x', zorder=10, c='black')
        ax.plot(z.real, z.imag, 'o', zorder=10, c='black', fillstyle='none')
        ax.set_aspect('equal')
        
    else:
        ax = fig.add_subplot(f'22{ix+1}', projection='3d')
        if params[stage]['cut'] is None:
            surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, vmin=vmin, vmax=vmax,
                                   zorder=-10)
        else: 
            select = (X<params[stage]['cut']).all(axis=0)
            surf = ax.plot_surface(X[:,select], Y[:,select], Z[:,select], 
                                   cmap=cm.coolwarm, vmin=vmin, vmax=vmax,
                                   zorder=-10)    
            cset = ax.contour(X, Y, Z, zdir='x', offset=0, levels=0, colors='red',linewidths=3)
        ax.set_zlabel('|H(s)|')
        ax.set_zlim(vmin, vmax)
        ax.azim = params[stage]['azim']
        ax.elev = params[stage]['elev']
        
        
    # Customize the z axis.
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel(r'$\Re(s)=\sigma$')
    ax.set_ylabel(r'$\Im(s)=\omega$')
    
    import sympy as sp
    from sympy.abc import s
    from sympy import latex
    sp.init_printing()
    ax.set_title('')
    # fig.update() 
fig.suptitle(r'$H(s)=%s$'%latex(sp.Poly(b,s)/sp.Poly(a,s)))