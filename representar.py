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

X = np.linspace(xmin, xmax, 2001) # real
Y = np.linspace(ymin, ymax, 2001) # imag
X, Y = np.meshgrid(X, Y)
s = X + 1j*Y # variable s

b = [1, 0, 1] # c_0s^2 + c_1s + c_2 en ese orden (inversa de polyval)
a = [1, 2*2*1, 1]

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


#%% RESONANCIAS
plt.close('all')
xi=0.5
w = np.linspace(0, 100, 4001)
s = 1j*w # variable s

b = [1,0, 50**2] # c_0s^2 + c_1s + c_2 en ese orden (inversa de polyval)

# Esta función define la función de transferencia
def h(s, b, a):
    z = polyval(s, b[::-1])/polyval(s, a[::-1])
    return z

fig, ax = plt.subplots(figsize=(6,5))
# for xi in [0,0.3, 0.5, 1, 1.5, 2, 5]:
for xi in [50, 100, 500]:
    a = [1, 0.5*xi, xi**2]
    Z = abs(h(s, b, a))
    ax.plot(w, Z, label=fr'$\xi={xi}$')
# ax.set_ylim(0, 2.5)
ax.set_ylabel(r'$|H(s)|$')
ax.set_xlabel(r'$\Im(s)=\omega$')
plt.legend()


#%% RESONANCIAS
import matplotlib
plt.close('all')
xi=0.5
w = np.linspace(0, 3, 4001)
s = 1j*w # variable s

b = [1,0,0.5] # c_0s^2 + c_1s + c_2 en ese orden (inversa de polyval)

# Esta función define la función de transferencia
def h(s, b, a):
    z = polyval(s, b[::-1])/polyval(s, a[::-1])
    return z

fig, ax = plt.subplots(figsize=(6,5))
for xi in [0,0.3, 0.5, 1, 1.5, 2, 5]:
    a = [1, 2*xi*1, 1]
    Z = np.angle(h(s, b, a))
    ax.plot(w, Z*180/np.pi, label=fr'$\xi={xi}$')
ax.set_ylim(-95, 185)
ax.set_ylabel(r'$\angle H(s) (degrees)$')
ax.set_xlabel(r'$\Im(s)=\omega$')
ax.set_yticks([0, -45, -90, 45, 90, 135, 180])
plt.legend()


#%% RESONANCIAS
plt.close('all')
xi=0.5
w = np.linspace(0, 3, 4001)
s = 1j*w # variable s

b = [1,0] # c_0s^2 + c_1s + c_2 en ese orden (inversa de polyval)

# Esta función define la función de transferencia
def h(s, b, a):
    z = polyval(s, b[::-1])/polyval(s, a[::-1])
    return z

fig, ax = plt.subplots(figsize=(4,3))
for xi in [1, 0.6]:
    a = [1, 2*xi*1, 1]
    Z = abs(h(s, b, a))
    ax.plot(w, Z, label=fr'$\xi={xi:.3f}$')
ax.set_ylim(0, 1.1)
ax.set_ylabel(r'$|H(s)|$')
ax.set_xlabel(r'$\Im(s)=\omega$')
plt.legend()
plt.tight_layout()

#%% DECIBELIOS
plt.close('all')
xi=0.5
w = np.linspace(0, 1000, 4001)
s = 1j*w # variable s

b = [1,0] # c_0s^2 + c_1s + c_2 en ese orden (inversa de polyval)

# Esta función define la función de transferencia
def h(s, b, a):
    z = polyval(s, b[::-1])/polyval(s, a[::-1])
    return z

fig, ax = plt.subplots(1, 2, figsize=(8,3))
xi=.6
wo = 100
a = [1, 2*xi*wo, wo**2]
Z = abs(h(s, b, a))
ax[0].plot(w, Z, label=fr'$\xi={xi:.3f}$')
ax[0].set_ylabel(r'$|H(s)|$')
ax[0].set_xlabel(r'$\Im(s)=\omega$')
ax[0].set_title('escala lineal')
ax[0].grid(True,which="both",ls="-")

Zdb = 20*np.log10(Z)
ax[1].semilogx(w, Zdb, label=fr'$\xi={xi:.3f}$')
ax[1].set_ylabel(r'$20\log|H(s)|$ (dB)')
ax[1].set_xlabel(r'$\Im(s)=\omega$')
ax[1].set_title('escala logaritmica')
ax[1].grid(True,which="both",ls="-")
plt.tight_layout()



#%% OP AMP
plt.close('all')
xi=0.5
w = np.linspace(0, 1e7, 4001)
s = 1j*w # variable s

# Esta función define la función de transferencia
def hs(s, Ao=2e5, R1=100, R2=200,wc=2*np.pi*5):
    z = - (Ao*wc*R2/(R2+R1))*(1/(s+wc*(1+R1*Ao/(R2+R1))))
    return z

fig, ax = plt.subplots(figsize=(6,5))
for R2 in [50, 100, 200, 500, 1000, 2000, 5000, 1e4]:
    Z = abs(hs(s, R2=R2))
    ax.semilogx(w, Z, label=fr'$R2={R2}$')
# ax.set_ylim(0, 2.5)
ax.set_ylabel(r'$|H(s)|$')
ax.set_xlabel(r'$\Im(s)=\omega$')
ax.set_title('Para R1=100')
plt.legend()