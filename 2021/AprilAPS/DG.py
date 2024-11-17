#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family']      = 'serif'
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

# Number of interpolation nodes
N = 3

# Interface points
xL  = -1.5
xLU = -0.5
xUL = +0.5
xU  = +1.5
Width = xU - xL

# Vertical extent
yL = -0.5
yU = +0.5
Height = yU - yL

# Font size
fs = 13

fig, ax =  plt.subplots( 1, 1, figsize = (6,2) )
plt.axis('off')

# Extent of element
ax.set_xlim( xL-0.1, xU+0.1 )
ax.set_ylim( yL, 1.1*yU )

# Interfaces
yScale = [ 0.45, 1.0 ]
y = [ yScale[0] * yL, yScale[1] * yU ]
ax.plot( [xL ,xL] , y, c = 'k', lw = 0.5 )
ax.plot( [xLU,xLU], y, c = 'k', lw = 0.5 )
ax.plot( [xUL,xUL], y, c = 'k', lw = 0.5 )
ax.plot( [xU ,xU] , y, c = 'k', lw = 0.5 )

# Interface labels
fsx = 10
ax.text( xLU-0.055, 1.35*yScale[0]*yL, r'$x_{L}$', fontsize = fsx )
ax.text( xUL-0.055, 1.35*yScale[0]*yL, r'$x_{U}$', fontsize = fsx )

# Polynomials

cp = 'red'
ck = 'black'
cn = 'blue'

# Assumes reference element defined on [ -1/2, +1/2 ]
if  ( N == 2 ):
    xx = np.array( [-1.0/np.sqrt(12.0),+1.0/np.sqrt(12.0)]    , np.float64 )
elif( N == 3 ):
    xx = np.array( [-np.sqrt(3.0/20.0),0.0,+np.sqrt(3.0/20.0)], np.float64 )

# x
Nx = 1000#N
xp = np.linspace( xL , xLU, Nx )
xk = np.linspace( xLU, xUL, Nx )
xn = np.linspace( xUL, xU , Nx )

# Plot interpolation points
ax.plot( xx-1.0, 0.5 * ( yL + yU ) * np.ones( N ), '.', c = cp )
ax.plot( xx    , 0.5 * ( yL + yU ) * np.ones( N ), '.', c = ck )
ax.plot( xx+1.0, 0.5 * ( yL + yU ) * np.ones( N ), '.', c = cn )

# Interpolation point labels
for i in range( N ):
    ax.text( xx[i]-0.05, -0.1, r'$x_{{{:}}}$'.format(i+1) )

def Lagrange( x, i ):
    L = 1.0
    for j in range( N ):
        if( j != i ):
            L *= ( x - xx[j] ) / ( xx[i] - xx[j] )
    return L

L = np.empty( (N,Nx), np.float64 )
xi = np.linspace( -0.5, +0.5, Nx )
for i in range( N ):
    for k in range( Nx ):
        L[i,k] = Lagrange( xi[k], i )

y0 = yL + 0.2
ax.text( xL  + 0.4 , y0, r'$K-1$' , color = cp, fontsize = fs )
ax.text( xLU + 0.45, y0, r'$K$'   , color = ck, fontsize = fs )
ax.text( xUL + 0.4 , y0 , r'$K+1$', color = cn, fontsize = fs )

# Functions to interpolate
if( N == 2 ):
    def fp(x):
        return 0.4 * x
    def fk(x):
        return 0.2 * x + 0.2
    def fn(x):
        return -0.2 * x + 0.15
if( N == 3 ):
    ys = 0.1
    def fp(x):
        return ys + 0.4 * x**2
    def fk(x):
        return ys + 0.2 * x**2 + 0.1
    def fn(x):
        return ys + -0.2 * x**2 + 0.25

# Previous element
fpp = fp(xx)
yp = np.zeros( Nx, np.float64 )
for i in range( N ):
    yp += fpp[i] * L[i]
ax.plot( xp, yp, c = cp )

# Current element
fkk = fk(xx)
yk = np.zeros( Nx, np.float64 )
for i in range( N ):
    yk += fkk[i] * L[i]
ax.plot( xk, yk, c = ck )

# Next element
fnn = fn(xx)
yn = np.zeros( Nx, np.float64 )
for i in range( N ):
    yn += fnn[i] * L[i]
ax.plot( xn, yn, c = cn )

# Plot nodal points
ax.plot( xx-1.0, fpp, 'x', c = cp )
ax.plot( xx    , fkk, 'x', c = ck )
ax.plot( xx+1.0, fnn, 'x', c = cn )
for i in range( N ):
    ax.text( xx[i]-0.05, fkk[i]+0.07, r'$U_{{{:}}}$'.format(i+1) )

# Display equation on top
ax.text( xL + 0.36 * Width, yL + 1.0*Height, \
         r'$U\left(x\right)=\sum_{i=1}^{N\left(=3\right)}U_{i}\,\ell_{i}\left(x\right)$' )

plt.show()
#plt.savefig( 'Poster/Images/DG.png', dpi = 300, bbox_inches = 'tight' )

