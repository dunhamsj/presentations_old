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
#ax.text( xLU-0.08, 1.35*yScale[0]*yL, r'$x_{L}$', fontsize = fs )
#ax.text( xUL-0.08, 1.35*yScale[0]*yL, r'$x_{U}$', fontsize = fs )

# Polynomials

cp = 'red'
ck = 'black'
cn = 'blue'

# Assumes reference element defined on [ -1/2, +1/2 ]
if  ( N == 2 ):
    xx = np.array( [-1.0/2.0,+1.0/2.0]    , np.float64 )
elif( N == 3 ):
    xx = np.array( [-1.0/2.0,0.0,+1.0/2.0], np.float64 )

# x
Nx = 1000#N
xp = np.linspace( xL , xLU, Nx )
xk = np.linspace( xLU, xUL, Nx )
xn = np.linspace( xUL, xU , Nx )

# Plot interpolation points
a0 = 0.3
a1 = 0.5
ax.plot( xx-1.0, a0 + a1 * ( yL + yU ) * np.ones( N ), '.', c = cp )
ax.plot( xx+1.0, a0 + a1 * ( yL + yU ) * np.ones( N ), '.', c = cn )
ax.plot( xx    , a0 + a1 * ( yL + yU ) * np.ones( N ), '.', c = 'magenta' )
ax.plot( xx[1] , a0 + a1 * ( yL + yU ) * np.ones( 1 ), '.', c = ck )

ax.text( xx[0]-0.25, 0.1, r'$x^{K-1}_{3}$', c = cp )
ax.text( xx[0]+0.05, 0.1, r'$x^{K}_{1}$'  , c = ck )
ax.arrow( xx[0], a0 + a1 * ( yL + yU ), -0.1, -0.1, color = cp )
ax.arrow( xx[0], a0 + a1 * ( yL + yU ), +0.1, -0.1, color = ck )

ax.text( xx[1]-0.07, 0.15, r'$x^{K}_{2}$', c = ck )

ax.text( xx[2]-0.25, 0.1, r'$x^{K}_{3}$'  , c = ck )
ax.text( xx[2]+0.05, 0.1, r'$x^{K+1}_{1}$', c = cn )
ax.arrow( xx[2], a0 + a1 * ( yL + yU ), -0.1, -0.1, color = ck )
ax.arrow( xx[2], a0 + a1 * ( yL + yU ), +0.1, -0.1, color = cn )

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

# Functions to interpolate
if( N == 2 ):
    def fp(x):
        return 0.4 * x
    def fk(x):
        return 0.2 * x + 0.2
    def fn(x):
        return -0.2 * x + 0.15
if( N == 3 ):
    a0 = -0.1
    a1 = 0.1
    a2 = 0.2
    def fp(x):
        x = x - 1.0
        return a0 + a1 * x + a2 * x**2
    def fk(x):
        return a0 + a1 * x + a2 * x**2
    def fn(x):
        x = x + 1.0
        return a0 + a1 * x + a2 * x**2

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
ax.plot( xx-1.0, fpp   , 'x', c = cp )
ax.plot( xx+1.0, fnn   , 'x', c = cn )
ax.plot( xx    , fkk   , 'x', c = 'magenta' )
ax.plot( xx[1] , fkk[1], 'x', c = ck )

y0 = yL + 0.22
ax.text( xx[0]-0.25, y0, r'$f^{K-1}_{3}$', c = cp )
ax.text( xx[0]+0.05, y0, r'$f^{K}_{1}$'  , c = ck )
ax.arrow( xx[0], fkk[0], -0.09, -0.07, color = cp )
ax.arrow( xx[0], fkk[0], +0.09, -0.07, color = ck )

y0 = yL + 0.25
ax.text( xx[1]-0.05, y0, r'$f^{K}_{2}$',c = ck )

y0 = yL + 0.3
ax.text( xx[2]-0.2, y0, r'$f^{K}_{3}$'  , c = ck )
ax.text( xx[2]+0.05, y0, r'$f^{K+1}_{1}$', c = cn )
ax.arrow( xx[2], fkk[2], -0.09, -0.07, color = ck )
ax.arrow( xx[2], fkk[2], +0.09, -0.07, color = cn )

# Element labels
fs = 10
y0 = yL + 0.1
ax.text( xL  + 0.4 , y0, r'$K-1$', fontsize = fs )
ax.text( xLU + 0.45, y0, r'$K$'  , fontsize = fs )
ax.text( xUL + 0.4 , y0 , r'$K+1$', fontsize = fs )

# Display equation on top
ax.text( xL + 0.385 * Width, yL + 1.0*Height, \
         r'$f\left(x\right)=\sum_{i=1}^{N\left(=3\right)}f_{i}\,\phi_{i}\left(x\right)$' )

#plt.show()
plt.savefig( 'CG.png', dpi = 300, bbox_inches = 'tight' )

