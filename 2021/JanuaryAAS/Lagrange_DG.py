#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family']      = 'serif'
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

# Number of interpolation nodes
N = 3

# Interface points
xL  = -0.5
xU  = +0.5
Width = xU - xL

# Vertical extent
yL = -1.0
yU = +2.0
Height = yU - yL

# Font size
fs = 13

fig, ax =  plt.subplots( 1, 1, figsize = (6,2) )
ax.set_xlim( xL, xU )
ax.set_ylim( yL, yU )

# Assumes reference element defined on [ -1/2, +1/2 ]
if  ( N == 2 ):
    xx = np.array( [-1.0/np.sqrt(12.0),+1.0/np.sqrt(12.0)]    , np.float64 )
elif( N == 3 ):
    xx = np.array( [-np.sqrt(3.0/20.0),0.0,+np.sqrt(3.0/20.0)], np.float64 )

# Plot grid points
plt.plot( xx, yL + 0.5 * Height * np.ones( xx.shape[0] ), 'k.' )

Nx = 1000
x = np.linspace( xL, xU, Nx )

def Lagrange( x, i ):
    L = 1.0
    for j in range( N ):
        if( j != i ):
            L *= ( x - xx[j] ) / ( xx[i] - xx[j] )
    return L

L = np.empty( (N,Nx), np.float64 )
for i in range( N ):
    for k in range( Nx ):
        L[i,k] = Lagrange( x[k], i )

# Plot polynomials
c = [ 'r', 'g', 'b' ]
for i in range( N ):
    plt.plot( x, L[i], c[i], label = r'$\ell_{{{:d}}}$'.format( i+1 ) )

plt.legend(loc=(0.7,0.64),labelspacing=0.01)
xticks = [ -0.5, -0.25, 0.0, +0.25, 0.5 ]
plt.xticks( xticks )
lw = 0.5

ax.text( -0.25, 1.5, r'$\ell_{i}\left(\xi_{j}\right)=\delta_{ij}$' )

for i in range( N ):
    plt.axvline( xx[i], ls = '--', color = 'black', lw = lw )

plt.axhline( 0.0, ls = '--', color = 'black', lw = lw )
plt.axhline( 1.0, ls = '--', color = 'black', lw = lw )

plt.xlabel( r'$\xi$' )
plt.title( 'DG Basis Polynomials' )

#plt.show()
plt.savefig( 'Lagrange_DG.png', dpi = 300, bbox_inches = 'tight' )

