#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family']      = 'serif'
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

# Number of interpolation nodes
N = 3

# Interface points
xL   = -1.0
xLU  =  1.0
xU   = +3.0
Width = xU - xL

# Vertical extent
yL = -1.0
yU = +2.0
Height = yU - yL

# Font size
fs = 13

fig, ax =  plt.subplots( 1, 1, figsize = (6,2) )
ax.set_xlim( xL-0.05, xU+0.05 )
ax.set_ylim( yL, yU )

# Assumes reference element defined on [ -1/2, +1/2 ]
if  ( N == 2 ):
    xx = np.array( [-1.0,+1.0]    , np.float64 )
elif( N == 3 ):
    xx = np.array( [-1.0,0.0,+1.0], np.float64 )

Nx = 1000
x1 = np.linspace( xL , xLU, Nx )
x2 = np.linspace( xLU, xU , Nx )

def Lagrange( x, i ):
    L = 1.0
    for j in range( N ):
        if( j != i ):
            L *= ( x - xx[j] ) / ( xx[i] - xx[j] )
    return L

L = np.empty( (N,Nx), np.float64 )
for i in range( N ):
    for k in range( Nx ):
        L[i,k] = Lagrange( x1[k], i )

# Plot polynomials
for i in range( N ):
    plt.plot( x1, L[i], 'r', label = r'$\ell_{{{:d}}}$'.format( i+1 ) )
    plt.plot( x2, L[i], 'b', label = r'$\ell_{{{:d}}}$'.format( i+1 ) )

#plt.legend(loc=(0.7,0.64),labelspacing=0.01)
xticks = [ -1.0, -0.5, 0.0, +0.5, +1.0, 1.5, 2.0, 2.5, 3.0 ]
plt.xticks( xticks )
lw = 0.5

plt.axvline( xx[0], ls = '--', color = 'black', lw = lw )
plt.axvline( xx[1], ls = '--', color = 'black', lw = lw )
plt.axvline( xx[1]+2.0, ls = '--', color = 'black', lw = lw )
plt.axvline( xx[2]+2.0, ls = '--', color = 'black', lw = lw )
plt.axvline( xx[-1], ls = '--', color = 'black', lw = 1.5 )

plt.axhline( 0.0, ls = '--', color = 'black', lw = lw )
plt.axhline( 1.0, ls = '--', color = 'black', lw = lw )

plt.text( 0.0, yL + 0.8*Height, 'K' )
plt.text( 2.0, yL + 0.8*Height, 'K+1' )

plt.xlabel( r'$\xi$' )
plt.title( 'CG Basis Polynomials' )
plt.xlabel( r'$\xi$' )

# Plot grid points
plt.plot( xx, yL + 0.5 * Height * np.ones( xx.shape[0] ), 'r.' )
plt.plot( xx+2.0, yL + 0.5 * Height * np.ones( xx.shape[0] ), 'b.' )
plt.plot( xx[-1], yL + 0.5 * Height, 'm.' )

#plt.show()
plt.savefig( 'Lagrange_CG.png', dpi = 300, bbox_inches = 'tight' )

