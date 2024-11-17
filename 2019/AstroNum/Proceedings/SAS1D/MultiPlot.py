#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure( figsize = (12,3) )
plt.rcParams['mathtext.fontset'] = 'stix'

lw = 1.5
ymin = 4.0e6
ymax = 3.0e10
xmin = 40.0
xmax = 540.0

nX = [ '032', '064', '128', '256' ]


for i in range( len( nX ) ):

    ax = fig.add_subplot( 1, 4, i+1 )
    rho = np.loadtxt( 'MovieData_PF_D_{:}.dat'.format( nX[i] ) )

    x = np.linspace( xmin, xmax, int( nX[i] ) )

    alpha = np.linspace( 0, 1, rho.shape[0] )
    for j in range( rho.shape[0] ):
        if( j % 10 == 0 ):
            ax.semilogy( x, rho[j], 'k-', alpha = alpha[j] )

    ax.semilogy( x, rho[0 ], 'r-',  label = 't = 0 ms'   )
    ax.semilogy( x, rho[-1], 'r--', label = 't = 300 ms' )

    ax.set_xlim( xmin, xmax )
    ax.set_ylim( ymin, ymax )

    ax.text( 200, 1e10, 'nX = {:}'.format( nX[i].lstrip( '0' ) ), fontsize = 18 )

    ax.spines['top'   ].set_linewidth( lw )
    ax.spines['bottom'].set_linewidth( lw )
    ax.spines['left'  ].set_linewidth( lw )
    ax.spines['right' ].set_linewidth( lw )

    ax.tick_params( which = 'both', axis = 'both', direction = 'in' )

    if( nX[i] == '032' ):
        ax.set_xlabel( 'Radial Coordinate [km]' )
        ax.set_ylabel( 'Rest-Mass--Density ' \
                         + r'$\left[\mathrm{g\,cm^{-3}}\right]$', \
                       fontsize = 12 )
    if( nX[i] == '064' ):
        ax.set_xlabel( 'Radial Coordinate [km]' )
        ax.set_yticks( [] )
    if( nX[i] == '128' ):
        ax.set_yticks( [] )
        ax.set_xlabel( 'Radial Coordinate [km]' )
    if( nX[i] == '256' ):
        ax.set_xlabel( 'Radial Coordinate [km]' )
        ax.set_yticks( [] )
        ax.legend( loc = (0.35,0.6) )

plt.subplots_adjust( hspace = 0.02, wspace = 0.02 )
#plt.show()
plt.savefig( 'SAS1D.png', bbox_inches = 'tight', pad_inches = 0.1, dpi = 300 )
