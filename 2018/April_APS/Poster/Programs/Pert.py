#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rc, rcParams
rcParams[ 'text.usetex'     ] = True
rcParams[ 'axes.labelsize'  ] = 15
rcParams[ 'xtick.labelsize' ] = 15
rcParams[ 'ytick.labelsize' ] = 15

variable = 'PF_D'

problem = np.array( [ [ '../Images/PertSine/TTT' , 'Reference, 2000 cells' , 'k-' ] , \
                      [ '../Images/PertSine/FFT' , 'No limiter'            , 'r-' ] , \
                      [ '../Images/PertSine/TFT' , 'Limiter, no detector'  , 'g-' ] , \
                      [ '../Images/PertSine/TTT' , 'Limiter and detector'  , 'b-' ] ] )

ax1 = plt.subplot( 111 )
axInset = plt.axes( [ 0.13 , 0.117 , 0.61 , 0.21 ] , facecolor='w')
xmin = 0.64
xmax = 0.83

for i in range( len( problem ) ):
    if ( i == 0 ):
        suffix = '_2000.dat'
        label = 'Reference, 2000 cells'
    else:
        suffix = '.dat'

    X , Yi , Yf = np.loadtxt( problem[i,0] + suffix , unpack = True )
    ax1.plot( X , Yf , problem[i,2] , label = problem[i,1] )
    arg = np.where( ( X > xmin ) & ( X < xmax ) )
    axInset.plot( X[arg] , Yf[arg] , problem[i,2] )

ax1.set_xlim( np.min( X ) , np.max( X ) )
ax1.set_xlabel( r'$x$'    , labelpad = -3 )
ax1.set_ylabel( r'$\rho$' )
ax1.legend( bbox_to_anchor= ( 0.45 , 1.0 ) , loc = 'upper center' )

axInset.set_xlim( np.min( X[arg] ) , np.max( X[arg] ) )
axInset.get_xaxis().set_ticks([])
axInset.get_yaxis().set_ticks([])

plt.suptitle( 'Perturtbed shock tube, SSP-RK3\nSecond-order accuracy, \
               200 cells, CFL=0.1, $\Gamma=5/3$' )
#plt.savefig( 'PertSine_' + variable + '.png' )
plt.show()

