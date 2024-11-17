#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

from matplotlib import rc , rcParams
rcParams[ 'text.usetex'     ] = True
rcParams[ 'axes.labelsize'  ] = 15
rcParams[ 'xtick.labelsize' ] = 15
rcParams[ 'ytick.labelsize' ] = 15

problem = [ '../Images/AdvectedSineTFF/' , \
            '../Images/AdvectedSineTTF/' , \
            '../Images/AdvectedSineFFF/' ]

variable = 'PF_D'

Np = np.array( [ '3' ] )#, '2' , '3' ] )
nX = np.array( [ '16' , '32' , '64' , '128' ] )

L1 = np.empty( ( len( problem ) , len( Np ) , len( nX ) )  , float )

for p in xrange( len( problem ) ):
    for N in xrange( len( Np ) ):
        for K in xrange( len( nX ) ):
            approx = np.loadtxt( problem[p] + 'Np' + Np[N] + '_' \
                                 + 'nX' + nX[K] + '_' + variable + '.dat' )
            L1[p,N,K] = np.sum( np.abs( approx[:,1] - approx[:,0] ) )

nXf = nX.astype( 'float' )
Npf = Np.astype( 'float' )

def f( x , m , b ):
    return m * x + b

fig = plt.figure()
fig.suptitle( 'Convergence Rate (CR) for Advected Sine Wave\nSSP-RK3, \
               Third-Order Accurate Method' )

ax1 = fig.add_subplot(111)

c = [ 'k^' , 'ko' , 'kx' ]
nDOF     = Npf[0] * nXf
logDOF   = np.log10( nDOF )
x        = logDOF
x_interp = np.linspace( np.min( x ) , np.max( x ) , 100 )

for i in range( len( problem ) ):

    y  = np.log10( L1[i,0,:] / nDOF )
    popt , pcov = curve_fit( f , x , y )
    yi = interp1d( x , popt[0] * x + popt[1] )

    ax1.plot( x_interp , yi(x_interp) , 'k-' )

    if   ( i == 0 ):
        label = 'Limiter, no detector: CR=%.2f' % ( popt[0] )
        ms    = 7
    elif ( i == 1 ):
        label = 'Limiter and detector: CR=%.2f' % ( popt[0] )
        ms    = 7
    elif ( i == 2 ):
        label = 'No limiter: CR=%.2f' % ( popt[0] )
        ms    = 15
    ax1.plot( x , y , c[i] , label = label , markersize = ms )

ax1.set_xlabel( r'$\ell og_{10}\left(nDOF\right)$' , labelpad = -0.05 )
ax1.set_ylabel( r'$\ell og_{10}\left(L^{1}error/nDOF\right)$' )

ax1.legend()
plt.savefig( '../Images/CR_' + variable + '.png' )
#plt.show()
plt.close()

