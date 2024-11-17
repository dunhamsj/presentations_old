#!/usr/bin/env python3

import yt
import numpy as np
import subprocess
from os import listdir
from os.path import isfile
from sys import argv, exit
import matplotlib.pyplot as plt

yt.funcs.mylog.setLevel(0) # Suppress initial yt output to screen

# --- Get user's HOME directory ---
HOME = subprocess.check_output( ["echo $HOME"], shell = True)
HOME = HOME[:-1].decode( "utf-8" ) + '/'

# Specify directory containing plotfiles

ProblemDirectory \
  = HOME + 'Research/DataFromPastRuns/AstroNum2019/SASI/'
ProblemName = 'SASI'

File = 'thornado_00000000'

ds = yt.load( '{:}'.format( ProblemDirectory + File ) )

MaxLevel = ds.index.max_level
Time     = ds.current_time
nX       = ds.domain_dimensions
xL       = ds.domain_left_edge
xH       = ds.domain_right_edge

CoveringGrid \
  = ds.covering_grid \
      ( level           = MaxLevel, \
        left_edge       = xL, \
        dims            = nX * 2**MaxLevel, \
        num_ghost_zones = nX[0] )

from scipy.integrate import simps
X2 = np.linspace( xL[1].to_ndarray(), xH[1].to_ndarray(), nX[1] )
def AngleAverage( integrand ):
    return 0.5 * simps( integrand * np.sin( X2 ), X2 )

def ComputeQuantities( CoveringGrid ):
    # Lorentz Factor
    c = 2.99792458e8

    # Convert Vi from km/s to SI units
    V1 = CoveringGrid['PF_V1'].to_ndarray()[:,:,0] * 1.0e3 # m/s
    V2 = CoveringGrid['PF_V2'].to_ndarray()[:,:,0] * 1.0e3 # rad/s
    V3 = CoveringGrid['PF_V3'].to_ndarray()[:,:,0] * 1.0e3 # rad/s

    # gamma_ii all in SI units
    g1 = CoveringGrid['GF_Gm_11'].to_ndarray()[:,:,0] # dimensionless
    g2 = CoveringGrid['GF_Gm_22'].to_ndarray()[:,:,0] # m^2
    g3 = CoveringGrid['GF_Gm_33'].to_ndarray()[:,:,0] # m^2

    W    = 1.0 / np.sqrt( 1.0 - ( g1*V1**2 + g2*V2**2 + g3*V3**2 ) / c**2 )
    W_av = np.empty( nX[0], float )

    # Specific Enthalpy
    c = 2.99792458e10

    rho = CoveringGrid['PF_D'].to_ndarray()[:,:,0]
    e   = CoveringGrid['PF_E'].to_ndarray()[:,:,0]
    p   = CoveringGrid['AF_P'].to_ndarray()[:,:,0]

    h = ( c**2 + ( e + p ) / rho ) / c**2
    h_av = np.empty( nX[0], float )

    for i in range( nX[0] ):
        W_av[i] = AngleAverage( W[i,:] )
        h_av[i] = AngleAverage( h[i,:] )

    psi   = CoveringGrid['GF_Psi'  ].to_ndarray()[:,0,0]
    alpha = CoveringGrid['GF_Alpha'].to_ndarray()[:,0,0]

    return W_av, h_av, psi, alpha

W_av, h_av, psi, alpha = ComputeQuantities( CoveringGrid )

x = np.linspace( xL[0].to_ndarray(), xH[0].to_ndarray(), nX[0] )

# Newtonian gravitational potential
M       = 2.8 * 2.0e30
G       = 6.67e-11
phi     = -G * M / ( x * 1.0e3 ) # convert x to meters
c       = 2.99792458e8
alpha_N = 1.0 + 2.0 * phi / c**2

# Plotting

plt.rcParams['text.usetex'] = True
fig = plt.figure()
ax1 = fig.add_subplot( 111 )

me = 10
ms = 3

ax1.plot( x, alpha / alpha_N, 'k-', \
          markevery = me, markersize = ms, label = r'$\alpha/\alpha_{N}$' )
ax1.plot( x, h_av,  'k--', \
          markevery = me, markersize = ms, label = r'$h/c^2$'  )
ax1.plot( x, psi,   'k-.', \
          markevery = me, markersize = ms, label = r'$\psi$'   )
ax1.plot( x, W_av,  'k:', \
          markevery = me, markersize = ms, label = r'$W$'      )
ax1.set_xlim( x.min(), x.max() )
ax1.set_xlabel( 'Radial Coordinate r [km]', fontsize = 13 )
ax1.tick_params( which = 'both', direction = 'in' )
ax1.legend( loc = 1 )#, prop = {'size': 8} )
lw = 1.5
ax1.spines['top'   ].set_linewidth( lw )
ax1.spines['bottom'].set_linewidth( lw )
ax1.spines['left'  ].set_linewidth( lw )
ax1.spines['right' ].set_linewidth( lw )

# Inset
vmin = 1.0e15

xL  = xL.to_ndarray()
xH  = xH.to_ndarray()
dX1 = ( xH[0] - xL[0] ) / nX[0]; dX1 = 0
dX2 = ( xH[1] - xL[1] ) / nX[1]; dX2 = 0
X1  = np.linspace( xL[0] + dX1, xH[0] - dX1, nX[0] )
X2  = np.linspace( xL[1] + dX2, xH[1] - dX2, nX[1] )
theta, r = np.meshgrid( X2, X1 )

from matplotlib.colors import LogNorm
norm = LogNorm()

rho   = CoveringGrid['PF_D' ].to_ndarray()[:,:,0]
P     = CoveringGrid['AF_P' ].to_ndarray()[:,:,0]
Gamma = CoveringGrid['AF_Gm'].to_ndarray()[:,:,0]

Data = P / rho**Gamma

ax2 = fig.add_axes( [0.35,0.4,0.4,0.4], polar = True )
im = ax2.pcolormesh( theta, r, Data, cmap = 'jet', norm = norm )
ax2.set_thetamin(180.0/np.pi*X2[0])
ax2.set_thetamax(180.0/np.pi*X2[-1])
ax2.set_theta_zero_location( 'N' )
ax2.set_theta_direction( -1 )
ax2.set_yticklabels( [] )
ax2.set_xticks( [ 0.0, np.pi/2, np.pi ] )

from matplotlib import ticker

cbaxes = fig.add_axes( [0.44, 0.4, 0.03, 0.4] )
cbar = fig.colorbar( im, \
                     cax = cbaxes, \
                     ticks = [ np.min(Data), np.max(Data) ] )

cbar.ax.yaxis.set_ticks_position( 'left' )

cbar.ax.set_ylabel( \
  'Polytropic Constant\n$[\mathrm{erg/cm^{3}/(g/cm^{3})^{4/3}}]$', \
  rotation = 90, labelpad = -50.0 )
cbar.ax.set_yticklabels( \
  [ '{:.0e}'.format( np.min( Data ) ), '{:.0e}'.format( np.max( Data ) ) ] )
#plt.show()
plt.savefig( 'GR1D.png' )
plt.close()
