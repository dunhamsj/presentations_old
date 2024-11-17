#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams[ 'text.usetex' ] = True

from NonRel_Bernoulli import *

# Lapse function, conformal factor, and Lorentz factor

'''
def alpha(r):
    return np.sqrt( 1.0 - 2.0 * G * M / ( c**(2.0) * r ) )

def psi(r):
    return alpha(r)**(-0.5)
'''

# Isotropic coordinates
def alpha(r):
    return np.abs( ( 1 - G * M / ( 2.0 * c**(2.0) * r ) ) \
                 / ( 1 + G * M / ( 2.0 * c**(2.0) * r ) ) )

def psi(r):
    return 1.0 + G * M / ( 2.0 * c**(2.0) * r )

def W(r,v):
    return ( 1.0 - psi(r)**(4.0) * v**(2.0) / c**(2.0) )**(-0.5)

########## Pre-shock fluid variables (pressure assumed to be zero)
v_pre   = -psi(r_pre)**(-2.0) * c * np.sqrt( 1.0 - alpha(r_pre)**(2.0) )
rho_pre = 1.0 / ( psi(r_pre)**(4.0) * np.sqrt( 1.0 - alpha(r_pre)**(2.0) ) ) \
            * Mdot / ( 4.0 * np.pi * r_pre**(2.0) * c )
h_pre   = c**(2.0)

########## Jump conditions
# Constants in jump conditions
C1 = rho_pre[-1] * alpha(r_pre[-1])**(-1.0) * v_pre[-1]
C2 = rho_pre[-1] * h_pre * alpha(r_pre[-1])**(-2.0) * ( v_pre[-1] / c )**(2.0)
C3 = rho_pre[-1] * h_pre * alpha(r_pre[-1])**(-2.0) * v_pre[-1]

# Find velocity just below shock with Newton's method
A = psi(Rs)**(8.0) * 1.0 /  ( gamma - 1.0 )**(2.0) * C3**(2.0) / c**(6.0)
B = -2.0 * psi(Rs)**(8.0) * gamma / ( gamma - 1.0 )**(2.0) * C2 * C3 / c**(4.0)
C = psi(Rs)**(4.0) / c**(2.0) * ( psi(Rs)**(4.0) * gamma*(2.0)         \
    / ( gamma - 1.0 )**(2.0) * C2**(2.0) + 2.0 * 1.0 / ( gamma - 1.0 ) \
    * C3**(2.0) / c**(2.0) + C1**(2.0) * c**(2.0) )
D = -2.0 * psi(Rs)**(4.0) * gamma / ( gamma - 1.0 ) * C2 * C3 / c**(2.0)
E = 1.0 / c**(2.0) * ( C3**(2.0) - C1**(2.0) * c**(4.0) )

# Newton's method

# Tolerance
epsilon = 1.e-16

# Max number of iterations
nMax = 100

def Newton( f , dfdv , v , r = 0.0 ):

    Delta_v = 1.

    nIter = 1

    while ( ( np.abs( Delta_v ) > epsilon ) and ( nIter < nMax ) ):

        Delta_v = -f( v , r ) / dfdv( v , r )

        v = Delta_v + v

        nIter += 1

    if ( nIter == nMax ):
        print 'Maxed out, dv = %.2e' , np.abs( Delta_v )

    return v

# Functions to get velocity just below the shock
def f(v,r):
    return A * v**(4.0) + B * v**(3.0) + C * v**(2.0) + D * v + E
def dfdv(v,r):
    return 4.0 * A * v**(3.0) + 3.0 * B * v**(2.0) + 2.0 * C * v + D

# Use Newtonian velocity just below the shock as initial guess
v2_N = -( gamma - 1. ) / ( gamma + 1. ) * np.sqrt( 2.0 * G * M / Rs )

v2 = Newton( f , dfdv , v2_N )

# Get density and pressure just below the shock
rho2 = np.abs(C1) * np.sqrt( 1.0 / v2**(2.0) - psi(Rs)**(4.0) / c**(2.0) )
p2 = ( gamma - 1.0 ) / gamma * ( C3 - rho2 * c**(2.0) * W(Rs,v2)**(2.0) * v2 ) \
                             / ( W(Rs,v2)**(2.0) * v2 )

# Get polytropic constant
K = p2 * rho2**(-gamma)

# CD, CE
def CD(r,v,rho):
    return psi(r)**(6.0) * alpha(r) * r**(2.0) * rho * W(r,v) * v

def CE(r,v,rho,h):
    return CD(r,v,rho) * alpha(r) * h * W(r,v)

CD_exact = CD(r_pre[-1],v_pre[-1],rho_pre[-1])

########## Post-shock fluid variables
def f(v,r):
    return gamma / ( gamma - 1.0 ) * K / c**(2.0) \
* ( CD_exact / ( psi(r)**(6.0) * alpha(r) * r**(2.0) * W(r,v) * v ) )**(gamma - 1.0) \
- 1.0 / ( alpha(r)*W(r,v) ) + 1.0

def dfdv(v,r):
    return -gamma * K / c**(2.0) \
* ( CD_exact / ( psi(r)**(6.0) * alpha(r) * r**(2.0) * W(r,v) * v ) )**(gamma-1.0) \
* ( psi(r)**(4.0) * v / c**(2.0) * W(r,v)**(2.0) + 1.0 / v ) \
+ psi(r)**(4.0) / alpha(r) * v / c**(2.0) * W(r,v)

v_post = np.empty(N, float)
for i in range( len( r_post ) ):
    v_post[i] = Newton(f,dfdv,v2,r_post[i])
    v2 = v_post[i]

rho_post = -1.0 / ( psi(r_post)**(6.0) * alpha(r_post) ) \
           * Mdot / ( 4.0 * np.pi * r_post**(2.0) * W(r_post,v_post) * v_post )
p_post = K * rho_post**(gamma)
h_post = c**(2.0) + gamma / ( gamma - 1.0 ) * p_post / rho_post

# Combine arrays
r   = np.hstack( ( r_post[::-1] , r_pre[::-1] ) )

h   = np.hstack( ( h_post[::-1] , c**(2.0) * np.ones( len( r_pre ) ) ) )

rho = np.hstack( ( rho_post[::-1] , rho_pre[::-1]            ) )
v   = np.hstack( ( v_post  [::-1] , v_pre  [::-1]            ) )
p   = np.hstack( ( p_post  [::-1] , np.zeros( len( v_pre ) ) ) )

### Specify pressure ahead of shock with constant Mach number M
M = 100.0

# Assume h = 1
p_pre = rho_pre * v_pre**2 / ( gamma * M**2 )
p = np.hstack( ( p_post[::-1] , p_pre[::-1] ) )

# Convert to cgs
Rs    *= 1.0e-3
Ri    *= 1.0e-3
Rf     = ( 2.0 * Rs )
r     *= 1.0e-3
R_PNS *= 1.0e-3

rho   *= 1.0e-3
p     *= 1.0e1
v     *= 1.0e-3

########## Plot fluid variables
plt.suptitle( 'Initial Conditions for SAS\n$M=1.4\,M_{\odot},\ \dot{M}=0.3\,M_{\odot}/s$' )
fs = 13

xlim = ( Ri , Rf )

ax_rho = plt.subplot(311)

ax_rho.semilogy( r , rho   , 'k-'   )
ax_rho.axvline( R_PNS , color = 'k' )
ax_rho.axvline( Rs , color = 'k' , ls = '--' )

ax_rho.text( R_PNS + 5 , 1.0e7  , '$R_{PNS}=40\,km$'        , fontsize = fs )
ax_rho.text( Rs    + 5 , 1.0e10 , '$R_{shock}=%i\,km$' % Rs , fontsize = fs )

ax_rho.set_xlim( xlim )
ax_rho.get_xaxis().set_ticks( [] )
ax_rho.set_ylabel( r'$\ell og_{10}\left(\frac{\rho}{g/cm^{3}}\right)$' , fontsize = fs )

ax_P = plt.subplot(312)
ax_P.semilogy( r , p   , 'k-'  )
ax_P.axvline( R_PNS , color = 'k' )
ax_P.axvline( Rs , color = 'k' , ls = '--' )

ax_P.set_xlim( xlim )
ax_P.get_xaxis().set_ticks( [] )
ax_P.set_ylabel( r'$\ell og_{10}\left(\frac{p}{dyne/cm^{2}}\right)$' , fontsize = fs )

ax_v = plt.subplot(313)
ax_v.plot( r , v , 'k-'  )
ax_v.axvline( R_PNS , color = 'k' )
ax_v.axvline( Rs , color = 'k' , ls = '--' )

ax_v.set_xlim( xlim )
ax_v.set_xlabel( r'$r\,\left[km\right]$'   , fontsize = fs )
ax_v.set_ylabel( r'$v\,\left[km/s\right]$' , fontsize = fs , labelpad = -2 )

plt.subplots_adjust( hspace = 0 )
#plt.show()
plt.savefig( 'SAS_IC.png' )
plt.close()
