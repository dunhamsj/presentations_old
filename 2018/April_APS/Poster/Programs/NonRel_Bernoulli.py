#!/usr/bin/python

import numpy as np
from sys import argv
np.seterr( all='ignore' )
##### Physical constants

# Adiabatic constant
gamma = 4.0 / 3.0 # [ dimensionless ]

Physical = True

if ( Physical == True ):
    # Real numbers
    c     = 2.99792458e8             # [ m / s ]
    G     = 6.673e-11                # [ m^3 / ( kg * s^2 ) ]
    Msun  = 1.99892e30               # [ kg ]
    M     = 1.4 * Msun               # [ kg ]
    Rs    = float( argv[1] ) * 1.0e3 # [ m ]
    Mdot  = float( argv[2] ) * Msun  # [ kg / s ]
    R_PNS = 40.0 * 1.0e3             # [ m ]
    R_sc  = 2.0 * G * M / c**(2.0)   # [ m ]
    Ri    = float( argv[3] ) * R_sc  # [ m ]
    Rf    = 10.0 * Rs

    print 'M/Msun:   ' , M / Msun
    print 'Ri/km:    ' , Ri / 1.0e3
    print 'R_PNS/km: ' , R_PNS / 1.0e3
    print 'Rs/km:    ' , Rs / 1.0e3

else:

    # Code Units (M~1.4Msun, Rs~180km, R_PNS~40km)
    c     = 1.0
    G     = 1.0
    Msun  = 1.0
    M     = 0.1
    Rs    = 1.0
    Mdot  = 4.0 * np.pi # ( Mdot at shock radius )
    R_PNS = 2.0 / 9.0
    R_sc  = 2.0 * G * M / c**(2.0) # Schwarzschild radius
    Ri    = 1.0 * R_sc
    Rf    = 10.0 * Rs

    print 'M:    ' , M / Msun
    print 'Ri:   ' , Ri
    print 'R_PNS:' , R_PNS
    print 'Rs:   ' , Rs

# Alpha from Blondin et. al., 2003, equation 6
Alpha = 4.0 * gamma / ( ( gamma - 1.0 ) * ( gamma + 1.0 ) ) \
          * ( ( gamma - 1.0 ) / ( gamma + 1.0 ) )**gamma # [ dimensionless ]

# Number of points
N = 10000

# Radial coordinate
r_pre  = np.linspace( Rf , Rs , N )
r_post = np.linspace( Rs , Ri , N )

# Define the function whose root we want to find
def f( u , r ):

    return r * u**(2.0) + Alpha * ( 2.0 * G * M )**( ( gamma + 1.0 ) / 2.0 ) \
             * Rs**( ( 3.0 * gamma - 5.0 ) / 2.0 ) * r**( 3.0 - 2.0 * gamma ) \
             *  u**( 1.0 - gamma ) - 2.0 * G * M

# Define the derivative of the function
def dfdu( u , r ):

    return 2.0 * r * u + Alpha * ( 2.0 * G * M )**( ( gamma + 1.0 ) / 2.0 ) \
             * Rs**( ( 3.0 * gamma - 5.0 ) / 2.0 ) * r**( 3.0 - 2.0 * gamma ) \
             * ( 1.0 - gamma ) * u**( -gamma )

# Initial value for u from RH conditions
u_Rs = -( gamma - 1.0 ) / ( gamma + 1.0 ) * np.sqrt( 2.0 * G * M / Rs )

u_post    = np.empty( N , float )
u_post[0] = np.abs( u_Rs )

def Bernoulli( r , i ):

    # Tolerance
    epsilon = 1.e-16
    Delta_u = 1.

    # Initial guess
    ui = u_post[i]

    # Max number of iterations
    nMax  = 100
    nIter = 1
    while ( ( np.abs( Delta_u ) > epsilon ) and ( nIter < nMax ) ):
        Delta_u = -f( ui , r ) / dfdu( ui , r )
        ui      = Delta_u + ui

        nIter += 1

    return ui

for i in range( N - 1 ):

    u_post[i+1] = np.abs( Bernoulli( r_post[i] , i ) )

# K constant from P = K * rho^gamma
K = 2.0 / ( gamma + 1.0 ) * ( ( gamma - 1.0 ) / ( gamma + 1.0 ) )**( gamma ) \
    * ( Mdot / ( 4.0 * np.pi ) )**( 1.0 - gamma ) \
    * ( 2.0 * G * M )**( ( gamma + 1.0 ) / 2.0 ) \
    * Rs**( ( 3.0 * gamma - 5.0 ) / 2.0 )

# Post-shock fluid variables
rho_post = ( ( gamma - 1.0 ) / ( 2.0 * gamma * K ) \
         * ( 2.0 * G * M / r_post - u_post**(2.0) ) )**( 1.0 / ( gamma - 1.0 ) )
p_post   = K * rho_post**( gamma )

# Pre-shock fluid variables
u_pre   = -np.sqrt( 2.0 * G * M / r_pre )
rho_pre = -Mdot / ( 4.0 * np.pi * r_pre**(2.0) * u_pre )

# Make copies of fluid variables for relativistic program
u_pre_N  = np.copy( u_pre )
u_post_N = np.copy( u_post )

rho_pre_N  = np.copy( rho_pre )
rho_post_N = np.copy( rho_post )

p_post_N = np.copy( p_post )

'''
import matplotlib.pyplot as plt

##### Plotting

### Density
plt.semilogy( r_pre / Rs ,  rho_pre * c**(2.0) ,  'r-'  , \
              label = r'$\ell og_{10}\left(\rho\,c^2\right)$' )
plt.semilogy( r_post / Rs , rho_post * c**(2.0) , 'r-' )

### Pressure
plt.semilogy( r_post / Rs , p_post , 'g-' , \
              label = r'$\ell og_{10}\left(p\right)$' )

### Velocity
#plt.semilogy( r_pre / Rs ,  np.abs( u_pre )   , 'b-' ,  \
#              label = r'$|u|$' )
#plt.semilogy( r_post / Rs , np.abs( u_post )  , 'b-' ,  \
#              label = r'$|u|$' )

plt.xlabel( r'$r/R_s$' )

plt.title( 'Fluid Variables' )
plt.legend()

plt.xlim( Ri / Rs , Rf / Rs )

plt.show()
'''

del u_pre , u_post , rho_pre , rho_post , p_post , f , dfdu , Bernoulli
