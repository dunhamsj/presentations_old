#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

# Define amplitude, wavenumber, range, and phase
AL, kL, xL, phiL = 1.0, 3.0 * np.pi, np.linspace( 0.0 , 1.0 , 100 ) , -0.5
AR, kR, xR, phiR = 1.0, 3.0 * np.pi, np.linspace( 1.0 , 2.0 , 100 ) , 0.5

ax = plt.figure( figsize = ( 10.25 , 2 ) )
plt.axis('off')

# Plot and label cell interface
plt.axvline( 1.0 , color = 'k' , linestyle = '--' )
plt.text( 1.05 , 1.0 , 'Cell interface' )

# Plot zero line
plt.axhline( 0.0 , color = 'k' )

# Plot nodal points
plt.plot( xL[25] , 0.0 , 'ko' )
plt.plot( xL[75] , 0.0 , 'ko' )
plt.plot( xR[25] , 0.0 , 'ko' )
plt.plot( xR[75] , 0.0 , 'ko' )
plt.text( xR[70] , -0.5 , 'Nodal\npoint' )

plt.text( 0.32 , 0.64 , 'Approximate\nsolution $u^{K}_{h}$' )
plt.text( 1.55 , 0.64 , 'Approximate\nsolution $u^{K+1}_{h}$' )

# Plot curves
plt.plot( xL , AL * np.sin( kL * xL + phiL ) , 'k-' )
plt.plot( xR , AR * np.sin( kR * xR + phiR ) , 'k-' )
plt.show()
