#!/usr/bin/python

import numpy as np 
import matplotlib
matplotlib.rcParams[ 'text.usetex' ] = False
import matplotlib.pyplot as plt 
import os 
from sys import exit 
np.seterr( divide = 'ignore' )

############# USER INPUT HERE
variable = 'PF_D'
xD       = 0.5
#############

names = { 'x'     : [ 'x'     , 0  ] ,
          'CF_D'  : [ 'CF_D'  , 1  ] ,
          'CF_S1' : [ 'CF_S1' , 2  ] ,
          'CF_S2' : [ 'CF_S2' , 3  ] ,
          'CF_S3' : [ 'CF_S3' , 4  ] ,
          'CF_E'  : [ 'CF_E'  , 5  ] ,
          'CF_Ne' : [ 'CF_Ne' , 6  ] ,
          'PF_D'  : [ 'PF_D'  , 7  ] ,
          'PF_V1' : [ 'PF_V1' , 8  ] ,
          'PF_V2' : [ 'PF_V2' , 9  ] ,
          'PF_V3' : [ 'PF_V3' , 10 ] ,
          'PF_E'  : [ 'PF_E'  , 11 ] ,
          'PF_Ne' : [ 'PF_Ne' , 12 ] ,
          'AF_P'  : [ 'AF_P'  , 13 ] ,
          'AF_T'  : [ 'AF_T'  , 14 ] ,
          'AF_Ye' : [ 'AF_Ye' , 15 ] ,
          'AF_S'  : [ 'AF_S'  , 16 ] ,
          'AF_E'  : [ 'AF_E'  , 17 ] ,
          'AF_Me' : [ 'AF_Me' , 18 ] ,
          'AF_Mp' : [ 'AF_Mp' , 19 ] ,
          'AF_Mn' : [ 'AF_Mn' , 20 ] ,
          'AF_Gm' : [ 'AF_Gm' , 21 ] ,
          'AF_Cs' : [ 'AF_Cs' , 22 ] }

# Get the number of degrees of freedom
nX = 400

x , p , rho , v , e = np.loadtxt( 'solution.dat' , unpack = True )

# Plotting
fig = plt.figure( num = variable , figsize = ( 11 , 2 ) )
plt.suptitle( variable, fontsize = 20 )

# Plot solution
plt.plot( x , rho , 'k-' , label='final'   )
plt.legend()
plt.ylabel( variable, fontsize = 15 )

# Plot discontinuity
#plt.axvline( x = xD, linestyle = '-', color = 'k' )

plt.xlim( np.min( x ), np.max( x ) )

plt.show()
