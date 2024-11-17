#!/usr/bin/env python3

import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
plt.style.use( './Publication.sty' )

from ReadFieldsYahil import ReadFields

print( '' )
print( 'Running MakeMovieHDF.py...' )
print( '---------------------------' )

# --- Get user's HOME directory ---
HOME = subprocess.check_output( ["echo $HOME"], shell = True)
HOME = HOME[:-1].decode( "utf-8" ) + '/'

# --- Get user's THORNADO_DIR directory ---
THORNADO_DIR = subprocess.check_output( ["echo $THORNADO_DIR"], shell = True)
THORNADO_DIR = THORNADO_DIR[:-1].decode( "utf-8" ) + '/'

############################ User Input ############################

# --- Define root path for data ---

RootPath1 = HOME + 'Desktop/JanuaryAAS/Data/Yahil_CFA/'
RootPath2 = HOME + 'Desktop/JanuaryAAS/Data/Yahil_Newtonian/'
suffix = ''

Problem1 = 'YahilCollapse'
Problem2 = 'YahilLattimerCollapse'

Fields             = [ 'PF_D', 'PF_V1' ]
TimeUnit          = 'ms'
LengthUnit        = 'km'
Dimension         = ['X1']
UseSemiLogXScale  = True
UseSemiLogYScale  = True
UseGeometryFields = True
iSSi = 1300
iSSf = 1468

FigTitle = 'Self-Similar Collapse'

N  = 3
nX = 256

############################

nSS = iSSf-iSSi+1

if N == 2:

    w = np.array( [0.5,0.5], np.float64 )

elif N == 3:

    w = np.array( [5.0/18.0,8.0/18.0,5.0/18.0], np.float64 )

else:

    exit( 'Invalid value for N. Exiting...' )

c = 2.99792458e5

Snapshots = np.linspace( iSSi, iSSf, nSS, dtype = np.int64 )

# --- Define where to look for data ---
PathToData1 = RootPath1 + suffix + Problem1
PathToData2 = RootPath2 + suffix + Problem2

Names1 = ReadFields( PathToData1, Snapshots, UseGeometryFields )
Names2 = ReadFields( PathToData2, Snapshots, UseGeometryFields )
Time  = Names1['Time']

x    = np.array( Names1[Dimension[0]] )
xMin = np.min(x)
xMax = np.max(x)
x1    = Names1['X1C']
x2    = Names2['X1C']

PF_D1  = Names1[Fields[0]]
PF_V11 = Names1[Fields[1]]
GF_Sg1 = Names1['GF_Sg']

PF_D2  = Names2[Fields[0]]
PF_V12 = Names2[Fields[1]]
GF_Sg2 = Names2['GF_Sg']

D1_K = np.empty( (nSS,nX), np.float64 )
D2_K = np.empty( (nSS,nX), np.float64 )

V1_K = np.empty( (nSS,nX), np.float64 )
V2_K = np.empty( (nSS,nX), np.float64 )

################ Plotting information

# Animation program adapted from
# https://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/
fig, axs = plt.subplots( 2, 1 )
plt.subplots_adjust( hspace = 0.0, wspace = 0.0 )
fig.suptitle( FigTitle, fontsize = 20 )

xticks = np.logspace( 0, 4, 5 )
xMin = 0.25
xMax = 1.0e4

# Compute cell averages
for t in range( nSS ):

    for i in range( nX ):

        d = np.linspace( i*N, (i+1)*N-1, N, dtype = int )

        V1 = np.sum( w * GF_Sg1[t][0,0,d] )
        V2 = np.sum( w * GF_Sg2[t][0,0,d] )

        D1_K[t][i] = np.sum( w * GF_Sg1[t][0,0,d] * PF_D1 [t][0,0,d] ) / V1
        V1_K[t][i] = np.sum( w * GF_Sg1[t][0,0,d] * PF_V11[t][0,0,d] ) / V1

        D2_K[t][i] = np.sum( w * GF_Sg2[t][0,0,d] * PF_D2 [t][0,0,d] ) / V2
        V2_K[t][i] = np.sum( w * GF_Sg2[t][0,0,d] * PF_V12[t][0,0,d] ) / V2

yMin = +np.inf
yMax = -np.inf
for t in range( nSS ):
    yMin = min( yMin, D1_K[t].min(), D2_K[t].min() )
    yMax = max( yMax, D1_K[t].max(), D2_K[t].max() )
axs[0].set_xlim( xMin, xMax )
axs[0].set_ylim( yMin, yMax )
axs[0].set_ylabel( r'$\mathrm{\rho\ [g\,cm^{-3}]}$' )
axs[0].set_yscale( 'log' )
axs[0].set_xscale( 'log' )
axs[0].xaxis.set_visible(False)

yMin = +np.inf
yMax = -np.inf
for t in range( nSS ):
    yMin = min( yMin, V1_K[t].min() / c, V2_K[t].min() / c )
    yMax = max( yMax, V1_K[t].max() / c, V2_K[t].max() / c )
axs[1].set_ylabel( r'$v/c$' )
axs[1].set_xlim( xMin, xMax )
axs[1].set_ylim( yMin, yMax )
axs[1].set_xticks( xticks )
axs[1].set_xscale( 'log' )

line11, = axs[0].plot( [], [], 'k-' , label = 'GR-CFA' )
line21, = axs[1].plot( [], [], 'k-' )
line12, = axs[0].plot( [], [], 'k--', label = 'Newtonian' )
line22, = axs[1].plot( [], [], 'k--' )

time_text = axs[0].text( 1.0e1, 2.0e15, '' )

def InitializeFrame():
    line11.set_data([],[])
    line21.set_data([],[])
    line12.set_data([],[])
    line22.set_data([],[])
    time_text.set_text('')
    return line11, line21, line12, line22, time_text

# Animation function
def UpdateFrame(t):

    line11.set_data( x1, D1_K[t] )

    line21.set_data( x1, V1_K[t]/c )

    line12.set_data( x2, D2_K[t] )

    line22.set_data( x2, V2_K[t]/c )

    time_text.set_text( 'Time = {:.3e} ms'.format( Time[t] ) )
    return line11, line21, line12, line22, time_text

axs[0].legend( loc = 1, prop = {'size':8} )
axs[1].set_xlabel( r'$\mathrm{Radial\ Coordinate\ [km]}$' )

# Call the animator
anim = animation.FuncAnimation \
         ( fig, \
           UpdateFrame, \
           init_func = InitializeFrame, \
           frames    = nSS, \
           blit      = True )

SaveFileAs = '../Movies/{:}_MultiPanel_DV_NewtonianVsGR.mp4'.format( Problem1 )

anim.save( SaveFileAs, dpi = 300, fps = int( nSS / 10.0 ) )

os.system( 'rm -f *.pyc' )
os.system( 'rm -rf __pycache__' )
