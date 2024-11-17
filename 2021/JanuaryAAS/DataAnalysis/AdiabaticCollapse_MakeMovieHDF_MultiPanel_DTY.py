#!/usr/bin/env python3

import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
plt.style.use( './Publication.sty' )

from AdiabaticCollapse_ReadFieldsHDF import ReadFields

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

RootPath1 = HOME + 'Desktop/JanuaryAAS/Data/AdiabaticCollapse_CFA/'
RootPath2 = HOME + 'Desktop/JanuaryAAS/Data/AdiabaticCollapse_Newtonian/'
suffix = ''

Problem = 'GravitationalCollapse'

Fields             = [ 'AF_P', 'AF_T', 'AF_Ye' ]
TimeUnit          = 'ms'
LengthUnit        = 'km'
Dimension         = ['X1']
UseSemiLogXScale  = True
UseSemiLogYScale  = True
UseGeometryFields = True
iSSi = 2000
iSSf = 2417

FigTitle = 'Adiabatic Collapse'

N  = 2
nX = 256

############################

nSS = iSSf-iSSi+1

if N == 2:

    w = np.array( [0.5,0.5], np.float64 )

elif N == 3:

    w = np.array( [5.0/18.0,8.0/18.0,5.0/18.0], np.float64 )

else:

    exit( 'Invalid value for N. Exiting...' )

Snapshots = np.linspace( iSSi, iSSf, nSS, dtype = np.int64 )

# --- Define where to look for data ---
PathToData1 = RootPath1 + suffix + Problem
PathToData2 = RootPath2 + suffix + Problem

print( 'Loading in data...\n' )

Names1 = ReadFields( PathToData1, Snapshots, UseGeometryFields )
Names2 = ReadFields( PathToData2, Snapshots, UseGeometryFields )
Time  = Names1['Time']

x    = np.array( Names1[Dimension[0]] )
xMin = np.min(x)
xMax = np.max(x)
x1   = Names1['X1C']
x2   = Names2['X1C']

AF_P1  = Names1[Fields[0]]
AF_T1  = Names1[Fields[1]]
AF_Ye1 = Names1[Fields[2]]
GF_Sg1 = Names1['GF_Sg']

AF_P2  = Names2[Fields[0]]
AF_T2  = Names2[Fields[1]]
AF_Ye2 = Names2[Fields[2]]
GF_Sg2 = Names2['GF_Sg']

P1_K = np.empty( (nSS,nX), np.float64 )
T1_K = np.empty( (nSS,nX), np.float64 )
Y1_K = np.empty( (nSS,nX), np.float64 )

P2_K = np.empty( (nSS,nX), np.float64 )
T2_K = np.empty( (nSS,nX), np.float64 )
Y2_K = np.empty( (nSS,nX), np.float64 )

# Compute cell averages
for t in range( nSS ):

    for i in range( nX ):

        d = np.linspace( i*N, (i+1)*N-1, N, dtype = int )

        V1 = np.sum( w * GF_Sg1[t][0,0,d] )
        V2 = np.sum( w * GF_Sg2[t][0,0,d] )

        P1_K[t][i] = np.sum( w * GF_Sg1[t][0,0,d] * AF_P1 [t][0,0,d] ) / V1
        T1_K[t][i] = np.sum( w * GF_Sg1[t][0,0,d] * AF_T1 [t][0,0,d] ) / V1
        Y1_K[t][i] = np.sum( w * GF_Sg1[t][0,0,d] * AF_Ye1[t][0,0,d] ) / V1

        P2_K[t][i] = np.sum( w * GF_Sg2[t][0,0,d] * AF_P2 [t][0,0,d] ) / V2
        T2_K[t][i] = np.sum( w * GF_Sg2[t][0,0,d] * AF_T2 [t][0,0,d] ) / V2
        Y2_K[t][i] = np.sum( w * GF_Sg2[t][0,0,d] * AF_Ye2[t][0,0,d] ) / V2

################ Plotting information

# Animation program adapted from
# https://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/
fig, axs = plt.subplots( 3 )
plt.subplots_adjust( hspace = 0.0, wspace = 0.0 )
fig.suptitle( FigTitle, fontsize = 20 )

xticks = np.logspace( 0, 3, 4 )
xMin = 0.25
xMax = 8.0e3

yMin = +np.inf
yMax = -np.inf
for t in range( nSS ):
    yMin = min( yMin, P1_K[t].min(), P2_K[t].min() )
    yMax = max( yMax, P1_K[t].max(), P2_K[t].max() )
axs[0].set_xlim( xMin, xMax )
axs[0].set_ylim( yMin, yMax )
axs[0].set_ylabel( r'$\mathrm{P [erg\,cm^{-3}]}$' )
axs[0].set_xscale( 'log' )
axs[0].set_yscale( 'log' )
axs[0].xaxis.set_visible(False)

yMin = +np.inf
yMax = -np.inf
for t in range( nSS ):
    yMin = min( yMin, T1_K[t].min(), T2_K[t].min() )
    yMax = max( yMax, T1_K[t].max(), T2_K[t].max() )
axs[1].set_xlim( xMin, xMax )
axs[1].set_ylim( yMin, yMax )
axs[1].set_ylabel( r'$\mathrm{T [K]}$' )
axs[1].xaxis.set_visible(False)
axs[1].set_xscale( 'log' )
axs[1].set_yscale( 'log' )

yMin = +np.inf
yMax = -np.inf
for t in range( nSS ):
    yMin = min( yMin, Y1_K[t].min(), Y2_K[t].min() )
    yMax = max( yMax, Y1_K[t].max(), Y2_K[t].max() )
axs[2].set_xlim( xMin, xMax )
axs[2].set_ylim( yMin, yMax )
axs[2].set_ylabel( r'$\mathrm{Y_{e}}$' )
axs[2].set_xscale( 'log' )
axs[2].set_xticks( xticks )
axs[2].set_xlabel( r'$\mathrm{Radial\ Coordinate\ [km]}$' )

line11, = axs[0].plot( [], [], 'k-', label = 'GR-CFA' )
line21, = axs[1].plot( [], [], 'k-' )
line31, = axs[2].plot( [], [], 'k-' )
line12, = axs[0].plot( [], [], 'k--', label = 'Newtonian' )
line22, = axs[1].plot( [], [], 'k--' )
line32, = axs[2].plot( [], [], 'k--' )

time_text = axs[0].text( 1.0e1, 5.0e34, '' )

def InitializeFrame():
    line11.set_data([],[])
    line21.set_data([],[])
    line31.set_data([],[])
    line12.set_data([],[])
    line22.set_data([],[])
    line32.set_data([],[])
    time_text.set_text('')
    return line11, line21, line31, line12, line22, line32, time_text

# Animation function
def UpdateFrame(t):
    if( t % 10 == 0 ):
        print( '{:d}/{:d}'.format( t, nSS ) )

    line11.set_data( x1, P1_K[t] )
    line21.set_data( x1, T1_K[t] )
    line31.set_data( x1, Y1_K[t] )

    line12.set_data( x2, P2_K[t] )
    line22.set_data( x2, T2_K[t] )
    line32.set_data( x2, Y2_K[t] )

    time_text.set_text( 'Time = {:.3e} ms'.format( Time[t] ) )
    return line11, line21, line31, line12, line22, line32, time_text

axs[0].legend()

# Call the animator
anim = animation.FuncAnimation \
         ( fig, \
           UpdateFrame, \
           init_func = InitializeFrame, \
           frames    = nSS, \
           blit      = True )

SaveFileAs = '../Movies/AdiabaticCollapse_MultiPanel_DTY_NewtonianVsGR.mp4'

anim.save( SaveFileAs, dpi = 300, fps = int( nSS / 10.0 ) )

os.system( 'rm -f *.pyc' )
os.system( 'rm -rf __pycache__' )
