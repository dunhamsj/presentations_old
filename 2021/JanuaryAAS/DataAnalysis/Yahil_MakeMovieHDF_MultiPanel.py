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

RootPath = HOME + 'Desktop/JanuaryAAS/Data/Yahil_CFA/'
suffix = ''

Problem = 'YahilCollapse'

Fields             = [ 'PF_D', 'PF_V1', 'GF_CF', 'GF_al', 'GF_b1' ]
TimeUnit          = 'ms'
LengthUnit        = 'km'
Dimension         = ['X1']
UseSemiLogXScale  = True
UseSemiLogYScale  = True
UseGeometryFields = True
iSSi = 1300
iSSf = 1468

FigTitle = 'GR-CFA Self-Similar Collapse'

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
PathToData = RootPath + suffix + Problem

Names = ReadFields( PathToData, Snapshots, UseGeometryFields )
Time  = Names['Time']

x    = np.array( Names[Dimension[0]] )
xMin = np.min(x)
xMax = np.max(x)
x    = Names['X1C']

PF_D  = Names[Fields[0]]
PF_V1 = Names[Fields[1]]
GF_CF = Names[Fields[2]]
GF_al = Names[Fields[3]]
GF_b1 = Names[Fields[4]]
GF_Sg = Names['GF_Sg']

D_K = np.empty( (nSS,nX), np.float64 )
V_K = np.empty( (nSS,nX), np.float64 )
C_K = np.empty( (nSS,nX), np.float64 )
A_K = np.empty( (nSS,nX), np.float64 )
B_K = np.empty( (nSS,nX), np.float64 )

# Compute cell averages
for t in range( nSS ):

    for i in range( nX ):

        d = np.linspace( i*N, (i+1)*N-1, N, dtype = int )

        V = np.sum( w * GF_Sg[t][0,0,d] )

        D_K[t][i] = np.sum( w * GF_Sg[t][0,0,d] * PF_D [t][0,0,d] ) / V
        V_K[t][i] = np.sum( w * GF_Sg[t][0,0,d] * PF_V1[t][0,0,d] ) / V
        C_K[t][i] = np.sum( w * GF_Sg[t][0,0,d] * GF_CF[t][0,0,d] ) / V
        A_K[t][i] = np.sum( w * GF_Sg[t][0,0,d] * GF_al[t][0,0,d] ) / V
        B_K[t][i] = np.sum( w * GF_Sg[t][0,0,d] * GF_b1[t][0,0,d] ) / V

################ Plotting information

# Animation program adapted from
# https://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/
fig, axs = plt.subplots( 2, 2 )
plt.subplots_adjust( hspace = 0.0, wspace = 0.0 )
fig.suptitle( FigTitle, fontsize = 20 )

xticks = np.logspace( 0, 4, 5 )
xMin = 0.25
xMax = 1.0e4

yMin = +np.inf
yMax = -np.inf
for t in range( nSS ):
    yMin = min( yMin, D_K[t].min() )
    yMax = max( yMax, D_K[t].max() )
axs[0,0].set_xlim( xMin, xMax )
axs[0,0].set_ylim( yMin, yMax )
axs[0,0].set_ylabel( r'$\mathrm{\rho\ [g\,cm^{-3}]}$' )
axs[0,0].set_yscale( 'log' )
axs[0,0].set_xscale( 'log' )
axs[0,0].xaxis.set_visible(False)

yMin = +np.inf
yMax = -np.inf
for t in range( nSS ):
    yMin = min( yMin, V_K[t].min() / c )
    yMax = max( yMax, V_K[t].max() / c )
axs[0,1].set_ylabel( r'$v/c$' )
axs[0,1].set_xlim( xMin, xMax )
axs[0,1].set_ylim( yMin, yMax )
axs[0,1].yaxis.set_label_position('right')
axs[0,1].yaxis.tick_right()
axs[0,1].xaxis.set_visible(False)
axs[0,1].set_xscale( 'log' )

yMin = +np.inf
yMax = -np.inf
for t in range( nSS ):
    yMin = min( yMin, A_K[t].min(), C_K[t].min() )
    yMax = max( yMax, A_K[t].max(), C_K[t].max() )
axs[1,0].set_xlim( xMin, xMax )
axs[1,0].set_ylim( yMin, yMax )
axs[1,0].set_xlabel( r'$\mathrm{Radial\ Coordinate\ [km]}$' )
axs[1,0].set_xscale( 'log' )
axs[1,0].set_xticks( xticks )

yMin = +np.inf
yMax = -np.inf
for t in range( nSS ):
    yMin = min( yMin, B_K[t].min() / c )
    yMax = max( yMax, B_K[t].max() / c )
axs[1,1].set_xlim( xMin, xMax )
axs[1,1].set_ylim( yMin, yMax )
axs[1,1].set_xscale( 'log' )
axs[1,1].set_ylabel( r'$\beta/c$' )
axs[1,1].yaxis.set_label_position('right')
axs[1,1].yaxis.tick_right()
axs[1,1].set_xlabel( r'$\mathrm{Radial\ Coordinate\ [km]}$' )
axs[1,1].set_xticks( xticks )

line1, = axs[0,0].plot( [], [], 'k-' )
line2, = axs[0,1].plot( [], [], 'k-' )
line3, = axs[1,0].plot( [], [], 'k-' , label = 'Conformal Factor' )
line4, = axs[1,0].plot( [], [], 'k--', label = 'Lapse Function' )
line5, = axs[1,1].plot( [], [], 'k-' )

time_text = axs[0,0].text( 8.0e2, 2.0e15, '' )

def InitializeFrame():
    line1.set_data([],[])
    line2.set_data([],[])
    line3.set_data([],[])
    line4.set_data([],[])
    line5.set_data([],[])
    time_text.set_text('')
    return line1, line2, line3, line4, line5, time_text

# Animation function
def UpdateFrame(t):

    line1.set_data( x, D_K[t] )

    line2.set_data( x, V_K[t]/c )

    line3.set_data( x, C_K[t] )

    line4.set_data( x, A_K[t] )

    line5.set_data( x, B_K[t]/c )

    time_text.set_text( 'Time = {:.3e} ms'.format( Time[t] ) )
    return line1, line2, line3, line4, line5, time_text

axs[1,0].legend( loc = 1, prop = {'size':8} )

# Call the animator
anim = animation.FuncAnimation \
         ( fig, \
           UpdateFrame, \
           init_func = InitializeFrame, \
           frames    = nSS, \
           blit      = True )

SaveFileAs = '../Movies/{:}_MultiPanel.mp4'.format( Problem )

anim.save( SaveFileAs, dpi = 300, fps = int( nSS / 10.0 ) )

os.system( 'rm -f *.pyc' )
os.system( 'rm -rf __pycache__' )
