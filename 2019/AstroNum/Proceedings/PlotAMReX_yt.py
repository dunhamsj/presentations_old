#!/usr/bin/env python3

import yt
import numpy as np
import subprocess
from os import listdir
from os.path import isfile
from sys import argv, exit
import matplotlib.pyplot as plt

"""
To-Do:
  - Make PNS circular
    - Fix plot-window limits
  - Add 2D movie-making capability
    - Implement movie-making capability for curvilinear coordinates
    - Allow custom variables, i.e. Entropy, to be made into movies
"""

yt.funcs.mylog.setLevel(0) # Suppress initial yt output to screen

# --- Get user's HOME directory ---
HOME = subprocess.check_output( ["echo $HOME"], shell = True)
HOME = HOME[:-1].decode( "utf-8" ) + '/'

# Specify directory containing plotfiles

K = 128
ProblemDirectory \
  = HOME + 'Research/DataFromPastRuns/AstroNum2019/SAS1D_nNodes2_HLL_NoTCI_{:}/'.format( str(K).zfill(4) )
ProblemName = 'SAS_R'
VariableToPlot = 'Entropy'
MakeDataFile = False
DataFileName = 'MovieData_{:}.dat'.format( str(K).zfill(4) )

UsePhysicalUnits = False
CoordinateSystem = 'spherical'

# Get last plotfile in directory
FileArray \
  = np.sort(np.array( [ file for file in listdir( ProblemDirectory ) ] ))
FileList = []
for iFile in range( FileArray.shape[0] ):
    sFile = FileArray[iFile]
    if( sFile[0:8] == 'thornado' ):
        FileList.append( sFile )
FileArray = np.array( FileList )
File = FileArray[-1]

# Remove "/" at end of filename, if present
if ( File[-1] == '/' ): File = File[:-1]

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

# XXX.to_ndarray() strips array of yt units

DataUnit = ''
if  ( VariableToPlot == 'PF_D'  ):
    Data = CoveringGrid['PF_D'].to_ndarray()
    if( UsePhysicalUnits ):
        DataUnit = 'g/cm**3'

x = np.linspace( xL[0].to_ndarray(), xH[0].to_ndarray(), nX[0] )

if( MakeDataFile ):

    Overwrite = True
    if( isfile( DataFileName ) ):
        Overwrite = input( 'File: "{:}" exists. overwrite? (Y/N): '.format \
                      ( DataFileName ) )
        if( not Overwrite == 'Y' ):
            print( 'Not overwriting file, using existing file for movie.' )
            Overwrite = False
        else:
            Overwrite = True

    if( Overwrite ):
        # Put all time-slices into one array to use for movie making
        Data = np.empty( (FileArray.shape[0],nX[0]), float )
        for i in range( FileArray.shape[0] ):
            print( '{:}/{:}'.format( i+1, FileArray.shape[0] ) )
            ds = yt.load( '{:}'.format( ProblemDirectory + FileArray[i] ) )

            CoveringGrid \
              = ds.covering_grid \
                  ( level           = MaxLevel, \
                    left_edge       = xL, \
                    dims            = nX * 2**MaxLevel, \
                    num_ghost_zones = nX[0] )

            Data[i] = CoveringGrid[VariableToPlot].to_ndarray()[:,0,0]

        np.savetxt( DataFileName, Data )

Data = np.loadtxt( 'MovieData_{:}.dat'.format( str(K).zfill(4) ) )
x = np.linspace( 40, 540, K )

'''
alpha = np.linspace( 0, 1, Data.shape[0] )
for i in range( Data.shape[0] ):
    if( i % 10 == 0 ):
        plt.semilogy( x, Data[i], 'k-', alpha = alpha[i] )
plt.semilogy( x, Data[0],  'r-',  label = 't = 0 ms' )
plt.semilogy( x, Data[-1], 'r--', label = 't = 300 ms' )

plt.xlabel( 'Radial Distance [km]' )
plt.ylabel( r'$\rho\ \left[g\,cm^{-3}\right]$' )

plt.xlim( 40, 540 )
plt.ylim( 4.0e6, 3.0e10 )
plt.text( 220, 1e10, 'nX = {:}'.format( K ), fontsize = 20 )

plt.legend()
#plt.show()
plt.savefig( 'SAS1D_{:}.png'.format( str(K).zfill(4) ) )
'''

h1 = CoveringGrid['GF_h_1'].to_ndarray()
h1Unit = ''
h2 = CoveringGrid['GF_h_2'].to_ndarray()
h2Unit = 'm'
h3 = CoveringGrid['GF_h_3'].to_ndarray()
h3Unit = 'm'

r     = np.linspace( 40.0e3, 540.0e3, nX[0] )
theta = np.linspace( 0.0, np.pi, nX[1] )

dr     = ( r    [-1] - r    [0] ) / nX[0]
dtheta = ( theta[-1] - theta[0] ) / nX[1]

print( dr )
