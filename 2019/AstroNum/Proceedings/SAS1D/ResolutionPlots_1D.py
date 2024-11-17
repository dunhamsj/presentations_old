#!/usr/bin/env python3

import yt
import numpy as np
import subprocess
from os import listdir
from os.path import isfile
from sys import argv, exit
import matplotlib.pyplot as plt
plt.rc( 'font', family = 'serif' )

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

K = int( argv[1] )
ProblemDirectory \
  = HOME + 'Research/DataFromPastRuns/AstroNum2019/SAS1D_HLL_NoTCI_nNodes3_{:}/'.format( str(K).zfill(3) )
VariableToPlot = 'PF_D'
MakeDataFile = True
DataFileName = 'MovieData_{:}_{:}.dat'.format( VariableToPlot, str(K).zfill(3) )

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

            if( VariableToPlot == 'Entropy' ):
                PF_D  = CoveringGrid['PF_D' ].to_ndarray()[:,0,0]
                AF_P  = CoveringGrid['AF_P' ].to_ndarray()[:,0,0]
                AF_Gm = CoveringGrid['AF_Gm'].to_ndarray()[:,0,0]
                Data[i] = AF_P / PF_D**AF_Gm
            else:
                Data[i] = CoveringGrid[VariableToPlot].to_ndarray()[:,0,0]

        np.savetxt( DataFileName, Data )

Data = np.loadtxt( 'MovieData_{:}_{:}.dat'.format( \
         VariableToPlot, str(K).zfill(3) ) )

x = np.linspace( x.min(), x.max(), K )

# Plotting

fig, ax = plt.subplots()
plt.rcParams['mathtext.fontset'] = 'stix'

alpha = np.linspace( 0, 1, Data.shape[0] )
for i in range( Data.shape[0] ):
    if( i % 10 == 0 ):
        ax.semilogy( x, Data[i], 'k-', alpha = alpha[i] )
ax.semilogy( x, Data[0 ], 'r-',  label = 't = 0 ms'   )
ax.semilogy( x, Data[-1], 'r--', label = 't = 300 ms' )

ax.set_xlabel( 'Radial Coordinate [km]' )
ax.set_ylabel( r'$\rho\ \left[\mathrm{g\,cm^{-3}}\right]$', fontsize = 13 )

ax.set_xlim( 40, 540 )
ax.set_ylim( 4.0e6, 3.0e10 )
ax.text( 220, 1e10, 'nX = {:}'.format( K ), fontsize = 20 )

lw = 1.5
ax.spines['top'   ].set_linewidth( lw )
ax.spines['bottom'].set_linewidth( lw )
ax.spines['left'  ].set_linewidth( lw )
ax.spines['right' ].set_linewidth( lw )

ax.tick_params( which = 'both', axis = 'both', direction = 'in' )

ax.legend()

#plt.show()
plt.savefig( 'SAS1D_{:}.png'.format( str(K).zfill(3) ) )
