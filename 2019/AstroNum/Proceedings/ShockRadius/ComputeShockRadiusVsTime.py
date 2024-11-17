#!/usr/bin/env python3

import yt
import numpy as np
import subprocess
from os import listdir

yt.funcs.mylog.setLevel(0) # Suppress initial yt output to screen

# --- Get user's HOME directory ---
HOME = subprocess.check_output( ["echo $HOME"], shell = True)
HOME = HOME[:-1].decode( "utf-8" ) + '/'

# Specify directory containing plotfiles

ProblemDirectory \
  = HOME + 'Research/DataFromPastRuns/AstroNum2019/SASI/'
MakeDataFile              = True
ComputeShockRadius        = True
PlotShockRadiusVsTime     = True
DataFileName              = 'Entropy' + '.dat'
TimeFileName              = 'Time' + '.dat'
ShockRadiusVsTimeFileName = 'ShockRadiusVsTime' + '.dat'
OutputFile                = 'Output' + '.dat'
SaveFigAs                 = 'ShockRadiusVsTime' + '.png'

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

with open( OutputFile, 'a+' ) as stdout:
    stdout.write( '# Output.dat\n' )

if( MakeDataFile ):
    # Put all time-slices into one array to use for movie making
    Data = np.empty( (FileArray.shape[0],nX[0],nX[1]), float )
    Time = np.empty( FileArray.shape[0], float )
    with open( OutputFile, 'a+' ) as stdout:
        stdout.write( '\nGenerating file {:}...\n'.format( DataFileName ) )
    for i in range( FileArray.shape[0] ):
        with open( OutputFile, 'a+' ) as stdout:
            stdout.write( '{:}/{:}\n'.format( i+1, FileArray.shape[0] ) )
        ds = yt.load( '{:}'.format( ProblemDirectory + FileArray[i] ) )

        CoveringGrid \
          = ds.covering_grid \
              ( level           = MaxLevel, \
                left_edge       = xL, \
                dims            = nX * 2**MaxLevel, \
                num_ghost_zones = nX[0] )

        PF_D  = CoveringGrid['PF_D' ].to_ndarray()[:,:,0]
        AF_P  = CoveringGrid['AF_P' ].to_ndarray()[:,:,0]
        AF_Gm = CoveringGrid['AF_Gm'].to_ndarray()[:,:,0]

        Data[i] = AF_P / PF_D**AF_Gm
        Time[i] = ds.current_time

    np.savetxt( TimeFileName, Time )

    # Save multi-D array with np.savetxt. Taken from:
    # https://stackoverflow.com/questions/3685265/how-to-write-a-multidimensional-array-to-a-text-file
    with open( DataFileName, 'w' ) as FileOut:
        FileOut.write( '# Array shape: {0}\n'.format( Data.shape ) )

        # Iterating through an n-dimensional array produces slices along
        # the last axis. This is equivalent to Data[i] in this case
        for TimeSlice in Data:
            np.savetxt( FileOut, TimeSlice )
            FileOut.write( '# New slice\n' )

if( ComputeShockRadius ):
    with open( OutputFile, 'a+' ) as stdout:
        stdout.write( 'Reading in the data file...\n' )
    Data = np.loadtxt( DataFileName ).reshape( \
             (FileArray.shape[0],nX[0],nX[1]) )
    Time = np.loadtxt( TimeFileName )

    DataUnit = 'erg/cm**3/(g/cm**3)**(4/3)'

    xL = xL.to_ndarray()
    xH = xH.to_ndarray()

    dX1 = ( xH[0] - xL[0] ) / nX[0]
    dX2 = ( xH[1] - xL[1] ) / nX[1]

    X1 = np.linspace( xL[0] + dX1, xH[0] - dX1, nX[0] )
    X2 = np.linspace( xL[1] + dX2, xH[1] + dX2, nX[1] )

    EntropyThreshold = 3.0e15

    Volume = np.zeros( FileArray.shape[0], float )

    with open( OutputFile, 'a+' ) as stdout:
        stdout.write( 'Computing volumes...\n' )

    for iT in range( Data.shape[0] ):
        with open( OutputFile, 'a+' ) as stdout:
            stdout.write( '{:}/{:}\n'.format( iT+1, Data.shape[0] ) )
        for iX1 in range( Data.shape[1] ):
            for iX2 in range( Data.shape[2] ):
                if( Data[iT,iX1,iX2] > EntropyThreshold ):
                    Volume[iT] \
                      += ( ( X1[iX1] + dX1 )**3 - X1[iX1]**3 ) \
                           * ( np.cos( X2[iX2] ) - np.cos( X2[iX2] + dX2 ) )

    ShockRadius = ( Volume / 2.0 )**( 1.0 / 3.0 )
    np.savetxt( ShockRadiusVsTimeFileName, np.vstack( ( Time, ShockRadius ) ) )

if( PlotShockRadiusVsTime ):
    import matplotlib.pyplot as plt
    plt.rc( 'font', family = 'serif' )
    t, r = np.loadtxt( ShockRadiusVsTimeFileName )
    plt.plot( t, r, 'k-' )
    plt.xlim( t[0], t[-1] )
    plt.xlabel( 'Time [ms]' )
    plt.ylabel( 'Average Shock Radius [km]' )
    plt.show()
    #plt.savefig( SaveFigAs )
    plt.close()
