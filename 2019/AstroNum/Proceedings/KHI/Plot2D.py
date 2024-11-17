#!/usr/bin/env python3

import yt
import numpy as np
import subprocess
from os import listdir
from os.path import isfile
from sys import argv, exit
import matplotlib.pyplot as plt

yt.funcs.mylog.setLevel(0) # Suppress initial yt output to screen

# --- Get user's HOME directory ---
HOME = subprocess.check_output( ["echo $HOME"], shell = True)
HOME = HOME[:-1].decode( "utf-8" ) + '/'

# Specify directory containing plotfiles

ProblemDirectory = HOME \
  + 'Research/DataFromPastRuns/AstroNum2019/KHI_20190904/'
ProblemName      = 'KHI'
VariableToPlot   = 'PF_D'

File = 'thornado_00006012' # 1.38
#File = 'thornado_00007889' # 1.80
#File = 'thornado_00010641' # 2.34
#File = 'thornado_00013499' # 3.00

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

Data = CoveringGrid[VariableToPlot].to_ndarray()
DataUnit = ''

data   = { VariableToPlot: (Data,DataUnit) }
fields = VariableToPlot

slc = yt.SlicePlot( ds, 'z', fields, aspect = 0.5 )

slc.set_log( fields, False )

vmin = 0.0
vmax = 2.0
slc.set_zlim( fields, vmin, vmax )

# First argument is location in figure, [0,0] is lower-left,
#                                       [1,1] is upper-right
slc.annotate_text( [0.02, 0.95], 't = {:.2f}'.format( Time.to_ndarray() ), \
                   coord_system = 'axis', \
                   text_args={'color':'black'} )

slc.set_cmap( field = fields, cmap = 'jet' )
slc.set_colorbar_label( fields, 'Primitive Rest-Mass-Density' )
slc.set_xlabel( 'x' )
slc.set_ylabel( 'y' )

lw = 1.5


ax = slc.plots[VariableToPlot].axes

ax.xaxis.set_ticks([])
ax.xaxis.set_ticklabels([])

ax.spines['top'   ].set_linewidth( lw )
ax.spines['bottom'].set_linewidth( lw )
ax.spines['left'  ].set_linewidth( lw )
ax.spines['right' ].set_linewidth( lw )

Time = float( Time.to_ndarray() * 100 )

slc.save( ProblemName + '_' + VariableToPlot \
            + '_{:d}.png'.format( int( round( Time ) ) ) )
