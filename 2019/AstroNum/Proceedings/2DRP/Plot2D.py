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
  + 'Research/DataFromPastRuns/AstroNum2019/2DRP_CompLimit/'
ProblemName      = '2DRP'
VariableToPlot   = 'AF_P'

#File = 'thornado_00049144' # CharLimit
File = 'thornado_00049151' # CompLimit

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

slc = yt.SlicePlot( ds, 'z', fields, origin = 'lower-left-window' )

vmin = 1.0e-2
vmax = 5.0e1
slc.set_zlim( fields, vmin, vmax )

slc.set_cmap( field = fields, cmap = 'Purples' )
slc.set_colorbar_label( fields, 'Pressure' )
slc.set_xlabel( 'x' )
slc.set_ylabel( 'y' )

if( File == 'thornado_00049144' ):
    Annotate = 'Characteristic Limiting'
    SaveAs = '2DRP_CharLimit.png'
else:
    Annotate = 'Component-Wise Limiting'
    SaveAs = '2DRP_CompLimit.png'

slc.annotate_text( [0.02, 0.95], Annotate, \
                   coord_system = 'axis', \
                   text_args={'color':'black','size':30} )

lw = 1.5

ax = slc.plots[VariableToPlot].axes

ax.xaxis.set_ticks([])
ax.xaxis.set_ticklabels([])

ax.spines['top'   ].set_linewidth( lw )
ax.spines['bottom'].set_linewidth( lw )
ax.spines['left'  ].set_linewidth( lw )
ax.spines['right' ].set_linewidth( lw )

if( File == 'thornado_00049144' ):
    SaveAs = '2DRP_CharLimit.png'
else:
    SaveAs = '2DRP_CompLimit.png'


slc.save( SaveAs )
