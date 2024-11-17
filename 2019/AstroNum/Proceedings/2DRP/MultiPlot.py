#!/usr/bin/env python3

import yt
import numpy as np
import subprocess
from os import listdir
from os.path import isfile
from sys import argv, exit
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

yt.funcs.mylog.setLevel(0) # Suppress initial yt output to screen

# --- Get user's HOME directory ---
HOME = subprocess.check_output( ["echo $HOME"], shell = True)
HOME = HOME[:-1].decode( "utf-8" ) + '/'

# Specify directory containing plotfiles

ProblemName    = '2DRP'
VariableToPlot = 'AF_P'

Root = HOME + 'Research/DataFromPastRuns/AstroNum2019/'

Files = [ Root + '2DRP_CharLimit/thornado_00049144', \
          Root + '2DRP_CompLimit/thornado_00049151' ]

Limiting = [ 'Characteristic Limiting', 'Component-Wise Limiting' ]

fig = plt.figure()

grid = AxesGrid( fig, \
                 rect          = (0.075,0.075,0.85,0.85), \
                 nrows_ncols   = (1,2), \
                 axes_pad      = 0.4, \
                 label_mode    = "L", \
                 share_all     = True, \
                 cbar_location = "right", \
                 cbar_mode     = "single", \
                 cbar_size     = "3%", \
                 cbar_pad      = "1%" )

lw = 1.5
vmin = 1.0e-2
vmax = 5.0e1

for i, File in enumerate( Files ):

    ds = yt.load( '{:}'.format( File ) )
    Time = ds.current_time

    p = yt.plot_2d( ds, VariableToPlot, origin = 'lower-left-window' )

    p.set_log( VariableToPlot, True )
    p.set_zlim( VariableToPlot, vmin, vmax )
    p.set_cmap( VariableToPlot, cmap = 'Purples' )
    p.annotate_text( [0.02,0.925], \
                     '{:}'.format( Limiting[i] ), \
                     coord_system = 'axis', \
                     text_args={'color':'black'} )
    p.set_colorbar_label( VariableToPlot, 'Pressure' )
    p.set_xlabel( 'x' )
    p.set_ylabel( 'y' )

    p.set_minorticks( 'all', 'off' )

    plot = p.plots[VariableToPlot]

    plot.figure = fig
    plot.axes   = grid[i].axes

    plot.axes.spines['top'   ].set_linewidth( lw )
    plot.axes.spines['bottom'].set_linewidth( lw )
    plot.axes.spines['left'  ].set_linewidth( lw )
    plot.axes.spines['right' ].set_linewidth( lw )

    #plot.axes.set_xlim( 0.0, 0.99 )

    plot.cax = grid.cbar_axes[i]

    p._setup_plots()

plt.savefig( '2DRP.png', bbox_inches = 'tight', pad_inches = 0.1, dpi = 300 )
