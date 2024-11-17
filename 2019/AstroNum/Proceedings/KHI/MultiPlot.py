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

ProblemDirectory = HOME \
  + 'Research/DataFromPastRuns/AstroNum2019/KHI_20190904/'
ProblemName      = 'KHI'
VariableToPlot   = 'PF_D'

Files = [ 'thornado_00006012', \
          'thornado_00007889', \
          'thornado_00010641', \
          'thornado_00013499' ]

fig = plt.figure()

grid = AxesGrid( fig, \
                 rect          = (0.075,0.075,0.85,0.85), \
                 nrows_ncols   = (2,2), \
                 axes_pad      = 0.2, \
                 label_mode    = "L", \
                 share_all     = True, \
                 cbar_location = "right", \
                 cbar_mode     = "single", \
                 cbar_size     = "5%", \
                 cbar_pad      = "1%" )

lw = 1.5
vmin = 0.0
vmax = 2.0

for i, File in enumerate( Files ):

    ds = yt.load( '{:}'.format( ProblemDirectory + File ) )
    Time = ds.current_time

    p = yt.plot_2d( ds, VariableToPlot, aspect = 0.5 )

    p.set_log( VariableToPlot, False )
    p.set_zlim( VariableToPlot, vmin, vmax )
    p.set_cmap( VariableToPlot, cmap = 'jet' )
    p.annotate_text( [0.02,0.9], \
                     't = {:.2f}'.format( Time.to_ndarray() ), \
                     coord_system = 'axis', \
                     text_args={'color':'black'} )
    p.set_colorbar_label( VariableToPlot, 'Rest-Mass--Density' )
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

    plot.cax = grid.cbar_axes[i]

    p._setup_plots()
plt.savefig( 'KHI.png', bbox_inches = 'tight', pad_inches = 0.1 )
