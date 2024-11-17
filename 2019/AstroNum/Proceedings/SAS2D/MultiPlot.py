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
  + 'Research/Data/AstroNum2019/SASI_MPNS2.8_Rs180_Mdot3.0/'
ProblemName      = 'SASI'
VariableToPlot   = 'Entropy'

Files = [ 'thornado_00118313', \
          'thornado_00174725' ]#, \
          #'thornado_00581839', \
          #'thornado_00709807' ]

fig = plt.figure()

grid = AxesGrid( fig, \
                 rect          = (0.075,0.075,0.85,0.85), \
                 nrows_ncols   = (1,4), \
                 axes_pad      = 0.05, \
                 label_mode    = "L", \
                 share_all     = True, \
                 cbar_location = "right", \
                 cbar_mode     = "single", \
                 cbar_size     = "5%", \
                 cbar_pad      = "3%" )

lw = 1.5
vmin = 1.0e15
vmax = 1.0e16

for i, File in enumerate( Files ):

    ds = yt.load( '{:}'.format( ProblemDirectory + File ) )
    MaxLevel = ds.index.max_level
    Time = ds.current_time
    nX       = ds.domain_dimensions
    xL       = ds.domain_left_edge
    xR       = ds.domain_right_edge
    print( nX );exit()

    CoveringGrid \
      = ds.covering_grid \
          ( level           = MaxLevel, \
            left_edge       = xL, \
            dims            = nX * 2**MaxLevel, \
            num_ghost_zones = nX[0] )

    PF_D  = CoveringGrid['PF_D' ].to_ndarray()
    AF_P  = CoveringGrid['AF_P' ].to_ndarray()
    AF_Gm = CoveringGrid['AF_Gm'].to_ndarray()
    Data  = AF_P / PF_D**AF_Gm
    DataUnit = 'erg/cm**3/(g/cm**3)**({:}/3)'.format( \
                 int( 3 * AF_Gm[0][0][0] ) )

    data = { VariableToPlot: (Data,DataUnit) }

    ds = yt.load_uniform_grid \
           ( data, \
             nX, \
             bbox = np.array( \
                      [ [xL[0],xR[0]], [xL[1],xR[1]], [xL[2],xR[2]] ] ), \
             length_unit = 'km', \
             geometry = 'spherical' )

    slc = yt.SlicePlot( ds, 'phi', VariableToPlot, \
                        axes_unit = 'km', \
                        origin = 'lower-left-window', \
                        width = ((1080,'km'),(1080,'km')) )
    slc.set_zlim( VariableToPlot, vmin, vmax )
    slc.set_log( VariableToPlot, True )

    slc.annotate_text( [0.4, 0.05], 't = {:d} ms'.format( int( Time ) ), \
                       coord_system = 'axis', \
                       text_args={'color':'black'} )

    slc.set_cmap( VariableToPlot, cmap = 'Purples' )
    slc.set_colorbar_label( VariableToPlot, r'$\mathrm{Polytropic\ Constant} \left(\mathrm{\frac{erg/cm^{3}}{\left(g/cm^{3}\right)^{4/3}}}\right)$' )

    plot = slc.plots[VariableToPlot]

    plot.figure = fig
    plot.axes   = grid[i].axes

    plot.axes.spines['top'   ].set_linewidth( lw )
    plot.axes.spines['bottom'].set_linewidth( lw )
    plot.axes.spines['left'  ].set_linewidth( lw )
    plot.axes.spines['right' ].set_linewidth( lw )

    plot.axes.set_xlim( 0,    540 )
    plot.axes.set_ylim( -540, 540 )

    plot.cax = grid.cbar_axes[i]

    slc._setup_plots()

plt.savefig( 'SAS2D.png', bbox_inches = 'tight', pad_inches = 0.1 )
