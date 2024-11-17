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
  + 'Research/DataFromPastRuns/AstroNum2019/SASI_nNodes3_HLL_NoTCI_200x128/'
ProblemName      = 'SASI'
VariableToPlot   = 'Entropy'
cmap             = 'Purples'
UseLogScale      = True
UseCustomTicks   = True

File = 'thornado_00119695' # 165 ms
File = 'thornado_00168201' # 218 ms
File = 'thornado_00194137' # 241 ms
File = 'thornado_00304527' # 333 ms (last slice before hitting boundary)

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

c = 2.99792458e10
if  ( VariableToPlot == 'PF_D'  ):
    Data = CoveringGrid['PF_D' ].to_ndarray()
    DataUnit = 'g/cm**3'
elif( VariableToPlot == 'PF_V1' ):
    Data = CoveringGrid['PF_V1'].to_ndarray() * 1.0e5
    DataUnit = 'km/s'
elif( VariableToPlot == 'PF_V2' ):
    Data = CoveringGrid['PF_V2'].to_ndarray() * 1.0e5
    DataUnit = 'km/s'
elif( VariableToPlot == 'PF_V3' ):
    Data = CoveringGrid['PF_V3'].to_ndarray() * 1.0e5
    DataUnit = 'km/s'
elif( VariableToPlot == 'PF_E'  ):
    Data = CoveringGrid['PF_E' ].to_ndarray()
    DataUnit = 'erg/cm**3'
elif( VariableToPlot == 'CF_D'  ):
    Data = CoveringGrid['CF_D' ].to_ndarray()
    DataUnit = 'g/cm**3'
elif( VariableToPlot == 'CF_S1' ):
    Data = CoveringGrid['CF_S1'].to_ndarray()
    DataUnit = 'g/cm**2/s'
elif( VariableToPlot == 'CF_S2' ):
    Data = CoveringGrid['CF_S2'].to_ndarray()
    DataUnit = 'g/cm**2/s'
elif( VariableToPlot == 'CF_S3' ):
    Data = CoveringGrid['CF_S3'].to_ndarray()
    DataUnit = 'g/cm**2/s'
elif( VariableToPlot == 'CF_E'  ):
    Data = CoveringGrid['CF_E' ].to_ndarray()
    DataUnit = 'erg/cm**3'
elif( VariableToPlot == 'AF_P'  ):
    Data = CoveringGrid['AF_P' ].to_ndarray()
    DataUnit = 'erg/cm**3'
elif( VariableToPlot == 'Entropy' ):
    PF_D  = CoveringGrid['PF_D' ].to_ndarray()
    AF_P  = CoveringGrid['AF_P' ].to_ndarray()
    AF_Gm = CoveringGrid['AF_Gm'].to_ndarray()
    Data  = AF_P / PF_D**AF_Gm
    DataUnit = 'erg/cm**3/(g/cm**3)**({:}/3)'.format( \
                 int( 3 * AF_Gm[0][0][0] ) )
elif( VariableToPlot == 'LorentzFactor' ):
    PF_V1 = CoveringGrid['PF_V1'   ].to_ndarray() * 1.0e5
    PF_V2 = CoveringGrid['PF_V2'   ].to_ndarray() * 1.0e5
    PF_V3 = CoveringGrid['PF_V3'   ].to_ndarray() * 1.0e5
    GF_g1 = CoveringGrid['GF_Gm_11'].to_ndarray()
    GF_g2 = CoveringGrid['GF_Gm_22'].to_ndarray()
    GF_g3 = CoveringGrid['GF_Gm_33'].to_ndarray()
    VSq = GF_g1 * PF_V1**2 + GF_g2 * PF_V2**2 + GF_g3 * PF_V3**2
    Data  = 1.0 / np.sqrt( 1.0 - VSq / c**2 )
    DataUnit = ''
elif( VariableToPlot == 'SpecificEnthalpy' ):
    PF_D = CoveringGrid['PF_D'].to_ndarray()
    PF_E = CoveringGrid['PF_E'].to_ndarray()
    AF_P = CoveringGrid['AF_P'].to_ndarray()
    Data  = ( c**2 + ( PF_E + AF_P ) / PF_D ) / c**2
    print( np.max( Data ) );exit()
    DataUnit = ''

data    = { VariableToPlot: (Data,DataUnit) }
fields  = VariableToPlot

ds = yt.load_uniform_grid \
       ( data, \
         nX, \
         bbox = np.array( \
                  [ [xL[0],xH[0]], [xL[1],xH[1]], [xL[2],xH[2]] ] ), \
         length_unit = 'km', \
         geometry = 'spherical' )

slc = yt.SlicePlot( ds, 'phi', fields, \
                    axes_unit = 'km', \
                    origin = 'lower-left-window', \
                    width = ((1080,'km'),(1080,'km')) )

if( UseLogScale ):
    slc.set_log( fields, True )
else:
    slc.set_log( fields, False )

if( UseCustomTicks ):
    vmin = 1.0e15
    vmax = 1.0e16
    slc.set_zlim( fields, vmin, vmax )

#slc.set_minorticks( 'all', 'off' )

# First argument is location in figure, [0,0] is lower-left,
#                                       [1,1] is upper-right
slc.annotate_text( [0.05, 0.9], 't = {:d} ms'.format( int( Time ) ), \
                   coord_system = 'axis', \
                   text_args={'color':'black'} )

slc.set_cmap( field = fields, cmap = cmap )
slc.set_colorbar_label( fields, r'$\mathrm{Polytropic\ Constant} \left(\mathrm{\frac{erg/cm^{3}}{\left(g/cm^{3}\right)^{4/3}}}\right)$' )

lw = 1.5

ax = slc.plots[VariableToPlot].axes

ax.spines['top'   ].set_linewidth( lw )
ax.spines['bottom'].set_linewidth( lw )
ax.spines['left'  ].set_linewidth( lw )
ax.spines['right' ].set_linewidth( lw )

slc.save( ProblemName + '_' + VariableToPlot \
            + '_{:}ms.png'.format( str( int( Time ) ).zfill(3) ) )
