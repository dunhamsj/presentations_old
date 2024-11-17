#!/usr/bin/env python3

# --- Import libraries ---
import warnings
warnings.simplefilter( action = 'ignore' , category = FutureWarning )

# --- Get rid of "Too many open files" error ---
import resource
resource.setrlimit( resource.RLIMIT_NOFILE, ( 10000, -1 ) )

import numpy as np
import h5py as h5

def ReadFields( PathToData, Snapshots, UseGeometryFields ):

    FF_root = PathToData + '_FluidFields_'
    GF_root = PathToData + '_GeometryFields_'

    nFiles = len( Snapshots )

    TimeSteps              = np.empty( nFiles, object )
    DataFileNames_Fluid    = np.empty( nFiles, object )
    DataFileNames_Geometry = np.empty( nFiles, object )

    # Arrays to hold data
    Data_FF = np.empty( nFiles, object )
    Data_GF = np.empty( nFiles, object )

    for i in range( nFiles ):

        TimeSteps[i] = str( Snapshots[i] ).zfill( 6 )

        DataFileNames_Fluid[i] = FF_root + str( TimeSteps[i] ) + '.h5'
        Data_FF[i]             = h5.File( DataFileNames_Fluid[i], 'r' )

        if( UseGeometryFields ):
            DataFileNames_Geometry[i] \
              = GF_root + str( TimeSteps[i] ) + '.h5'
            Data_GF[i] \
              = h5.File( DataFileNames_Geometry[i], 'r' )

    # Get the spatial grid
    X  = Data_FF[0][ 'Spatial Grid' ]
    X1 = np.array( X[ 'X1' ] )
    X2 = np.array( X[ 'X2' ] )
    X3 = np.array( X[ 'X3' ] )
    X1C = np.array( X[ 'X1_C' ] )
    X2C = np.array( X[ 'X2_C' ] )
    X3C = np.array( X[ 'X3_C' ] )

    T  = np.empty( nFiles, float  )
    FF = np.empty( nFiles, object )
    PF = np.empty( nFiles, object )
    CF = np.empty( nFiles, object )
    AF = np.empty( nFiles, object )
    DF = np.empty( nFiles, object )
    GF = np.empty( nFiles, object )

    # Conserved fields
    CF_D  = np.empty( nFiles, object )
    CF_S1 = np.empty( nFiles, object )
    CF_S2 = np.empty( nFiles, object )
    CF_S3 = np.empty( nFiles, object )
    CF_E  = np.empty( nFiles, object )
    CF_Ne = np.empty( nFiles, object )

    # Primitive fields
    PF_D  = np.empty( nFiles, object )
    PF_V1 = np.empty( nFiles, object )
    PF_V2 = np.empty( nFiles, object )
    PF_V3 = np.empty( nFiles, object )
    PF_E  = np.empty( nFiles, object )
    PF_Ne = np.empty( nFiles, object )

    # Auxiliary fields
    AF_P  = np.empty( nFiles, object )
    AF_T  = np.empty( nFiles, object )
    AF_Ye = np.empty( nFiles, object )
    AF_Cs = np.empty( nFiles, object )
    AF_Gm = np.empty( nFiles, object )

    # Diagnostic fields
    DF_TCI   = np.empty( nFiles, object )
    DF_Sh_X1 = np.empty( nFiles, object )
    DF_Sh_X2 = np.empty( nFiles, object )
    DF_Sh_X3 = np.empty( nFiles, object )

    # Geometry fields
    GF_CF = np.empty( nFiles, object )
    GF_al = np.empty( nFiles, object )
    GF_NP = np.empty( nFiles, object )
    GF_b1 = np.empty( nFiles, object )
    GF_b2 = np.empty( nFiles, object )
    GF_b3 = np.empty( nFiles, object )
    GF_g1 = np.empty( nFiles, object )
    GF_g2 = np.empty( nFiles, object )
    GF_g3 = np.empty( nFiles, object )
    GF_h1 = np.empty( nFiles, object )
    GF_h2 = np.empty( nFiles, object )
    GF_h3 = np.empty( nFiles, object )
    GF_Sg = np.empty( nFiles, object )

    # Derived fields
    PolytropicConstant = np.empty( nFiles, object )
    BernoulliConstant  = np.empty( nFiles, object )
    MassConstant       = np.empty( nFiles, object )

    for i in range( nFiles ):

        # First level groups

        FF[i] = Data_FF[i][ 'Fluid Fields'    ]
        if( UseGeometryFields ): GF[i] = Data_GF[i][ 'Geometry Fields' ]

        # Second level groups

        PF[i] = FF[i][ 'Primitive' ]
        CF[i] = FF[i][ 'Conserved' ]
        AF[i] = FF[i][ 'Auxiliary' ]
        DF[i] = FF[i][ 'Diagnostic' ]

        # Third level groups

        CF_D [i] = CF[i][ 'Conserved Baryon Density'       ][:][:][:]
        CF_S1[i] = CF[i][ 'Conserved Momentum Density (1)' ][:][:][:]
        CF_S2[i] = CF[i][ 'Conserved Momentum Density (2)' ][:][:][:]
        CF_S3[i] = CF[i][ 'Conserved Momentum Density (3)' ][:][:][:]
        CF_E [i] = CF[i][ 'Conserved Energy Density'       ][:][:][:]
        CF_Ne[i] = CF[i][ 'Conserved Electron Density'     ][:][:][:]

        PF_D [i] = PF[i][ 'Comoving Baryon Density'   ][:][:][:]
        PF_V1[i] = PF[i][ 'Three-Velocity (1)'        ][:][:][:]
        PF_V2[i] = PF[i][ 'Three-Velocity (2)'        ][:][:][:]
        PF_V3[i] = PF[i][ 'Three-Velocity (3)'        ][:][:][:]
        PF_E [i] = PF[i][ 'Internal Energy Density'   ][:][:][:]
        PF_Ne[i] = PF[i][ 'Comoving Electron Density' ][:][:][:]

        AF_P [i] = AF[i][ 'Pressure'                        ][:][:][:]
        AF_T [i] = AF[i][ 'Temperature'                     ][:][:][:]
        AF_Ye[i] = AF[i][ 'Electron Fraction'               ][:][:][:]
        AF_Cs[i] = AF[i][ 'Sound Speed'                     ][:][:][:]
        AF_Gm[i] = AF[i][ 'Ratio of Specific Heats (Gamma)' ][:][:][:]

        DF_TCI  [i] = DF[i][ 'TCI'        ][:][:][:]
        DF_Sh_X1[i] = DF[i][ 'Shock (X1)' ][:][:][:]
        DF_Sh_X2[i] = DF[i][ 'Shock (X2)' ][:][:][:]
        DF_Sh_X3[i] = DF[i][ 'Shock (X3)' ][:][:][:]

        if( UseGeometryFields ):

            GF_CF[i] = GF[i][ 'Conformal Factor'                 ][:][:][:]
            GF_al[i] = GF[i][ 'Lapse Function'                   ][:][:][:]
            GF_NP[i] = GF[i][ 'Newtonian Potential'              ][:][:][:]
            GF_b1[i] = GF[i][ 'Shift Vector (1)'                 ][:][:][:]
            GF_b2[i] = GF[i][ 'Shift Vector (2)'                 ][:][:][:]
            GF_b3[i] = GF[i][ 'Shift Vector (3)'                 ][:][:][:]
            GF_g1[i] = GF[i][ 'Spatial Metric Component (11)'    ][:][:][:]
            GF_g2[i] = GF[i][ 'Spatial Metric Component (22)'    ][:][:][:]
            GF_g3[i] = GF[i][ 'Spatial Metric Component (33)'    ][:][:][:]
            GF_h1[i] = GF[i][ 'Spatial Scale Factor (1)'         ][:][:][:]
            GF_h2[i] = GF[i][ 'Spatial Scale Factor (2)'         ][:][:][:]
            GF_h3[i] = GF[i][ 'Spatial Scale Factor (3)'         ][:][:][:]
            GF_Sg[i] = GF[i][ 'Sqrt of Spatial Metric Determina' ][:][:][:]

        # Derived fields

        c = 2.99792458e10
        h = c**2 + ( PF_E[i] + AF_P[i] ) / PF_D[i]
        W = 1.0 / np.sqrt( 1.0 - GF_CF[i]**4 * (PF_V1[i]*1.0e5)**2 / c**2 )

        PolytropicConstant[i] = AF_P[i] / PF_D[i]**AF_Gm[i]
        BernoulliConstant [i] = GF_al[i] * h * W
#        MassConstant      [i] = 4.0 * np.pi * GF_al[i] * (GF_Sg[i]*1.0e10) \
#                                  * PF_D[i] * W * (PF_V1[i]*1.0e5)
        MassConstant      [i] = 4.0 * np.pi * GF_al[i] * GF_CF[i]**6 * X1**2 \
                                  * PF_D[i] * W * (PF_V1[i]*1.0e5)

        # Time

        T[i] = Data_FF[i][ 'Time' ][0]

    names =                             \
    {                                   \
      'Time'               : T        , \
      'CF_D'               : CF_D     , \
      'CF_S1'              : CF_S1    , \
      'CF_S2'              : CF_S2    , \
      'CF_S3'              : CF_S3    , \
      'CF_E'               : CF_E     , \
      'CF_Ne'              : CF_Ne    , \
      'PF_D'               : PF_D     , \
      'PF_V1'              : PF_V1    , \
      'PF_V2'              : PF_V2    , \
      'PF_V3'              : PF_V3    , \
      'PF_E'               : PF_E     , \
      'PF_Ne'              : PF_Ne    , \
      'AF_P'               : AF_P     , \
      'AF_T'               : AF_T     , \
      'AF_Ye'              : AF_Ye    , \
      'AF_Cs'              : AF_Cs    , \
      'AF_Gm'              : AF_Gm    , \
      'GF_CF'              : GF_CF    , \
      'GF_al'              : GF_al    , \
      'GF_NP'              : GF_NP    , \
      'GF_b1'              : GF_b1    , \
      'GF_b2'              : GF_b2    , \
      'GF_b2'              : GF_b3    , \
      'GF_g1'              : GF_g1    , \
      'GF_g2'              : GF_g2    , \
      'GF_g3'              : GF_g3    , \
      'GF_h1'              : GF_h1    , \
      'GF_h2'              : GF_h2    , \
      'GF_h3'              : GF_h3    , \
      'GF_Sg'              : GF_Sg    , \
      'X1C'                : X1C      , \
      'X2C'                : X2C      , \
      'X3C'                : X3C      , \
      'X1'                 : X1       , \
      'X2'                 : X2       , \
      'X3'                 : X3       , \
      'DF_TCI'             : DF_TCI   , \
      'DF_Sh_X1'           : DF_Sh_X1 , \
      'DF_Sh_X2'           : DF_Sh_X2 , \
      'DF_Sh_X3'           : DF_Sh_X3 , \
      'PolytropicConstant' : PolytropicConstant, \
      'BernoulliConstant'  : BernoulliConstant, \
      'MassConstant'       : MassConstant
    }

    return names
