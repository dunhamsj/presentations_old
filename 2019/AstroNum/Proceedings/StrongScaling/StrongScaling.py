import numpy as np
import matplotlib.pyplot as plt
plt.rc('text',usetex=True)

n  = np.array( [ 1, 2, 4, 8, 16 ] )
tt = 1/n

t = np.array( [ \
1.243875E+003
,4.419989E+002
,3.266438E+002
,9.427181E+001
,4.848672E+001] )

lw = 1.5

# Log
f, ax = plt.subplots( figsize = (8,2) )
ax.plot( [n[0],n[-1]], [tt[0]/tt[0],tt[0]/tt[-1]], 'k:', \
         label = 'Linear Speed-Up' )
ax.plot( n, t[0]/t, 'k.-' )
ax.set_xlabel( r'$N$' )
ax.set_ylabel( 'Speed-Up' )
ax.legend()
ax.tick_params( axis = 'both', direction = 'in' )
ax.grid(True)
ax.spines['top'   ].set_linewidth( lw )
ax.spines['bottom'].set_linewidth( lw )
ax.spines['left'  ].set_linewidth( lw )
ax.spines['right' ].set_linewidth( lw )

#plt.subplots_adjust( hspace = 0.3 )
#plt.show()
plt.savefig( 'StrongScaling.png', bbox_inches = 'tight', pad_inches = 0.1 )
