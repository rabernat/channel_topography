from MITgcmutils.mds import rdmds
import bump_channel
import os
import mycolors
from pylab import *

base_dir = os.path.join(os.environ['D'],'projects','drag_strat')

grid_dir = '/Users/rpa/Data/DATASTORE.RPA/projects/drag_strat/grid_bump'
data_dir_flat = '/Volumes/scratch/tmp/channel/taux2000_rb0110_flat/mds/'
data_dir_bump = '/Volumes/scratch/tmp/channel/taux2000_rb0110_bump/mds_orig/'
output_dir = '/Volumes/scratch/tmp/channel/taux2000_rb0110_bump/vtk'

#mod_flat = MITgcmmodel.ModelInstance(data_dir_flat, grid_dir=grid_dir)
#mod_bump = MITgcmmodel.ModelInstance(data_dir_bump, grid_dir=grid_dir)
bflat = bump_channel.BumpChannel(base_dir, 'taux2000_rb0110_flat')
bflat.load_default_fields()
bbump = bump_channel.BumpChannel(base_dir, 'taux2000_rb0110_bump')
bbump.load_default_fields()

iters = arange(0,131520+1,192)

nb = 50
nf = 145

Tf = rdmds(data_dir_flat + 'THETA', iters[nf], lev=0)
Tb = rdmds(data_dir_bump + 'THETA', iters[nb], lev=0)

clevs = arange(0,8.1,0.5)
X = linspace(0,2000,bflat.Nx)
Y = X
myticks = [500,1000,1500]

close('all')
rc('figure.subplot', left=0.08, right=0.88, bottom=0.05, top=0.9, wspace=0.12, hspace=0.22)
rcParams['legend.fontsize'] = 7
rcParams['font.size'] = 7.
rcParams['lines.markersize'] = 4

figure(figsize=(7,6.0))
clf()
subplot(221);
#contourf(X,Y,Tf, clevs, extend='both', cmap=get_cmap('sst'))
pcolormesh(X,Y,Tf, cmap=get_cmap('sst'), rasterized=True); clim([0,8])
contour(X,Y,Tf, clevs, extend='both', colors='k', linewidths=0.3)
xticks(myticks); yticks(myticks)
#xlabel('x (km)');
ylabel('y (km)')
title(r'$\theta_s$ - Flat @ %g days' % (nf*2))

subplot(222);
#CF=contourf(X,Y,Tb, clevs, extend='both', cmap=get_cmap('sst'))
PC=pcolormesh(X,Y,Tb, cmap=get_cmap('sst'), rasterized=True); clim([0,8])
contour(X,Y,Tb, clevs, extend='both', colors='k', linewidths=0.3)
xticks(myticks); yticks(myticks, [])
#xlabel('x (km)'); ylabel('y (km)')
title(r'$\theta_s$ - Ridge @ %g days' % (nb*2))

y0 = gca().get_position().y0; dy = gca().get_position().ymax - y0
cb=colorbar(PC, cax=axes([0.92,y0,0.02,dy]), ticks=clevs)
gcf().text(0.92, 0.92, r'$^\circ$C')

wblevs = arange(-6,6) + 0.5
subplot(223)
contourf(X,Y, bbump.g * bbump.tAlpha * bflat.integrate_vertical(bflat.WpTp) * 1e5,
    wblevs, cmap=get_cmap('posneg'), extend='both')
xticks(myticks); yticks(myticks)
xlabel('x (km)'); ylabel('y (km)')
title(r"$\int \overline{w'b'} dz$ - Flat")
subplot(224)
CF=contourf(X,Y, bbump.g * bbump.tAlpha * bbump.integrate_vertical(bbump.WpTp) * 1e5,
    wblevs, cmap=get_cmap('posneg'), extend='both')
xticks(myticks); yticks(myticks, [])
xlabel('x (km)'); #ylabel('y (km)')
title(r"$\int \overline{w'b'} dz$ - Ridge")

y0 = gca().get_position().y0; dy = gca().get_position().ymax - y0
cb=colorbar(CF, cax=axes([0.92,y0,0.02,dy]), ticks=wblevs[1:] - 0.5)
gcf().text(0.90, 0.45, r'10$^{-5}$ m$^3$ s$^{-3}$')

show()

savefig('figures/paper/theta_spinup_new.pdf')

            