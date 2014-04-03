from MITgcmdata import MITgcmmodel
import bump_channel
import os
import mycolors
from pylab import *

grid_dir = '/Users/rpa/Data/DATASTORE.RPA/projects/drag_strat/grid_bump'
data_dir_flat = '/Volumes/scratch/tmp/channel/taux2000_rb0110_flat/mds'
data_dir_bump = '/Volumes/scratch/tmp/channel/taux2000_rb0110_bump/mds_orig'
output_dir = '/Volumes/scratch/tmp/channel/taux2000_rb0110_bump/vtk'

mod_flat = MITgcmmodel.ModelInstance(data_dir_flat, grid_dir=grid_dir)
mod_bump = MITgcmmodel.ModelInstance(data_dir_bump, grid_dir=grid_dir)

iters = arange(0,131520+1,192)

clevs = arange(0,8.1,0.5)
X = linspace(0,2000,mod_flat.Nx)
Y = X
myticks = [500,1000,1500]

close('all')
rc('figure.subplot', left=0.08, right=0.88, bottom=0.15, top=0.9, wspace=0.12, hspace=0.12)
rcParams['legend.fontsize'] = 7
rcParams['font.size'] = 7.
rcParams['lines.markersize'] = 4

figure(figsize=(7,3.0))

for n in arange(len(iters)):
    
    Tf = mod_flat.rdmds('THETA', iters[n], lev=0)
    Tb = mod_bump.rdmds('THETA', iters[n], lev=0)

    clf()
    subplot(121);
    #contourf(X,Y,Tf, clevs, extend='both', cmap=get_cmap('sst'))
    pcolormesh(X,Y,Tf, cmap=get_cmap('sst')); clim([0,8])
    contour(X,Y,Tf, clevs, extend='both', colors='k', linewidths=0.5)
    xticks(myticks); yticks(myticks)
    xlabel('x (km)'); ylabel('y (km)')
    title('Flat - %g days' % (n*2))
    subplot(122);

    #CF=contourf(X,Y,Tb, clevs, extend='both', cmap=get_cmap('sst'))
    PC=pcolormesh(X,Y,Tb, cmap=get_cmap('sst')); clim([0,8])
    contour(X,Y,Tb, clevs, extend='both', colors='k', linewidths=0.5)
    xticks(myticks); yticks(myticks, [])
    xlabel('x (km)'); #ylabel('y (km)')
    title('Ridge - %g days' % (n*2))

    cb=colorbar(PC, cax=axes([0.91,0.15,0.03,0.75]), ticks=clevs)
    gcf().text(0.91, 0.92, r'$\theta$ ($^\circ$C)')
    show()

    savefig('movies/surface_theta_spinup/frame_%04d.png' % n)

            