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
js = [200,250,300] # the latitudes where to make hovmuller diagrams

Nt = len(iters)
Nx = mod_flat.Nx
Nj = len(js)

hov_vort_flat = zeros((Nj,Nt,Nx))
hov_vort_bump = zeros((Nj,Nt,Nx))
hov_theta_flat = zeros((Nj,Nt,Nx))
hov_theta_bump = zeros((Nj,Nt,Nx))

hovs_vort = {'flat': hov_vort_flat, 'bump': hov_vort_bump}
hovs_theta = {'flat': hov_theta_flat, 'bump': hov_theta_bump}
mods = {'flat': mod_flat, 'bump': mod_bump}

DX = mod_flat.dxc[0]
for n in arange(Nt):
    for mt in ['flat', 'bump']:
        m = mods[mt]
        h_vort = hovs_vort[mt]
        h_theta = hovs_theta[mt]
        UV = m.rdmds('UV', iters[n], useMask=False, lev=0)
        T = m.rdmds('THETA', iters[n], useMask=False, lev=0)
        U,V = UV[0], UV[1]
        vort = ( -(U - roll(U,1,axis=0)) + (V - roll(V,1,axis=1)))/DX
        for nj in arange(Nj):
            h_vort[nj,n] = vort[js[nj]]
            h_theta[nj,n] = T[js[nj]]
            
# for mean fields
base_dir = os.path.join(os.environ['D'],'projects','drag_strat')
b_bump = bump_channel.BumpChannel(base_dir, 'taux2000_rb0110_bump')
b_flat = bump_channel.BumpChannel(base_dir, 'taux2000_rb0110_flat')

U_bump = b_bump.rdmds('Uveltave')
U_flat = b_flat.rdmds('Uveltave')
T_bump = b_bump.rdmds('Ttave', useMask=False, lev=0)
T_flat = b_flat.rdmds('Ttave', useMask=False, lev=0)

Us_bump = U_bump[0].mean(axis=1)
Us_flat = U_flat[0].mean(axis=1)
Ubt_bump = b_bump.average_vertical(U_bump).mean(axis=1)
Ubt_flat = b_flat.average_vertical(U_flat).mean(axis=1)


T = iters*900 / (24*60*60.)
X = arange(0,Nx)*5
Tlims = [500,880]
Tstart = 600
deltaT = 1000
vortlevs = (arange(-6,6)+0.5)*1e-5
theta_levs = (arange(-1,1,0.2) + 0.1)*2
for nj in arange(Nj):
    figure(nj, figsize=(10,6))
    clf()
    subplot(2,1,1)
    contourf(T,X,hov_vort_flat[nj].T , vortlevs,
        cmap=get_cmap('bwr'), extend='both' , rasterized=True)
    plot([Tstart, Tstart+deltaT],[0, deltaT*(24*60*60)/1000*Ubt_flat[js[nj]]], 'k-')
    plot([Tstart, Tstart+deltaT],[0, deltaT*(24*60*60)/1000*Us_flat[js[nj]]], 'k--')
    xlim(Tlims); ylim(X[r_[0,-1]])
    ylabel('X (km)'); xlabel('Time (days)')
    title('Flat, Y = %g km' % (5*js[nj]))
    subplot(2,1,2)
    contourf(T,X,hov_vort_bump[nj].T , vortlevs,
        cmap=get_cmap('bwr'), extend='both' , rasterized=True)
    plot([Tstart, Tstart+deltaT],[0, deltaT*(24*60*60)/1000*Ubt_bump[js[nj]]], 'k-')
    plot([Tstart, Tstart+deltaT],[0, deltaT*(24*60*60)/1000*Us_bump[js[nj]]],'k--')
    xlim(Tlims); ylim(X[r_[0,-1]])
    ylabel('X (km)'); xlabel('Time (days)')
    title('Bump, Y = %g km' % (5*js[nj]))
    tight_layout()
    #colorbar()
    savefig('figures/hov/vort_Y%g.png' % (5*js[nj]))

    figure(nj+3, figsize=(10,6))
    clf()
    subplot(2,1,1)
    contourf(T,X,hov_theta_flat[nj].T - T_flat[js[nj]][:,newaxis], theta_levs,
        cmap=get_cmap('posneg'), extend='both' , rasterized=True)
    plot([Tstart, Tstart+deltaT],[0, deltaT*(24*60*60)/1000*Ubt_flat[js[nj]]], 'k-')
    plot([Tstart, Tstart+deltaT],[0, deltaT*(24*60*60)/1000*Us_flat[js[nj]]], 'k--')
    xlim(Tlims); ylim(X[r_[0,-1]])
    ylabel('X (km)'); xlabel('Time (days)')
    title('Flat, Y = %g km' % (5*js[nj]))
    subplot(2,1,2)
    contourf(T,X,hov_theta_bump[nj].T - T_bump[js[nj]][:,newaxis], theta_levs,
        cmap=get_cmap('posneg'), extend='both' , rasterized=True)
    plot([Tstart, Tstart+deltaT],[0, deltaT*(24*60*60)/1000*Ubt_bump[js[nj]]], 'k-')
    plot([Tstart, Tstart+deltaT],[0, deltaT*(24*60*60)/1000*Us_bump[js[nj]]],'k--')
    xlim(Tlims); ylim(X[r_[0,-1]])
    ylabel('X (km)'); xlabel('Time (days)')
    title('Bump, Y = %g km' % (5*js[nj]))
    tight_layout()
    savefig('figures/hov/theta_Y%g.png' % (5*js[nj]))    
    

            