from pylab import *
import os
import bump_channel

run_name = 'taux2000_rb0110_bump'
base_dir = os.path.join(os.environ['D'],'projects','drag_strat')
b = bump_channel.BumpChannel(base_dir, run_name)
outdir = '/Volumes/scratch/tmp/output_for_paraview/' + run_name

T = b.rdmds('Ttave')
V = b.rdmds('Vveltave')
U = b.rdmds('Uveltave')
W = b.rdmds('Wveltave')
P = b.rdmds('Phhytave')

UT = b.rdmds('UTtave')
VT = b.rdmds('VTtave')
WT = b.rdmds('WTtave')
Tdif = b.rdmds('Tdiftave')

Tbar = T.mean(axis=2)[:,:,newaxis]
Vbar = V.mean(axis=2)[:,:,newaxis]
Ubar = U.mean(axis=2)[:,:,newaxis]
Wbar = W.mean(axis=2)[:,:,newaxis]

Panom = P - P[:,:,0][:,:,newaxis]

# terms in heat budget
dx_UT = b.ddx_ugrid_to_cgrid(UT)
dy_VT = b.ddy_vgrid_to_cgrid(VT)
dz_WT = b.ddz_wgrid_to_cgrid(WT)

dz_Tdif = b.ddz_wgrid_to_cgrid(Tdif)

# integrated in vertical already
month = 24*60*60*30.
forc = - (T[0] - linspace(0,8,b.Ny)[:,newaxis]) / month * b.dzf[0]

# mixed layer heat budget
krange = r_[:11]
ml_dx_UT = b.average_vertical(dx_UT, krange=krange)
ml_dy_VT = b.average_vertical(dy_VT, krange=krange)
ml_dz_WT = b.average_vertical(dz_WT, krange=krange)
ml_dz_Tdif=b.average_vertical(dz_Tdif, krange=krange)
ml_forc = forc / b.dzf[krange].sum() 



x = b.xc / 1000
y = b.yc / 1000
mld = b.dzf[krange].sum()
mylim = array([-5,5])
clf()
subplot(231);
pcolormesh(x, y, month*(ml_dx_UT + ml_dy_VT), cmap=get_cmap('bwr'))
title(r'$\nabla_h u\theta$')
clim(mylim)
colorbar()
subplot(232);
pcolormesh(x, y, month*(ml_dz_WT), cmap=get_cmap('bwr'))
clim(mylim)
title(r'$\partial_z w\theta$')
colorbar()
subplot(233);
pcolormesh(x, y, month*(ml_dx_UT + ml_dy_VT+ml_dz_WT), cmap=get_cmap('bwr'))
clim(mylim/5)
title(r'$\nabla v\theta$')
colorbar()
subplot(234);
pcolormesh(x, y, month*(ml_dz_Tdif), cmap=get_cmap('bwr'))
clim(mylim/5)
title(r'$\partial_z \kappa \partial_z\theta$')
colorbar()
subplot(235);
pcolormesh(x, y, month*(ml_dz_Tdif+ml_dx_UT + ml_dy_VT+ml_dz_WT), cmap=get_cmap('bwr'))
clim(mylim/5)
title(r'$\nabla v\theta + \partial_z \kappa \partial_z\theta$')
colorbar()
subplot(236);
pcolormesh(x, y, month*(ml_forc), cmap=get_cmap('bwr'))
clim(mylim/5)
title('forcing')
colorbar()
t = gcf().text(0.5,0.95, 'Heat Budget (deg / month) Upper %g m' % mld, fontsize='large', ha='center')
savefig('figures/heat_budget_upper_%gm.png' % mld)


# eddy flux convergence
UpTp = UT - U*T
VpTp = VT - V*T
WpTp = WT - W*T
Vwave = V - Vbar
Uwave = U - Ubar
Wwave = W - Wbar
Twave = T - Tbar
VTwave = Vwave * Twave
UTwave = Uwave * Twave
WTwave = Wwave * Twave

T0 = tile(b.average_vertical(T)[newaxis,:,:], [b.Nz, 1,1])
Psi = tile( (cumsum(b.integrate_vertical(U), axis=0) * b.dyc[0] / 1e6)[newaxis,:,:], [b.Nz, 1,1])

# output VTK
outdir = '/Volumes/scratch/tmp/channel/wave/'
b.output_VTK(outdir + 'wave',
            {'T': T, 'U': U, 'V': V, 'P': P, 'T0': T0, 'PSI': Psi,
              'UTwave': UTwave, 'VTwave': VTwave, 'WTwave': WTwave,
              'UpTp': UpTp, 'VpTp': VpTp, 'WpTp': WpTp} )

eddy_flux_conv = b.ddx_ugrid_to_cgrid(UpTp) + b.ddy_vgrid_to_cgrid(VpTp) + b.ddz_wgrid_to_cgrid(WpTp)

# now do the quasilinear decomposition
# linearized mixed layer budget
ml_Ubar = b.ugrid_to_cgrid(b.average_vertical(Ubar, krange=krange))
ml_Tybar = b.ddy_cgrid_centered(b.average_vertical(Tbar, krange=krange))
ml_Vwave = b.vgrid_to_cgrid(b.average_vertical(Vwave, krange=krange))
ml_Txwave = b.ddx_cgrid_centered(b.average_vertical(Twave, krange=krange))
ml_forcwave = - Twave[0] / month * b.dzf[0] / b.dzf[krange].sum()
ml_eddy_flux = b.average_vertical(eddy_flux_conv, krange=krange)
#ml_dx_UT = b.average_vertical


