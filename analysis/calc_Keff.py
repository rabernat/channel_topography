from pylab import *
import bump_channel
from effdiff import effdiff
from fnmatch import fnmatch
import os

base_dir = os.path.join(os.environ['D'],'DATASTORE.RPA','projects','drag_strat')
Nr = len(bump_channel.wind_runs)

N = 50 # resolution of Keff calculation

# set up result variables
Le2 = {'surface_mean': zeros((Nr, N)),
       'barotropic_mean': zeros((Nr, N)),
       'surface_full': zeros((Nr, N)),
       'barotropic_full': zeros((Nr, N))}

# set up effective diffusivity engine
b = bump_channel.BumpChannel(base_dir, 'taux2000_rb0110_flat'),      
rac = ma.masked_array(b[0].rac, b[0].mask[0])
dx = tile(b[0].dxc,[b[0].Ny,1])
dy = tile(b[0].dyc,[b[0].Nx,1]).T
eng = effdiff.EffDiffEngine(rac,dx,dy,N)

nr = -1
for run in bump_channel.wind_runs:
    nr += 1
    B = bump_channel.BumpChannel(base_dir, run)
    Tbar = B.rdmds('Ttave')
    Tbar_b = B.average_vertical(Tbar)
    Le2['surface_mean'][nr] = eng.calc_Le2(Tbar[0])['Le2']
    Le2['barotropic_mean'][nr] = eng.calc_Le2(Tbar_b)['Le2']
    
    # find matching snapshots
    match_iters = []
    for f in os.listdir(B.output_dir):
        if fnmatch(f,'THETA_2month.*.data'):
            match_iters.append(int(f[-15:-5]))
    Nt = len(match_iters)
    
    for n in range(Nt):
        T = B.rdmds('THETA_2month', match_iters[n])
        Tb = B.average_vertical(T)
        Le2['surface_full'][nr] += eng.calc_Le2(T[0])['Le2']
        Le2['barotropic_full'][nr] += eng.calc_Le2(Tb)['Le2']
    if Nt > 0:
        Le2['surface_full'][nr] /= Nt
        Le2['barotropic_full'][nr] /= Nt

Y = eng.A / B.Lx / 1e3
savez('data/Le2.npz', Y=Y, **Le2)

yrange = r_[25:45]
idx_flat,idx_bump = r_[:7],r_[7:14]

tau = logspace(0,6,7,base=2) * 0.0125

close('all')
figure(figsize=(6.5,5))
clf()
subplot(221)
semilogx(tau,ma.masked_invalid(Le2['barotropic_full'])[:,yrange].mean(axis=1)[idx_flat] / B.Lx**2, 'k.-')
semilogx(tau,ma.masked_invalid(Le2['barotropic_full'])[:,yrange].mean(axis=1)[idx_bump] / B.Lx**2, 'k^--')
xlabel(r'$\tau_0$')
ylabel(r'$L_{eq}^2$')

subplot(222)
plot(tau,ma.masked_invalid(Le2['barotropic_mean'])[:,yrange].mean(axis=1)[idx_bump] / B.Lx**2, 'ko-')
plot(tau,ma.masked_invalid(Le2['barotropic_full'])[:,yrange].mean(axis=1)[idx_bump] /
     ma.masked_invalid(Le2['barotropic_mean'])[:,yrange].mean(axis=1)[idx_bump] , 'k*--')
xlabel(r'$\tau_0$')
ylabel(r'$L_{eq}^2$')
title('Barotropic')

subplot(223)
plot(tau,ma.masked_invalid(Le2['surface_full'])[:,yrange].mean(axis=1)[idx_flat] / B.Lx**2, 'k.-')
plot(tau,ma.masked_invalid(Le2['surface_full'])[:,yrange].mean(axis=1)[idx_bump] / B.Lx**2, 'k^--')
xlabel(r'$\tau_0$')
ylabel(r'$L_{eq}^2$')

subplot(224)
plot(tau,ma.masked_invalid(Le2['surface_mean'])[:,yrange].mean(axis=1)[idx_bump] / B.Lx**2, 'ko-')
plot(tau,ma.masked_invalid(Le2['surface_full'])[:,yrange].mean(axis=1)[idx_bump] /
     ma.masked_invalid(Le2['barotropic_mean'])[:,yrange].mean(axis=1)[idx_bump] , 'k*--')
xlabel(r'$\tau_0$')
ylabel(r'$L_{eq}^2 / L_x^2$')
title('Surface')
tight_layout()

# one panel figures
rcParams['font.size'] = 8
figure(figsize=(3.25,5))
clf()
semilogx(tau,ma.masked_invalid(Le2['barotropic_full'])[:,yrange].mean(axis=1)[idx_flat] / B.Lx**2, 'k.-')
semilogx(tau,ma.masked_invalid(Le2['barotropic_full'])[:,yrange].mean(axis=1)[idx_bump] / B.Lx**2, 'k^--')
semilogx(tau,ma.masked_invalid(Le2['barotropic_mean'])[:,yrange].mean(axis=1)[idx_bump] / B.Lx**2, 'r^-')
semilogx(tau,ma.masked_invalid(Le2['barotropic_full'])[:,yrange].mean(axis=1)[idx_bump] /
     ma.masked_invalid(Le2['barotropic_mean'])[:,yrange].mean(axis=1)[idx_bump] , 'b^-')
xlabel(r'$\tau_0$')
ylabel(r'$L_{eq}^2$')
title(r'$\Theta$ - Equivalent Length')
tight_layout()
grid()
legend(['flat','ridge','ridge (mean)','ridge (full/mean)'], loc='upper left')
savefig('figures/Leq_vs_tau.pdf')

# figure(figsize=(6.5,3.5))
# subplot(121)
# plot(Y, ma.masked_invalid(Le2_bt_m[0])/B.Lx**2, 'k-')
# plot(Y, ma.masked_invalid(Le2_bt_m[1])/B.Lx**2, 'r-')
# xlabel(r'$Y_{eq}$ (km)')
# ylabel(r'$\overline{L}_{eq}^2 / L_x^2$ (km$^2$)')
# grid()
# title('Mean Flow Equivalent Lenth')
# subplot(122)
# plot(Y, ma.masked_invalid(Le2_bt[0].mean(axis=0))/B.Lx**2, 'k-')
# plot(Y, ma.masked_invalid(Le2_bt[1].mean(axis=0))/B.Lx**2, 'r-')
# plot(Y, ma.masked_invalid(
#         Le2_bt_m[1]/B.Lx**2 * Le2_bt[0].mean(axis=0))/B.Lx**2, 'b-')
# xlabel(r'$Y_{eq}$ (km)')
# ylabel(r'$L_{eq}^2 / L_x^2$')
# grid()
# legend(['flat','ridge','flat (scaled)'])
# title('Full Equivalent Lenth')
 
