from pylab import *
d = load('data/GAV_data.npz')

#idx_flat = arange(13,N)
idx_wind = arange(7)
idx_flat = arange(13,13+7)

tau = d['GAV_tau0'][idx_wind]
Sv = 1e6

rcParams['font.size'] = 8
close('all')
figure(figsize=(3.1,6.2))
clf()
s1=subplot(211)
plot( tau, d['GAV_Utrans'][idx_flat] / Sv, 'ko-')
plot( tau, d['GAV_Utrans'][idx_wind] / Sv, 'k^--')
plot( tau, d['GAV_Utrans_BT'][idx_flat] / Sv, 'bo-')
plot( tau, d['GAV_Utrans_BT'][idx_wind] / Sv, 'b^--')
plot( tau, d['GAV_Utrans_TW'][idx_flat] / Sv, 'ro-')
plot( tau, d['GAV_Utrans_TW'][idx_wind] / Sv, 'r^--')
legend(['Total (flat)','Total (ridge)',
        'BT (flat)','BT (ridge)',
        'TW (flat)','TW (ridge)',], loc='upper left')
ylim([200,2000])
title('Zonal Transport')
ylabel('Sv')

s2=subplot(212)
plot( tau, d['GAV_Utrans'][idx_flat] / Sv, 'ko-')
plot( tau, d['GAV_Utrans'][idx_wind] / Sv, 'k^--')
plot( tau, d['GAV_Utrans_BT'][idx_flat] / Sv, 'bo-')
plot( tau, d['GAV_Utrans_BT'][idx_wind] / Sv, 'b^--')
plot( tau, d['GAV_Utrans_TW'][idx_flat] / Sv, 'ro-')
plot( tau, d['GAV_Utrans_TW'][idx_wind] / Sv, 'r^--')
ylim([-20,200])
ylabel('Sv')
xlabel(r'$\tau_0$ (N m$^{-2}$)')
s1.set_xticks(s2.get_xticks())
s1.set_xticklabels([])

tight_layout()
savefig('figures/transport.pdf')

rcParams['legend.fontsize'] = 7
figure(figsize=(3.2,3.0))
loglog( tau, d['GAV_Utrans'][idx_flat] / Sv, 'ko-')
loglog( tau, d['GAV_Utrans'][idx_wind] / Sv, 'k^--')
#loglog( tau, d['GAV_Utrans_BT'][idx_flat] / Sv, 'bo-')
#loglog( tau, d['GAV_Utrans_BT'][idx_wind] / Sv, 'b^--')
loglog( tau, d['GAV_Utrans_TW'][idx_flat] / Sv, 'ro-')
loglog( tau, d['GAV_Utrans_TW'][idx_wind] / Sv, 'r^--')
ylim([3e1,3e3])
grid()
grid(which='minor')
legend(['Total (flat)','Total (ridge)',
        'TW (flat)','TW (ridge)',], loc='upper left')
xlabel(r'$\tau_0$ (N m$^{-2}$)')
ylabel('Transport (Sv)')
title('Transport Components vs. Wind')
tight_layout()
savefig('figures/transport_loglog.pdf')
