from pylab import *
import bump_channel
#import mycolors
import os

base_dir = os.path.join(os.environ['D'],'DATASTORE.RPA','projects','drag_strat')

tau0 = 0
suf = 'flatshort'
b = bump_channel.BumpChannel(base_dir, 'taux%04d_rb0110_%s' % (tau0*10000, suf))
T = b.rdmds('Ttave')
Tbar = T.mean(axis=2)
Psi_bar = b.get_psi_bar()
Psi_iso = b.get_psi_iso()
Psi_iso_z = b.get_psi_iso_z()
G = b.get_layers_computer().G.squeeze()
Psi_eddy = Psi_iso_z - Psi_bar

Psi_levs = (arange(-10,10)+0.5)/100.

rcParams['font.size'] = 8
rc('figure.subplot', left=0.07, right=0.85, bottom=0.15, top=0.85, wspace=0.3, hspace=0.4)
close('all')
figure(figsize=(6.5,2.7))
#ax1 = subplot2grid([1,20],[0,0],colspan=8)
ax1 = subplot(121)
contourf(b.yc / 1e3, G, Psi_iso/1e6, Psi_levs, cmap=get_cmap('bwr'))
xlabel('Y (km)'); ylabel(r'$\theta$ ($^\circ$ C)')
title(r'$\Psi$ (Sv) - Isopycnal')
#ax2 = subplot2grid([1,20],[0,9],colspan=8)
ax2 = subplot(122)
cf=contourf(b.yc / 1e3, b.zc, Psi_iso_z/1e6, Psi_levs, cmap=get_cmap('bwr'))
c=contour(b.yc / 1e3, b.zc, Tbar, arange(0,8,0.5), colors='k')
clabel(c, fmt='%2.1f')
xlabel('Y (km)'); ylabel('Z (m)')
title(r'$\Psi$ (Sv) - Depth')
#cb_ax = subplot2grid([1,20],[0,19])
cb_ax = axes([0.9,0.15,0.02,0.7])
colorbar(cf, ticks=arange(-0.08,0.09,0.02), cax=cb_ax)
#cb_ax.set_visible(False)
savefig('figures/nowind/Psi.pdf')

figure(figsize=(6.5,2.7))
ax1 = subplot(121)
cf=contourf(b.yc / 1e3, b.zc, Psi_bar/1e6, Psi_levs, cmap=get_cmap('bwr'))
c=contour(b.yc / 1e3, b.zc, Tbar, arange(0,8,0.5), colors='k')
clabel(c, fmt='%2.1f')
xlabel('Y (km)'); ylabel('Z (m)')
title(r'$\overline{\Psi}$ (Sv)')
ax1 = subplot(122)
cf=contourf(b.yc / 1e3, b.zc, Psi_eddy/1e6, Psi_levs, cmap=get_cmap('bwr'))
c=contour(b.yc / 1e3, b.zc, Tbar, arange(0,8,0.5), colors='k')
clabel(c, fmt='%2.1f')
xlabel('Y (km)'); ylabel('Z (m)')
title(r'$\Psi^\ast$ (Sv)')
#cb_ax = subplot2grid([1,20],[0,19])
cb_ax = axes([0.9,0.15,0.02,0.7])
colorbar(cf, ticks=arange(-0.08,0.09,0.02), cax=cb_ax)
#cb_ax.set_visible(False)
savefig('figures/nowind/Psi_eddy.pdf')

show()