from pylab import *
import os
import bump_channel

run_name = 'taux2000_rb0110_bump'
base_dir = os.path.join(os.environ['D'],'projects','drag_strat')
b = bump_channel.BumpChannel(base_dir, run_name)

B = b.g * b.tAlpha * b.rdmds('Ttave')
P = b.rdmds('Phhytave') * b.rho0

X = linspace(0,2000.,b.Nx)
Y = linspace(0,2000.,b.Ny)

Pref = P.mean(axis=2).mean(axis=1)
Pa = P - Pref[:,newaxis,newaxis]
Pfull = P - b.rho0 * b.g * b.zc[:,newaxis,newaxis]

Pb = ma.masked_array(b.value_at_bottom(P), P.mask[0])
Pab = ma.masked_array(b.value_at_bottom(Pa), P.mask[0])
Pfb = ma.masked_array(b.value_at_bottom(Pfull), P.mask[0])
H = ma.masked_array(b.depth, P.mask[0])
dHdx = b.ddx_cgrid_centered(H)
Fd = (Pab * dHdx).mean(axis=1)

dbar = 1e-4
rcParams['font.size'] = 8
close('all')
figure(figsize=(6.5,5))
subplot(221)
pc1=pcolormesh(X,Y,H, rasterized=True)
title('H (m)')
xticks([]); yticks([]);
subplot(222)
pc2=pcolormesh(X,Y,Pfb * dbar, rasterized=True)
#contour(X,Y,H, 3, colors='k')
title(r'$p_b$ (dBar)')
xticks([]); yticks([]);
subplot(223)
pc3=pcolormesh(X,Y,Pb * dbar, rasterized=True)
#contour(X,Y,H, 3, colors='k')
title(r'$p_b - g\rho_0 z$ (dBar)')
xticks([]); yticks([]);
subplot(224)
pc4=pcolormesh(X,Y,Pab * dbar, rasterized=True)
#contour(X,Y,H, 3, colors='k')
title(r'$p_b - p_{0b}$ (dBar)')
xticks([]); yticks([]);
tight_layout()
subplot(221)
colorbar(pc1)# orientation='horizontal')
subplot(222)
colorbar(pc2)#, orientation='horizontal')
subplot(223)
colorbar(pc3)#, orientation='horizontal')
subplot(224)
colorbar(pc4)#, orientation='horizontal')

show()

savefig('figures/bottom_pressure.pdf')

j = 200

fig=figure(figsize=(6.5,2.2))

#ax1 = fig.add_subplot(111)
ax1 = fig.add_axes([0.1,0.2,0.8,0.7])
ax1.plot(Y, -H[j] + b.H, 'b')
ax1.set_xlabel('X (km)')
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel(r"$h'$ (m)", color='b')
for tl in ax1.get_yticklabels():
    tl.set_color('b')

ax2 = ax1.twinx()
ax2.plot(Y, Pab[j] * dbar, 'k')
ax2.set_ylabel(r"$p_b'$ (dBar)", color='k')
for tl in ax2.get_yticklabels():
    tl.set_color('k')
ax2.grid()

ir = r_[j-20,j+20]
ax2.plot(Y[ir], Pab[j,ir] * dbar, 'ok')
ir = r_[j-40,j+40]
ax2.plot(Y[ir], Pab[j,ir] * dbar, '^k')
#tight_layout()

savefig('figures/pressure_section.pdf')

figure(figsize=(3.25,2.5))
for di in arange(50):
    ir = r_[j-di,j+di]
    DX = diff(Y[ir])*1000.
    hprime = H[j,ir[0]] - H[j,j]
    plot(2*di*5, -diff(Pab[j,ir]) * hprime / DX , 'k.' )
plot([0,2*di*5], -ones(2)*Fd[j], 'k--')
xlabel('$\Delta$x (km)')
ylabel(r'(N m$^{-2}$)')
title('Estimated Form Drag')
tight_layout()
savefig('figures/form_drag_estimate.pdf')

