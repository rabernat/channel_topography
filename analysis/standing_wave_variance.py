from pylab import *
import os
import bump_channel
import mycolors

rcParams['figure.figsize'] = [6, 3.5]
base_dir = os.path.join(os.environ['D'],'projects','drag_strat')
figdir = './figures_variance/'

r = 'taux2000_rb0110_bump'
b = bump_channel.BumpChannel(base_dir, r)

T = b.rdmds('Ttave')
Tw = b.cgrid_to_ugrid(T) # on the U grid
Ts = b.cgrid_to_vgrid(T) # on the V grid
Tt = b.cgrid_to_wgrid(T) # on the W grid
U = b.rdmds('Uveltave')
V = b.rdmds('Vveltave')
W = b.rdmds('Wveltave')
P = b.rdmds('Phhytave')
Tbar = T.mean(axis=2)
Ubar = U.mean(axis=2)
Vbar = V.mean(axis=2)
Wbar = W.mean(axis=2)

Twave = T - Tbar[:,:,newaxis]
Uwave = U - Ubar[:,:,newaxis]
Vwave = V - Vbar[:,:,newaxis]
Wwave = W - Wbar[:,:,newaxis]
Twaves = b.cgrid_to_vgrid(Twave)
Twavet = b.cgrid_to_wgrid(Twave)

UT = b.rdmds('UTtave')
VT = b.rdmds('VTtave')
WT = b.rdmds('WTtave')
UpTp = UT - U*Tw
VpTp = VT - V*Ts
WpTp = WT - W*Tt

### zonal variance budget ###
# 0: <vbar> * <ddy( T^2 / 2)>
# 1: <wbar> * <ddz( T^2 / 2)>
# 2: <vdag Tdag> ddy(Tbar)
# 3: <wdag Tdag> ddz(Tbar)
# 4: <Tdag ddx( UpTp ) >
# 5: <Tdag ddy( VpTp ) >
# 6: <Tdag ddz( WpTp ) >
leg = [ r'$\partial_y \langle \bar{v} \theta^{\dagger^2} / 2 \rangle$',
        r'$\partial_z \langle \bar{w} \theta^{\dagger^2} / 2 \rangle$',
        r'$\langle v^\dagger \theta^\dagger \rangle \partial_y \langle \theta \rangle$',
        r'$\langle w^\dagger \theta^\dagger \rangle \partial_z \langle \theta \rangle$',
        r"$\langle \theta^\dagger \rangle \partial_x \langle u'\theta' \rangle$",
        r"$\langle \theta^\dagger \rangle \partial_y \langle v'\theta' \rangle$",
        r"$\langle \theta^\dagger \rangle \partial_z \langle w'\theta' \rangle$",]

TVB = empty((7,b.Nz,b.Ny))

TVB[0] = b.ddy_vgrid_to_cgrid( (V * Twaves**2).mean(axis=2)/2 )
TVB[1] = b.ddz_wgrid_to_cgrid( (W * Twavet**2).mean(axis=2)/2 )
TVB[2] = b.vgrid_to_cgrid(Vwave * Twaves).mean(axis=2) * b.ddy_cgrid_centered(Tbar)
TVB[3] = b.wgrid_to_cgrid(Wwave * Twavet).mean(axis=2) * b.ddz_cgrid_centered(Tbar)
TVB[4] = ( Twave * b.ddx_ugrid_to_cgrid(UpTp) ).mean(axis=2)
TVB[5] = ( Twave * b.ddy_vgrid_to_cgrid(VpTp) ).mean(axis=2)
TVB[6] = ( Twave * b.ddz_wgrid_to_cgrid(WpTp) ).mean(axis=2)

Tybar = b.ddy_cgrid_centered(Tbar)
Tzbar = b.ddz_cgrid_centered(Tbar)

# transients
TVB_T = empty((2,b.Nz,b.Ny))
TVB_T[0] = VpTp.mean(axis=2) * Tybar
TVB_T[1] = WpTp.mean(axis=2) * Tzbar

# have to mask near the boundaries
TVB = ma.masked_array(TVB,zeros(TVB.shape))
TVB.mask[:,:,:5] = True
TVB.mask[:,:,-5:] = True

TVB_T = ma.masked_array(TVB_T,zeros(TVB_T.shape))
TVB_T.mask[:,:,:5] = True
TVB_T.mask[:,:,-5:] = True

#Tybar = ma.masked_array(Tybar, abs(Tybar)<1e-7)

div_UpTp = b.ddx_ugrid_to_cgrid(UpTp) + b.ddy_vgrid_to_cgrid(VpTp) + b.ddz_wgrid_to_cgrid(WpTp)
T0 = b.average_vertical(T)

close('all')
mylim = array([-1,1])*5e-8
Tlev = arange(0,8)
figure(figsize=(6,7))
myx = [500,1000,1500]
myy = [-1000,-500,0]
for n in arange(len(leg)):
    subplot(4,2,n+1)
    pcolormesh(b.yc/1e3, b.zc, TVB[n], cmap=get_cmap('posneg'), rasterized=True)
    clim(mylim)
    #colorbar()
    contour(b.yc/1e3, b.zc, Tbar, Tlev, colors='k')
    ylim([-1500,0])
    if n==6:
        xticks(myx); yticks(myy);
        xlabel('Y (km)'); ylabel('Z (m)')
    else:
        xticks(myx,[]); yticks(myy,[])
    title(leg[n])

subplot(4,2,8)
pcolormesh(b.yc/1e3, b.zc, TVB.sum(axis=0), cmap=get_cmap('posneg'), rasterized=True)
clim(mylim)
#colorbar()
contour(b.yc/1e3, b.zc, Tbar, Tlev, colors='k')
ylim([-1500,0])
xticks(myx,[]); yticks(myy,[])
title('Sum LHS')
tight_layout()
savefig(figdir + 'wave_variance_%s.pdf' % r)

figure(figsize=(6,2))
subplot(1,2,1)
pcolormesh(b.yc/1e3, b.zc, TVB[2:4].sum(axis=0), cmap=get_cmap('posneg'), rasterized=True)
clim(mylim)
contour(b.yc/1e3, b.zc, Tbar, Tlev, colors='k')
ylim([-1500,0])
xticks(myx); yticks(myy);
xlabel('Y (km)'); ylabel('Z (m)')
title(r'$\langle \mathbf{u}^\dagger \theta^\dagger \rangle \cdot \nabla \langle \theta \rangle$')
subplot(1,2,2)
pcolormesh(b.yc/1e3, b.zc, TVB[4:-1].sum(axis=0), cmap=get_cmap('posneg'), rasterized=True)
clim(mylim)
ylim([-1500,0])
contour(b.yc/1e3, b.zc, Tbar, Tlev, colors='k')
xticks(myx,[]); yticks(myy),[];
title(r"$\langle \theta^\dagger \rangle \nabla \cdot \langle \mathbf{u}'\theta' \rangle$")
tight_layout()
savefig(figdir + 'wave_variance_approx_%s.pdf' % r)


figure(figsize=(6,2.5))
subplot(1,2,1)
pcolormesh(b.xc, b.yc, b.average_vertical(Twave * div_UpTp), cmap=get_cmap('posneg'), rasterized=True )
clim(mylim)
contour(b.xc, b.yc, T0, 8, colors='k')
xticks([]); yticks([])
xlabel('X'); ylabel('Y')
title(r"$H^{-1}\int \theta^\dagger \nabla \cdot (\mathbf{u}'\theta') dz$")
subplot(1,2,2)
pcolormesh(b.xc, b.yc, 1e-7*b.average_vertical(Twave**2), cmap=get_cmap('posneg'), rasterized=True )
clim(mylim)
contour(b.xc, b.yc, T0, 8, colors='k')
xticks([]); yticks([])
xlabel('X'); ylabel('Y')
title(r"$H^{-1}\int \gamma \theta^{\dagger 2} dz$")
tight_layout()
savefig(figdir + 'wave_damping_%s.pdf' % r)

y = b.yc / 1e3
Qfac = b.rho0 * 3994 * b.Lx / 1e12 # gives PW
rcParams['font.size'] = 8
rcParams['legend.fontsize'] = 7
figure(figsize=(3.5,5))
subplot(211)
myleg = [r'$\mathcal{H}_{SE}$', r'$\mathcal{H}_{SE}^w$',
         r'$\mathcal{H}_{SE}^{tc}$', r'$\mathcal{H}_{SE}^{te}$' ]
# plot(y, Qfac * b.integrate_vertical(TVB[2] / Tybar), 'k', linewidth=2)
# plot(y, Qfac * -b.integrate_vertical(TVB[3] / Tybar), 'b', )
# plot(y, Qfac * -b.integrate_vertical(TVB[0:2].sum(axis=0) / Tybar), 'g', )
# plot(y, Qfac * -b.integrate_vertical(TVB[4:7].sum(axis=0) / Tybar), 'r', )
# plot(y, Qfac * -b.integrate_vertical(TVB[r_[0:2,3,4:7]].sum(axis=0) / Tybar), 'k--', )
plot(y, b.average_vertical(TVB[2]) * 1e8, 'k', linewidth=2)
plot(y, -b.average_vertical(TVB[3]) * 1e8, 'b', )
plot(y, -b.average_vertical(TVB[0:2].sum(axis=0)) * 1e8, 'g', )
plot(y, -b.average_vertical(TVB[4:7].sum(axis=0)) * 1e8, 'r', )
#plot(y, -b.average_vertical(TVB[r_[0:2,3,4:7]].sum(axis=0)) * 1e8, 'k--', )
plot(y, -b.average_vertical(TVB[r_[0:7]].sum(axis=0)) * 1e8, '0.7', )
#legend(myleg, loc='upper left')
legend([leg[2], r'$-$'+leg[3], r"$-\nabla \cdot \langle \mathbf{u}^\dagger \theta^{\dagger 2} /2 \rangle$",
        r"$-\langle \theta^\dagger \rangle \nabla \cdot \langle \mathbf{u}'\theta' \rangle$"],
        loc='upper left')
xlim([0,2000]); #ylim([-120,120])
ylim([-1,1])
#title(r'$\mathcal{H}_{SE}$ Components')
title('Standing Eddy Variance Budget')
#xlabel('y (km)')
ylabel(r'10$^{-8}$ K$^{2}$ s$^{-1}$')
#tight_layout()
grid()
#savefig('figures/paper/standing_variance_components.pdf')

subplot(212)
myleg = [r'$\mathcal{H}_{SE}$', r'$\mathcal{H}_{SE}^w$',
         r'$\mathcal{H}_{SE}^{tc}$', r'$\mathcal{H}_{SE}^{te}$' ]

plot(y, b.average_vertical(TVB_T[0]) * 1e8, 'k', linewidth=2)
plot(y, -b.average_vertical(TVB_T[1]) * 1e8, 'b' )
plot(y, b.average_vertical(TVB[4:7].sum(axis=0)) * 1e8, 'r' )
plot(y, (-b.average_vertical(TVB[4:7].sum(axis=0)) + b.average_vertical(TVB_T.sum(axis=0))) * 1e8, '0.7' )
#legend(myleg, loc='upper left')
legend([r"$\bar{v'\theta'} \langle \theta \rangle_y$", r"$-\bar{w'\theta'} \langle \theta \rangle_z$",
        r"$\langle \theta^\dagger \rangle \nabla \cdot \langle \mathbf{u}'\theta' \rangle$"],
        loc='upper left')
xlim([0,2000]); #ylim([-120,120])
ylim([-1,1])
#title(r'$\mathcal{H}_{SE}$ Components')
title('Transient Eddy Variance Budget')
xlabel('y (km)')
ylabel(r'10$^{-8}$ K$^{2}$ s$^{-1}$')
tight_layout()
grid()
savefig('figures/paper/variance_components_both.pdf')