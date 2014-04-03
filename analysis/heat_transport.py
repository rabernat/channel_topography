from pylab import *
import bump_channel
#import mycolors
import os

base_dir = os.path.join(os.environ['D'],'DATASTORE.RPA','projects','drag_strat')

b_bump = bump_channel.BumpChannel(base_dir, 'taux2000_rb0110_bump')
b_flat = bump_channel.BumpChannel(base_dir, 'taux2000_rb0110_flat')

# k indices
kr = r_[0,4,8,12,15,17:b_flat.Nz:2]

tits = ['flat', 'ridge']

for b in [b_bump, b_flat]:
    b.load_default_fields()
    b.VTwave = (b.Vwave * b.cgrid_to_vgrid(b.Twave)).mean(axis=2)
    b.WTwave = (b.Wwave * b.cgrid_to_wgrid(b.Twave)).mean(axis=2)
    b.VTg = b.VpTpbar + b.VTwave
    b.WTg = b.WpTpbar + b.WTwave
    
deltaTheta = 8.
Qfac = b.rho0 * 3994 * b.Lx
for b in [b_bump, b_flat]:
    b.D = -2 * b.integrate_vertical( b.zc[:,newaxis] * b.Tbar) /  b.integrate_vertical(b.Tbar)[-2]
    b.Kg = - b.integrate_vertical(b.VTg) * b.Ly / (b.D * deltaTheta) 
    b.Hek = Qfac * b.integrate_vertical(b.Tbar * b.Vbar)
    b.Hg = Qfac * b.integrate_vertical(b.VTg)
    b.H_TE = Qfac * b.integrate_vertical(b.VpTpbar)
    b.H_SE = Qfac * b.integrate_vertical(b.VTwave)
Hek_aprx = Qfac * b.tau0 * sin(pi*b.yc/b.Ly) / abs(b.f) / b.rho0 * deltaTheta * b.yc/b.Ly



close('all')
rcParams['font.size'] = 8.
rcParams['figure.figsize'] = [7, 2.5]
rc('figure.subplot', left=0.1, right=0.89, bottom=0.17, top=0.88, wspace=0.12)

Tlevs = arange(0,7,0.5)
Ylabs = [0,500,1000,1500,2000]
VTlevs = 1e-3*(arange(-20,20,2)+1)
VTticks = 1e-3*(arange(-16,17,4))
figure(1)
ascale = 1e-7
aspace = 20
clf()
n=0
for b in [b_flat, b_bump]:
    subplot(1,2,n+1)
    cf=contourf(b.yg , b.zc, b.VTg, VTlevs, cmap=get_cmap('posneg'), extend='both')
    clim([-0.02,0.02])
    quiver(b.yc[::aspace] , b.zc[kr], b.VTg[kr][:,::aspace], b.WTg[kr][:,::aspace],
        angles='xy', scale_units='xy', scale=ascale)
    c=contour(b.yc , b.zc, b.Tbar, Tlevs, colors='0.5')
    clabel(c,[0.5,], fmt='%1.1f')
    title(r"$\langle v_g \theta \rangle$ (Kms$^{-1}$) : %s" % tits[n])
    xticks(array(Ylabs)*1e3, Ylabs)
    xlabel('Y (km)');
    ylim([-2000,0])
    if n==0:
        ylabel('Z (m)')
        myyticks = yticks()
    else:
        yticks(myyticks[0],[])
        gca().fill_between(b.yc[r_[0,-1]], -array([2058.0,2058.0]),-array([b.H,b.H]),
            edgecolor='none', facecolor='0.7', alpha=0.3)       
    n += 1
cb=colorbar(cf, cax=axes([0.91,0.2,0.01,0.7]), ticks=VTticks)
savefig('figures/paper/VTg.pdf')  

rcParams['legend.fontsize'] = 7
figure(3, figsize=(7,2))
clf()
n=0
for b in [b_flat, b_bump]:
    subplot(1,2,n+1)
    plot(b.yc/1e3, b.Hek/1e12, 'k')
    plot(b.yc/1e3, b.Hg/1e12, 'c')
    plot(b.yc/1e3, (b.Hek + b.Hg)/1e12, '0.5')
    plot(b.yc/1e3, Hek_aprx/1e12, 'k:')
    xlabel('Y (km)')
    title('Heat Transport : %s' % tits[n])
    if n==0:
        myyticks = yticks()
        legend([r'$H_{Ek}$', r'$H_g$', r'$\mathcal{H}$'], loc='upper left')
        ylabel(r'$\mathcal{H}$ (TW)')
    else:
        yticks(myyticks[0],[])
    ylim([-120,120])
    n+=1
#tight_layout()
savefig('figures/paper/H.pdf')

rc('figure.subplot', left=0.2, right=0.93, bottom=0.12, top=0.95,
                    hspace=0.25, wspace=0.12)
figure(2, figsize=(3.5,6))
clf()
b = b_bump
for n in arange(2):
    if n==0:
        vt = b_bump.VTwave
        wt = b_bump.WTwave
        myname = r'$\langle v^\dagger \theta^\dagger \rangle$'
    else:
        vt = b_bump.VpTpbar
        wt = b_bump.WpTpbar
        myname =  r"$\langle \bar{v' \theta'} \rangle$"
    subplot(3,1,n+2)
    cf=contourf(b.yg , b.zc, vt, VTlevs, cmap=get_cmap('posneg'), extend='both')
    clim([-0.02,0.02])
    #quiver(b.yc[::aspace] , b.zc[kr], vt[kr][:,::aspace], wt[kr][:,::aspace],
    #    angles='xy', scale_units='xy', scale=ascale)
    c=contour(b.yc , b.zc, b.Tbar, Tlevs, colors='0.5')
    clabel(c,[0.5,], fmt='%1.1f')
    title(r"%s (Kms$^{-1}$) : ridge" % myname)
    ylabel('Z (m)')
    xlim([0,2000e3])
    ylim([-2000,0])
    if n==0:
        xticks(array(Ylabs)*1e3, [])
    else:
        xticks(array(Ylabs)*1e3, Ylabs)
        xlabel('Y (km)')
    gca().fill_between(b.yc[r_[0,-1]], -array([2058.0,2058.0]),-array([b.H,b.H]),
        edgecolor='none', facecolor='0.7', alpha=0.3)       
cb=colorbar(cf, cax=axes([0.05,0.04,0.9,0.01]), ticks=VTticks[::2], orientation='horizontal')
#savefig('figures/paper/VTridge.pdf')

rcParams['legend.fontsize'] = 7
#figure(4, figsize=(3.25,2))
#clf()
subplot(3,1,1)
b = b_bump
plot(b.yc, b.Hek/1e12, 'k')
#plot(b.yc, b.Hg/1e12, 'c')
plot(b.yc, b.H_SE/1e12, 'c')
plot(b.yc, b.H_TE/1e12, 'm')
plot(b.yc, (b.Hek + b.Hg)/1e12, '0.5')
plot(b.yc, Hek_aprx/1e12, 'k:')
#xlabel('Y (km)')
xticks(array(Ylabs)*1e3, [])
xlim([0,2000e3])
title('Meridional Heat Transport : %s' % tits[n])
#legend([r'$H_{Ek}$', r'$H_g$', r'$H_{SE}$', r'$H_{TE}$', r'$\mathcal{H}$'], loc='upper left')
legend([r'$H_{Ek}$', r'$H_{SE}$', r'$H_{TE}$'], loc='upper left')
ylabel(r'$\mathcal{H}$ (TW)')
ylim([-120,120])
#tight_layout()
savefig('figures/paper/VTridge.pdf')


