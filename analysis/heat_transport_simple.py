from pylab import *
import bump_channel
import mycolors
import os

base_dir = os.path.join(os.environ['D'],'DATASTORE.RPA','projects','drag_strat')

b = bump_channel.BumpChannel(base_dir, 'taux2000_rb0110_flat')

b.load_default_fields()
b.VTwave = (b.Vwave * b.cgrid_to_vgrid(b.Twave)).mean(axis=2)
b.WTwave = (b.Wwave * b.cgrid_to_wgrid(b.Twave)).mean(axis=2)
b.VTg = b.VpTpbar + b.VTwave
b.WTg = b.WpTpbar + b.WTwave
    
deltaTheta = 8.
Qfac = b.rho0 * 3994 * b.Lx

b.D = -2 * b.integrate_vertical( b.zc[:,newaxis] * b.Tbar) /  b.integrate_vertical(b.Tbar)[-2]
b.Kg = - b.integrate_vertical(b.VTg) * b.Ly / (b.D * deltaTheta) 
b.Hek = Qfac * b.integrate_vertical(b.Tbar * b.Vbar)
b.Hg = Qfac * b.integrate_vertical(b.VTg)
b.H_TE = Qfac * b.integrate_vertical(b.VpTpbar)
b.H_SE = Qfac * b.integrate_vertical(b.VTwave)
Hek_aprx = Qfac * b.tau0 * sin(pi*b.yc/b.Ly) / abs(b.f) / b.rho0 * deltaTheta * b.yc/b.Ly



Tlevs = arange(0,7,0.5)
Ylabs = [0,500,1000,1500,2000]
VTlevs = 1e-3*(arange(-20,20,2)+1)
VTticks = 1e-3*(arange(-16,17,4))
ascale = 1e-7
aspace = 20

close('all')
rcParams['font.size'] = 8.
rcParams['legend.fontsize'] = 7
#rc('figure.subplot', left=0.1, right=0.89, bottom=0.17, top=0.88, wspace=0.12)
figure(figsize=[3.5, 4.5])

ax1 = subplot2grid((3,1),(0, 0))
plot(b.yc/1e3, (b.Hek + b.Hg)/1e12, 'k')
plot(b.yc/1e3, b.Hek/1e12, 'c')
plot(b.yc/1e3, b.Hg/1e12, 'm')
plot(b.yc/1e3, Hek_aprx/1e12, 'c:')
#xlabel('Y (km)')
title('Meridional Heat Transport : Flat')
legend([r'$\mathcal{H}$', r'$\mathcal{H}_{Ek}$', r'$\mathcal{H}_{eddy}$'], loc='upper left')
ylabel(r'(TW)')
ylim([-120,120])
xticks(Ylabs,[])

ax2 = plt.subplot2grid((3,1),(1, 0),rowspan=2)
cf=contourf(b.yg/1e3 , b.zc, b.VTg * 1e3, VTlevs*1e3, cmap=get_cmap('posneg'), extend='both')
#clim([-0.02,0.02])
#quiver(b.yc[::aspace] , b.zc[kr], b.VTg[kr][:,::aspace], b.WTg[kr][:,::aspace],
#    angles='xy', scale_units='xy', scale=ascale)
c=contour(b.yc /1e3 , b.zc, b.Tbar, Tlevs, colors='0.5')
clabel(c,[0.5,], fmt='%1.1f')
title(r"$\langle v_g \theta \rangle$ (10$^3$ Kms$^{-1}$)")
xticks(array(Ylabs), Ylabs)
xlabel('Y (km)');
ylim([-2000,0])
ylabel('Z (m)')
#tight_layout()
#cb=colorbar(cf, cax=axes([0.91,0.2,0.01,0.7]), ticks=VTticks)
cb = colorbar(cf, orientation='horizontal', ticks=VTticks*1e3)
tight_layout()
savefig('figures/paper/VT_flat.pdf')  

