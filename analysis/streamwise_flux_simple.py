from pylab import *
import os
import bump_channel
import mycolors

rcParams['figure.figsize'] = [6, 3.5]
rcParams['legend.fontsize'] = 7
rcParams['font.size'] = 8

base_dir = os.path.join(os.environ['D'],'DATASTORE.RPA','projects','drag_strat')

b = bump_channel.BumpChannel(base_dir, 'taux2000_rb0110_bump')
#b = bump_channel.BumpChannel(base_dir, 'taux2000_rb0110_flat')

T = b.rdmds('Ttave')
Tw = b.cgrid_to_ugrid(T) # on the U grid
Ts = b.cgrid_to_vgrid(T) # on the V grid
V = b.rdmds('Vveltave')
U = b.rdmds('Uveltave')
UT = b.rdmds('UTtave')
VT = b.rdmds('VTtave')

# vertically integrated fluxes
THETA = b.average_vertical(T)
THETAw = b.cgrid_to_ugrid(THETA)
THETAs = b.cgrid_to_ugrid(THETA)

# full time-mean flux
divTflux_m = b.ddx_ugrid_to_cgrid(UT) + b.ddy_vgrid_to_cgrid(VT)
# eddy flux
divTflux_e = b.ddx_ugrid_to_cgrid(UT - U*Tw) + b.ddy_vgrid_to_cgrid(VT - V*Ts)
# vertical integrals
# total
divTflux_m_vint = b.integrate_vertical(divTflux_m)
# eddy flux
divTflux_e_vint = b.integrate_vertical(divTflux_e)

# streamwise coord
Nt = 200
avg='THETA'
T0 = THETA
tmean = T0.mean(axis=1)
t = linspace(tmean.min(), tmean.max(), Nt)

# output arrays
# vertical structure for eddy flux
Tflux_e_t_z = zeros((b.Nz,Nt))
# depth integrated
Tflux_e_t = zeros(Nt)
Tflux_m_t = zeros(Nt)

grad_T0 = (b.ddx_cgrid_centered(T0)**2 + b.ddy_cgrid_centered(T0)**2)**0.5
area = zeros(Nt)
aint_grad_T0 = zeros(Nt)
aint_T_gradT0 = zeros((b.Nz,Nt))

# divergence method
tbnds = hstack([-inf, 0.5*(t[1:]+t[:-1]), inf])
for n in arange(Nt):
    idx = T0<=tbnds[n+1]
    area[n] = b.rac[idx].sum()
    aint_grad_T0[n] = (b.rac*grad_T0)[idx].sum()
    Tflux_m_t[n] = (b.rac*divTflux_m_vint)[idx].sum()
    Tflux_e_t[n] = (b.rac*divTflux_e_vint)[idx].sum()
    for k in arange(b.Nz):
        aint_T_gradT0[k,n] = (b.rac*T[k]*grad_T0)[idx].sum()
        Tflux_e_t_z[k,n] = (b.rac*divTflux_e[k])[idx].sum()

Y0 = area / b.Lx
Y00 = 0.5*(Y0[1:] + Y0[:-1])
#Y0 = interp(t,T0.mean(axis=1), b.yc)
# I can't believe this worked...
Tbar_T0 = -diff(aint_T_gradT0) / (t[1]-t[2]) / (-diff(aint_grad_T0)/(t[1]-t[2]))[newaxis,:]

# vertical integrals
Qfac = b.rho0 * 3994

H = Qfac*Tflux_m_t
H_e = Qfac*Tflux_e_t
H_m = H - H_e

deltaTheta = 8.
#Hek_aprx = Qfac * b.tau0 * sin(pi*Y0/b.Ly) / abs(b.f) / b.rho0 * deltaTheta * b.yc/b.Ly


# figures
# copied from heat_transport_simple

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
plot(Y0/1e3, H/1e12, 'k')
plot(Y0/1e3, H_m/1e12, 'c')
plot(Y0/1e3, H_e/1e12, 'm')
#plot(Y0/1e3, Hek_aprx/1e12, 'c:')
#xlabel('Y (km)')
title(r'Cross-$\Theta$ Heat Transport : Ridge')
legend([r'$\mathcal{H}^\Theta$', r'$\mathcal{H}_{Ek}^\Theta$', r'$\mathcal{H}_{eddy}^\Theta$'], loc='upper left')
ylabel(r'(TW)')
ylim([-120,120])
xticks(Ylabs,[])

ax2 = plt.subplot2grid((3,1),(1, 0),rowspan=2)
cf=contourf(Y0/1e3 , b.zc, Tflux_e_t_z*1e3/b.Lx, VTlevs*1e3, cmap=get_cmap('posneg'), extend='both')
c=contour(Y00 /1e3 , b.zc, Tbar_T0, Tlevs, colors='0.5')
clabel(c,[0.5,], fmt='%1.1f')
title(r"$\langle v_g \theta \rangle$ (10$^3$ Kms$^{-1}$)")
xticks(array(Ylabs), Ylabs)
xlabel(r'Y$_{eq}$ (km)');
ylim([-2000,0])
ylabel('Z (m)')
#tight_layout()
#cb=colorbar(cf, cax=axes([0.91,0.2,0.01,0.7]), ticks=VTticks)
cb = colorbar(cf, orientation='horizontal', ticks=VTticks*1e3)
tight_layout()
savefig('figures/paper/VT_flat-test.pdf')  

