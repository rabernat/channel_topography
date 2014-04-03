from pylab import *
import os
import bump_channel

# parameters
H = 2985 # depth
Lx = 2000e3 # zonal width
dx = 5e3 # grid spacing

rho0 = 1000.
N2s = 1e-5 # near-surface stratification
f0 = -0.8e-4 # coriolis param
d = 500 # strat depth
alpha = exp(-H/d)

U = 0.01 # barotropic velocity
Gam = 0.1 # baroclinic velocity
lam = 2592000.**-1 # relaxation rate

x = dx*arange(-200,200)
h0 = 927. # bump height
sig= 15.*dx
h = h0 * exp(-x**2/ (sig**2))
hx= -2*x/(sig**2)*h
hxx = ( 4*x**2 / sig**4 - 2 / sig**2 ) * h

# amplitude of B
B0 = N2s * d * (U + alpha*Gam) / (f0*lam)
B = B0 * hx
A = Gam**-1 * B0 * ( U*hx + lam*h)

# 2d grid
z = linspace(-H,0,100)
xx,zz = meshgrid(x, z)
# streamfunction
Psi = A[newaxis,:] + B[newaxis,:]*exp(zz/d)
# pressure (in decibars)
p = f0 * Psi
# buoyancy
buoy = f0 * B[newaxis,:] / d * exp(zz/d)
t = buoy / ( 2e-4 * 9.8)
# merid. velocity
v = B0 * ( Gam**-1*(U*hxx + lam*hx)[newaxis,:] + hxx[newaxis,:]*exp(zz/d))

# contour levels
plev = arange(-1.7,1.8,0.2); ptick = arange(-1.5,1.6,0.5)
Tlev = arange(-1.3,1.4,0.2); Ttick = arange(-1.,1.1,0.5)
vlev = arange(-42.5,42.6,5.); vtick = arange(-40,41,20)

# settings
close('all')
rc('figure.subplot', left=0.1, right=0.95, bottom=0.17, top=0.88, wspace=0.12)
rc('font', size=7)
figsize=(6.5,2.2)

# figures
X = x/1000.
figure(1, figsize=figsize); clf()
subplot(131)
contourf(X, z, p, plev, extend='both', cmap=get_cmap('bwr'))
xlabel('X (km)'); ylabel('depth (m)'); title(r'p/$\rho_0$ (m$^2$/s$^2$)');
gca().fill_between(X,-H+h,-H, color='k')
colorbar(ticks=ptick)
subplot(132)
contourf(X, z, t, Tlev, extend='both', cmap=get_cmap('bwr'))
xlabel('X (km)');  title('T (deg.)')
gca().fill_between(X,-H+h,-H, color='k')
gca().set_yticklabels([])
colorbar(ticks=Ttick)
subplot(133)
contourf(X, z, v*100., vlev, extend='both', cmap=get_cmap('bwr'))
xlabel('X (km)');  title('v (cm/s)')
gca().fill_between(X,-H+h,-H, color='k')
gca().set_yticklabels([])
colorbar(ticks=vtick)
savefig('figures/wave/wave_theory.pdf')



# simulation stuff

run_name = 'taux2000_rb0110_bump'
base_dir = os.path.join(os.environ['D'],'projects','drag_strat')
b = bump_channel.BumpChannel(base_dir, run_name)
outdir = '/Volumes/scratch/tmp/output_for_paraview/' + run_name

T = b.rdmds('Ttave')
V = b.rdmds('Vveltave')
P = b.rdmds('Phhytave')

Twave = T - T.mean(axis=2)[:,:,newaxis]
Vwave = V - V.mean(axis=2)[:,:,newaxis]
Panom = P - P[:,:,0][:,:,newaxis]

# figures
j = 200
figure(2, figsize=figsize); clf()
subplot(131)
contourf(X, b.zc, Panom[:,j,:], plev, extend='both', cmap=get_cmap('bwr'))
xlabel('X (km)'); ylabel('depth (m)'); title(r'p/$\rho_0$ (m$^2$/s$^2$)');
gca().fill_between(X,-H+h,-H, color='k')
colorbar(ticks=ptick)
subplot(132)
contourf(X, b.zc, Twave[:,j,:], Tlev, extend='both', cmap=get_cmap('bwr'))
xlabel('X (km)');  title('T (deg.)')
gca().fill_between(X,-H+h,-H, color='k')
gca().set_yticklabels([])
colorbar(ticks=Ttick)
subplot(133)
contourf(X, b.zc, Vwave[:,j,:]*100., vlev, extend='both', cmap=get_cmap('bwr'))
xlabel('X (km)');  title('v (cm/s)')
gca().fill_between(X,-H+h,-H, color='k')
gca().set_yticklabels([])
colorbar(ticks=vtick)
savefig('figures/wave/wave_simulation.pdf')


