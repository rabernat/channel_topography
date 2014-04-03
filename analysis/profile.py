from pylab import *
import os
import bump_channel

run_name = 'taux2000_rb0110_bump'
base_dir = os.path.join(os.environ['D'],'projects','drag_strat')
b = bump_channel.BumpChannel(base_dir, run_name)

B = b.g * b.tAlpha * b.rdmds('Ttave')
P = b.rdmds('Phhytave')

Bbar = B.mean(axis=2)
Bz = b.ddz_cgrid_centered(Bbar)
By = b.ddy_cgrid_centered(Bbar)
Qy = b.beta - b.f0 * b.ddz_cgrid_centered(By/Bz)
#Blevs = arange(0,8.1,0.5)*b.g*b.tAlpha
#clf()
#contourf(b.yc,b.zc,Qy/b.beta, arange(-10,10)+0.5, cmap=get_cmap('bwr'), extend='both')
#colorbar()
#contour(b.yc,b.zc,Bbar,Blevs,colors='k')
#show()

j = 200
figure(figsize=(6.5,2.7))
rcParams['font.size'] = 8
subplot(131)
plot(Bz[:,j]/1e-5,b.zc, 'k-')
xlabel(r'$N^2$ (10$^{-5}$ s$^{-1}$)')
ylabel('z (m)')
locs,labels = yticks()
grid()
subplot(132)
plot(-b.f0**-1*By[:,j]/1e-5,b.zc, 'k-')
xlabel(r'$U_z$ (10$^{-5}$ s$^{-1}$)')
yticks(locs,[])
grid()
subplot(133)
plot(Qy[:,j]/b.beta,b.zc, 'k-')
xlim([-10,10])
xlabel(r'$Q_y / \beta$')
yticks(locs,[])
grid()
tight_layout()

savefig('figures/runs/profile_%s.pdf' % run_name)

# side view
X = linspace(0,2000,400)
eta = b.rdmds('Etatave')
clf()
subplot(121)
contourf(X,X,B[0],arange(0,8.1,0.5)*b.g*b.tAlpha, extend='both')
contour(X,X, eta+1e6, 9, colors='k')
xlabel('x (km)'); ylabel('y (km)')
title("Surface Buoyancy and Pressure")
subplot(122)
contourf(X,b.zc,B[:,j,:].filled(0.), arange(0,5.1,0.5)*b.g*b.tAlpha, extend='both')
clim([0,5*b.g*b.tAlpha])
#colorbar()
contour(X,b.zc,P[:,j,:] - P.mean(axis=2)[:,j][:,newaxis], 11, colors='k')
fill_between(X,-b.depth[j,:],-b.H, color='w')
xlabel('x (km)')
ylabel('z (m)')
title('Cross-Section (y=1000 km)')
tight_layout()

savefig('figures/runs/wave_%s.pdf' % run_name)

show()