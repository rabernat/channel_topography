from pylab import *
import bump_channel
import os
import mycolors
import MITgcmdata.poisson

# plot settings
figdir = './figures/'

# the list of runs
mytau = array([0.0125, 0.0250, 0.05, 0.1, 0.2, 0.4, 0.8])
runs = ['flat', 'bump', 'bumplong']
run_taus = {'flat': mytau, 'bump': mytau, 'bumplong': array([0.2,])}

base_dir = os.path.join(os.environ['D'],'DATASTORE.RPA','projects','drag_strat')

Nx = 400
N = len(mytau)+1

cross_stream_flux_all = zeros((N,Nx))
cross_stream_grad_all = zeros((N,Nx))
EKE_all = zeros((N,Nx))

# length of streamwise coord
#Ns = 800
# use same streamwise coordinate for all runs
# some will have longer streamlines than others
Ns = 800
S0 = linspace(0,5.5e6,Ns)

fields = ['flux','grad','EKE','h']
# big data array
data_S = dict()
data_X = dict()
for k in fields:
    data_S[k] = ma.masked_array(zeros((N,Ns)),zeros((N,Ns)))
    data_X[k] = zeros((N,Nx*2))

def avg_cross_stream(q,s,S0):
    # s should be a masked array
    Ns = len(S0)
    qsum = zeros(Ns)
    area = zeros(Ns)
    for n in arange(Ns):
        fullmask = s.mask | (s>S0[n])
        area[n] = (1-fullmask).sum()
        qsum[n] = sum(ma.masked_array(q,fullmask))
    #subplot(121); plot(area); subplot(122); plot(qsum)
    Dqsum = diff(qsum)
    Darea = diff(area)
    qS = Dqsum/Darea
    # interpolate nans
    mask = isnan(qS)
    qS[mask] = interp(flatnonzero(mask), flatnonzero(~mask), qS[~mask])
    return ma.masked_array(hstack([qS, 0]), S0>s.max())


n=-1
for tau in mytau:
    if tau==0.2:
        myruns = ['bump','bumplong']
    else:
        myruns = ['bump',]
    for run in myruns:
        n+=1

        r = 'taux%04d_rb0110_%s' % (tau*10000, run)
        try:
            b.reinit(r)
        except NameError:
            b = bump_channel.BumpChannel(base_dir, r)
        sol = MITgcmdata.poisson.solver(b)
        
        T = b.rdmds('Ttave')
        U = b.rdmds('Uveltave')
        V = b.rdmds('Vveltave')
        UT = b.rdmds('UTtave')
        VT = b.rdmds('VTtave')
        EKE = 0.5*b.average_vertical((b.rdmds('UUtave')-U**2+b.rdmds('VVtave')-V**2))
        X,Y = meshgrid(b.xc,b.yc)
        
        UpTp = UT - U*b.cgrid_to_ugrid(T)
        VpTp = VT - V*b.cgrid_to_vgrid(T)
        
        div_eddy_flux = b.integrate_vertical(
                     b.ddx_ugrid_to_cgrid(UpTp) +
                     b.ddy_vgrid_to_cgrid(VpTp) )

        #while (err2**0.5 > 2):
        sol = MITgcmdata.poisson.solver(b)
        sol.solve(div_eddy_flux.copy())
        UpTp_div, VpTp_div = sol.get_vectors()
    
        err2 = mean(( (b.ddx_ugrid_to_cgrid(UpTp_div)
               +b.ddy_vgrid_to_cgrid(VpTp_div))[2:]
               -div_eddy_flux[2:] )**2) / mean(div_eddy_flux**2)
           
        print 'Normalized error: %g' % err2**0.5
        print 'Solvers error %g' % sol.mynorm
        print 'Congergence ratio %g' % (sol.residuals[-1]/sol.residuals[0])**(1.0/len(sol.residuals))

        T0 = b.average_vertical(T)
        T0x = b.ddx_cgrid_centered(T0)
        T0y = b.ddy_cgrid_centered(T0)
        gradT = (T0x**2 + T0y**2)**(0.5)

        ### FIND MAXIMUM CONTOUR!
        Nt = 100
        t = linspace(T0.mean(axis=1)[1],T0.mean(axis=1)[-1],Nt)
        transport = zeros(Nt)
        tbnds = hstack([-inf, 0.5*(t[1:]+t[:-1]), inf])
        for m in arange(Nt):
            idx = T0<=tbnds[m+1]
            transport[m] = (b.rac*UpTp_div)[idx].sum() 
        mmax = argmax(abs(transport))
        
        useMax = False
        if useMax:
            T00 = t[mmax]
            Tbnds = [T00-0.1,T00+0.1]
        else:
            Tbnds = T0.mean(axis=1)[r_[b.Ny/2, b.Ny*2/3]]
            T00 = Tbnds.mean()        
        mask = (T0<Tbnds[0])|(T0>Tbnds[1])
        
        # project everything normal
        # divergent component only
        Tflux_div = b.ugrid_to_cgrid(UpTp_div) * T0x / gradT + b.vgrid_to_cgrid(VpTp_div) * T0y / gradT
                
        # along-stream coordinate
        c=contour(b.xc, b.yc,T0, [T00], colors='k')
        p = c.collections[0].get_paths()[0].vertices
        Np = p.shape[0]
        dx,dy = diff(p[:,0]),diff(p[:,1])
        s = hstack([0,cumsum((dx**2 + dy**2)**0.5)])
        S = empty((b.Ny,b.Nx)).flatten()
        for i in find(~mask):
            k = argmin( (X.ravel()[i]-p[:,0])**2 + (Y.ravel()[i]-p[:,1])**2 )
            S[i] = s[k]
        S.shape = X.shape
        S = ma.masked_array(S, mask)
        #S0 = linspace(0,s[-1],Ns)
        
        # average along stream
        Smask = (S0>S.max())
        data_S['h'][n] = avg_cross_stream(b.H-b.depth,S,S0)
        data_S['flux'][n] = avg_cross_stream(Tflux_div,S,S0)
        data_S['grad'][n] = avg_cross_stream(gradT,S,S0)
        data_S['EKE'][n] = avg_cross_stream(EKE,S,S0)
        
        #data_X['flux'][n] = ma.masked_array(Tflux_div,mask).mean(axis=0)
        #data_X['grad'][n] = ma.masked_array(gradT,mask).mean(axis=0)
        #data_X['EKE'][n] = ma.masked_array(EKE,mask).mean(axis=0)
        #cross_stream_flux_all[n] = Tflux_div.mean(axis=0)
        #cross_stream_grad_all[n] = gradT.mean(axis=0)
        #EKE_all[n] = EKE.mean(axis=0)
            
Xmoments = dict()    
X2moments = dict()    

colors = ['r','b','g','c','m','y','k']
for varname in ['cross_stream_flux','cross_stream_grad','EKE']:
    #figure()
    for vartype in ['all','max']:
        f= varname+'_'+vartype
        v = locals()[f]
        Xmoments[f] = zeros(N)
        X2moments[f] = zeros(N)
        for n in arange(N):
            if vartype=='all':
                ls=colors[n]+'-'
            else:
                ls=colors[n]+'--'
            #plot(b.xc/1e3, v[n]/max(abs(v[n])),ls)
            #X = sum(v[n]*b.xc)/sum(v[n])
            i0 = argmax(abs(v)[n])
            X = sum(roll(v[n],b.Nx/2-i0)*b.xc)/sum(v[n])+b.xc[i0]-b.Lx/2
            i = argmin(abs(b.xc - X))
            X2 = sum( roll(v[n],b.Nx/2-i)*(b.xc-b.Lx/2)**2) / sum(v[n])
            Xmoments[f][n] = X
            X2moments[f][n] = X2

np.savez('data/Xmoments.npz', **Xmoments)
np.savez('data/X2moments.npz', **X2moments)
        

for f in fields:
    data_S[f] = ma.masked_invalid(data_S[f])
tau_labels = hstack([mytau[:5], '0.2 (double)',mytau[5:]])

fluxfac = mytau / b.rho0 / b.f0 * 4

close('all')
rcParams['font.size']=7

x = arange(N+1)
xm = 0.5*(x[1:]+x[:-1])
clevs = arange(0,101,5)/100.
figure(figsize=(6.5,5.0))
subplot(221)
pcolormesh(S0 / 1e3, x, data_S['grad'], rasterized=True)
contour(S0 / 1e3, xm, ma.masked_invalid(data_S['h']), [250,500,750], colors='k')
xlabel('S (km)')
yticks(xm, tau_labels)
ylabel(r'$\tau_0$ (N m$^{-2}$)')
title(r'$|\nabla \Theta|$')
grid()

subplot(222)
pcolormesh(S0 / 1e3, x, -data_S['flux'], rasterized=True)
#clim([0,2])
contour(S0 / 1e3, xm, ma.masked_invalid(data_S['h']), [250,500,750], colors='k')
xlabel('S (km)')
yticks(xm,tau_labels)
ylabel(r'$\tau_0$ (N m$^{-2}$)')
title(r"$\mathbf{F} \cdot \hat n$")
grid()

subplot(223)
pcolormesh(S0 / 1e3, x, data_S['grad']/data_S['grad'].max(axis=1)[:,newaxis], rasterized=True)
clim([0,1])
contour(S0 / 1e3, xm, ma.masked_invalid(data_S['h']), [250,500,750], colors='k')
xlabel('S (km)')
yticks(xm, tau_labels)
ylabel(r'$\tau_0$ (N m$^{-2}$)')
title(r'Normalized $|\nabla \Theta|$')
grid()

subplot(224)
pcolormesh(S0 / 1e3, x, data_S['flux']/data_S['flux'].min(axis=1)[:,newaxis], rasterized=True)
#contourf(S0 / 1e3, x, data_S['flux']/fluxfac[:,newaxis], clevs*2, extend='both')
clim([0,1])
contour(S0 / 1e3, xm, ma.masked_invalid(data_S['h']), [250,500,750], colors='k')
xlabel('S (km)')
yticks(xm,tau_labels)
ylabel(r'$\tau_0$ (N m$^{-2}$)')
title(r"Normalized $\mathbf{F} \cdot \hat n$")
grid()
tight_layout()
savefig('figures/normalized_cross_stream_structure_S_new.pdf')

figure(figsize=(6.5,2.5))
subplot(121)
loglog(mytau, (Xmoments['cross_stream_flux_max']-b.Lx/2) , 'k.-')
loglog(mytau, (Xmoments['cross_stream_grad_max']-b.Lx/2) , 'b.-')
loglog(mytau, (Xmoments['EKE_max']-b.Lx/2) , 'r.-')
loglog(mytau, (Xmoments['cross_stream_flux_all']-b.Lx/2) , 'k--')
loglog(mytau, (Xmoments['cross_stream_grad_all']-b.Lx/2) , 'b--')
loglog(mytau, (Xmoments['EKE_all']-b.Lx/2) , 'r--')
ylim([2e4,7e5])
xlabel(r'$\tau_0$ (N m$^2$)')
title('Location of Peak Flux')
ylabel(r'X (m)')
legend(['flux','grad','EKE'],loc='upper left')
subplot(122)
loglog(mytau, X2moments['cross_stream_flux_max']**0.5 , 'k.-')
loglog(mytau, X2moments['cross_stream_grad_max']**0.5 , 'b.-')
loglog(mytau, X2moments['EKE_max']**0.5 , 'r.-')
loglog(mytau, X2moments['cross_stream_flux_all']**0.5 , 'k--')
loglog(mytau, X2moments['cross_stream_grad_all']**0.5 , 'b--')
loglog(mytau, X2moments['EKE_all']**0.5 , 'r--')
title('Width of Flux Region')
#ylabel('(km)')
ylim([3e5,8e5])
xlabel(r'$\tau_0$ (N m$^2$)')
#ylim([0,250])
tight_layout()
savefig('figures/cross_stream_flux_location.pdf')
